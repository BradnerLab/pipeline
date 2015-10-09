#include "score_matrix.h"
#include "fimo_style_printer.h"
#include "fasta_reader.h"
#include "bamliquidator_regions.h"

#include <samtools/bam.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <iostream>
#include <fstream>

using namespace liquidator;

enum InputType
{
    bam_input_type,
    fasta_input_type,
    invalid_input_type
};

int process_command_line(int argc,
                         char** argv,
                         std::string& input_file_path,
                         InputType& input_type,
                         std::ifstream& motif,
                         std::array<double, AlphabetSize>& background_array,
                         std::string& region_file_path,
                         std::string& ouput_file_path)
{
    namespace po = boost::program_options;

    std::string motif_file_path, background_file_path;

    // todo: add more info to help, like this:
    //   meme style position weight matrix (pwm) file
    //   .fasta file to search for motifs
    po::options_description options("Usage: motif_liquidator [options] motif fasta|bam\noptions");
    options.add_options()
        ("help,h", "produce help message")
        ("background,b", po::value(&background_file_path), "meme style background frequency file")
        ("region,r", po::value(&region_file_path), ".bed region file for filtering bam input")
        ("output,o", po::value(&ouput_file_path), "file to write matches to.  output is fasta style for fimo input, and output is a .bam for bam input.")
    ;

    // todo: manually check if a positional argument is omitted,
    // since the po exception message describes it as a non-positional argument, which could be confusing.
    // aside: boost program_options handling of positional arguments is disappointing, but I don't expect another C++ argument parser to be much better while being as easy to package
    po::options_description hidden;
    hidden.add_options()
        ("motif", po::value(&motif_file_path)->required())
        ("fasta_or_bam", po::value(&input_file_path)->required())
    ;

    po::options_description combined;
    combined.add(options).add(hidden);

    po::positional_options_description positional;
    positional.add("motif", 1);
    positional.add("fasta_or_bam", 1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(combined).positional(positional).run(), vm);

        if (vm.count("help"))
        {
            std::cerr << options << std::endl;
            return 1;
        }

        po::notify(vm);

        const std::string input_extension = boost::filesystem::extension(input_file_path);
        if (input_extension == ".bam")
        {
            input_type = bam_input_type;
        }
        else if (input_extension == ".fasta")
        {
            input_type = fasta_input_type;
        }
        else
        {
            std::cerr << "only .bam and .fasta extensions are supported at this time" << std::endl;
            return 1;
        }

        if (vm.count("region_file_path") && input_type != bam_input_type)
        {
            std::cerr << "only .bam input files support region filtering" << std::endl;
            return 1;
        }

        motif.open(motif_file_path);
        if (!motif)
        {
            std::cerr << "failed to open motif file " << motif_file_path << std::endl;
            return 1;
        }

        if (vm.count("background"))
        {
            std::ifstream background(background_file_path);
            if (!background)
            {
                std::cerr << "failed to open background file " << background_file_path << std::endl;
                return 1;
            }
            background_array = ScoreMatrix::read_background(background);
        }
        else
        {
            background_array = {.25, .25, .25, .25};
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << options << "\n";
        return 1;
    }

    return 0;
}

void process_fasta(const std::vector<ScoreMatrix>& matrices, const std::string& fasta_file_path, const std::string& output_file_path)
{
    std::ifstream fasta_input(fasta_file_path);
    if (!fasta_input)
    {
        throw std::runtime_error("failed to open " + fasta_file_path);
    }

    std::ofstream output(output_file_path);
    FimoStylePrinter printer(output);

    FastaReader fasta_reader(fasta_input);
    std::string sequence;
    std::string sequence_name;
    while (fasta_reader.next_read(sequence, sequence_name))
    {
        for (const auto& matrix : matrices)
        {
            matrix.score(sequence, sequence_name, printer);
        }
    }
}

class BamScorer
{
public:
    BamScorer(const std::string& bam_input_file_path,
              const std::vector<ScoreMatrix>& matrices,
              const std::string& bam_output_file_path,
              const std::string& region_file_path = "")
    :
        m_input(bam_open(bam_input_file_path.c_str(), "r")),
        m_output(0),
        m_header(bam_header_read(m_input)),
        m_index(bam_index_load(bam_input_file_path.c_str())),
        m_matrices(matrices),
        m_read(0),
        m_read_count(0),
        m_read_hit_count(0),
        m_total_hit_count(0)
    {
        if (m_input == 0 || m_header == 0 || m_index == 0)
        {
            throw std::runtime_error("failed to open " + bam_input_file_path);
        }

        if (!bam_output_file_path.empty())
        {
            m_output = bam_open(bam_output_file_path.c_str(), "w");
            bam_header_write(m_output, m_header);
        }

        if (!region_file_path.empty())
        {
            score_regions(region_file_path);
        }
        else
        {
            score_all_reads();
        }
    }

    ~BamScorer()
    {
        const auto percentage_printer = [&](const size_t hits) {
            std::cout << hits << "/" << m_read_count << " = " << 100*(double(hits)/m_read_count) << "%" << std::endl;
        };
        percentage_printer(m_read_hit_count);
        percentage_printer(m_total_hit_count);

        bam_index_destroy(m_index);
        bam_header_destroy(m_header);
        bam_close(m_input);

        if (m_output)
        {
            bam_close(m_output);
        }
    }

    void operator()(const std::string& motif_name,
                    const std::string& sequence_name,
                    size_t start,
                    size_t stop,
                    const ScoreMatrix::Score& score)
    {
        if (score.pvalue() < 0.0001)
        {
            ++m_total_hit_count;
            if (m_output)
            {
                bam_write1(m_output, m_read);
            }
        }
    }

private:
    void score_all_reads()
    {
        bam1_t* read = bam_init1();
        try
        {
            while (bam_read1(m_input, read) >= 0)
            {
                m_read = read;
                score_read();
            }
        }
        catch(...)
        {
            bam_destroy1(read);
            throw;
        }
        bam_destroy1(read);
    }

    void score_regions(const std::string& region_file_path)
    {
        for (const Region& region : parse_regions(region_file_path, "bed", 0))
        {
            // todo: don't I just need to parse_region once per chromosome to get the tid? perhaps there is a faster way to do this without parsing a whole region string?
            // todo: consider adding a util function to do this and remove duplicate code in bamliquidator.cpp
             std::stringstream coord;
            coord << region.chromosome << ':' << region.start << '-' << region.stop;

            int ref,beg,end;
            const int region_parse_rc = bam_parse_region(m_header,coord.str().c_str(), &ref, &beg, &end);
            if (region_parse_rc != 0)
            {
                std::stringstream error_msg;
                error_msg << "bam_parse_region failed with return code " << region_parse_rc;
                throw std::runtime_error(error_msg.str());
            }
            if(ref<0)
            {
                // this bam doesn't have this chromosome
                continue;
            }

            const int fetch_rc = bam_fetch(m_input, m_index, ref, beg, end, this, bam_fetch_func);
            if (fetch_rc != 0)
            {
                std::stringstream error_msg;
                error_msg << "bam_fetch failed with return code " << fetch_rc;
                throw std::runtime_error(error_msg.str());
            }
        }
    }

    void score_read()
    {
        const bam1_core_t *c = &m_read->core;
        uint8_t *s = bam1_seq(m_read);

        // [s, s+c->l_qseq) is the sequence, with two bases packed into each byte.
        // I bet we could directly search that instead of first copying into a string
        // but lets get something simple working first. An intermediate step could be
        // to search integers without using bam_nt16_rev_table (and I wouldn't have
        // to worry about the packing complexity).

        if (m_sequence.size() != size_t(c->l_qseq))
        {
            // assuming that all reads are uniform length, this will only happen once
            m_sequence = std::string(c->l_qseq, ' ');
        }
        for (int i = 0; i < c->l_qseq; ++i)
        {
            m_sequence[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
        }
        ++m_read_count;

        const size_t hit_count_before_this_read = m_total_hit_count;
        const static std::string no_sequence_name;
        for (const auto& matrix : m_matrices)
        {
            matrix.score(m_sequence, no_sequence_name, *this);
        }
        if (m_total_hit_count > hit_count_before_this_read)
        {
            ++m_read_hit_count;
        }
    }

    static int bam_fetch_func(const bam1_t* read, void* handle)
    {
        BamScorer& scorer = *static_cast<BamScorer*>(handle);
        scorer.m_read = read;
        scorer.score_read();
        return 0;
    }

private:
    bamFile m_input;
    bamFile m_output;
    bam_header_t* m_header;
    bam_index_t* m_index;
    const std::vector<ScoreMatrix>& m_matrices;
    const bam1_t* m_read;
    size_t m_read_count;
    size_t m_read_hit_count;
    size_t m_total_hit_count;
    std::string m_sequence;
};

int main(int argc, char** argv)
{
    try
    {
        std::string input_file_path, region_file_path, ouput_file_path;
        std::ifstream motif;
        InputType input_type = invalid_input_type;
        std::array<double, AlphabetSize> background;

        const int rc = process_command_line(argc, argv, input_file_path, input_type, motif, background, region_file_path, ouput_file_path);
        if ( rc ) return rc;

        std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif, background);

        if (input_type == bam_input_type)
        {
            BamScorer(input_file_path, matrices, ouput_file_path, region_file_path);
        }
        else if (input_type == fasta_input_type)
        {
            process_fasta(matrices, input_file_path, ouput_file_path);
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

/* The MIT License (MIT) 

   Copyright (c) 2015 John DiMatteo (jdimatteo@gmail.com)

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE. 
 */

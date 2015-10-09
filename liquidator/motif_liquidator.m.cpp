#include "score_matrix.h"
#include "fimo_style_printer.h"
#include "fasta_reader.h"

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
                         std::array<double, AlphabetSize>& background_array)
{
    namespace po = boost::program_options;

    std::string motif_file_path, background_file_path;

    // todo: add more info to help, like this:
    //   meme style position weight matrix (pwm) file
    //   .fasta file to search for motifs
    po::options_description options("Usage: motif_liquidator [options] motif fasta|bam\noptions");
    options.add_options()
        ("help,h", "produce help message")
        ("background,b", po::value(&background_file_path), "meme style background frequency file");

    // todo: manually check if a positional argument is omitted,
    // since the po exception message describes it as a non-positional argument, which could be confusing.
    // aside: boost program_options handling of positional arguments is disappointing, but I don't expect another C++ argument parser to be much better while being as easy to package
    po::options_description hidden;
    hidden.add_options()
        ("motif", po::value(&motif_file_path)->required())
        ("fasta_or_bam", po::value(&input_file_path)->required());

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

void process_fasta(const std::vector<ScoreMatrix>& matrices, const std::string& fasta_file_path)
{
    std::ifstream fasta_input(fasta_file_path);
    if (!fasta_input)
    {
        throw std::runtime_error("failed to open " + fasta_file_path);
    }

    FimoStylePrinter printer(std::cout);

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

void process_bam(const std::vector<ScoreMatrix>& matrices, const std::string& bam_file_path)
{
    bam1_t* read = bam_init1();
    bamFile input = bam_open(bam_file_path.c_str(), "r");
    bam_header_t* header = bam_header_read(input);

    if (read == 0 || input == 0 || header == 0)
    {
        throw std::runtime_error("failed to open " + bam_file_path);
    }

    std::string sequence;
    size_t read_count = 0;
    size_t hit_count = 0;

    auto hit_counter = [&hit_count](const std::string& motif_name,
                                    const std::string& sequence_name,
                                    size_t start,
                                    size_t stop,
                                    const ScoreMatrix::Score& score) {
        if (score.pvalue() < 0.0001)
        {
            ++hit_count;
        }
    };

    while (bam_read1(input, read) >= 0)
    {
        const bam1_core_t *c = &read->core;
        uint8_t *s = bam1_seq(read);

        // [s, s+c->l_qseq) is the sequence, with two bases packed into each byte.
        // I bet we could directly search that instead of first copying into a string
        // but lets get something simple working first. An intermediate step could be
        // to search integers without using bam_nt16_rev_table (and I wouldn't have
        // to worry about the packing complexity).

        if (sequence.size() != size_t(c->l_qseq))
        {
            // assuming that all reads are uniform length, this will only happen once
            sequence = std::string(c->l_qseq, ' ');
        }

        for (int i = 0; i < c->l_qseq; ++i)
        {
            sequence[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
        }
        ++read_count;

        const static std::string sequence_name = "???";
        for (const auto& matrix : matrices)
        {
            matrix.score(sequence, sequence_name, hit_counter);
        }
    }

    bam_header_destroy(header);
    bam_close(input);
    bam_destroy1(read);

    std::cout << hit_count << "/" << read_count << " = " << (double(hit_count)/read_count) << std::endl;
}

int main(int argc, char** argv)
{
    try
    {
        std::string input_file_path;
        std::ifstream motif;
        InputType input_type = invalid_input_type;
        std::array<double, AlphabetSize> background;

        const int rc = process_command_line(argc, argv, input_file_path, input_type, motif, background);
        if ( rc ) return rc;

        std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif, background);

        if (input_type == bam_input_type)
        {
            process_bam(matrices, input_file_path);
        }
        else if (input_type == fasta_input_type)
        {
            process_fasta(matrices, input_file_path);
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

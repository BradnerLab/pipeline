#include "bam_scorer.h"
#include "fasta_scorer.h"
#include "score_matrix.h"
#include "version.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

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
                         std::ifstream& motif_file,
                         std::array<double, AlphabetSize>& background_array,
                         std::string& region_file_path,
                         std::string& ouput_file_path,
                         std::string& motif_name,
                         BamScorer::PrintStyle& bam_print_style,
                         bool& unmapped_only)
{
    namespace po = boost::program_options;

    std::string motif_file_path, background_file_path, print_argument;

    po::options_description options(std::string("usage: motif_liquidator [options] motif fasta|bam")
                                  + "\nversion " + std::string(version) +
                                  + "\n\nScore DNA sequences against a motif, with exact same scores as MEME's FIMO."
                                  + "\nFor more info, see https://github.com/BradnerLab/pipeline/wiki/motif_liquidator"
                                  + "\n\npositional arguments:"
                                  + "\n  motif                   A MEME style position weight matrix file."
                                  + "\n  fasta|bam               A fasta file or (sorted and indexed) bam file."
                                  + "\n\noptional arguments");
    options.add_options()
        ("motif-name,m", po::value(&motif_name), "Motif name (default is all motifs in the motif file)")
        ("background,b", po::value(&background_file_path), "MEME style background frequency file.  Note that only 0-order background (single nucleotide) frequenceis are currently used (just like FIMO).  Backgrounds specified in the motif file are never used (just like default FIMO behavior).")
        ("help,h", "Display this help and exit.")
        ("output,o", po::value(&ouput_file_path), "File to write matches to. Output is fimo style for fasta input, and output is a (sorted/indexed) .bam for bam input.")
        ("region,r", po::value(&region_file_path), ".bed or .gff region file for filtering bam input.")
        ("unmapped-only,u", "Only scores unmapped reads from bam.")
        ("print,p", po::value(&print_argument), "For bams, additionally prints detailed fimo style output to stdout.  Specify '-p fimo' for fimo style output or '-p mapped-fimo' for the sequence name to include the chromosome and the for start/stop values to be the mapped positions.")
    ;

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

        // manually check if a positional argument is omitted before calling notify,
        // since the po exception message describes it as a non-positional argument, which is confusing to a user.
        if (vm.count("motif") != 1 || vm.count("fasta_or_bam") != 1)
        {
            std::cerr << "invalid positional arguments" << std::endl;
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

        motif_file.open(motif_file_path);
        if (!motif_file)
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

        if (vm.count("print"))
        {
            if (print_argument == "fimo")
            {
                bam_print_style = BamScorer::Fimo;
            } else if (print_argument == "mapped-fimo")
            {
                bam_print_style = BamScorer::MappedFimo;
            }
            else
            {
                std::cerr << "Invalid print argument '" << print_argument << "'; please specify either 'fimo' or 'mapped-fimo'." << std::endl;
                return 1;
            }
        }
        else
        {
            bam_print_style = BamScorer::None;
        }

        unmapped_only = vm.count("unmapped-only") > 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << options << "\n";
        return 1;
    }

    return 0;
}

int main(int argc, char** argv)
{
    try
    {
        std::string input_file_path, region_file_path, ouput_file_path, motif_name;
        std::ifstream motif_file;
        InputType input_type = invalid_input_type;
        BamScorer::PrintStyle bam_print_style = BamScorer::None;
        bool unmapped_only = false;
        std::array<double, AlphabetSize> background = ScoreMatrix::default_acgt_background;

        const int rc = process_command_line(argc, argv, input_file_path, input_type, motif_file, background, region_file_path, ouput_file_path, motif_name, bam_print_style, unmapped_only);
        if (rc) return rc;

        std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif_file, background, motif_name);

        if (input_type == bam_input_type)
        {
            BamScorer(input_file_path, matrices, bam_print_style, unmapped_only, ouput_file_path, region_file_path);
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

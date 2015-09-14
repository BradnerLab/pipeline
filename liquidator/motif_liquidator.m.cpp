#include "score_matrix.h"
#include "fimo_style_printer.h"
#include "fasta_reader.h"

#include <boost/program_options.hpp>

#include <string>
#include <iostream>
#include <fstream>

using namespace liquidator;

int process_command_line(int argc,
                         char** argv,
                         std::ifstream& fasta,
                         std::ifstream& motif)
{
    namespace po = boost::program_options;

    std::string fasta_file_path, motif_file_path, background_file_path;

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
        ("fasta_or_bam", po::value(&fasta_file_path)->required());

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

        fasta.open(fasta_file_path);
        if ( !fasta )
        {
            std::cerr << "failed to open fasta file " << fasta_file_path << std::endl;
            return 1;
        }

        motif.open(motif_file_path);
        if ( !motif )
        {
            std::cerr << "failed to open motif file " << motif_file_path << std::endl;
            return 1;
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

int main(int argc, char** argv)
{
    std::ifstream fasta;
    std::ifstream motif;

    // todo: read in background from an optional file, defaulting to {.25, .25, .25, .25}
    const std::array<double, AlphabetSize> background = {.256, .244, .244, .256};

    const int rc = process_command_line(argc, argv, fasta, motif);
    if ( rc ) return rc;

    std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif, background);

    FimoStylePrinter printer(std::cout);

    FastaReader fasta_reader(fasta);
    std::string sequence;
    std::string sequence_name;
    while (fasta_reader.next_read(sequence, sequence_name))
    {
        for (const auto& matrix : matrices)
        {
            matrix.score(sequence, sequence_name, printer);
        }
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

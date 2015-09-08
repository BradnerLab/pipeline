#include "score_matrix.h"
#include "fimo_style_printer.h"

#include <boost/program_options.hpp>

#include <string>
#include <iostream>
#include <fstream>

int process_command_line(int argc,
                         char** argv,
                         std::ifstream& fasta,
                         std::ifstream& motif,
                         std::ifstream& background)
{
    namespace po = boost::program_options;

    std::string fasta_file_path, motif_file_path, background_file_path;

    po::options_description desc("usage");
    desc.add_options()
        ("help,h", "produce help message")
        ("fasta,f", po::value(&fasta_file_path)->required(), ".fasta file to search for motifs")
        ("motif,m", po::value(&motif_file_path)->required(), "meme style position weight matrix (pwm) file")
        ("background,bg", po::value(&background_file_path)->required(), "meme style background frequency file");

    po::variables_map vm;

    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cerr << desc << std::endl;
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

        background.open(background_file_path);
        if ( !background )
        {
            std::cerr << "failed to open background file " << background_file_path << std::endl;
            return 1;
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << desc << "\n";
        return 1;
    }

    return 0;
}

int main(int argc, char** argv)
{
    using namespace liquidator;

    std::ifstream fasta;
    std::ifstream motif;
    std::ifstream background;

    const int rc = process_command_line(argc, argv, fasta, motif, background);
    if ( rc ) return rc;

    // todo: use background file
    std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif, {.256, .244, .244, .256});

    FimoStylePrinter printer(std::cout);

    // todo: FastaReader fasta_reader(fasta);
    std::string sequence;
    std::string sequence_name;
    char strand = '.';
    //while (fasta_reader.next_read(sequence, sequence_name, strand))
    {
        for (const auto& matrix : matrices)
        {
            // todo: consider making score function take '+', '-', '.' arg instead of bool
            if (strand == '+' || strand == '.')
            {
                matrix.score(sequence, sequence_name, true, printer);
            }
            if (strand == '-' || strand == '.')
            {
                matrix.score(sequence, sequence_name, false, printer);
            }
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

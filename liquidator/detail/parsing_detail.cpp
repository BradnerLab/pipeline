// parsing is separated into a separate .cpp file because compiling boost spirit can be slow

#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3

#include "parsing_detail.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

namespace liquidator { namespace detail {

// Input format described at http://meme.ebi.edu.au/meme/doc/meme-format.html .
std::vector<PWM> read_pwm(std::istream& input)
{
    // todo: use spirit::istream_iterator instead of copying to a string
    std::stringstream input_buffer;
    input_buffer << input.rdbuf();
    const std::string input_str = input_buffer.str();
    auto begin = input_str.begin();
    auto end = input_str.end();

    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phoenix = boost::phoenix;

    uint meme_version = 0;
    const std::string default_alphabet = "ACGT";
    std::string alphabet = default_alphabet;
    std::string strands = "+ -";

    std::vector<std::string> motif_names;

    // todo: learn how to use spirit properly so we don't need this current and vector non-sense
    std::vector<unsigned> number_of_sites_vector;
    const unsigned default_number_of_sites = 20; // match meme's default of 20
    unsigned number_of_sites = default_number_of_sites;

    std::vector<std::vector<std::vector<double>>> matrices;
    std::vector<std::vector<double>> current_matrix;
    std::vector<double> current_line;

    bool success = qi::parse(begin, end,
          qi::lit("MEME version ")
          >> qi::uint_[ phoenix::ref(meme_version) = qi::_1 ] >> *ascii::space
          >> -(qi::lit("ALPHABET= ")
               >> qi::as_string[ +qi::upper ]
                               [ phoenix::ref(alphabet) = qi::_1 ]
              )
          >> *ascii::space
          >> -(qi::lit("strands: ")
               >> qi::as_string[ +(qi::print - qi::eol) ]
                               [ phoenix::ref(strands) = qi::_1 ]
              )
          >> *ascii::space

          // we ignore motif file specified background frequencies, and I am pretty sure fimo ignores these too
          >> -(qi::lit("Background letter frequencies")
               >> *ascii::space
               >> *(qi::print - qi::eol)
               >> *ascii::space
              )
          >> *ascii::space

          >> +(qi::lit("MOTIF ")
               >> qi::as_string[ +(qi::print - ascii::space) ]
                               [ phoenix::bind([&] (const std::string& motif_name)
                                 {
                                     motif_names.push_back(motif_name);
                                 },
                                 qi::_1)
                               ]

               // ignore any alternative motif name
               >> -(ascii::blank >> +(qi::print - ascii::space))

               >> *ascii::space
               >> qi::lit("letter-probability matrix:")

               // Ignore any provided alphabet and motif lengthes
               // -- these aren't required and can be calculated from matrix anyway.
               >> -(ascii::blank >> qi::lit("alength= ") >> qi::uint_)
               >> -(ascii::blank >> qi::lit("w= ") >> qi::uint_)

               >> -(ascii::blank
                    >> qi::lit("nsites= ") >> qi::uint_[ phoenix::ref(number_of_sites) = qi::_1 ])

               // ignore any E-value
               >> -(ascii::blank >> qi::lit("E= ") >> qi::double_)

               >> *ascii::space

               // finally we parse the matrix
               >> +( +(qi::double_[ phoenix::bind([&] (double& value)
                                    {
                                        current_line.push_back(value);
                                    },
                                    qi::_1)
                                  ]

                       >> *ascii::blank
                      )
                     >> *ascii::space
                   )[ phoenix::bind([&] ()
                      {
                        current_matrix.push_back(current_line);
                        current_line.clear();
                      })
                    ]
              )[ phoenix::bind([&] ()
                 {
                    matrices.push_back(current_matrix);
                    current_matrix.clear();
                    number_of_sites_vector.push_back(number_of_sites);
                    number_of_sites = default_number_of_sites;
                 })
               ]
        );

    if (!success || end != begin || matrices.empty())
    {
        throw std::runtime_error("motif pwm parse failed, stopped at: [" + std::string(begin, end) + "]");
    }

    if (alphabet != default_alphabet)
    {
        throw std::runtime_error("only " + default_alphabet + " motif pwm alphabet supported, but found " + alphabet);
    }

    if (strands != "+ -")
    {
        throw std::runtime_error("only + - motif pwm strands supported, but found " + strands);
    }

    if (matrices.size() != motif_names.size()
     || matrices.size() != number_of_sites_vector.size())
    {
        throw std::runtime_error("motif pwm parsing error: mismatch matrices and motif name sizes");
    }

    std::vector<PWM> pwms;
    for (size_t i=0; i < matrices.size(); ++i)
    {
        const auto& matrix = matrices
                [i];

        PWM pwm { number_of_sites_vector[i], motif_names[i] };
        for (const auto& row: matrix)
        {
            if (row.size() != AlphabetSize)
            {
                throw std::runtime_error("motif pwm matrix row has wrong number of entries");
            }
            std::array<double, AlphabetSize> pwm_row;
            for (size_t k=0; k < AlphabetSize; ++k)
            {
                pwm_row[k] = row[k];
            }
            pwm.matrix.push_back(pwm_row);
        }
        pwms.push_back(pwm);
    }

    return pwms;
}

}}

/* The MIT License (MIT)

   Copyright (c) 2016 Boulder Labs (jdimatteo@boulderlabs.com)

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

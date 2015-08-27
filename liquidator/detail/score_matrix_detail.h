#ifndef LIQUIDATOR_DETAIL_SCORE_MATRIX_DETAIL_H_INCLUDED
#define LIQUIDATOR_DETAIL_SCORE_MATRIX_DETAIL_H_INCLUDED

#include "liquidator_util.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>

namespace liquidator { namespace detail {

// Input format described at http://meme.ebi.edu.au/meme/doc/meme-format.html .
struct PWM
{
    const size_t number_of_sites;
    std::string name;
    std::vector<std::array<double, AlphabetSize>> matrix;
};

struct ScaledPWM
{
    const size_t number_of_sites;
    const int min;
    const int scale;
    const int range;
    std::string name;
    std::vector<std::array<int, AlphabetSize>> matrix;
};

// Transforms PWM probability values into log pseudo-site-adjusted likelihood ratio values.
// Returns min/max ratio values.
inline std::pair<double, double>
log_adjusted_likelihood_ratio(PWM& pwm, 
                              const std::array<double, AlphabetSize>& background, 
                              double number_of_pseudo_sites=.1)
{
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    for (auto& row : pwm.matrix)
    {
        for (int i=0; i < AlphabetSize; ++i)
        {
            const double adjusted = (row[i]*pwm.number_of_sites + number_of_pseudo_sites*background[i])/(pwm.number_of_sites+number_of_pseudo_sites);
            const double ratio = adjusted / background[i];
            const double log_ratio = std::log2(ratio);
            min = std::min(min, log_ratio);
            max = std::max(max, log_ratio);
            row[i] = log_ratio;
        }
    }
    return std::make_pair(min, max);
}

inline ScaledPWM
scale(const PWM& pwm, std::pair<double, double> min_max, int range)
{
    double min = min_max.first;
    const double max = min_max.second;
    if (min == max)
    {
        min = max - 1;        
    }
    min = std::floor(min);

    const int scale = std::floor(range/(max-min));
    ScaledPWM scaled_pwm { /*number_of_sites=*/ pwm.number_of_sites,
                           /*min=*/ int(min),
                           /*scale=*/ scale,
                           /*range=*/ range,
                           /*name=*/ pwm.name };
    scaled_pwm.matrix.reserve(pwm.matrix.size());
    for (auto& row : pwm.matrix)
    {
        std::array<int, AlphabetSize> scaled_row;
        for (int alphabet_index=0; alphabet_index < AlphabetSize; ++alphabet_index)
        {
            scaled_row[alphabet_index] = std::round((row[alphabet_index] - min) * scale);  
        }
        scaled_pwm.matrix.push_back(scaled_row);
    }
    return scaled_pwm;
}

// precondition: (end-begin) <= matrix.size() && sequence.size() <= end && end >= begin
// postcondition: returns score or -1 if not scorable due to invalid alphabet char
int score(const std::vector<std::array<int, AlphabetSize>>& matrix, const std::string& sequence, size_t begin, size_t end)
{
    int score = 0;
    for (size_t position=begin, row=0; position < end; ++position, ++row)
    {
        auto column = alphabet_index(sequence[position]);
        if (column >= AlphabetSize)
        {
            return -1;
        }

        score += matrix[row][column];
    }
    return score;
}

inline std::vector<PWM> read_pwm(std::istream& input)
{
    // todo: replace below hack with spirit parsing,
    //       should throw if alphabet size isn't ACGT,
    //       should read nsites and maybe E,
    //       should read multiple matrices and return an empty vector if none (not throw),
    //       and add tests for all these changes 
    PWM pwm { /*number_of_sites =*/ 18 };
    // todo: read number_of_sites in if provided, else default to 20 like FIMO

    bool matrixFirstLineSeen = false;
    for(std::string line; getline(input, line); )
    {
        boost::trim(line);

        std::vector<std::string> split;
        boost::split(split, line, boost::is_any_of(" "), boost::token_compress_on);
        if (pwm.name.empty())
        {
            if (split.size() >= 2 && split[0] == "MOTIF")
            {
                pwm.name = split[1];
            }
        }
        else
        {
            if (matrixFirstLineSeen == false)
            {
                if (!split.empty() && split.front() == "letter-probability")
                {
                    matrixFirstLineSeen = true;
                }
            }
            else if (split.size() == AlphabetSize)
            {
                std::array<double, AlphabetSize> row;
                for (int i=0; i < AlphabetSize; ++i)
                {
                    row[i] = boost::lexical_cast<double>(split[i]);
                }
                pwm.matrix.push_back(row);
            }
            else
            {
                break;
            }
        }
    }
    if (pwm.name.empty() || pwm.matrix.empty())
    {
        throw std::runtime_error("Failed to read a motif from input");
    }
    std::vector<PWM> pwms;
    pwms.push_back(pwm);
    return pwms;
}

} }

#endif

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

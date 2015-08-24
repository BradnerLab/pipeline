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
    std::string name;
    std::vector<std::array<double, AlphabetSize>> matrix;
    size_t number_of_sites;
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

inline std::vector<PWM> read_pwm(std::istream& input)
{
    // todo: replace below hack with spirit parsing,
    //       should throw if alphabet size isn't ACGT,
    //       should read nsites and maybe E,
    //       should read multiple matrices and return an empty vector if none (not throw),
    //       and add tests for all these changes 
    PWM pwm;
    pwm.number_of_sites = 18; // todo: read this in if provided, else default to 20 like FIMO

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

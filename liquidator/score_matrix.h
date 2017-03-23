#ifndef LIQUIDATOR_SCORE_MATRIX_H_INCLUDED
#define LIQUIDATOR_SCORE_MATRIX_H_INCLUDED

#include "liquidator_util.h"

#include <array>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

namespace liquidator
{

// motif position weight matrix (pwm) for scoring sequences 
class ScoreMatrix
{
public:
    static constexpr std::array<double, AlphabetSize> default_acgt_background = {{0.281774, 0.222020, 0.228876, 0.267330}};

    // Input format described at http://meme.ebi.edu.au/meme/doc/meme-format.html .
    // Psuedo count logic described at http://meme-suite.org/doc/general-faq.html .
    static std::vector<ScoreMatrix> read(std::istream& meme_style_pwm,
                                         const std::array<double, AlphabetSize>& acgt_background = default_acgt_background,
                                         const std::string &motif_name = "",
                                         bool include_reverse_complement = true,
                                         double pseudo_sites = 0.1);

    // Background format described at http://meme.ebi.edu.au/meme/doc/bfile-format.html .
    // Note that only order 0 values are used.
    static std::array<double, AlphabetSize> read_background(std::istream& background);

    ScoreMatrix(const std::string& name,
                const std::array<double, AlphabetSize>& background,
                bool average_background_for_reverse,
                const std::vector<std::array<double, AlphabetSize>>& pwm,
                unsigned number_of_sites,
                bool is_reverse_complement = false,
                double pseudo_sites = 0.1);
 
    // Scores reference a sequence string so are intended to be used only
    // in the scope of a ScoreConsumer operator.
    class Score
    {
        public:
            // writes the matched sequence
            friend std::ostream& operator<<(std::ostream& out, const Score& score);
            
            // returns a copy of the matched sequence
            inline std::string matched_sequence() const;

            // The pvalue, or NAN if sequence was not scorable.
            // Note that NaN < x is false for any double x.
            double pvalue() const { return m_pvalue; }

            // The score, or 0 if sequence was not scorable.
            double score() const { return m_score; }

            bool is_reverse_complement() const { return m_is_reverse_complement; }

            Score(const std::string& sequence, bool is_reverse_complement, size_t begin, size_t end, double pvalue, double score);

        private:
            const std::string& m_sequence;
            const bool m_is_reverse_complement;
            const size_t m_begin;
            const size_t m_end;
            const double m_pvalue;
            const double m_score;
    };

    // See fimo_style_printer.h for example of a ScoreConsumer.
    template <typename ScoreConsumer>
    void score(const std::string& sequence, ScoreConsumer& consumer) const
    {
        for (size_t start = 1, stop = m_matrix.size(); stop <= sequence.size(); ++start, ++stop)
        {
            const Score score = score_sequence(sequence, start-1, stop);
            consumer(m_name, start, stop, score);
        }
    }

    std::string name() { return m_name; }
    size_t length() { return m_matrix.size(); }

    // Matrix value for the sequence position (row) and base letter (column).
    // Value is a log likelihood ratio, adjusted with a psuedo count and scaled.
    // Base should be ACGT/acgt -- else throws exception.
    // position is 0 based and must be < length() -- else undefined behavior.
    int value(size_t position, char base)
    {
        const auto column = alphabet_index(base);
        if ( column >= AlphabetSize ) throw std::runtime_error("Invalid base " + std::string(1, base));
        return m_matrix[position][column];
    }

    const std::string& name() const
    {
        return m_name;
    }

    bool is_reverse_complement() const
    {
        return m_is_reverse_complement;
    }

    const std::vector<std::array<unsigned, AlphabetSize>>& matrix() const
    {
        return m_matrix;
    }

    double scale() const
    {
        return m_scale;
    }

    double min_before_scaling() const
    {
        return m_min_before_scaling;
    }

    const std::vector<double>& pvalues() const
    {
        return m_pvalues;
    }

private:
    const std::string m_name;
    const bool m_is_reverse_complement;
    std::vector<std::array<unsigned, AlphabetSize>> m_matrix;
    double m_scale;
    double m_min_before_scaling;
    std::vector<double> m_pvalues;

    Score score_sequence(const std::string& sequence, size_t begin, size_t end) const;
};

inline std::ostream& operator<<(std::ostream& out, const ScoreMatrix::Score& score)
{
    if (score.m_is_reverse_complement)
    {
        if (score.m_end > score.m_begin)
        {
            for (size_t i=score.m_end-1; ; --i)
            {
                out << char(std::toupper(complement(score.m_sequence[i])));
                if ( i == score.m_begin )
                {
                    break;
                }
            }
        }
    }
    else
    {
        for (size_t i=score.m_begin; i < score.m_end; ++i)
        {
            out << char(std::toupper(score.m_sequence[i]));
        }
    }
    return out;
}

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
#endif

#ifndef LIQUIDATOR_SCORE_MATRIX_H_INCLUDED
#define LIQUIDATOR_SCORE_MATRIX_H_INCLUDED

#include "liquidator_util.h"

#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace liquidator
{

// motif position weight matrix (pwm) for scoring sequences 
class ScoreMatrix
{
public:
    // Input format described at http://meme.ebi.edu.au/meme/doc/meme-format.html .
    // Psuedo count logic described at http://meme-suite.org/doc/general-faq.html .
    static std::vector<ScoreMatrix> read(std::istream& meme_style_pwm,
                                         std::array<double, AlphabetSize> acgt_background = { .25, .25, .25, .25},
                                         double pseudo_sites = 0.1);
 
    // Scores reference the sequence string so are intended to be used only 
    // in the scope of a ScoreConsumer operator.
    class Score
    {
        public:
            // writes the matched sequence
            inline std::ostream& operator()(std::ostream& out) const;
            
            // returns a copy of the matched sequence
            inline std::string matched_sequence() const;

            // The pvalue, or NAN if sequence was not scorable.
            // Note that NaN < x is false for any double x.
            const double pvalue() { return m_pvalue; }

            // The score, or -1 if sequence was not scorable. 
            const double score() { return m_score; }

        private:
            Score(const std::string& sequence, size_t begin, size_t end, float pvalue, float score);

            Score(const Score&) = delete;

            const std::string& m_sequence;
            const size_t m_begin;
            const size_t m_end;
            const double m_pvalue;
            const double m_score;
    };

    // See fimo_style_printer.h for example of a ScoreConsumer.
    template <typename ScoreConsumer>
    void score(const std::string& sequence, const std::string& sequence_name, bool forward_strand, ScoreConsumer& consumer)
    {
        for (size_t start = 1, stop = m_matrix.size(); stop <= sequence.size(); ++start, ++stop)
        {
            const Score score = score_sequence(sequence, start-1, stop);
            consumer(sequence_name, forward_strand, start, stop, score);
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
        const int column = alphabet_index(base);
        if ( column >= AlphabetSize ) throw std::runtime_error("Invalid base " + std::string(1, base));
        return m_matrix[position][column];
    }

private:
    ScoreMatrix();

    std::string m_name;
    std::vector<std::array<int, AlphabetSize>> m_matrix;
    double m_min;

    Score score_sequence(const std::string& sequence, size_t begin, size_t end);
};

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

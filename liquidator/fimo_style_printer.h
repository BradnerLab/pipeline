#ifndef LIQUIDATOR_FIMO_STYLE_PRINTER_H_INCLUDED 
#define LIQUIDATOR_FIMO_STYLE_PRINTER_H_INCLUDED 

#include <iostream>
#include "score_matrix.h"

namespace liquidator
{

// Example of a ScoreConsumer 
class FimoStylePrinter
{
public:
    FimoStylePrinter(std::ostream& out, bool include_header = true, double threshold = 0.0001)
    :
        sequence_name(0),
        m_out(out),
        m_threshold(threshold)
    {
        if (include_header)
        {
            m_out << "#pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence" << std::endl;
        }
    }

    void operator()(const std::string& motif_name,
                    size_t start,
                    size_t stop,
                    const ScoreMatrix::Score& score)
    {
        if (score.pvalue() < m_threshold)
        {
            m_out << motif_name << '\t' 
                  << (sequence_name ? *sequence_name : std::string()) << '\t'
                  << start << '\t'
                  << stop << '\t'
                  << (score.is_reverse_complement() ? '-' : '+') << '\t';

            m_out.precision(6);
            m_out << score.score() << '\t';
            m_out.precision(3);

            m_out << score.pvalue() << '\t'
                  << '\t' // omit q-value for now
                  << score << std::endl; 
        }
    }

    const std::string* sequence_name;

private:
    std::ostream& m_out;
    const double m_threshold;
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

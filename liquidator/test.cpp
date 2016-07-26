#include "gtest/gtest.h"

#include "score_matrix.h"
#include "detail/score_matrix_detail.h"

using namespace liquidator;

const std::array<double, AlphabetSize> uniform_bg = {.25, .25, .25, .25};

TEST(ScoreMatrix, read_pwm_matrix)
{
    const std::string input_str = R"(MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.29 C 0.21 G 0.21 T 0.29

MOTIF JASPAR2014.MA0107.1 RELA

letter-probability matrix: alength= 4 w= 10 nsites= 18 E= 0
  0.000000        0.222222        0.611111        0.166667
  0.000000        0.000000        0.944444        0.055556
  0.000000        0.000000        1.000000        0.000000
  0.611111        0.000000        0.388889        0.000000
  0.555556        0.166667        0.222222        0.055556
  0.111111        0.000000        0.000000        0.888889
  0.000000        0.000000        0.000000        1.000000
  0.000000        0.111111        0.000000        0.888889
  0.000000        1.000000        0.000000        0.000000
  0.000000        1.000000        0.000000        0.000000)";

    std::istringstream ss(input_str);
    std::vector<detail::PWM> pwms = detail::read_pwm(ss);
    ASSERT_EQ(1, pwms.size());
    
    const auto& name = pwms[0].name;
    EXPECT_EQ("JASPAR2014.MA0107.1", name);

    const auto& matrix = pwms[0].matrix;
    ASSERT_EQ(10, matrix.size());
    const auto& row = matrix[0];
    ASSERT_EQ(4, row.size());

    EXPECT_FLOAT_EQ(0, matrix[0][0]);
    EXPECT_FLOAT_EQ(0.222222, matrix[0][1]);
    EXPECT_FLOAT_EQ(0.388889, matrix[3][2]);
    EXPECT_FLOAT_EQ(1, matrix[6][3]);
    EXPECT_FLOAT_EQ(1, matrix[9][1]);
}

TEST(ScoreMatrix, read_multiple_pwm)
{
    const std::string input_str = R"(MEME version 4

            ALPHABET= ACGT

            strands: + -

            Background letter frequencies
            A 0.303 C 0.183 G 0.209 T 0.306

            MOTIF crp
            letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009
             0.000000  0.176471  0.000000  0.823529
             0.000000  0.058824  0.647059  0.294118
             0.000000  0.058824  0.000000  0.941176
             0.176471  0.000000  0.764706  0.058824
             0.823529  0.058824  0.000000  0.117647
             0.294118  0.176471  0.176471  0.352941
             0.294118  0.352941  0.235294  0.117647
             0.117647  0.235294  0.352941  0.294118
             0.529412  0.000000  0.176471  0.294118
             0.058824  0.235294  0.588235  0.117647
             0.176471  0.235294  0.294118  0.294118
             0.000000  0.058824  0.117647  0.823529
             0.058824  0.882353  0.000000  0.058824
             0.764706  0.000000  0.176471  0.058824
             0.058824  0.882353  0.000000  0.058824
             0.823529  0.058824  0.058824  0.058824
             0.176471  0.411765  0.058824  0.352941
             0.411765  0.000000  0.000000  0.588235
             0.352941  0.058824  0.000000  0.588235

            MOTIF lexA
            letter-probability matrix: alength= 4 w= 18 nsites= 14 E= 3.2e-035
             0.214286  0.000000  0.000000  0.785714
             0.857143  0.000000  0.071429  0.071429
             0.000000  1.000000  0.000000  0.000000
             0.000000  0.000000  0.000000  1.000000
             0.000000  0.000000  1.000000  0.000000
             0.000000  0.000000  0.000000  1.000000
             0.857143  0.000000  0.071429  0.071429
             0.000000  0.071429  0.000000  0.928571
             0.857143  0.000000  0.071429  0.071429
             0.142857  0.000000  0.000000  0.857143
             0.571429  0.071429  0.214286  0.142857
             0.285714  0.285714  0.000000  0.428571
             1.000000  0.000000  0.000000  0.000000
             0.285714  0.214286  0.000000  0.500000
             0.428571  0.500000  0.000000  0.071429
             0.000000  1.000000  0.000000  0.000000
             1.000000  0.000000  0.000000  0.000000
             0.000000  0.000000  0.785714  0.214286
)";

    std::istringstream ss(input_str);
    std::vector<detail::PWM> pwms = detail::read_pwm(ss);
    ASSERT_EQ(2, pwms.size());

    EXPECT_EQ("crp", pwms[0].name);
    EXPECT_EQ("lexA", pwms[1].name);

    EXPECT_EQ(17, pwms[0].number_of_sites);
    EXPECT_EQ(14, pwms[1].number_of_sites);

    {
        const auto& matrix = pwms[0].matrix;
        ASSERT_EQ(19, matrix.size());

        EXPECT_FLOAT_EQ(0, matrix[0][0]);
        EXPECT_FLOAT_EQ(0.176471, matrix[0][1]);
        EXPECT_FLOAT_EQ(0.764706, matrix[3][2]);
        EXPECT_FLOAT_EQ(0.117647, matrix[6][3]);
        EXPECT_FLOAT_EQ(0.058824, matrix[18][1]);
    }

    {
        const auto& matrix = pwms[1].matrix;
        ASSERT_EQ(18, matrix.size());

        EXPECT_FLOAT_EQ(0.214286, matrix[0][0]);
        EXPECT_FLOAT_EQ(0, matrix[17][0]);
        EXPECT_FLOAT_EQ(0.214286, matrix[17][3]);
    }
}

TEST(ScoreMatrix, log_adjusted_likelihood_ratio)
{
    const unsigned number_of_sites = 18;
    detail::PWM pwm { /*number_of_sites =*/ number_of_sites };
    pwm.matrix = { {.25, .25, .25, .25},
                   {0, 0, 1, 0} };
    const double number_of_pseudo_sites=.1;

    const auto min_max = detail::log_adjusted_likelihood_ratio(pwm, uniform_bg, uniform_bg);
    EXPECT_EQ(number_of_sites, pwm.number_of_sites);
    ASSERT_EQ(2, pwm.matrix.size());
    ASSERT_EQ(4, pwm.matrix[0].size());

    const double zero = std::log2(number_of_pseudo_sites * uniform_bg[0] / ( number_of_sites + number_of_pseudo_sites)/uniform_bg[0]);
    const double one = std::log2((number_of_sites+number_of_pseudo_sites*uniform_bg[0]) / (number_of_sites + number_of_pseudo_sites)/uniform_bg[0]);
    const double quarter = 0; // matching a base for a position where all bases are equally likely scores zero points
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][0]);
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][1]);
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][2]);
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][3]);
    EXPECT_FLOAT_EQ(zero,    pwm.matrix[1][0]);
    EXPECT_FLOAT_EQ(zero,    pwm.matrix[1][1]);
    EXPECT_FLOAT_EQ(one,     pwm.matrix[1][2]);
    EXPECT_FLOAT_EQ(zero,    pwm.matrix[1][3]);

    EXPECT_FLOAT_EQ(min_max.first,  zero);
    EXPECT_FLOAT_EQ(min_max.second, one);
}

TEST(ScoreMatrix, scale)
{
    detail::PWM pwm { /*number_of_sites =*/ 10 };
    pwm.matrix = { {0, 0, 0, 0},
                   {-8, -8, 2, -8} };
    const detail::ScaledPWM scaled = detail::scale(pwm, std::make_pair(-8.0, 2.0), 30);

    const auto& matrix = scaled.matrix;
    ASSERT_EQ(2, matrix.size());

    EXPECT_EQ(10, scaled.number_of_sites);
    EXPECT_EQ(-8, scaled.min_before_scaling);
    EXPECT_EQ(3, scaled.scale); // max - min = 10, 10*3 = 30, so scale is 3
    EXPECT_EQ(30, scaled.range);

    EXPECT_EQ(24, matrix[0][0]); // (0 - -8) * 3 = 24
    EXPECT_EQ(24, matrix[0][1]);
    EXPECT_EQ(24, matrix[0][2]);
    EXPECT_EQ(24, matrix[0][3]);

    EXPECT_EQ(0, matrix[1][0]); // (-8 - -8) * 3
    EXPECT_EQ(0, matrix[1][1]);
    EXPECT_EQ(30, matrix[1][2]); // (2 - -8) * 3
    EXPECT_EQ(0, matrix[1][3]);
}

TEST(ScoreMatrix, scaled_score)
{
    const std::vector<std::array<unsigned, AlphabetSize>> matrix =
    //  A   C   G   T
    { { 24, 24, 24, 0 },
      { 0,  0,  30, 0 } };

    EXPECT_EQ(0, detail::score(matrix, "", 0, 0));
    EXPECT_EQ(0, detail::score(matrix, "AA", 0, 0));
    EXPECT_EQ(0, detail::score(matrix, "AG", 2, 2));

    EXPECT_EQ(24, detail::score(matrix, "A", 0, 1));
    EXPECT_EQ(0,  detail::score(matrix, "T", 0, 1));
    EXPECT_EQ(0, detail::score(matrix, "N", 0, 1));
    EXPECT_EQ(0, detail::score(matrix, "Z", 0, 1));

    EXPECT_EQ(24, detail::score(matrix, "AA", 0, 2));
    EXPECT_EQ(24, detail::score(matrix, "AA", 1, 2));
    EXPECT_EQ(54, detail::score(matrix, "AG", 0, 2));
    EXPECT_EQ(54, detail::score(matrix, "ag", 0, 2));
    EXPECT_EQ(54, detail::score(matrix, "AGN", 0, 2));
    EXPECT_EQ(54, detail::score(matrix, "NAGN", 1, 3));
}

TEST(ScoreMatrix, probability_distribution)
{
    using namespace detail;

    // score of 0 is 100% probable when empty matrix
    const std::vector<std::array<unsigned, AlphabetSize>> empty_matrix;
    std::vector<double> probabilities = probability_distribution(empty_matrix, uniform_bg);
    ASSERT_EQ(1, probabilities.size());
    EXPECT_FLOAT_EQ(1, probabilities[0]);

    // score of 0 is 100% probable for a zero matrix
    const std::vector<std::array<unsigned, AlphabetSize>> zero_matrix =
    //  A   C   G   T
    { { 0,  0,  0, 0 },
      { 0,  0,  0, 0 } };
    probabilities = probability_distribution(zero_matrix, uniform_bg);
    ASSERT_EQ(1, probabilities.size());
    EXPECT_FLOAT_EQ(1, probabilities[0]);

    const std::vector<std::array<unsigned, AlphabetSize>> length_one_matrix =
    //  A   C   G   T
    { { 0,  0,  1, 0 } };

    // sequence length 1 with max 1 per base can have values 0 or 1
    // value of 0 with 75% probability
    // value of 1 with 25% probability
    probabilities = probability_distribution(length_one_matrix, uniform_bg);
    ASSERT_EQ(2, probabilities.size());
    EXPECT_FLOAT_EQ(.75, probabilities[0]);
    EXPECT_FLOAT_EQ(.25, probabilities[1]);

    const std::vector<std::array<unsigned, AlphabetSize>> length_two_matrix =
    //  A   C   G   T
    { { 0,  0,  1, 1 },
      { 1,  0,  1, 0 } };
    // Scores for every possible sequence:
    // AA: 1, AC: 0, AG: 1, AT: 0
    // CA: 1, CC: 0, CG: 1, CT: 0
    // GA: 2, GC: 1, GG: 2, GT: 1
    // TA: 2, TC: 1, TG: 2, TT: 1
    // 16 sequences
    // 4 ways to score 0: 25%
    // 8 ways to score 1: 50%
    // 4 ways to score 2: 25%
    probabilities = probability_distribution(length_two_matrix, uniform_bg);
    ASSERT_EQ(3, probabilities.size());
    EXPECT_FLOAT_EQ(.25, probabilities[0]);
    EXPECT_FLOAT_EQ(.50, probabilities[1]);
    EXPECT_FLOAT_EQ(.25, probabilities[2]);
}

TEST(ScoreMatrix, pvalues)
{
    std::vector<double> a;
    detail::pdf_to_pvalues(a);
    EXPECT_TRUE(a.empty());

    a = { .1 };
    detail::pdf_to_pvalues(a);
    ASSERT_EQ(1, a.size());
    EXPECT_FLOAT_EQ(.1, a[0]);

    a = { .1, .2 };
    detail::pdf_to_pvalues(a);
    ASSERT_EQ(2, a.size());
    EXPECT_FLOAT_EQ(.2, a[1]);
    EXPECT_FLOAT_EQ(.1+.2, a[0]);

    a = { .1, .2, .3 };
    detail::pdf_to_pvalues(a);
    ASSERT_EQ(3, a.size());
    EXPECT_FLOAT_EQ(.3, a[2]);
    EXPECT_FLOAT_EQ(.3+.2, a[1]);
    EXPECT_FLOAT_EQ(.3+.2+.1, a[0]);

    a = { .9, .2, .3 };
    detail::pdf_to_pvalues(a);
    ASSERT_EQ(3, a.size());
    EXPECT_FLOAT_EQ(.3, a[2]);
    EXPECT_FLOAT_EQ(.3+.2, a[1]);
    EXPECT_FLOAT_EQ(1, a[0]);
}

TEST(ScoreMatrix, reverse_complement)
{
    std::vector<std::array<double, AlphabetSize>> matrix =
    { {1, 0, 0, 0},
      {0, 2, 0, 0},
      {0, 0, 3, 0},
      {0, 0, 0, 4},
      {1, 2, 3, 4} };
    const std::vector<std::array<double, AlphabetSize>> reverse_complement =
    { {4, 3, 2, 1},
      {4, 0, 0, 0},
      {0, 3, 0, 0},
      {0, 0, 2, 0},
      {0, 0, 0, 1} };
    detail::reverse_complement(matrix);
    EXPECT_EQ(reverse_complement, matrix);
}

TEST(ScoreMatrix, compress_sequence)
{
    std::string sequence = "AAAA";
    // Compressing ACGT sequences requires 2 bits per base pair,
    // so we can fit 4 base pairs into each uint8_t.

    // Since we need to quickly be able to get the compressed value at any start position,
    // we compress the sequence 4 times, corresponding to offset 0, 1, 2, and 3.
    // We don't need an offset 4 (or greater), since offset 4 starts at the second byte of offset 0 (and so on for larger offsets).
    std::array<std::vector<uint8_t>, 4> compressed_by_offset = detail::compress_sequence(sequence);
    for (const auto& compressed : compressed_by_offset)
    {
        ASSERT_EQ(1, compressed.size());
        EXPECT_EQ(0, compressed[0]);
    }

    sequence = "ACGTC";
    compressed_by_offset = detail::compress_sequence(sequence);

    // with no offset, ACGTA is compressed into two bytes
    ASSERT_EQ(2, compressed_by_offset[0].size());

    // first byte compresses ACGT, bits 0-1 storing A, 2-3 storing C, 4-5 storing G, 6-7 storing T
    // A is 00, C is 01, G is 10, T is 11, so
    //                                          11 10 01 00
    // high order bits are the last base pair T -^       ^- low order bits are first base pair (A)
    // 11100100 is decimal 228
    EXPECT_EQ(228, compressed_by_offset[0][0]);

    // next byte compresses the remaining C
    // 00 00 00 01
    EXPECT_EQ(1, compressed_by_offset[0][1]);

    // offset 1 starts at the position 1, just single byte needed to compress CGTC
    ASSERT_EQ(1, compressed_by_offset[1].size());

    // CGTC -> C T G C -> 01 11 10 01 -> decimal 121
    EXPECT_EQ(121, compressed_by_offset[1][0]);

    // offset 2 starts at position 2, just single byte needed to compress GTC
    ASSERT_EQ(1, compressed_by_offset[2].size());

    // GTC_ -> _ C T G -> 00 01 11 10 -> decimal 57
    //                    ^- unspecified base pair left as zero value
    EXPECT_EQ(30, compressed_by_offset[2][0]);

    // offset by 3 starts at position 3, just single byte needed to compress TC
    ASSERT_EQ(1, compressed_by_offset[3].size());

    // TC__ -> _ _ C T -> 00 00 01 11 -> 7

    // no bytes needed to compress empty string
    sequence = "";
    compressed_by_offset = detail::compress_sequence(sequence);
    for (const auto& compressed : compressed_by_offset)
    {
        ASSERT_TRUE(compressed.empty());
    }

    sequence = "ACTAGAACG";
    compressed_by_offset = detail::compress_sequence(sequence);

    // offset 0:
    // ACTA GAAC G___
    // A T C A  C A A G  _ _ _ G
    // 00110100 01000010 00000010
    // 52       66       2
    ASSERT_EQ(3, compressed_by_offset[0].size());
    EXPECT_EQ(52, compressed_by_offset[0][0]);
    EXPECT_EQ(66, compressed_by_offset[0][1]);
    EXPECT_EQ(2, compressed_by_offset[0][2]);

    // offset 1:
    // CTAG AACG
    // G A T C  G C A A
    // 10001101 10010000
    // 141      144
    ASSERT_EQ(2, compressed_by_offset[1].size());
    EXPECT_EQ(141, compressed_by_offset[1][0]);
    EXPECT_EQ(144, compressed_by_offset[1][1]);

    // offset 2:
    // TAGA ACG_
    // A G A T  _ G C A
    // 00100011 00100100
    // 35       36
    ASSERT_EQ(2, compressed_by_offset[2].size());
    EXPECT_EQ(35, compressed_by_offset[2][0]);
    EXPECT_EQ(36, compressed_by_offset[2][1]);

    // offset 3:
    // AGAA CG__
    // A A G A  _ _ G C
    // 00001000 00001001
    // 8        9
    ASSERT_EQ(2, compressed_by_offset[3].size());
    EXPECT_EQ(8, compressed_by_offset[3][0]);
    EXPECT_EQ(9, compressed_by_offset[3][1]);
}

TEST(ScoreMatrix, binary_matrix)
{
    std::vector<std::array<unsigned, AlphabetSize>> uncompressed_matrix =
    // A  C  G  T
    { {1, 3, 5, 9},   // 0
      {0, 2, 6, 0},   // 1
      {9, 0, 3, 0},   // 2
      {6, 0, 3, 4},   // 3
      {1, 5, 3, 4} }; // 4

    auto binary_matrix = detail::compress_matrix(uncompressed_matrix);

    // 5 long sequence, 5 long matrix
    std::string sequence = "ACTAG";
    auto compressed_by_offset = detail::compress_sequence(sequence);
    size_t start = 0;
    EXPECT_EQ(detail::score(uncompressed_matrix, sequence, start, start+uncompressed_matrix.size()),
              detail::score(binary_matrix, compressed_by_offset, start));

    // 6 long sequence, 5 long matrix
    sequence = "ACTAGA";
    compressed_by_offset = detail::compress_sequence(sequence);
    EXPECT_EQ(detail::score(uncompressed_matrix, sequence, start, start+uncompressed_matrix.size()),
              detail::score(binary_matrix, compressed_by_offset, start));

    // 9 long sequence, 5 long matrix
    sequence = "ACTAGAACG";
    //          1+2+0+6+3=12
    ASSERT_EQ(12, detail::score(uncompressed_matrix, sequence, start, start+uncompressed_matrix.size()));
    compressed_by_offset = detail::compress_sequence(sequence);
    //
    EXPECT_EQ(12, detail::score(binary_matrix, compressed_by_offset, start));

    // 18 long sequence, 5 long matrix
    sequence = "ACTAGAACGCGCGTCACC";
    compressed_by_offset = detail::compress_sequence(sequence);
    for (start = 0; start <= (sequence.size() - uncompressed_matrix.size()); ++start)
    {
        EXPECT_EQ(detail::score(uncompressed_matrix, sequence, start, start+uncompressed_matrix.size()),
                  detail::score(binary_matrix, compressed_by_offset, start));
    }
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
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

#include "gtest/gtest.h"

#include "score_matrix.h"
#include "detail/score_matrix_detail.h"
#include "fimo_style_printer.h"

using namespace liquidator;

const std::array<double, AlphabetSize> uniform_bg { {.25, .25, .25, .25} };

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

TEST(ScoreMatrix, read_pwm_misc)
{
    // adopted from .meme files in http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.15.tgz
    const std::string input_str = R"(MEME version 4.4
                                  MEME BigFoot DMMPMM version MAY-26-2009

                                  ALPHABET= ACGT

                                  strands: + -

                                  Background letter frequencies (from uniform background):
                                  A 0.25000 C 0.25000 G 0.25000 T 0.25000

                                  MOTIF MYB52 MYB52

                                  log-odds matrix: alength= 4 w= 8 E= 0
                                      98	    27	   -77	  -208
                                      -7	   106	  -129	   -83
                                    -127	  -123	  -283	   159
                                     189	  -422	  -297	  -306
                                     179	  -137	  -417	  -343
                                    -389	   189	  -300	  -315
                                    -149	  -224	   166	  -186
                                     -95	  -347	   166	  -205

                                  letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 4.9e-949
                                    0.493200	  0.300600	  0.147100	  0.059100
                                    0.237500	  0.519800	  0.102500	  0.140200
                                    0.103890	  0.106689	  0.035196	  0.754225
                                    0.924900	  0.013400	  0.031800	  0.029900
                                    0.865913	  0.096990	  0.013899	  0.023198
                                    0.016900	  0.923700	  0.031200	  0.028200
                                    0.089100	  0.053000	  0.788800	  0.069100
                                    0.129400	  0.022500	  0.787800	  0.060300

                                  URL http://floresta.eead.csic.es/footprintdb/index.php?db=ArabidopsisPBM:20140210&motif=MYB52

                                  MOTIF abd-A
                                  letter-probability matrix: alength= 4 w= 6 nsites= 45 E= 0
                                   0.609756097560976  0.121951219512195  0.0975609756097561  0.170731707317073
                                   0.0  0.0731707317073171  0.024390243902439  0.902439024390244
                                   0.902439024390244  0.024390243902439  0.0731707317073171  0.0
                                   0.902439024390244  0.0  0.0  0.0975609756097561
                                   0.317073170731707  0.0  0.0  0.682926829268293
                                   0.0  0.0731707317073171  0.146341463414634  0.780487804878049
                                  URL http://autosome.ru/DMMPMM/DMMPMM_18_MAR_RELEASE/comparison-DMMPMM/ABD_A.html
)";

    std::istringstream ss(input_str);
    std::vector<detail::PWM> pwms = detail::read_pwm(ss);
    ASSERT_EQ(2, pwms.size());
}

TEST(ScoreMatrix, log_adjusted_likelihood_ratio)
{
    const unsigned number_of_sites = 18;
    detail::PWM pwm { /*number_of_sites =*/ number_of_sites };
    pwm.matrix = { { { {.25, .25, .25, .25} },
                     { {0, 0, 1, 0        } } } };
    const double number_of_pseudo_sites=.1;

    const auto min_max = detail::log_adjusted_likelihood_ratio(pwm, uniform_bg, uniform_bg, number_of_pseudo_sites);
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
    pwm.matrix = { { { {0, 0, 0, 0   } },
                     { {-8, -8, 2, -8} } } };
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
    { { { { 24, 24, 24, 0 } },
        { { 0,  0,  30, 0 } } } };

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
    {{ {{ 0,  0,  0, 0 }},
       {{ 0,  0,  0, 0 }} }};
    probabilities = probability_distribution(zero_matrix, uniform_bg);
    ASSERT_EQ(1, probabilities.size());
    EXPECT_FLOAT_EQ(1, probabilities[0]);

    const std::vector<std::array<unsigned, AlphabetSize>> length_one_matrix =
    //  A   C   G   T
    {{ {{ 0,  0,  1, 0 }} }};

    // sequence length 1 with max 1 per base can have values 0 or 1
    // value of 0 with 75% probability
    // value of 1 with 25% probability
    probabilities = probability_distribution(length_one_matrix, uniform_bg);
    ASSERT_EQ(2, probabilities.size());
    EXPECT_FLOAT_EQ(.75, probabilities[0]);
    EXPECT_FLOAT_EQ(.25, probabilities[1]);

    const std::vector<std::array<unsigned, AlphabetSize>> length_two_matrix =
    //  A   C   G   T
    {{ {{ 0,  0,  1, 1 }},
       {{ 1,  0,  1, 0 }} }};
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
    {{ {{1, 0, 0, 0}},
       {{0, 2, 0, 0}},
       {{0, 0, 3, 0}},
       {{0, 0, 0, 4}},
       {{1, 2, 3, 4}} }};
    const std::vector<std::array<double, AlphabetSize>> reverse_complement =
    {{ {{4, 3, 2, 1}},
       {{4, 0, 0, 0}},
       {{0, 3, 0, 0}},
       {{0, 0, 2, 0}},
       {{0, 0, 0, 1}} }};
    detail::reverse_complement(matrix);
    EXPECT_EQ(reverse_complement, matrix);
}

TEST(ScoreMatrix, read_arguments)
{
    const std::string motif = R"(MEME version 4

        ALPHABET= ACGT

        strands: + -

        Background letter frequencies (from file `../HOCOMOCO/bkg.txt'):
        A 0.29182 C 0.20818 G 0.20818 T 0.29182
        MOTIF AP2D_f1 AP2D

        letter-probability matrix: alength= 4 w= 14 nsites= 5 E= 0
          0.200000	  0.400000	  0.200000	  0.200000
          0.000000	  0.000000	  1.000000	  0.000000
          0.000000	  0.600000	  0.200000	  0.200000
          0.000000	  1.000000	  0.000000	  0.000000
          0.000000	  0.600000	  0.000000	  0.400000
          0.000000	  0.000000	  1.000000	  0.000000
          0.400000	  0.200000	  0.400000	  0.000000
          0.000000	  0.000000	  1.000000	  0.000000
          0.000000	  0.000000	  1.000000	  0.000000
          0.000000	  0.800000	  0.200000	  0.000000
          0.200000	  0.200000	  0.400000	  0.200000
          0.000000	  0.800000	  0.200000	  0.000000
          0.000000	  0.200000	  0.800000	  0.000000
          0.200000	  0.200000	  0.000000	  0.600000

        MOTIF AHR_si AHR

        letter-probability matrix: alength= 4 w= 9 nsites= 154 E= 0
          0.266234	  0.116883	  0.363636	  0.253247
          0.071429	  0.077922	  0.227273	  0.623377
          0.142857	  0.285714	  0.136364	  0.435065
          0.019481	  0.006494	  0.948052	  0.025974
          0.006494	  0.974026	  0.006494	  0.012987
          0.019481	  0.006494	  0.967532	  0.006494
          0.000000	  0.019481	  0.006494	  0.974026
          0.000000	  0.000000	  1.000000	  0.000000
          0.279221	  0.435065	  0.103896	  0.181818
    )";
    {
        std::stringstream motif_ss(motif);
        std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif_ss);

        // 2 motifs, each with matrix in both the forward and reverse direction
        ASSERT_EQ(4, matrices.size());

        ASSERT_EQ("AP2D_f1", matrices[0].name());
        ASSERT_EQ(false, matrices[0].is_reverse_complement());
        ASSERT_EQ("AP2D_f1", matrices[1].name());
        ASSERT_EQ(true, matrices[1].is_reverse_complement());

        ASSERT_EQ("AHR_si", matrices[2].name());
        ASSERT_EQ(false, matrices[2].is_reverse_complement());
        ASSERT_EQ("AHR_si", matrices[3].name());
        ASSERT_EQ(true, matrices[3].is_reverse_complement());
    }

    {
        std::stringstream motif_ss(motif);
        std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif_ss,
                                                              ScoreMatrix::default_acgt_background,
                                                              "", // all motifs
                                                              false); // don't include reverse

        // 2 motifs, each with matrix in just the forward direction
        ASSERT_EQ(2, matrices.size());
        ASSERT_EQ("AP2D_f1", matrices[0].name());
        ASSERT_EQ(false, matrices[0].is_reverse_complement());
        ASSERT_EQ("AHR_si", matrices[1].name());
        ASSERT_EQ(false, matrices[1].is_reverse_complement());
    }

    {
        std::stringstream motif_ss(motif);
        std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif_ss,
                                                              ScoreMatrix::default_acgt_background,
                                                              "AHR_si");

        // 1 named motif
        ASSERT_EQ(2, matrices.size());
        ASSERT_EQ("AHR_si", matrices[0].name());
        ASSERT_EQ(false, matrices[0].is_reverse_complement());
        ASSERT_EQ("AHR_si", matrices[1].name());
        ASSERT_EQ(true, matrices[1].is_reverse_complement());
    }

    {
        std::stringstream motif_ss(motif);
        std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif_ss,
                                                              ScoreMatrix::default_acgt_background,
                                                              "AHR_si",
                                                              false); // don't include reverse

        // 1 named motif, only forward
        ASSERT_EQ(1, matrices.size());
        ASSERT_EQ("AHR_si", matrices[0].name());
        ASSERT_EQ(false, matrices[0].is_reverse_complement());
    }
}

std::string fimo_style_line(const ScoreMatrix& matrix, const std::string& sequence)
{
    std::stringstream ss;
    FimoStylePrinter printer(ss, false);
    matrix.score(sequence, printer);
    return ss.str();
}

TEST(ScoreMatrix, score_misc_fimo_style)
{
    const std::string motif = R"(MEME version 4
    ALPHABET= ACGT

    strands: + -

    Background letter frequencies (from file `../HOCOMOCO/bkg.txt'):
    A 0.29182 C 0.20818 G 0.20818 T 0.29182
    MOTIF AP2D_f1 AP2D

    letter-probability matrix: alength= 4 w= 14 nsites= 5 E= 0
      0.200000	  0.400000	  0.200000	  0.200000
      0.000000	  0.000000	  1.000000	  0.000000
      0.000000	  0.600000	  0.200000	  0.200000
      0.000000	  1.000000	  0.000000	  0.000000
      0.000000	  0.600000	  0.000000	  0.400000
      0.000000	  0.000000	  1.000000	  0.000000
      0.400000	  0.200000	  0.400000	  0.000000
      0.000000	  0.000000	  1.000000	  0.000000
      0.000000	  0.000000	  1.000000	  0.000000
      0.000000	  0.800000	  0.200000	  0.000000
      0.200000	  0.200000	  0.400000	  0.200000
      0.000000	  0.800000	  0.200000	  0.000000
      0.000000	  0.200000	  0.800000	  0.000000
      0.200000	  0.200000	  0.000000	  0.600000

    MOTIF AIRE_f2 AIRE

    letter-probability matrix: alength= 4 w= 18 nsites= 41 E= 0
    0.390244	  0.268293	  0.121951	  0.219512
    0.195122	  0.195122	  0.146341	  0.463415
    0.146341	  0.146341	  0.195122	  0.512195
    0.048780	  0.000000	  0.878049	  0.073171
    0.000000	  0.000000	  0.804878	  0.195122
    0.317073	  0.048780	  0.024390	  0.609756
    0.390244	  0.097561	  0.024390	  0.487805
    0.365854	  0.146341	  0.146341	  0.341463
    0.341463	  0.073171	  0.146341	  0.439024
    0.512195	  0.097561	  0.146341	  0.243902
    0.390244	  0.048780	  0.073171	  0.487805
    0.219512	  0.146341	  0.073171	  0.560976
    0.000000	  0.024390	  0.878049	  0.097561
    0.000000	  0.000000	  0.975610	  0.024390
    0.219512	  0.097561	  0.243902	  0.439024
    0.073171	  0.195122	  0.195122	  0.536585
    0.439024	  0.024390	  0.170732	  0.365854
    0.414634	  0.268293	  0.097561	  0.219512

    )";
    std::stringstream motif_ss(motif);
    std::vector<ScoreMatrix> matrices = ScoreMatrix::read(motif_ss);
    ASSERT_EQ(4, matrices.size());
    auto& AP2D_f1 = matrices[0];
    ASSERT_EQ("AP2D_f1", AP2D_f1.name());
    ASSERT_EQ(false, AP2D_f1.is_reverse_complement());
    auto& AIRE_f2 = matrices[2];
    ASSERT_EQ("AIRE_f2", AIRE_f2.name());
    ASSERT_EQ(false, AIRE_f2.is_reverse_complement());

    ASSERT_EQ("AP2D_f1\t\t1\t14\t+\t15.626\t3.66e-06\t\tGGCCTGCGGGGGGT\n",
              fimo_style_line(AP2D_f1, "GGCCTGCGGGGGGT"));
    ASSERT_EQ("AP2D_f1\t\t1\t14\t+\t15.626\t3.66e-06\t\tGGCCTGCGGGGGGT\n",
              fimo_style_line(AP2D_f1, "GGTGTAACTGGAC")); // CCCCCGCAGGCCCAAGGGAAGGCTGTN
    //const auto score = AP2D_f1.score_sequence("GGTGTAACTGGACCCCCCGCAGGCCCAAGGGAAGGCTGTN", 1, 14);
    //std::cout << score.score() << "\t" << score.pvalue() << std::endl;

    // In meme version 4.11.3, following pvalue changes from 7.75e-07 to 7.74e-07.
    // More recent meme versions behavior changes seem somewhat arbitrary,
    // so sticking with meme version 4.10.1 as expected values for now.
    ASSERT_EQ("AIRE_f2\t\t1\t18\t+\t15.0111\t7.75e-07\t\tCTTGGATTTATTGGACTA\n",
              fimo_style_line(AIRE_f2, "CTTGGATTTATTGGACTA"));
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

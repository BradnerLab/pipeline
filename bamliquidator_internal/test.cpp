#include "gtest/gtest.h"

#include "bamliquidator_util.h"

TEST(count, misc)
{
  EXPECT_EQ(0, count("", ""));
  EXPECT_EQ(0, count("", "A"));
  EXPECT_EQ(0, count("A", "B"));

  EXPECT_EQ(1, count("A", "A"));
  EXPECT_EQ(1, count("abc", "a"));
  EXPECT_EQ(1, count("abc", "bc"));

  EXPECT_EQ(1, count("bcbcbcbcdbcbc", "bcd"));
  //                     1. bcd

  EXPECT_EQ(0, count("bcbcbcbcdbcbc", "bccd"));
  EXPECT_EQ(0, count("abc", "d"));

  EXPECT_EQ(1, count("N", "N"));
  EXPECT_EQ(0, count("N", "A"));
  EXPECT_EQ(1, count("A", "N"));

  EXPECT_EQ(2, count("AB", "N"));
  EXPECT_EQ(2, count("NB", "N"));
  EXPECT_EQ(2, count("BN", "N"));
  
  EXPECT_EQ(0, count("BB", "AN"));
  EXPECT_EQ(0, count("BN", "AN"));
  EXPECT_EQ(1, count("BN", "BN"));
  EXPECT_EQ(1, count("ABB", "ANB"));

  EXPECT_EQ(1, count("GGGGGTAGAAGAGGAAGAGAGGAGGGGGGAAATCCCCTTT", "GGGAAATCCCCT"));
  //                                         1. GGGAAATCCCCT

  EXPECT_EQ(1, count("TTGGAAGGTACTCCTTTTTAGTAAGGGAAATCCCCTCTTC", "GGGAAATCCCCT"));
  //                                       1. GGGAAATCCCCT

  EXPECT_EQ(1, count("CATGGACACGGGACAGGTATTCAGCGGAAATTCCTTCCAG", "GGNNNTTCC"));
  //                                        1. GGNNNTTCC

  EXPECT_EQ(0, count("GATGGATCACAGGTCTATCACCCTATTAACCACTCACGGG", "AGGGGATTTCCC"));
  EXPECT_EQ(0, count("GAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTG", "GGNNNTTCC"));

  EXPECT_EQ(2, count("TTGGAAGATACTCCTTTTAAGAAG", "AAGA"));
  //               1.     AAGA
  //               2.                   AAGA

  EXPECT_EQ(3, count("TTGGAAGATACTCCTTTTAAGAAGAGGAAATCCCCTCTTC", "AAGA"));
  //               1.     AAGA
  //               2.                   AAGA
  //               3.                      AAGA
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

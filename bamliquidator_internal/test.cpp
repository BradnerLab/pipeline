#include "gtest/gtest.h"

#include "bamliquidator_util.h"

TEST(contains, misc)
{
  EXPECT_FALSE(contains("", ""));
  EXPECT_FALSE(contains("", "A"));
  EXPECT_FALSE(contains("A", "B"));

  EXPECT_TRUE(contains("A", "A"));
  EXPECT_TRUE(contains("abc", "a"));
  EXPECT_TRUE(contains("abc", "bc"));
  EXPECT_TRUE(contains("bcbcbcbcdbcbc", "bcd"));
  EXPECT_FALSE(contains("bcbcbcbcdbcbc", "bccd"));
  EXPECT_FALSE(contains("abc", "d"));

  EXPECT_TRUE(contains("N", "N"));
  EXPECT_FALSE(contains("N", "A"));
  EXPECT_TRUE(contains("A", "N"));

  EXPECT_TRUE(contains("AB", "N"));
  EXPECT_TRUE(contains("NB", "N"));
  EXPECT_TRUE(contains("BN", "N"));
  
  EXPECT_FALSE(contains("BB", "AN"));
  EXPECT_FALSE(contains("BN", "AN"));
  EXPECT_TRUE(contains("BN", "BN"));
  EXPECT_TRUE(contains("ABB", "ANB"));

  EXPECT_TRUE(contains("GGGGGTAGAAGAGGAAGAGAGGAGGGGGGAAATCCCCTTT", "GGGAAATCCCCT"));
  EXPECT_TRUE(contains("TTGGAAGGTACTCCTTTTTAGTAAGGGAAATCCCCTCTTC", "GGGAAATCCCCT"));
  EXPECT_TRUE(contains("CATGGACACGGGACAGGTATTCAGCGGAAATTCCTTCCAG", "GGNNNTTCC"));

  EXPECT_FALSE(contains("GATGGATCACAGGTCTATCACCCTATTAACCACTCACGGG", "AGGGGATTTCCC"));
  EXPECT_FALSE(contains("GAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTG", "GGNNNTTCC"));
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

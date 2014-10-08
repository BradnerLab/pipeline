#include "gtest/gtest.h"

TEST(initial, passing)
{
  ASSERT_TRUE(true);
}

TEST(initial, failing)
{
  ASSERT_TRUE(false);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

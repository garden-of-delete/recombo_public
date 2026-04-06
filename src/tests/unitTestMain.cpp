#include <gtest/gtest.h>

#include "argprocessorTests.h"
#include "autocorrTests.h"
#include "clkTests.h"
#include "randomTests.h"
#include "recomboCriteriaTests.h"

TEST(AutocorrTest, Autocorr) { EXPECT_TRUE(testAutocorr()); }

TEST(ClkTest, Data) { EXPECT_TRUE(testData()); }
TEST(ClkTest, CopyScalarToVector) { EXPECT_TRUE(testCopyScalarToVector()); }
TEST(ClkTest, GenericOutput) { EXPECT_TRUE(testGenericOutput()); }
TEST(ClkTest, Rog) { EXPECT_TRUE(testRog()); }
TEST(ClkTest, WrAcn) { EXPECT_TRUE(testWrAcn()); }
TEST(ClkTest, ReadFromText) { EXPECT_TRUE(testReadFromText()); }
TEST(ClkTest, ReadFromCoords) { EXPECT_TRUE(testReadFromCoords()); }
TEST(ClkTest, RandomReset) { EXPECT_TRUE(testRandomReset()); }
TEST(ClkTest, Bfacf3) { EXPECT_TRUE(testBfacf3()); }
TEST(ClkTest, Bfacf3SetZ) { EXPECT_TRUE(testBfacf3SetZ()); }
TEST(ClkTest, Bfacf3Run) { EXPECT_TRUE(testBfacf3Run()); }

TEST(RandomTest, RawSequence) { EXPECT_TRUE(testRandom()); }
TEST(RandomTest, Integers) { EXPECT_TRUE(testRandomInteger()); }

TEST(BfacfProbsTest, Precomputed) { EXPECT_TRUE(testPrecomputedBfacf3Probs()); }
TEST(BfacfProbsTest, HandComputed) { EXPECT_TRUE(testBfacf3ProbsHandComputed()); }

TEST(RecomboTest, ParallelRecombination) { EXPECT_TRUE(testParallelRecombination()); }

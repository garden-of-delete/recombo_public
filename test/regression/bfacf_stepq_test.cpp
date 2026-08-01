// Data-driven regression tests for stepQ()/init_Q() -- a structurally
// different code path from step() (used by bfacf_critical_z_test.cpp and
// bfacf_topology_invariance_test.cpp), which has its own q-dependent
// precomputed probability table (probMap) rather than step()'s direct
// perform_move(). None of the other regression suites in this project
// exercise it. These runs are deliberately shorter than the critical-z/
// topology suites -- the goal here is code-path coverage across several q
// values and knot/link types, not equilibrium or critical behavior.
//
// Cases and expected_newsud values are generated, not hand-typed -- see
// test/data/bfacf_stepq_cases.json. z-values reused from
// bfacf_critical_z_test.cpp (same provenance: recovered historical
// z-value table for 3_1/4_1, CRITICAL_Z for 0_1, DEFAULT_Z for 2_2_1).
//
// Observation: figure8_q1_seed1 and figure8_q2_seed1 produced IDENTICAL
// output despite q=1 and q=2 using different probabilities (verified
// directly via probMap: e.g. at length 34, p_plus2 is 0.0318 at q=1 vs
// 0.0335 at q=2 -- genuinely different, not a bug). The two probability
// sets are close enough that this particular 300-step, seed-1 sequence
// happened not to cross a differing accept/reject decision. Left in as a
// real, verified regression case rather than re-seeded to hide it -- it's
// a noteworthy result in its own right (q can fail to visibly matter for a
// given short trajectory even though the underlying probabilities differ).
//
// To confirm this is a coincidence specific to that seed rather than a
// sign q is silently ignored, the same comparison was tried for seeds
// 1..100: seed=1 matched (as above); seed=2 diverged immediately (q=1 ->
// length 38, q=2 -> length 30, no relation between the two strings at
// all). Both seed=1 and seed=2 are kept as separate cases below
// (figure8_q1/q2_seed1, figure8_q1/q2_seed2) so the suite documents both
// outcomes rather than just the more "expected" one.

#include <gtest/gtest.h>
#include "test_data_util.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"

using namespace std;

struct BfacfStepQCase {
    string name;
    string knot_type;
    string source_file;
    int component_count;
    int seed;
    double z;
    int q;
    int steps;
    vector<string> expected_newsud;
};

static vector<BfacfStepQCase> loadBfacfStepQCases()
{
    json11::Json data = load_test_data("bfacf_stepq_cases.json");
    vector<BfacfStepQCase> cases;

    for (const json11::Json& item : data["cases"].array_items())
    {
        BfacfStepQCase c;
        c.name = item["name"].string_value();
        c.knot_type = item["knot_type"].string_value();
        c.source_file = item["source_file"].string_value();
        c.component_count = item["component_count"].int_value();
        c.seed = item["seed"].int_value();
        c.z = item["z"].number_value();
        c.q = item["q"].int_value();
        c.steps = item["steps"].int_value();
        for (const json11::Json& s : item["expected_newsud"].array_items())
            c.expected_newsud.push_back(s.string_value());
        cases.push_back(c);
    }
    return cases;
}

class BfacfStepQTest : public ::testing::TestWithParam<BfacfStepQCase>
{
};

TEST_P(BfacfStepQTest, MatchesRecordedResult)
{
    const BfacfStepQCase& c = GetParam();

    string fullPath = string(REPO_ROOT_DIR) + "/" + c.source_file;
    ifstream file(fullPath.c_str());
    ASSERT_TRUE(file.is_open()) << "Could not open " << fullPath;

    clkConformationAsList comp0, comp1;
    ASSERT_TRUE(comp0.readFromCoords(file)) << "Could not read component 0 from " << c.source_file;

    unique_ptr<clkConformationBfacf3> knot;
    if (c.component_count == 2)
    {
        ASSERT_TRUE(comp1.readFromCoords(file)) << "Could not read component 1 from " << c.source_file;
        knot.reset(new clkConformationBfacf3(comp0, comp1));
    }
    else
    {
        knot.reset(new clkConformationBfacf3(comp0));
    }

    knot->setSeed(c.seed);
    knot->init_Q(c.z, c.q);
    for (int i = 0; i < c.steps; i++)
        knot->stepQ(c.q, c.z);

    ASSERT_EQ((int)c.expected_newsud.size(), c.component_count) << c.name;
    for (int i = 0; i < c.component_count; i++)
    {
        clkConformationAsList result(knot->getComponent(i));
        EXPECT_EQ(result.writeAsNewsud(), c.expected_newsud[i])
            << c.name << ": component " << i << " mismatch";
    }
}

INSTANTIATE_TEST_SUITE_P(
    RecordedCases,
    BfacfStepQTest,
    ::testing::ValuesIn(loadBfacfStepQCases()),
    [](const ::testing::TestParamInfo<BfacfStepQCase>& info) {
        return info.param.name;
    });

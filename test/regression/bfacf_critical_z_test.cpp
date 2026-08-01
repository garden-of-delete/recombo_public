// Data-driven BFACF regression tests: run BFACF at a specific z-value for a
// known number of steps, from a real initial/* starting conformation, and
// check the resulting NEWSUD encoding against a previously-recorded value.
//
// Cases and expected_newsud values live in test/data/bfacf_critical_z_cases.json
// and are NOT hand-typed -- they were produced by actually running this same
// simulation once and recording its output. Adding a new case is a data
// change, not a code change: the TEST_P body below handles any number of
// cases from the JSON file.
//
// z-value provenance per case, briefly:
//  - 3_1, 4_1: a per-knot-type, per-target-length z-value table used to
//    exist in src/legacyBfacf.cpp (added 2014, removed in "Remove legacy
//    z-value tables", commit a54c754). These values are recovered from
//    that table (`git show a54c754^:src/legacyBfacf.cpp`) -- the target
//    length used for each case is in its name, e.g. trefoil_z60_seed1 used
//    the entry calibrated for average length 60.
//  - 0_1: CRITICAL_Z = 1/4.6852 (Schmirler 2012), the same constant already
//    used in bfacf_regression_test.cpp. No entry exists for 0_1 in the
//    recovered table -- it was UNKNOWN_Z at every target length there, and
//    it also resisted a live zAnalyzer run (persistent non-convergence, see
//    the commit history/PR discussion around this test for details).
//  - 2_2_1: no calibrated value exists anywhere for this link; DEFAULT_Z
//    (0.20815) is used as the codebase's general-purpose operating default,
//    not a value verified for this specific topology.

#include <gtest/gtest.h>
#include "test_data_util.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"

using namespace std;

struct BfacfCriticalZCase {
    string name;
    string knot_type;
    string source_file;
    int component_count;
    int seed;
    double z;
    int steps;
    vector<string> expected_newsud;
};

static vector<BfacfCriticalZCase> loadBfacfCriticalZCases()
{
    json11::Json data = load_test_data("bfacf_critical_z_cases.json");
    vector<BfacfCriticalZCase> cases;

    for (const json11::Json& item : data["cases"].array_items())
    {
        BfacfCriticalZCase c;
        c.name = item["name"].string_value();
        c.knot_type = item["knot_type"].string_value();
        c.source_file = item["source_file"].string_value();
        c.component_count = item["component_count"].int_value();
        c.seed = item["seed"].int_value();
        c.z = item["z"].number_value();
        c.steps = item["steps"].int_value();
        for (const json11::Json& s : item["expected_newsud"].array_items())
            c.expected_newsud.push_back(s.string_value());
        cases.push_back(c);
    }
    return cases;
}

class BfacfCriticalZTest : public ::testing::TestWithParam<BfacfCriticalZCase>
{
};

TEST_P(BfacfCriticalZTest, MatchesRecordedResult)
{
    const BfacfCriticalZCase& c = GetParam();

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
    knot->setZ(c.z);
    for (int i = 0; i < c.steps; i++)
        knot->step();

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
    BfacfCriticalZTest,
    ::testing::ValuesIn(loadBfacfCriticalZCases()),
    [](const ::testing::TestParamInfo<BfacfCriticalZCase>& info) {
        return info.param.name;
    });

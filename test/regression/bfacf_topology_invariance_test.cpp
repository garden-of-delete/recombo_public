// Topology-invariance regression tests: grow a conformation at a relatively
// high z, then shrink it at a very low z, and check the resulting length.
//
// The premise: a 4-edge closed self-avoiding polygon on the cubic lattice
// can only ever be the unknot. BFACF's elementary moves are
// topology-preserving by construction, so shrinking at low z can only reach
// that trivial 4-edge square if the conformation genuinely is the unknot.
// A conformation that started (or remained, after growth) genuinely knotted
// can never reach 4 edges -- it will stabilize at some higher length
// specific to its own knot type. Reaching exactly 4 is therefore only
// expected for 0_1; any other case reaching 4 would indicate topology was
// NOT preserved during the run.
//
// Caveat: this specific "did length reach 4" check is only valid evidence
// for the topology of an individual, single-component conformation. It is
// NOT valid evidence that a link's components remain linked to each other --
// many nontrivial links (e.g. the Hopf link, which 2_2_1 in this project's
// naming likely is) have components that are each individually an unknot;
// only the pair is topologically nontrivial. So the 2_2_1 case here only
// pins down the observed regression value -- it does not assert length > 4,
// since a component reaching a short length near 4 would not by itself be
// evidence of a bug.
//
// Cases and expected_final_newsud values are generated, not hand-typed --
// see test/data/bfacf_topology_invariance_cases.json.

#include <gtest/gtest.h>
#include "test_data_util.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"

using namespace std;

struct BfacfTopologyCase {
    string name;
    string knot_type;
    string source_file;
    int component_count;
    int seed;
    double high_z;
    int high_steps;
    double low_z;
    int low_steps;
    vector<string> expected_final_newsud;
    bool expect_trivial_square;
};

static vector<BfacfTopologyCase> loadBfacfTopologyCases()
{
    json11::Json data = load_test_data("bfacf_topology_invariance_cases.json");
    vector<BfacfTopologyCase> cases;

    for (const json11::Json& item : data["cases"].array_items())
    {
        BfacfTopologyCase c;
        c.name = item["name"].string_value();
        c.knot_type = item["knot_type"].string_value();
        c.source_file = item["source_file"].string_value();
        c.component_count = item["component_count"].int_value();
        c.seed = item["seed"].int_value();
        c.high_z = item["high_z"].number_value();
        c.high_steps = item["high_steps"].int_value();
        c.low_z = item["low_z"].number_value();
        c.low_steps = item["low_steps"].int_value();
        for (const json11::Json& s : item["expected_final_newsud"].array_items())
            c.expected_final_newsud.push_back(s.string_value());
        c.expect_trivial_square = item["expect_trivial_square"].bool_value();
        cases.push_back(c);
    }
    return cases;
}

class BfacfTopologyInvarianceTest : public ::testing::TestWithParam<BfacfTopologyCase>
{
};

TEST_P(BfacfTopologyInvarianceTest, ShrinkAtLowZAfterGrowthAtHighZ)
{
    const BfacfTopologyCase& c = GetParam();

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

    // Grow at relatively high z, to move well away from the starting
    // conformation before testing whether shrinking can undo it.
    knot->setZ(c.high_z);
    for (int i = 0; i < c.high_steps; i++)
        knot->step();

    // Shrink at very low z.
    knot->setZ(c.low_z);
    for (int i = 0; i < c.low_steps; i++)
        knot->step();

    ASSERT_EQ((int)c.expected_final_newsud.size(), c.component_count) << c.name;

    vector<string> resultNewsud;
    for (int i = 0; i < c.component_count; i++)
    {
        clkConformationAsList result(knot->getComponent(i));
        resultNewsud.push_back(result.writeAsNewsud());
        EXPECT_EQ(resultNewsud[i], c.expected_final_newsud[i])
            << c.name << ": component " << i << " mismatch";
    }

    // The topology-invariance signal itself: only a genuine unknot should
    // ever be able to reach the trivial 4-edge square. See the file header
    // for why this check is skipped for multi-component (link) cases.
    if (c.component_count == 1)
    {
        if (c.expect_trivial_square)
            EXPECT_EQ(resultNewsud[0].size(), 4u)
                << c.name << ": expected the unknot to reach the trivial 4-edge square";
        else
            EXPECT_GT(resultNewsud[0].size(), 4u)
                << c.name << ": a genuinely knotted conformation reached the trivial "
                << "4-edge square, which is only possible for the unknot -- "
                << "this indicates topology was NOT preserved";
    }
}

INSTANTIATE_TEST_SUITE_P(
    RecordedCases,
    BfacfTopologyInvarianceTest,
    ::testing::ValuesIn(loadBfacfTopologyCases()),
    [](const ::testing::TestParamInfo<BfacfTopologyCase>& info) {
        return info.param.name;
    });

#include <gtest/gtest.h>
#include "test_data_util.h"

#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include <algorithm>
#include <deque>
#include <list>
#include <cmath>

#include <threevector.h>
#include <genericConformation.h>
#include <clkCigar.h>
#include <clkConformationAsList.h>
#include <clkConformationBfacf3.h>
#include <bfacfProbabilities.h>
#include <bfacfProbabilitiesFromZ.h>
#include <bfacfProbabilitiesFromZFixed.h>
#include <pseudorandom.h>

using namespace std;

class RecomboCriteriaTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        json11::Json data = load_test_data("recombo_conformations.json");
        int conflen = data["conformation_length"].int_value();
        const json11::Json::array& preArr = data["pre_recombo_unknot"].array_items();
        const json11::Json::array& postArr = data["post_recombo_unknot"].array_items();

        // Each vertex is three ints, and clkConformationAsList reads exactly
        // 3 * conflen of them. A shorter array would be read past the end.
        ASSERT_EQ(3 * conflen, (int)preArr.size()) << "pre_recombo_unknot length";
        ASSERT_EQ(3 * conflen, (int)postArr.size()) << "post_recombo_unknot length";

        vector<int> preCoords(preArr.size());
        vector<int> postCoords(postArr.size());
        for (size_t i = 0; i < preArr.size(); i++)
            preCoords[i] = preArr[i].int_value();
        for (size_t i = 0; i < postArr.size(); i++)
            postCoords[i] = postArr[i].int_value();

        clkConformationAsList preList(preCoords.data(), conflen);
        clkConformationAsList postList(postCoords.data(), conflen);
        preConformation.reset(new clkConformationBfacf3(preList));
        postConformation.reset(new clkConformationBfacf3(postList));
    }

    unique_ptr<clkConformationBfacf3> preConformation;
    unique_ptr<clkConformationBfacf3> postConformation;
};

TEST_F(RecomboCriteriaTest, ParallelRecombination)
{
    int sites = 0;
    int sitechoice = 0;
    int min_arc = 18, max_arc = 22;
    int sequence_type = 1;
    int recombo_type = 1;
    int n_components = 1;
    int total_para_site, total_anti_site, Para_site, Anti_site;
    sites=preConformation->countRecomboSites(min_arc, max_arc, sequence_type, recombo_type, total_para_site, total_anti_site, Para_site, Anti_site);
    preConformation->performRecombination(sitechoice);
    ASSERT_EQ(preConformation->size(), 1);
    ASSERT_TRUE(preConformation->getComponent(0)==postConformation->getComponent(0));
}

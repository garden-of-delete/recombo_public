#include <gtest/gtest.h>
#include "test_data_util.h"

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <deque>
#include <list>
#include <cmath>

#include "mmchain.h"
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

class recomboCriteriaTestClass
{
    public:
    clkConformationBfacf3* preConformation;
    clkConformationBfacf3* postConformation;
    mmchain mmctest;
    recomboCriteriaTestClass();
    ~recomboCriteriaTestClass();
};

recomboCriteriaTestClass::recomboCriteriaTestClass()
{
    json11::Json data = load_test_data("recombo_conformations.json");
    int conflen = data["conformation_length"].int_value();
    const json11::Json::array& preArr = data["pre_recombo_unknot"].array_items();
    const json11::Json::array& postArr = data["post_recombo_unknot"].array_items();

    vector<int> preCoords(preArr.size());
    vector<int> postCoords(postArr.size());
    for (size_t i = 0; i < preArr.size(); i++)
        preCoords[i] = preArr[i].int_value();
    for (size_t i = 0; i < postArr.size(); i++)
        postCoords[i] = postArr[i].int_value();

    clkConformationAsList preList(preCoords.data(), conflen);
    clkConformationAsList postList(postCoords.data(), conflen);
    preConformation = new clkConformationBfacf3(preList);
    postConformation = new clkConformationBfacf3(postList);
}

recomboCriteriaTestClass::~recomboCriteriaTestClass()
{
    delete preConformation;
    delete postConformation;
}

TEST(RecomboTest, ParallelRecombination)
{
    recomboCriteriaTestClass rctc;

    int sites = 0;
    int sitechoice = 0;
    int min_arc = 18, max_arc = 22;
    int sequence_type = 1;
    int recombo_type = 1;
    int n_components = 1;
    int total_para_site, total_anti_site, Para_site, Anti_site;
    sites=rctc.preConformation->countRecomboSites(min_arc, max_arc, sequence_type, recombo_type, total_para_site, total_anti_site, Para_site, Anti_site);
    rctc.preConformation->performRecombination(sitechoice);
    ASSERT_EQ(rctc.preConformation->size(), 1);
    ASSERT_TRUE(rctc.preConformation->getComponent(0)==rctc.postConformation->getComponent(0));
}

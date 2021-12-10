#include "testFramework.h"

#include <vector>
#include <string>

#include "recomboCriteriaTests.h"
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

#include <sstream>
#include <iostream>

#include <algorithm>
#include <deque>
#include <list>
#include <cmath>

using namespace std;

class recomboCriteriaTestClass 
{
    private:
    static const int conflen= 120;
    static const int preRecomboUnknot[conflen];
    static const int postRecomboUnknot[conflen];
    static const clkConformationAsList preConformationList;
    static const clkConformationAsList postConformationList;
    
    public:
    clkConformationBfacf3 preConformation;
    clkConformationBfacf3 postConformation;
    mmchain mmctest;
    recomboCriteriaTestClass(); 
    bool testParallelRecombination();
};

const int recomboCriteriaTestClass::preRecomboUnknot[120]={0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0, 5, 0, 0, 6, 0, 1, 6, 0, 2, 6, 0, 3, 6, 0, 4, 6, 0, 4, 5, 0, 5, 5, 0, 6, 5, 0, 7, 5, 0, 7, 4, 0, 7, 3, 0, 7, 2, 0, 7, 1, 0, 7, 0, 0, 7, -1, 0, 6, -1, 0, 5, -1, 0, 4, -1, 0, 3, -1, 0, 3, 0, 0, 3, 0, 1, 3, 1, 1, 3, 2, 1, 3, 2, 0, 3, 3, 0, 3, 4, 0, 4, 4, 0, 5, 4, 0, 5, 3, 0, 5, 2, 0, 5, 1, 0, 4, 1, 0, 3, 1, 0, 2, 1, 0, 1, 1, 0};
const int recomboCriteriaTestClass::postRecomboUnknot[120]={0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0, 5, 0, 0, 6, 0, 1, 6, 0, 2, 6, 0, 3, 6, 0, 4, 6, 0, 4, 5, 0, 4, 4, 0, 3, 4, 0, 3, 3, 0, 3, 2, 0, 3, 2, 1, 3, 1, 1, 3, 0, 1, 3, 0, 0, 3, -1, 0, 4, -1, 0, 5, -1, 0, 6, -1, 0, 7, -1, 0, 7, 0, 0, 7, 1, 0, 7, 2, 0, 7, 3, 0, 7, 4, 0, 7, 5, 0, 6, 5, 0, 5, 5, 0, 5, 4, 0, 5, 3, 0, 5, 2, 0, 5, 1, 0, 4, 1, 0, 3, 1, 0, 2, 1, 0, 1, 1};
const clkConformationAsList recomboCriteriaTestClass::preConformationList(preRecomboUnknot,40);
const clkConformationAsList recomboCriteriaTestClass::postConformationList(postRecomboUnknot,40);

recomboCriteriaTestClass::recomboCriteriaTestClass() : preConformation(preConformationList), postConformation(postConformationList)
{}

bool testParallelRecombination()
{
    recomboCriteriaTestClass rctc;

    int sites = 0;
    int sitechoice = 0;
    int min_arc = 18, max_arc = 22;
    int sequence_type = 1;
    int recombo_type = 1;
    int n_components = 1;
    int total_para_site, total_anti_site, Para_site, Anti_site;
    sites=rctc.preConformation.countRecomboSites(min_arc, max_arc, sequence_type, recombo_type, total_para_site, total_anti_site, Para_site, Anti_site);
    rctc.preConformation.performRecombination(cout, sequence_type, recombo_type, n_components, sitechoice);
    ASSERT_EQUAL(rctc.preConformation.size(),1);
    ASSERT(rctc.preConformation.getComponent(0)==rctc.postConformation.getComponent(0));
    
    return true;
}
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
{



}

bool testParallelRecombination()
{
    recomboCriteriaTestClass rctc;
    //The choice between either as the before recombination is arbitrary; 
    //they recombine into eachother without outside influence

    int sites = 0;
    int sitechoice = 0;
    //result.getVertex(0,v);
    //v*=-1;
    //result.translate(v);
    //clkConformationBfacf3 postConformation(reclass.postConformationList);
    char orientation = 'p';
    sites=rctc.preConformation.countRecomboSites(18,22,orientation);
    rctc.preConformation.performRecombination(sitechoice);
    ASSERT_EQUAL(rctc.preConformation.size(),1);
    ASSERT(rctc.preConformation.getComponent(0)==rctc.postConformation.getComponent(0));
    //clkConformationAsList result(rctc.preConformation.getComponent(0));
    
    return true;
}
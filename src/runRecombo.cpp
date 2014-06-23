/* 
 * File:   runRecombo.cpp
 * Author: kmo
 * 
 * Created on February 28, 2013, 9:11 AM
 */

#include "runRecombo.h"

using namespace std;

runRecombo::runRecombo(const runRecombo& orig) : args(orig.args) { }

runRecombo::runRecombo(recomboArgprocessor& args) : args(args),
conformation(0), siteSelector(args.getSeed()*17 + 53),
savedComponent0(0), savedComponent1(0),
numberOfComponents(0)
{
   // record the time
   time(&starttime);

}

runRecombo::~runRecombo()
{
   //   cout << "~runRecombo" << endl;
   if (conformation) delete conformation;
   if (savedComponent0) delete savedComponent0;
   if (savedComponent1) delete savedComponent1;
}

bool runRecombo::finished(long conformationsRecorded, long iterationsDone) const
{
   if (args.useNumberOfIterations())
      return conformationsRecorded >= args.getNumberOfConformations();
   return iterationsDone >= args.getNumberOfIterations();
}

void runRecombo::warmup()
{
   conformation->step(args.getWarmup());
}

void runRecombo::timestamp(ostream& os, time_t time) const
{
   char buffer[256];
   strftime(buffer, 256, "%c", localtime(&time));
   os << buffer;
}

void runRecombo::recordConformation(ostream& os, clkConformationBfacf3 *conformation)
{
   // TODO: check for exceptions
   list<clkConformationAsList> components;
   conformation->getComponents(components);
   list<clkConformationAsList>::const_iterator i;
   for (i = components.begin(); i != components.end(); i++)
   {
      i->writeAsCube(os);
      statsfile << i->size() << " ";
   }
}

void runRecombo::writeLog(ostream& os) const
{
   os << "Initial coords taken from file " << args.getInitialConformationFileName();
   if (getNumberOfComponents() == 1)
      os << ", which contains a knot." << endl;
   else if (getNumberOfComponents() == 2)
      os << ", which contains a two component link." << endl;
   else
      os << ", which contains a " << getNumberOfComponents()
      << " component link." << endl;
   //   const string& outfilename = args.outputFileName();
   os << "Output data to files " << args.getBeforeFileName() << ", "
           << args.getAfterFileName() << ", " << args.getStatsFileName()
           << " and " << args.getLogFileName()
           << "." << endl;
   os << "Set z to " << args.getZ() << "." << endl;
   os << "Random number seed set to " << args.getSeed();

   if (args.seedWasSetUsingSystemTimer())
      os << " from the system timer." << endl;
   else
      os << "." << endl;

   if (args.getWarmup() > 0)
      os << "Begin with a warmup of " << args.getWarmup() << " iterations." << endl;
   else
      os << "No warmup." << endl;

   if (args.useNumberOfIterations())
      os << "Run the simulation until " << args.getNumberOfConformations() << " conformations output." << endl;
   else
      os << "Run the simulation for " << args.getNumberOfIterations() << " BFACF iterations performed." << endl;

   if (args.getMinarc() == args.getMaxarc())
      os << "Recombo criterion is arclength exactly " << args.getMinarc() << "." << endl;
   else
      os << "Recombo criterion is arclength from " << args.getMinarc() << " to " << args.getMaxarc() << "." << endl;

   os << "Sample at most every " << args.getSampleRate() << " iterations using protocol ";
   if (args.aProtocol()) os << "A." << endl;
   else if (args.hProtocol()) os << "H." << endl;
   else os << "B." << endl;

   os << "Experiment begun ";
   timestamp(os, starttime);
   os << "." << endl;
}

bool runRecombo::openfiles()
{
   logfile.open(args.getLogFileName().c_str(), ios::out);
   if (!logfile)
      return false;

   statsfile.open(args.getStatsFileName().c_str(), ios::out | ios::binary);
   if (!statsfile)
   {
      logfile.close();
      return false;
   }

   before.open(args.getBeforeFileName().c_str(), ios::out | ios::binary);
   if (!before)
   {
      logfile.close();
      statsfile.close();
      return false;
   }

   after.open(args.getAfterFileName().c_str(), ios::out | ios::binary);
   if (!after)
   {
      logfile.close();
      statsfile.close();
      before.close();
      return false;
   }

   // record log
   writeLog(logfile);

   return true;
}

void runRecombo::closefiles()
{
   // log time
   time_t endtime;
   time(&endtime);
   logfile << "Experiment ended ";
   timestamp(logfile, endtime);
   logfile << "." << endl << "Experiment duration was "
           << difftime(endtime, starttime) << " seconds." << endl;

   logfile.close();
   statsfile.close();
   before.close();
   after.close();
}

void runRecombo::markConformation()
{
   if (savedComponent0) delete savedComponent0;
   if (savedComponent1) delete savedComponent1;

   if (conformation->size() == 1)
   {
      savedComponent0 = new clkConformationAsList(conformation->getComponent(0));
      delete conformation;
      conformation = new clkConformationBfacf3(*savedComponent0);

      //      cout << "One component marked." << endl;
      //      cout << *savedComponent0 << endl; 
   }
   else
   {
      savedComponent0 = new clkConformationAsList(conformation->getComponent(0));
      savedComponent1 = new clkConformationAsList(conformation->getComponent(1));
      delete conformation;
      conformation = new clkConformationBfacf3(*savedComponent0, *savedComponent1);

      //      cout << "Two components marked." << endl;
      //      cout << *savedComponent0 << endl;
      //      cout << *savedComponent1 << endl;
   }
   conformation->setZ(args.getZ());
}

void runRecombo::restoreConformation()
{
   delete conformation;
   if (savedComponent1)
      conformation = new clkConformationBfacf3(*savedComponent0, *savedComponent1);
   else
      conformation = new clkConformationBfacf3(*savedComponent0);
   conformation->setZ(args.getZ());
}

runRecomboA::runRecomboA(recomboArgprocessor& args) : runRecombo(args) { }

runRecomboA::~runRecomboA() { }

void runRecomboA::dorecombo()
{
   cout << "Protocol A." << endl;
   cout << getInitialComp0() << endl;
   if (getNumberOfComponents() > 1)
      cout << getInitialComp1() << endl;

   long iteration, recorded;
   iteration = recorded = 0;
   int rejected = 0;

   warmup();

   while (!finished(recorded, iteration))
   {
      markConformation();
      //      pseudorandom s = saveRandomState();
      conformation->step(args.getSampleRate());
      //      cout << "?" << conformation->getComponent(0) << endl;
      cout << conformation->getComponent(0).size() << " " << conformation->getZ() << endl;
      iteration += args.getSampleRate();
      int numberOfSites = conformation->countRecomboSites(args.getMinarc(), args.getMaxarc());
      if (numberOfSites > 0)
      {
         int site = siteSelector.rand_integer(0, numberOfSites);
         statsfile << iteration << " " << rejected << " " << numberOfSites
                 << " " << site << " ";
         recordConformation(before, conformation);
         conformation->performRecombination(site);
         recordConformation(after, conformation);
         conformation->undoRecombination();
         statsfile << endl;
         recorded++;
         rejected = 0;

         cout << "Conformation " << recorded << " recorded." << endl;
      }
      else
      {
         restoreConformation();
         rejected++;

         cout << "Conformation rejected." << endl;
      }
   }
}

runRecomboH::runRecomboH(recomboArgprocessor& args) : runRecombo(args) { }

runRecomboH::~runRecomboH() {
   //   cout << "~runRecomboH" << endl;
}

void runRecomboH::dorecombo()
{
   cout << "Protocol H." << endl;
   cout << getInitialComp0() << endl;
   if (getNumberOfComponents() > 1)
      cout << getInitialComp1() << endl;

   long iteration, recorded;
   iteration = recorded = 0;
   int rejected = 0;

   warmup();

   while (!finished(recorded, iteration))
   {
      conformation->step(args.getSampleRate());
      iteration += args.getSampleRate();
      int numberOfSites = conformation->countRecomboSites(args.getMinarc(), args.getMaxarc());
      if (numberOfSites > 0)
      {
         int site = siteSelector.rand_integer(0, numberOfSites);
         statsfile << iteration << " " << rejected << " " << numberOfSites
                 << " " << site << " ";
         recordConformation(before, conformation);
         conformation->performRecombination(site);
         recordConformation(after, conformation);
         conformation->undoRecombination();
         statsfile << endl;
         recorded++;
         rejected = 0;

         cout << "Conformation " << recorded << " recorded." << endl;
      }
      else
      {
         rejected++;

         cout << "Conformation rejected." << endl;
      }
   }
}

runRecomboB::runRecomboB(recomboArgprocessor& args) : runRecombo(args) { }

runRecomboB::~runRecomboB() { }

void runRecomboB::dorecombo()
{
   // record log
   //   cout << "Protocol B." << endl;
   //   cout << getInitialComp0() << endl;
   //   if (getNumberOfComponents() > 1)
   //      cout << getInitialComp1() << endl;

   long iteration, recorded;
   iteration = recorded = 0;
   int rejected = 0;
   int conformationsFound;
   int stepsize = 100;
   bool found;
   int numberOfSites;

   warmup();

   pseudorandom r, s;

   while (!finished(recorded, iteration))
   {
      iteration += 2 * args.getSampleRate();

      conformation->step(args.getSampleRate());
      markConformation();
      //      cout << conformation->getComponent(0) << endl;
      s = saveRandomState();
      conformationsFound = 0;
      for (int i = 0; i < args.getSampleRate(); i += stepsize)
      {
         conformation->step(stepsize);
         if (conformation->countRecomboSites(args.getMinarc(), args.getMaxarc()) > 0)
            conformationsFound++;
      }

      if (conformationsFound > 0)
      {
         // restore the marked conformation and the state of the random number generator
         restoreConformation();
         //         cout << conformation->getComponent(0) << endl;
         copyRandomState(s);
         r = saveRandomState();

         found = false;
         int chosenConformation = siteSelector.rand_integer(1, conformationsFound);
         for (int i = 0; i < args.getSampleRate(); i += stepsize)
         {
            conformation->step(stepsize);
            numberOfSites = conformation->countRecomboSites(args.getMinarc(), args.getMaxarc());
            if (numberOfSites > 0)
            {
               chosenConformation--;
               if (chosenConformation == 0)
               {
                  found = true;
                  int site = siteSelector.rand_integer(0, numberOfSites);
                  statsfile << iteration << " " << rejected << " " << numberOfSites
                          << " " << site << " ";
                  recordConformation(before, conformation);
                  conformation->performRecombination(site);
                  recordConformation(after, conformation);
                  conformation->undoRecombination();
                  statsfile << endl;
                  recorded++;
                  rejected = 0;

                  cout << "Conformation " << recorded << " recorded." << endl;
               }
            }
         }

         if (!found)
         {
            cerr << "Problem with protocol B." << endl;
            return;
         }

      }
      else
      {
         //         restoreConformation();
         rejected++;

         cout << "All conformations rejected." << endl;
      }
   }
}

bool runRecombo::readInitialCoords(istream& is) //IMPORTANT!!!
{
   if (!initialComp0.readFromCoords(is))
      return false;
   //   cout << initialComp0 << endl;
   numberOfComponents = 1;
   if (initialComp1.readFromCoords(is))
      numberOfComponents = 2;

   if (getNumberOfComponents() == 2)
      conformation = new clkConformationBfacf3(getInitialComp0(), getInitialComp1());
   else if (getNumberOfComponents() == 1)
      conformation = new clkConformationBfacf3(getInitialComp0());

   // set z and seed the random number generator
   conformation->setZ(args.getZ());
   conformation->setSeed(args.getSeed());


   return true;
}

bool runRecombo::readInitialCoords()
{
   ifstream init(args.getInitialConformationFileName().c_str());
   if (!init)
   {
      cerr << "Could not open " << args.getInitialConformationFileName() << " for input." << endl;
      return false;
   }
   if (!readInitialCoords(init))
   {
      cerr << "File " << args.getInitialConformationFileName() << " does not contain a conformation." << endl;
   }
   init.close();
   return numberOfComponents > 0;
}

int runRecombo::getNumberOfComponents() const
{
   return numberOfComponents;
}

const clk& runRecombo::getInitialComp0() const
{
   return initialComp0;
}

const clk& runRecombo::getInitialComp1() const
{
   return initialComp1;
}

void runRerecombo::timestamp(std::ostream& os, time_t time) const
{
   char buffer[256];
   strftime(buffer, 256, "%c", localtime(&time));
   os << buffer;
}

void runRerecombo::recordConformation(std::ostream& os, clkConformationBfacf3 *conformation)
{
   // TODO: check for exceptions
   list<clkConformationAsList> components;
   conformation->getComponents(components);
   list<clkConformationAsList>::const_iterator i;
   for (i = components.begin(); i != components.end(); i++)
   {
      i->writeAsCube(os);
   }
}

clkConformationBfacf3* runRerecombo::readConformation(std::istream& is)
{
   if(is)
   {
      clkConformationAsList comp0;
      if(!comp0.readFromCube(is)) return 0;
      
      if(args.twocomonentlinks())
      {
         clkConformationAsList comp1;
         if(!comp1.readFromCube(is))
         {
            cerr << "Unexpected end of file." << endl;
            return 0;
         }
         return new clkConformationBfacf3(comp0, comp1);
      }
      
      return new clkConformationBfacf3(comp0);
   }
   
   return 0;
}

runRerecombo::runRerecombo(rerecomboArgprocessor& args) : args(args),
        siteSelector(args.getSeed())
{
   // record the time
   time(&starttime);
}

runRerecombo::~runRerecombo() { }

void runRerecombo::writeLog(ostream& os) const
{
   os << "Conformation sequence loaded from file " << args.getBeforeFileName();
   if (args.twocomonentlinks())
      os << " as two-component links." << endl;
   else
      os << " as knots." << endl;
   os << "Output data to files " << args.getBeforeFileName()
           << " and " << args.getLogFileName()
           << "." << endl;
   os << "Random number seed set to " << args.getSeed();

   if (args.seedWasSetUsingSystemTimer())
      os << " from the system timer." << endl;
   else
      os << "." << endl;

   if (args.getMinarc() == args.getMaxarc())
      os << "Recombo criterion is arclength exactly " << args.getMinarc() << "." << endl;
   else
      os << "Recombo criterion is arclength from " << args.getMinarc() << " to " << args.getMaxarc() << "." << endl;

   os << "Experiment begun ";
   timestamp(os, starttime);
   os << "." << endl;
}

bool runRerecombo::openfiles() 
{ 
   before.open(args.getBeforeFileName().c_str(), ios::in | ios::binary);
   if (!before)
   {
      cerr << "Could not open " << args.getBeforeFileName() << " for input." << endl;
      return false;
   }

   logfile.open(args.getLogFileName().c_str(), ios::out);
   if (!logfile)
   {
      cerr << "Could not open " << args.getLogFileName() << " for output." << endl;
      before.close();
      return false;
   }

   statsfile.open(args.getStatsFileName().c_str(), ios::out | ios::binary);
   if (!statsfile)
   {
      cerr << "Could not open " << args.getStatsFileName() << " for output." << endl;
      logfile.close();
      before.close();
      return false;
   }

   after.open(args.getAfterFileName().c_str(), ios::out | ios::binary);
   if (!after)
   {
      cerr << "Could not open " << args.getAfterFileName() << " for output." << endl;
      logfile.close();
      statsfile.close();
      before.close();
      return false;
   }

   // record log
   writeLog(logfile);

   return true;
}

void runRerecombo::closefiles()
{
   // log time
   time_t endtime;
   time(&endtime);
   logfile << "Experiment ended ";
   timestamp(logfile, endtime);
   logfile << "." << endl << "Experiment duration was "
           << difftime(endtime, starttime) << " seconds." << endl;

   logfile.close();
   statsfile.close();
   before.close();
   after.close();
}

void runRerecombo::dorecombo()
{
   clkConformationBfacf3* conformation;
   int count = 0, sites, choice;
   while(conformation = readConformation(before))
   {
      count++;
      sites = conformation->countRecomboSites(args.getMinarc(), args.getMaxarc());
      if(sites > 0)
      {
         choice = siteSelector.rand_integer(0, sites);
         conformation->performRecombination(choice);
         recordConformation(after, conformation);
      }
      statsfile << count << " " << sites << " " << choice << endl;
      delete conformation;
   }
}


/* 
 * File:   argprocessor.h
 * Author: kmo
 *
 * Created on September 17, 2012, 1:25 PM
 */

#ifndef ARGPROCESSOR_H
#define	ARGPROCESSOR_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <sstream>

//using namespace std;

class ArgumentProcessorArgument
{
protected:
   virtual std::ostream& display(std::ostream& os) const = 0;

public:
   std::string description;
   std::string shortDescription;
   bool found;
   std::string content;

   ArgumentProcessorArgument(const std::string& desc);

   ArgumentProcessorArgument(const std::string& desc, const std::string& shortDesc);

   friend std::ostream &operator<<(std::ostream&, const ArgumentProcessorArgument&);
};

class ArgumentProcessorBaseOption : public ArgumentProcessorArgument
{
protected:

   std::ostream& display(std::ostream& os) const;

public:
   std::string preargument;

   ArgumentProcessorBaseOption();
   ArgumentProcessorBaseOption(const std::string& desc, const std::string& prearg);
   ArgumentProcessorBaseOption(const std::string& desc, const std::string& prearg, const std::string& shortDsc);

   ArgumentProcessorBaseOption& operator=(const ArgumentProcessorBaseOption &rhs);
   bool operator==(const ArgumentProcessorBaseOption &rhs) const;
   bool operator<(const ArgumentProcessorBaseOption &rhs) const;
};

class ArgumentProcessorFlag : public ArgumentProcessorBaseOption
{
public:

   ArgumentProcessorFlag();

   ArgumentProcessorFlag(const std::string& desc, const std::string& prearg);

};

class ArgumentProcessorOption : public ArgumentProcessorBaseOption
{
public:

   ArgumentProcessorOption();

   ArgumentProcessorOption(const std::string& desc, const std::string& prearg, const std::string& shortDesc);

};

class ArgumentProcessorOrderedArgument : public ArgumentProcessorArgument
{
protected:

   std::ostream& display(std::ostream& os) const;

public:
   unsigned order;

   ArgumentProcessorOrderedArgument(const std::string& desc, unsigned ord, const std::string& shortDesc);

   ArgumentProcessorOrderedArgument& operator=(const ArgumentProcessorOrderedArgument &rhs);
   bool operator==(const ArgumentProcessorOrderedArgument &rhs) const;
   bool operator<(const ArgumentProcessorOrderedArgument &rhs) const;
};

class ArgumentProcessor
{
protected:
   std::map<std::string, ArgumentProcessorFlag> flags;
   std::map<std::string, ArgumentProcessorOption> options;
   std::vector<ArgumentProcessorOrderedArgument> ordered;

   virtual void displayOptions(std::ostream& os) const;

   virtual void displayFlagPattern(std::ostream& os) const;

   virtual void displayOrdered(std::ostream& os) const;

   virtual bool processVector(std::vector<std::string> &arg);

   std::string tostr(double d) const;
   std::string tostr(int i) const;
   std::string tostr(long i) const;

public:
   const std::string commandName;

   ArgumentProcessor(const std::string& name);

   //   virtual void addFlag(const char* flag, const char* description);
   virtual void addFlag(const std::string& flag, const std::string& description);

   virtual void addOption(const std::string& option, const std::string& shortDescription, const std::string& description);

   virtual void addOrdered(const std::string& shortDesc, const std::string& description);

   virtual bool processOrdered(const std::string& arg, int n);

   virtual bool processFlag(const std::string& flag);

   virtual bool processMulticharFlag(const std::string& arg);

   virtual bool processSinglecharFlag(char f);

   virtual bool processSinglecharFlags(const std::string& arg);

   virtual bool processFlags(const std::string& arg);

   virtual bool processOption(const std::string& arg1, const std::string& arg2);

   bool process(std::vector<std::string> &arg);

   bool process(int argc, char** argv);
      
   void splitForArgv(const std::string &argString, std::vector<std::string> &argv);

   friend std::ostream &operator<<(std::ostream&, const ArgumentProcessor&);
   //   friend void argumentProcessorTest();
};


#endif	/* ARGPROCESSOR_H */


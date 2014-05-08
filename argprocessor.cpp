/* 
 * File:   argprocessor.cpp
 * Author: kmo
 * 
 * Created on September 17, 2012, 1:25 PM
 */

#include "argprocessor.h"

#define DISPLAY_GAP 16

using namespace std;

ArgumentProcessorArgument::ArgumentProcessorArgument(const string& desc) : description(desc), found(false) { }

ArgumentProcessorArgument::ArgumentProcessorArgument(const string& desc, const string& shortDesc) : description(desc), shortDescription(shortDesc), found(false) { }

ostream &operator<<(ostream& os, const ArgumentProcessorArgument& t)
{
   return t.display(os);
}

ostream& ArgumentProcessorBaseOption::display(ostream& os) const
{
   int i = 0;
   os << "-";
   if (preargument.length() != 1)
   {
      os << "-";
      i = preargument.length();
   }
   os << preargument;

   if (shortDescription.length() > 0)
   {
      os << " " << shortDescription;
      i += 1 + shortDescription.length();
   }

   while (i < DISPLAY_GAP)
   {
      os << " ";
      i++;
   }
   os << description;
   return os;
}

ArgumentProcessorBaseOption::ArgumentProcessorBaseOption() : ArgumentProcessorArgument("") { }

ArgumentProcessorBaseOption::ArgumentProcessorBaseOption(const string& desc, const string& prearg) : ArgumentProcessorArgument(desc), preargument(prearg) { }

ArgumentProcessorBaseOption::ArgumentProcessorBaseOption(const string& desc, const string& prearg, const string& shortDsc) : ArgumentProcessorArgument(desc, shortDsc), preargument(prearg) { }

ArgumentProcessorBaseOption& ArgumentProcessorBaseOption::operator=(const ArgumentProcessorBaseOption &rhs)
{
   preargument = rhs.preargument;
   description = rhs.description;
   shortDescription = rhs.shortDescription;
   return *this;
}

bool ArgumentProcessorBaseOption::operator==(const ArgumentProcessorBaseOption &rhs) const
{
   return preargument == rhs.preargument
           && description == rhs.description
           && shortDescription == rhs.shortDescription;
}

bool ArgumentProcessorBaseOption::operator<(const ArgumentProcessorBaseOption &rhs) const
{
   if (preargument == rhs.preargument)
      return description < rhs.description;
   return preargument < rhs.preargument;
}

ArgumentProcessorFlag::ArgumentProcessorFlag() : ArgumentProcessorBaseOption() { }

ArgumentProcessorFlag::ArgumentProcessorFlag(const string& desc, const string& prearg) : ArgumentProcessorBaseOption(desc, prearg) { }

ArgumentProcessorOption::ArgumentProcessorOption() : ArgumentProcessorBaseOption() { }

ArgumentProcessorOption::ArgumentProcessorOption(const string& desc, const string& prearg, const string& shortDesc) : ArgumentProcessorBaseOption(desc, prearg, shortDesc) { }

ostream& ArgumentProcessorOrderedArgument::display(ostream& os) const
{
   int i = -2;
   if (shortDescription.length() > 0)
   {
      i += shortDescription.length();
      os << shortDescription;
   }
   while (i < DISPLAY_GAP)
   {
      os << " ";
      i++;
   }
   return os << description;
}

ArgumentProcessorOrderedArgument::ArgumentProcessorOrderedArgument(const string& desc, unsigned ord, const string& shortDesc) : ArgumentProcessorArgument(desc, shortDesc), order(ord) { }

ArgumentProcessorOrderedArgument& ArgumentProcessorOrderedArgument::operator=(const ArgumentProcessorOrderedArgument &rhs)
{
   description = rhs.description;
   shortDescription = rhs.shortDescription;
   order = rhs.order;
   return *this;
}

bool ArgumentProcessorOrderedArgument::operator==(const ArgumentProcessorOrderedArgument &rhs) const
{
   return order == rhs.order && description == rhs.description;
}

bool ArgumentProcessorOrderedArgument::operator<(const ArgumentProcessorOrderedArgument &rhs) const
{
   if (order == rhs.order)
      return description < rhs.description;
   return order < rhs.order;
}

void ArgumentProcessor::displayOptions(ostream& os) const
{
   map<string, ArgumentProcessorOption>::const_iterator i;
   for (i = options.begin(); i != options.end(); i++)
   {
      ArgumentProcessorOption f = i->second;
      os << " -";
      if (f.preargument.length() != 1)
         os << '-';
      os << f.preargument;
      if (f.shortDescription.length() > 0)
         os << ' ' << f.shortDescription;
   }
}

void ArgumentProcessor::displayFlagPattern(ostream& os) const
{
   if (flags.empty()) return;

   // first the one character flags
   os << " ";
   //      os << " [";
   bool foundFlag = false;
   map<string, ArgumentProcessorFlag>::const_iterator i;
   for (i = flags.begin(); i != flags.end(); i++)
   {
      ArgumentProcessorFlag f = i->second;
      if (f.preargument.length() == 1)
      {
         if (!foundFlag)
         {
            os << "-";
            foundFlag = true;
         }
         os << f.preargument;
      }
   }

   // next display the two character flags
   for (i = flags.begin(); i != flags.end(); i++)
   {
      ArgumentProcessorFlag f = (*i).second;
      if (f.preargument.length() != 1)
      {
         if (foundFlag) os << " ";
         os << "--" << f.preargument;
      }
   }

   //      os << "]";
}

void ArgumentProcessor::displayOrdered(ostream& os) const
{
   vector<ArgumentProcessorOrderedArgument>::const_iterator i;
   for (i = ordered.begin(); i != ordered.end(); i++)
   {
      os << ' ' << i->shortDescription;
   }
}

ArgumentProcessor::ArgumentProcessor(const string& name) : commandName(name) { }

void ArgumentProcessor::addFlag(const string& flag, const string& description)
{
   flags[flag] = ArgumentProcessorFlag(description, flag);
}

void ArgumentProcessor::addOption(const string& option, const string& shortDescription, const string& description)
{
   options[option] = ArgumentProcessorOption(description, option, shortDescription);
}

void ArgumentProcessor::addOrdered(const string& shortDesc, const string& description)
{
   int n = ordered.size();
   ordered.push_back(ArgumentProcessorOrderedArgument(description, n, shortDesc));
}

bool ArgumentProcessor::processOrdered(const string& arg, int n)
{
   if (n >= ordered.size())
   {
      //      cout << n << ", no such argument." << endl;
      return false;
   }
   ordered[n].found = true;
   ordered[n].content = arg;
   return true;
}

bool ArgumentProcessor::processFlag(const string& flag)
{
   map<string, ArgumentProcessorFlag>::iterator f = flags.find(flag.c_str());
   if (f == flags.end())
   {
      // no such flag found
      return false;
   }

   f->second.found = true;
   return true;
}

bool ArgumentProcessor::processMulticharFlag(const string& arg)
{
   return processFlag(arg.substr(2, arg.length() - 2));
}

bool ArgumentProcessor::processSinglecharFlag(char f)
{
   // there must be another way to do this
   // initialize string by creating a null terminated c-string
   char tocstring[2];
   tocstring[0] = f;
   tocstring[1] = 0;
   return processFlag(tocstring);
}

bool ArgumentProcessor::processSinglecharFlags(const string& arg)
{
   // the first character is '-', which is not a flag
   for (int i = 1; i < arg.length(); i++)
   {
      if (!processSinglecharFlag(arg[i])) return false;
   }
   return true;
}

bool ArgumentProcessor::processFlags(const string& arg)
{
   // assume it is already checked that arg[0] == '-'
   // first see if this might be a multi char flag
   if (arg.length() > 1 && arg[1] == '-')
      return processMulticharFlag(arg);

   // these must be single char flags
   return processSinglecharFlags(arg);
}

bool ArgumentProcessor::processOption(const string& arg1, const string& arg2)
{
   int n = 0;
   // strip the leading dashes
   while (arg1[n] == '-')
   {
      if (n == arg1.length() - 1)
      {
         // string consists entirely of dashes
         return false;
      }
      n++;
   }
   string option = arg1.substr(n, arg1.length() - n); // strip leading dashes
   //   cout << "search for option " << option << endl;

   map<string, ArgumentProcessorOption>::iterator f = options.find(option.c_str());
   if (f == options.end())
   {
      // no such option found
      return false;
   }

   f->second.found = true;
   f->second.content = arg2;
   return true;
}

bool ArgumentProcessor::processVector(vector<string> &arg)
{
   int n = 0; // to keep track of which ordered argument we are on
   // discard the first argument, which contains the command as called

   for (int i = 0; i < arg.size(); i++)
   {
      // first check if arg[i] could be a flag or option
      if (arg[i].length() > 0 && arg[i][0] == '-')
      {
         if (!processFlags(arg[i]))
         {
            // these are not flags, is this an option?
            // for this to be possible, arg[i+1] must exist
            if (i >= arg.size() - 1 || !processOption(arg[i], arg[i + 1]))
               return false;

            // successfully processed two arguments as an option
            i++;
         }
      }
      else
      {
         if (!processOrdered(arg[i], n))
            return false;
         n++;
      }
   }

   return true;
}

bool ArgumentProcessor::process(vector<string> &arg)
{
   processVector(arg);

   return true; //(added for windows compatibility)
}

bool ArgumentProcessor::process(int argc, char** argv)
{
   vector<string> arg;
   // want to skip the first string in argv, which contains the command name
   for (int i = 1; i < argc; i++)
      arg.push_back(argv[i]);
   return process(arg);
}

ostream &operator<<(ostream& os, const ArgumentProcessor& t)
{
   os << t.commandName;
   if (t.flags.empty() && t.options.empty() && t.ordered.empty())
      return os;

   t.displayFlagPattern(os);
   t.displayOptions(os);
   t.displayOrdered(os);

   os << endl << "Flags:" << endl;
   map<string, ArgumentProcessorFlag>::const_iterator f;
   for (f = t.flags.begin(); f != t.flags.end(); f++)
      if (f->second.preargument.length() == 1)
         os << "   " << f->second << endl;
   for (f = t.flags.begin(); f != t.flags.end(); f++)
      if (f->second.preargument.length() != 1)
         os << "   " << f->second << endl;

   os << endl << "Options:" << endl;
   map<string, ArgumentProcessorOption>::const_iterator o;
   for (o = t.options.begin(); o != t.options.end(); o++)
      os << "   " << o->second << endl;

   os << endl << "Ordered arguments:" << endl;
   vector<ArgumentProcessorOrderedArgument>::const_iterator a;
   for (a = t.ordered.begin(); a != t.ordered.end(); a++)
      os << "   " << *a << endl;

   return os;
}

string ArgumentProcessor::tostr(double d) const
{
   stringstream str;
   str << d;
   return str.str();
}

string ArgumentProcessor::tostr(int i) const
{
   stringstream str;
   str << i;
   return str.str();
}

string ArgumentProcessor::tostr(long i) const
{
   stringstream str;
   str << i;
   return str.str();
}

void ArgumentProcessor::splitForArgv(const string &argString, vector<string> &argv)
{
   argv.clear();
   int p = 0;
   int q = argString.find_first_of(" ");
   while (q != string::npos)
   {
      argv.push_back(argString.substr(p, q - p));
      p = q + 1;
      q = argString.find(" ", p);
   }
   argv.push_back(argString.substr(p, argString.length() - p));
}

#include <testFramework.h>

#include <argprocessor.h>

using namespace std;

class testArgumentProcessor : public ArgumentProcessor
{
public:
   testArgumentProcessor(const string& name) : ArgumentProcessor(name) { }
   friend bool testArgprocessor();
};

testArgumentProcessor* getArgumentProcessor()
{
   testArgumentProcessor *b = new testArgumentProcessor("asdf");

   b->addFlag("f", "flag flag");
   b->addFlag("g", "gee flag");
   b->addFlag("hold", "hold on flag");
   b->addFlag("axe", "axe flag");

   b->addOption("arc", "len", "arc length len");
   b->addOption("e", "en", "energy en");

   b->addOrdered("file1", "first file");
   b->addOrdered("binary", "second file, binary");
   b->addOrdered("richText3", "third file, rich text");

   return b;
}

bool testArgprocessor()
{
   testArgumentProcessor a("asdf");
   a.addFlag("f", "flag flag");
   a.addFlag("g", "gee flag");
   a.addFlag("hold", "hold on flag");
   a.addFlag("axe", "axe flag");

   a.addOption("arc", "len", "arc length len");
   a.addOption("e", "en", "energy en");

   a.addOrdered("file1", "first file");
   a.addOrdered("binary", "second file, binary");
   a.addOrdered("richText3", "third file, rich text");

   ASSERT(a.flags.size() == 4);
   ASSERT(a.options.size() == 2);
   ASSERT(a.ordered.size() == 3);

   ArgumentProcessorFlag f1 = a.flags["f"];
   ArgumentProcessorFlag f2 = a.flags["g"];
   ArgumentProcessorFlag f3 = a.flags["hold"];
   ArgumentProcessorFlag f4 = a.flags["axe"];

   ASSERT(f1.preargument == "f");
   ASSERT(f1.shortDescription == "");
   ASSERT(f1.description == "flag flag");

   ASSERT(f2.preargument == "g");
   ASSERT(f2.shortDescription == "");
   ASSERT(f2.description == "gee flag");

   ASSERT(f3.preargument == "hold");
   ASSERT(f3.shortDescription == "");
   ASSERT(f3.description == "hold on flag");

   ASSERT(f4.preargument == "axe");
   ASSERT(f4.shortDescription == "");
   ASSERT(f4.description == "axe flag");

   ArgumentProcessorOption o1 = a.options["arc"];
   ArgumentProcessorOption o2 = a.options["e"];

   ASSERT(o1.preargument == "arc");
   ASSERT(o1.shortDescription == "len");
   ASSERT(o1.description == "arc length len");

   ASSERT(o2.preargument == "e");
   ASSERT(o2.shortDescription == "en");
   ASSERT(o2.description == "energy en");

   ArgumentProcessorOrderedArgument a1 = a.ordered[0];
   ArgumentProcessorOrderedArgument a2 = a.ordered[1];
   ArgumentProcessorOrderedArgument a3 = a.ordered[2];

   ASSERT(a1.order == 0);
   ASSERT(a1.shortDescription == "file1");
   ASSERT(a1.description == "first file");

   ASSERT(a2.order == 1);
   ASSERT(a2.shortDescription == "binary");
   ASSERT(a2.description == "second file, binary");

   ASSERT(a3.order == 2);
   ASSERT(a3.shortDescription == "richText3");
   ASSERT(a3.description == "third file, rich text");

   // create a new argument set with similar arguments as a
   testArgumentProcessor *b = getArgumentProcessor();

   ASSERT(b->flags.size() == 4);
   ASSERT(b->options.size() == 2);
   ASSERT(b->ordered.size() == 3);

   map<string, ArgumentProcessorFlag>::iterator bflag;
   for (bflag = b->flags.begin(); bflag != b->flags.end(); bflag++)
      ASSERT(!bflag->second.found);

   ASSERT(!b->processFlag("no such flag"));
   for (bflag = b->flags.begin(); bflag != b->flags.end(); bflag++)
      ASSERT(!bflag->second.found);

   ASSERT(b->processFlag("hold"));
   ASSERT(b->flags["hold"].found);

   ASSERT(b->processMulticharFlag("--axe"));
   ASSERT(b->flags["axe"].found);

   ASSERT(b->processSinglecharFlag('f'));
   ASSERT(b->flags["f"].found);

   ASSERT(b->processSinglecharFlags("-fg"));
   ASSERT(b->flags["g"].found);
   for (bflag = b->flags.begin(); bflag != b->flags.end(); bflag++)
      ASSERT(bflag->second.found);

   map<string, ArgumentProcessorOption>::iterator boption;
   for (boption = b->options.begin(); boption != b->options.end(); boption++)
      ASSERT(!boption->second.found);

   ASSERT(!b->processOption("--nonesuch", "asdf"));

   ASSERT(b->processOption("-e", "20"));
   ASSERT(b->options["e"].found);
   ASSERT(b->options["e"].content == "20");

   ASSERT(b->processOption("--e", "30"));
   ASSERT(b->options["e"].found);
   ASSERT(b->options["e"].content == "30");

   ASSERT(b->processOption("--arc", "73"));
   ASSERT(b->options["arc"].found);
   ASSERT(b->options["arc"].content == "73");

   for (boption = b->options.begin(); boption != b->options.end(); boption++)
      ASSERT(boption->second.found);

   for (int i = 0; i < b->ordered.size(); i++)
      ASSERT(!b->ordered[i].found);

   for (int i = 0; i < b->ordered.size(); i++)
   {
      stringstream ss;

      ss << "file" << i;
      string number;
      ss >> number;
      ASSERT(b->processOrdered(number, i));
      ASSERT(b->ordered[i].found);
      ASSERT(b->ordered[i].content == number);
   }

   delete b;

   b = getArgumentProcessor();

   vector<string> arg;

   ASSERT(b->process(arg));

   arg.push_back("asdf");
   ASSERT(b->process(arg));

   arg.push_back("qwer");
   ASSERT(b->process(arg));

   arg.push_back("zcxv");
   ASSERT(b->process(arg));

   arg.push_back("asdf");
   ASSERT(!b->process(arg));

   delete b;

   b = getArgumentProcessor();
   string argstring("asdf qwer zxcv");
   b->splitForArgv(argstring, arg);

   ASSERT(b->process(arg));

   ASSERT(b->ordered[0].found && b->ordered[0].content == "asdf");

   ASSERT(b->ordered[1].found && b->ordered[1].content == "qwer");

   ASSERT(b->ordered[2].found && b->ordered[2].content == "zxcv");

   delete b;

   b = getArgumentProcessor();
   argstring = "--hold asdf -fg qwer --axe zxcv";
   b->splitForArgv(argstring, arg);

   ASSERT(b->process(arg));

   ASSERT(b->ordered[0].found && b->ordered[0].content == "asdf");

   ASSERT(b->ordered[1].found && b->ordered[1].content == "qwer");

   ASSERT(b->ordered[2].found && b->ordered[2].content == "zxcv");

   ASSERT(b->flags["axe"].found);

   ASSERT(b->flags["hold"].found);

   ASSERT(b->flags["f"].found);

   ASSERT(b->flags["g"].found);

   delete b;

   b = getArgumentProcessor();
   argstring = "--hold -e asdf -fg qwer --axe --arc zxcv";
   b->splitForArgv(argstring, arg);

   ASSERT(b->process(arg));

   ASSERT(b->ordered[0].found && b->ordered[0].content == "qwer");

   ASSERT(!b->ordered[1].found);

   ASSERT(!b->ordered[2].found);

   ASSERT(b->flags["axe"].found);

   ASSERT(b->flags["hold"].found);

   ASSERT(b->flags["f"].found);

   ASSERT(b->flags["g"].found);

   ASSERT(b->options["arc"].found && b->options["arc"].content == "zxcv");

   ASSERT(b->options["e"].found && b->options["e"].content == "asdf");

   delete b;
}

/****************************************************************
 phase-trio.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Essex.H"
#include "BOOM/VcfReader.H"
using namespace std;
using namespace BOOM;

class Application {
  void phase(const Genotype &mother,const Genotype &father,Genotype &child);
  Genotype getEssexGT(Essex::Node *siteGenotype,String label);
  int getEssexNumericChild(Essex::CompositeNode *,int whichChild);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("phase-trio <in.essex> <out.essex>");
  const String infile=cmd.arg(0);
  const String outfile=cmd.arg(1);

  // Process each gene in the input file
  Essex::Parser parser(infile);
  Essex::Node *root;
  while(root=parser.nextElem()) {
    Vector<Essex::Node*> sites;
    root->findDescendents("genotypes",sites);
    for(Vector<Essex::Node*>::iterator cur=sites.begin(), end=sites.end() ;
	cur!=end ; ++cur) {
      Essex::Node *site=*cur;
      Genotype mother=getEssexGT(site,"mother");
      Genotype father=getEssexGT(site,"father");
      Genotype child=getEssexGT(site,"child");
      phase(mother,father,child);
    }
  }

  return 0;
}



int Application::getEssexNumericChild(Essex::CompositeNode *node,int which)
{
  Essex::NumericNode *child=
    dynamic_cast<Essex::NumericNode*>(node->getIthChild(which));
  if(!child) throw "Child is not numeric in getEssexNumericChild";
  return child->getValue();
}


Genotype Application::getEssexGT(Essex::Node *site,String label)
{
  Essex::CompositeNode *child=
    dynamic_cast<Essex::CompositeNode*>(site->findChild(label));
  if(!child) throw label+" not found in getEssexGT";
  if(child->getNumChildren()!=2) throw label+" has wrong number of children";
  Genotype g;
  g.addAllele(getEssexNumericChild(child,0));
  g.addAllele(getEssexNumericChild(child,1));
  return g;
}



void Application::phase(const Genotype &mother,const Genotype &father,
			Genotype &child)
{
  cout<<mother<<"\t"<<father<<"\t"<<child<<endl;
}





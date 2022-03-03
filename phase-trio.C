/****************************************************************
 phase-trio.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Essex.H"
#include "BOOM/VcfReader.H"
#include "BOOM/Map.H"
using namespace std;
using namespace BOOM;

class Application {
  Map<String,String> phasingMap;
  void initMap(Map<String,String> &);
  bool phase(Genotype &mother,Genotype &father,Genotype &child);
  Genotype getEssexGT(Essex::Node *siteGenotype,String label);
  void installGT(Essex::Node *,const String &label,const Genotype &);
  int getEssexNumericChild(Essex::CompositeNode *,int whichChild);
  void setEssexNumericChild(Essex::CompositeNode *,int which,int value);
  String compactString(const Genotype &);
  void install(const String &encoded,int first,int second,Genotype &);
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

  initMap(phasingMap);
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("phase-trio <in.essex> <out.essex>");
  const String infile=cmd.arg(0);
  const String outfile=cmd.arg(1);

  // Create output file
  ofstream os(outfile);
  
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
      cout<<mother<<" "<<father<<" "<<child<<"  =>  ";
      bool success=phase(mother,father,child);
      cout<<mother<<" "<<father<<" "<<child<<endl;
      installGT(site,"child",child);
      installGT(site,"mother",mother);
      installGT(site,"father",father);
    }
    root->printOn(os); os<<endl;
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



void Application::setEssexNumericChild(Essex::CompositeNode *node,int which,
				      int value)
{
  Essex::NumericNode *child=
    dynamic_cast<Essex::NumericNode*>(node->getIthChild(which));
  if(!child) throw "Child is not numeric in getEssexNumericChild";
  child->setValue(value);
}



void Application::installGT(Essex::Node *site,const String &label,
			    const Genotype &G)
{
  Essex::CompositeNode *child=
    dynamic_cast<Essex::CompositeNode*>(site->findChild(label));
  if(!child) throw label+" not found in getEssexGT";
  if(child->getNumChildren()!=2) throw label+" has wrong number of children";
  setEssexNumericChild(child,0,G[0]);
  setEssexNumericChild(child,1,G[1]);
}



Genotype Application::getEssexGT(Essex::Node *site,String label)
{
  Essex::CompositeNode *child=
    dynamic_cast<Essex::CompositeNode*>(site->findChild(label));
  if(!child) throw label+" not found in getEssexGT";
  if(child->getNumChildren()!=2) throw label+" has wrong number of children";
  Genotype g;
  int allele1=getEssexNumericChild(child,0);
  int allele2=getEssexNumericChild(child,1);
  g.addAllele(allele1);
  g.addAllele(allele2);
  return g;
}



String Application::compactString(const Genotype &G)
{
  String s="";
  const int numAlleles=G.numAlleles();
  for(int i=0 ; i<numAlleles ; ++i) s+=String(G[i]);
  if(s=="10") s="01"; // For unphased genotypes, these are equivalent
  return s;
}



void Application::initMap(Map<String,String> &M)
{
  // PRECONDITION: each individual has been normalized, so that 10
  // doesn't occur in the inputs (for an individual)

  // The first two columns are mother, second pair is father, third is child
  
  M["000000"]="000000";
  M["000100"]="000100";
  M["000101"]="001001";
  M["001101"]="001101";
  M["010000"]="010000";
  M["010001"]="100010";
  M["010100"]="010100";
  M["010111"]="101011";
  M["011101"]="011101";
  M["011111"]="101111";
  M["110001"]="110010";
  M["110101"]="110110";
  M["110111"]="111011";
  M["111111"]="111111";

  //M["010101"]=""; TRIP HET, CAN'T PHASE
  //M["110000"]=""; DE NOVO
  //M["000001"]=""; DE NOVO, CAN'T PHASE
  //M["000010"]=""; DE NOVO, CAN'T PHASE
  //M["000011"]=""; TWO DE NOVOS
  //M["000110"]=""; UNNORMALIZED
  //M["000111"]=""; DE NOVO
  //M["001100"]=""; DE NOVO
  //M["001000"]=""; UNNORMALIZED
  //M["001001"]=""; UNNORMALIZED
  //M["001010"]=""; UNNORMALIZED
  //M["001011"]=""; UNNORMALIZED
  //M["001110"]=""; UNNORMALIZED
  //M["001111"]=""; DE NOVO
  //M["010010"]=""; UNNORMALIZED
  //M["010011"]=""; DE NOVO
  //M["010110"]=""; UNNORMALIZED
  //M["011000"]=""; UNNORMALIZED
  //M["011001"]=""; UNNORMALIZED
  //M["011010"]=""; UNNORMALIZED
  //M["011011"]=""; UNNORMALIZED
  //M["011100"]=""; DE NOVO
  //M["011110"]=""; UNNORMALIZED
  //M["100000"]=""; UNNORMALIZED
  //M["100001"]=""; UNNORMALIZED
  //M["100010"]=""; UNNORMALIZED
  //M["100011"]=""; UNNORMALIZED
  //M["100100"]=""; UNNORMALIZED
  //M["100101"]=""; UNNORMALIZED
  //M["100110"]=""; UNNORMALIZED
  //M["100111"]=""; UNNORMALIZED
  //M["101000"]=""; UNNORMALIZED
  //M["101001"]=""; UNNORMALIZED
  //M["101010"]=""; UNNORMALIZED
  //M["101011"]=""; UNNORMALIZED
  //M["101100"]=""; UNNORMALIZED
  //M["101101"]=""; UNNORMALIZED
  //M["101110"]=""; UNNORMALIZED
  //M["101111"]=""; UNNORMALIZED
  //M["110010"]=""; UNNORMALIZED
  //M["110011"]=""; DE NOVO
  //M["110100"]=""; DE NOVO
  //M["110110"]=""; UNNORMALIZED
  //M["111000"]=""; UNNORMALIZED
  //M["111001"]=""; UNNORMALIZED
  //M["111010"]=""; UNNORMALIZED
  //M["111011"]=""; UNNORMALIZED
  //M["111100"]=""; DE NOVO
  //M["111101"]=""; DE NOVO
  //M["111110"]=""; UNNORMALIZED
}



bool Application::phase(Genotype &mother,Genotype &father,Genotype &child)
{
  // PRECONDITION: all three genotypes have two alleles each

  //cout<<mother<<"\t"<<father<<"\t"<<child<<endl;
  const String encoded=
    compactString(mother)+compactString(father)+compactString(child);
  if(encoded=="010101") return false; // triple het, can't phase
  if(!phasingMap.isDefined(encoded))
    throw String("Genotype encoding is not defined: ")+encoded;
  String phased=phasingMap[encoded];
  //cout<<encoded<<" => "<<phased<<endl;
  install(phased,0,1,mother);
  install(phased,2,3,father);
  install(phased,4,5,child);
  return true;
}



void Application::install(const String &encoded,int first,int second,
			  Genotype &G)
{
  Vector<int> &v=G.asVector();
  const char *p=encoded.c_str();
  v[0]=p[first]-'0';
  v[1]=p[second]-'0';
}




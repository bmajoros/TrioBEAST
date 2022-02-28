/****************************************************************
 sim1.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/VcfReader.H"
#include "BOOM/GSL/Random.H"
using namespace std;
using namespace BOOM;

class Application {
  int motherIndex, fatherIndex; // Indices in VCF #CHROM line
  bool motherMaternal, fatherMaternal; // Which copy is inherited
  VariantAndGenotypes simNext(VcfReader &);
  Genotype simGenotype(const Genotype &mother,const Genotype &father);
  bool sameHomozygotes(const Genotype &mother,const Genotype &father);
  void chooseInheritedCopies();
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

  GSL::Random::randomize();
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=6)
    throw String("sim1 <in.vcf> <mother-ID> <father-ID> <#genes> <variants-per-gene> <theta>");
  const String VCF_FILE=cmd.arg(0);
  const String MOTHER_ID=cmd.arg(1);
  const String FATHER_ID=cmd.arg(2);
  const int NUM_GENES=cmd.arg(3).asInt();
  const int VARIANTS_PER_GENE=cmd.arg(4).asInt();
  const float THETA=cmd.arg(5).asFloat();
  
  VcfReader reader(VCF_FILE);
  reader.hashSampleIDs();
  motherIndex=reader.getSampleIndex(MOTHER_ID);
  fatherIndex=reader.getSampleIndex(FATHER_ID);
  for(int geneNum=0 ; geneNum<NUM_GENES ; ++geneNum) {
    chooseInheritedCopies();
    for(int varNum=0 ; varNum<VARIANTS_PER_GENE ; ++varNum) {
      VariantAndGenotypes vg=simNext(reader);
      cout<<vg.variant<<"\t"<<vg.genotypes[0]<<"\t"<<vg.genotypes[1]
	  <<"\t"<<vg.genotypes[2]<<endl;
    }
  }

  return 0;
}



void Application::chooseInheritedCopies()
{
  // This function decides whether the child inherits his mother's
  // maternal vs. paternal copy, and similarly for his father's maternal/
  // paternal copy.
  
  motherMaternal=GSL::Random::randomBool();
  fatherMaternal=GSL::Random::randomBool();
}



VariantAndGenotypes Application::simNext(VcfReader &reader)
{
  VariantAndGenotypes vg;
  while(true) {
    if(!reader.nextVariant(vg)) {
      reader.rewind();
      if(!reader.nextVariant(vg)) throw "Cannot rewind in simNext()";
    }
    const Variant v=vg.variant;
    if(v.containsNonstandardAlleles() || v.isIndel() ||
       v.numAlleles()!=2) continue;
    const Vector<Genotype> &genotypes=vg.genotypes;
    const Genotype motherGT=genotypes[motherIndex];
    const Genotype fatherGT=genotypes[fatherIndex];
    if(sameHomozygotes(motherGT,fatherGT)) continue;
    const Genotype childGT=simGenotype(motherGT,fatherGT);
    VariantAndGenotypes r;
    r.variant=v;
    r.genotypes.push_back(motherGT);
    r.genotypes.push_back(fatherGT);
    r.genotypes.push_back(childGT);
    return r;
  }
}



bool Application::sameHomozygotes(const Genotype &mother,
				  const Genotype &father)
{
  return !mother.isHet() && !father.isHet() && mother[0]==father[0];
}



Genotype Application::simGenotype(const Genotype &mother,
				  const Genotype &father)
{
  Genotype child;
  child.addAllele(motherMaternal ? mother[0] : mother[1]);
  child.addAllele(fatherMaternal ? father[0] : father[1]);
  return child;
}





/****************************************************************
 sim1.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/VcfReader.H"
#include "BOOM/GSL/Random.H"
#include "BOOM/GSL/GslBinomial.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
using namespace std;
using namespace BOOM;

enum Individual { MOTHER=0, FATHER=1, CHILD=2 };
enum Allele { REF=0, ALT=1 };
enum MaternalPaternal { MAT=0, PAT=1 };

class Application {
  int motherIndex, fatherIndex; // Indices in VCF #CHROM line
  Array1D<bool> parentMaternal; // For this parent, which copy is passed down
  Array2D<int> V; // Affected status; indexed as: V[individual][mat/pat]
  Array3D<int> counts; // Indexed as V[site][individual][ref/alt]
  float RECOMB; // Recombination rate (between gene and causal variant)
  Array1D<bool> recombined; // indexed by Individual
  VariantAndGenotypes simNext(VcfReader &);
  Genotype simGenotype(const Genotype &mother,const Genotype &father);
  bool sameHomozygotes(const Genotype &mother,const Genotype &father);
  void chooseInheritedCopies();
  void simAffectedStatus();
  void unphaseGenotypes(Vector<Genotype> &); // mother/father/child
  void swapGenotype(Genotype &);
  void printGenotypes(const Vector<Genotype> &,ostream &);
  void writeTruth(int geneNum,const Vector<VariantAndGenotypes> &,
		  float theta,int total,ostream &);
  //void writeVCF(int geneNum,const Vector<VariantAndGenotypes> &,ostream &);
  void writeData(int geneNum,const Vector<VariantAndGenotypes> &,
		 const Array3D<int> &counts,const int total,ostream &);
  void writePhasedSites(const Vector<VariantAndGenotypes> &,const float theta,
			const int total,ostream &);
  void writeUnphasedSites(const Vector<VariantAndGenotypes> &,
			  const int total,ostream &);
  int inheritAffected(Individual parent);
  void swapCopy(int &c);
  String phasedV(BOOM::Array2D<int>::RowIn2DArray<int> row);
  String unphasedV(BOOM::Array2D<int>::RowIn2DArray<int> row);
  String phased(const Genotype &);
  String unphased(const Genotype &);
  void simCounts(const float theta,const Vector<VariantAndGenotypes> &,
		 const int total,Array3D<int> &counts);
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
  : V(3,2), parentMaternal(2), recombined(2)
{
  // ctor

  GSL::Random::randomize();
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=10)
    throw String("sim1 <in.vcf> <mother-ID> <father-ID> <#genes> <variants-per-gene> <reads-per-site> <recombination-rate> <theta> <out-truth.essex> <out-data.essex>");
  const String VCF_FILE=cmd.arg(0);
  const String MOTHER_ID=cmd.arg(1);
  const String FATHER_ID=cmd.arg(2);
  const int NUM_GENES=cmd.arg(3).asInt();
  const int VARIANTS_PER_GENE=cmd.arg(4).asInt();
  const int READS_PER_SITE=cmd.arg(5).asInt();
  RECOMB=cmd.arg(6).asFloat();
  const float THETA=cmd.arg(7).asFloat();
  const String truthFileName=cmd.arg(8);
  const String dataFileName=cmd.arg(9);

  counts.resize(VARIANTS_PER_GENE,3,2); // sites, individuals, alleles
  ofstream truthFile(truthFileName), dataFile(dataFileName);
  VcfReader reader(VCF_FILE);
  reader.hashSampleIDs();
  motherIndex=reader.getSampleIndex(MOTHER_ID);
  fatherIndex=reader.getSampleIndex(FATHER_ID);
  for(int geneNum=0 ; geneNum<NUM_GENES ; ++geneNum) {
    simAffectedStatus();
    chooseInheritedCopies();
    Vector<VariantAndGenotypes> variants;
    for(int varNum=0 ; varNum<VARIANTS_PER_GENE ; ++varNum)
      variants.push_back(simNext(reader));
    simCounts(THETA,variants,READS_PER_SITE,counts);
    writeTruth(geneNum,variants,THETA,READS_PER_SITE,truthFile);
    writeData(geneNum,variants,counts,READS_PER_SITE,dataFile);
  }

  return 0;
}



String Application::phasedV(BOOM::Array2D<int>::RowIn2DArray<int> row)
{
  // This has been revised to omit the "|" since Essex can't handle that
  
  return String(row[0])+" "+String(row[1]);
}



String Application::unphasedV(BOOM::Array2D<int>::RowIn2DArray<int> row)
{
  // This has been revised to omit the "/" since Essex can't handle that

  return String(row[0])+" "+String(row[1]);
}



void Application::writePhasedSites(const Vector<VariantAndGenotypes> &variants,
				  const float theta,const int total,
				  ostream &os)
{
  const int numSites=variants.size();
  for(int i=0 ; i<numSites ; ++i) {
    const VariantAndGenotypes &vg=variants[i];
    const Variant &v=vg.variant;
    const Vector<Genotype> &genotypes=vg.genotypes;
    os<<"\t(site "<<i<<" ";
    os<<"(genotypes "
      <<"(mother "<<phased(genotypes[MOTHER])<<") "
      <<"(father "<<phased(genotypes[FATHER])<<") "
      <<"(child "<<phased(genotypes[CHILD])<<"))\n";
    os<<"\t\t(counts "
      <<"(mother "<<counts[i][MOTHER][REF]<<" "<<counts[i][MOTHER][ALT]<<") "
      <<"(father "<<counts[i][FATHER][REF]<<" "<<counts[i][FATHER][ALT]<<") "
      <<"(child "<<counts[i][CHILD][REF]<<" "<<counts[i][CHILD][ALT]<<")))"
      <<endl;
  }
}



void Application::writeUnphasedSites(const Vector<VariantAndGenotypes> &
				    const_variants,
				    const int total,ostream &os)
{
  Vector<VariantAndGenotypes> variants=const_variants;
  const int numSites=variants.size();
  for(int i=0 ; i<numSites ; ++i) {
    VariantAndGenotypes &vg=variants[i];
    //    const Variant &v=vg.variant;
    Vector<Genotype> &genotypes=vg.genotypes;
    unphaseGenotypes(genotypes);
    os<<"\t(site "<<i<<" ";
    os<<"(genotypes "
      <<"(mother "<<unphased(genotypes[MOTHER])<<") "
      <<"(father "<<unphased(genotypes[FATHER])<<") "
      <<"(child "<<unphased(genotypes[CHILD])<<"))\n";
    os<<"\t\t(counts "
      <<"(mother "<<counts[i][MOTHER][REF]<<" "<<counts[i][MOTHER][ALT]<<") "
      <<"(father "<<counts[i][FATHER][REF]<<" "<<counts[i][FATHER][ALT]<<") "
      <<"(child "<<counts[i][CHILD][REF]<<" "<<counts[i][CHILD][ALT]<<")))"
      <<endl;
  }
}



void Application::writeTruth(const int geneNum,
			     const Vector<VariantAndGenotypes> &variants,
			     float theta,int total,ostream &os)
{
  os<<"(gene GENE"<<geneNum<<endl;
  os<<"\t(theta "<<theta<<")"<<endl;
  os<<"\t(affected (mother "<<phasedV(V[MOTHER])<<") (father "
    <<phasedV(V[FATHER])<<") (child "<<phasedV(V[CHILD])<<"))"<<endl;
  os<<"\t(inherited_copy (mother "<<(parentMaternal[MOTHER]?0:1)
    <<") (father "<<(parentMaternal[FATHER]?0:1)<<"))"<<endl;
  os<<"\t(recombined (mother "<<(recombined[MOTHER]?1:0)
    <<") (father "<<(recombined[FATHER]?1:0)<<"))"<<endl;
  writePhasedSites(variants,theta,total,os);
  os<<")"<<endl;
}



void Application::writeData(int geneNum,
			    const Vector<VariantAndGenotypes> &variants,
			    const Array3D<int> &counts,const int total,
			    ostream &os)
{
  os<<"(gene GENE"<<geneNum<<endl;
  writeUnphasedSites(variants,total,os);
  os<<")"<<endl;
}



/*void Application::writeVCF(const int geneNum,
			   const Vector<VariantAndGenotypes> &const_variants,
			   ostream &os)
{
  os<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMOTHER\tFATHER\tCHILD"<<endl;
  Vector<VariantAndGenotypes> variants=const_variants;
  const int numSites=variants.size();
  for(int site=0 ; site<numSites ; ++site) {
    VariantAndGenotypes &vg=variants[site];
    Variant &v=vg.variant;
    unphaseGenotypes(vg.genotypes);
    os<<v.getChr()<<"\t"<<v.getPos()<<"\t"<<v.getID()<<"\t"
      <<v.getAllele(0)<<"\t"<<v.getAllele(1)<<"\t100\tPASS\tVT=SNP\tGT";
    for(int indiv=0 ; indiv<3 ; ++indiv) {
      const Genotype &g=vg.genotypes[indiv];
      os<<"\t"<<g[0]<<"/"<<g[1];
    }
    os<<endl;
  }
  }*/



void Application::chooseInheritedCopies()
{
  // This function decides whether the child inherits his mother's
  // maternal vs. paternal copy, and similarly for his father's maternal/
  // paternal copy.

  //parentMaternal[MOTHER]=GSL::Random::randomBool();
  //parentMaternal[FATHER]=GSL::Random::randomBool();

  // UPDATE: Without loss of generality, and to simplify the representation,
  // we always assume each parent passes its maternal gene copy to the child
  parentMaternal[MOTHER]=true;
  parentMaternal[FATHER]=true;
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
    if(sameHomozygotes(motherGT,fatherGT)) continue; // triple homozygote
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
  child.addAllele(parentMaternal[MOTHER] ? mother[0] : mother[1]);
  child.addAllele(parentMaternal[FATHER] ? father[0] : father[1]);
  return child;
}



void Application::simAffectedStatus()
{
  // First, simulate that exactly one parent has ASE:
  switch(GSL::Random::randomInt(1,4)) {
  case 1: // 10 00
    V[MOTHER][0]=1; V[MOTHER][1]=0; V[FATHER][0]=0; V[FATHER][1]=0; break;
  case 2: // 01 00
    V[MOTHER][0]=0; V[MOTHER][1]=1; V[FATHER][0]=0; V[FATHER][1]=0; break;
  case 3: // 00 10
    V[MOTHER][0]=0; V[MOTHER][1]=0; V[FATHER][0]=1; V[FATHER][1]=0; break;
  case 4: // 00 01
    V[MOTHER][0]=0; V[MOTHER][1]=0; V[FATHER][0]=0; V[FATHER][1]=1; break;
  default: throw "Error in simAffectedStatus()";
  }

  // Now simulate what the child inherits:
  V[CHILD][0]=inheritAffected(MOTHER);
  V[CHILD][1]=inheritAffected(FATHER);
}



void Application::swapCopy(int &c)
{
  switch(c) {
  case 0: c=1; break;
  case 1: c=0; break;
  default: throw "Error in swapCopy()";
  }
}



int Application::inheritAffected(Individual parent)
{
  int parentsCopy=parentMaternal[parent] ? 0 : 1;
  if(GSL::Random::randomFloat(0,1)<RECOMB) {
    swapCopy(parentsCopy);
    recombined[parent]=true;
  }
  else recombined[parent]=false;
  return V[parent][parentsCopy];
}



void Application::swapGenotype(Genotype &g)
{
  Vector<int> &v=g.asVector();
  const int temp=v[0];
  v[0]=v[1];
  v[1]=temp;
}


void Application::unphaseGenotypes(Vector<Genotype> &G)
{
  // This function randomizes the phasing of the genotypes (M/F/C)

  for(int i=0 ; i<3 ; ++i)
    if(GSL::Random::randomBool()) swapGenotype(G[i]);
}



void Application::printGenotypes(const Vector<Genotype> &G,ostream &os)
{
  for(int i=0 ; i<3 ; ++i) {
    const Genotype &g=G[i];
    os<<g[0]<<'/'<<g[1];
    if(i+1<3) os<<'\t';
  }
}



String Application::phased(const Genotype &g)
{
  return String(g[0])+" "+String(g[1]);
}




String Application::unphased(const Genotype &g)
{
  return String(g[0])+" "+String(g[1]);
}



void Application::simCounts(const float theta,
			    const Vector<VariantAndGenotypes> &variants,
			    const int N,Array3D<int> &counts)
{
  const int numSites=variants.size();
  for(int site=0 ; site<numSites ; ++site) {
    const VariantAndGenotypes &vg=variants[site];
    const Vector<Genotype> &genotypes=vg.genotypes;
    for(int indiv=0 ; indiv<3 ; ++indiv) {
      const Genotype &g=genotypes[indiv];
      if(!g.isHet()) {
	const int allele=g[0]; const int otherAllele=1-allele;
	counts[site][indiv][allele]=N; counts[site][indiv][otherAllele]=0;
	continue; }
      const bool hasASE=V[indiv][0]!=V[indiv][1];
      const float myTheta=hasASE ? theta : 1.0;
      const float p=myTheta/(myTheta+1);
      GSL::GslBinomial binom(p);
      const int maternal=binom.random(N);
      const int alt=g[MAT]==ALT ? maternal : N-maternal;
      const int ref=N-alt;
      counts[site][indiv][REF]=ref;
      counts[site][indiv][ALT]=alt;
    }
  }
}





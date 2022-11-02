/****************************************************************
 sim-ped-genotypes.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include <queue>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/VcfReader.H"
#include "BOOM/GSL/Random.H"
#include "BOOM/GSL/GslBinomial.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
#include "BOOM/Set.H"
using namespace std;
using namespace BOOM;

/****************************************************************
                       enums and utilities
 ****************************************************************/
enum Allele { REF=0, ALT=1 };
enum MaternalPaternal { MAT=0, PAT=1, MATPAT_UNKNOWN=2 };
enum MotherFather { MOTHER=0, FATHER=1, MF_UNKNOWN=2 };
enum Sex { FEMALE=0, MALE=1, SEX_UNKNOWN };
Sex stringToSex(const String &);
ostream &operator<<(ostream &,Sex);

/****************************************************************
                         class Individual
 ****************************************************************/
class Individual {
public:
  Individual(const String &id,Sex,const String &motherID,
	     const String &fatherID);
  const String &getID() const;
  Individual *getParent(MotherFather) const;
  void setParent(MotherFather,Individual *);
  const String &getParentID(MotherFather) const;
  Sex getSex() const;
  void addChild(Individual *);
  Vector<Individual*> &getChildren();
  virtual bool isRoot() const;
  bool isLeaf() const;
  MaternalPaternal &getInherit(MotherFather); // which copy I inherit
  void resizeGenotypes(int);
  Array1D<Genotype> &getGenotypes();
private:
  String ID;
  Sex sex;
  Array1D<String> parentID;
  Array1D<Individual*> parents;
  Vector<Individual*> children;
  Array1D<MaternalPaternal> inherit; // which copy I inherit from Mom and Dad
  Array1D<Genotype> genotypes;
};

/****************************************************************
                           class Root
 ****************************************************************/
class Root : public Individual {
public:
  Root(const String &id,Sex);
  void setVcfID(const String &);
  const String &getVcfID() const;
  virtual bool isRoot() const;
private:
  String vcfID;
};

/****************************************************************
                         class Pedigree
 ****************************************************************/
class Pedigree {
public:
  Pedigree();
  virtual ~Pedigree();
  static Pedigree *loadFromTextFile(const String &filename);
  void addIndividual(Individual *);
  int size() const;
  Individual *operator[](int);
  void installPointers();
  Individual *findIndiv(const String ID) const;
  void printOn(ostream &) const;
  void topologicalSort(Vector<Individual*> &into);
  void simInherit();
  void getRoots(Vector<Root*> &into);
  void resizeGenotypes(int);
private:
  Vector<Individual*> individuals;
  void dfs(Individual *,Set<Individual*> &seen,Vector<Individual*> &stack);
};
ostream &operator<<(ostream &,const Pedigree &);

/****************************************************************
                         class Application
 ****************************************************************/
class Application {
  Vector<String> vcfIDs;
  void setRootGenotypes(Vector<Root*> &,VcfReader &,
			const Vector<VariantAndGenotypes> &variants);
  void inherit(Vector<Individual*> &topsort);
public:
  Application();
  int main(int argc,char *argv[]);
};

/****************************************************************
                            main()
 ****************************************************************/
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



/****************************************************************
                          Application methods
 ****************************************************************/
Application::Application()
{
  // ctor

  GSL::Random::randomize();
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=7)
    throw String("sim1 <*.pedigree> <in.vcf> <ID1,ID2,ID3,...> <#genes> <variants-per-gene> <out-phased.txt> <out-unphased.txt>");
  const String PED_FILE=cmd.arg(0);
  const String VCF_FILE=cmd.arg(1);
  const String ID_LIST=cmd.arg(2);
  const int NUM_GENES=cmd.arg(3).asInt();
  const int VARIANTS_PER_GENE=cmd.arg(4).asInt();
  const String phasedFilename=cmd.arg(5);
  const String unphasedFilename=cmd.arg(6);

  // Parse VCF identifier list
  ID_LIST.getFields(vcfIDs,",");

  // Load pedigree and simulate meioses
  Pedigree *pedigree=Pedigree::loadFromTextFile(PED_FILE);
  Vector<Individual*> topsort;
  pedigree->topologicalSort(topsort);
  pedigree->simInherit();
  pedigree->resizeGenotypes(VARIANTS_PER_GENE);
  Vector<Root*> roots; pedigree->getRoots(roots);
  cout<<roots.size()<<" roots found"<<endl;
  if(roots.size()!=vcfIDs.size())
    throw RootException("Number of VCF IDs does not match number of roots");
  for(int i=0 ; i<roots.size() ; ++i)
    roots[i]->setVcfID(vcfIDs[i]);
  /*cout<<"TOPOLOGICAL SORT:"<<endl;
  for(Vector<Individual*>::iterator cur=topsort.begin(), end=topsort.end() ;
      cur!=end ; ++cur)
      cout<<(*cur)->getID()<<endl;*/
  
  // Process VCF file
  ofstream phasedFile(phasedFilename), unphasedFile(unphasedFilename);
  VcfReader reader(VCF_FILE);
  reader.hashSampleIDs();
  for(int geneNum=0 ; geneNum<NUM_GENES ; ++geneNum) {
    cout<<"Simulating gene "<<(geneNum+1)<<endl;
    Vector<VariantAndGenotypes> variants;
    VariantAndGenotypes vg;
    for(int varNum=0 ; varNum<VARIANTS_PER_GENE ; ++varNum)
      { reader.nextVariant(vg); variants.push_back(vg); }
    //cout<<variants.size()<<" variants loaded"<<endl;
    
    setRootGenotypes(roots,reader,variants);
    inherit(topsort);
  }

  return 0;
}



void Application::setRootGenotypes(Vector<Root*> &roots,
				   VcfReader &reader,
				   const Vector<VariantAndGenotypes> &variants)
{
  const int N_VAR=variants.size();
  for(Vector<Root*>::iterator cur=roots.begin(), end=roots.end() ;
      cur!=end ; ++cur) {
    Root *root=*cur;
    const int sampleIndex=reader.getSampleIndex(root->getVcfID());
    Array1D<Genotype> &genotypes=root->getGenotypes();
    for(int i=0 ; i<N_VAR ; ++i) {
      genotypes[i]=variants[i].genotypes[sampleIndex];
    }
  }
}



void Application::inherit(Vector<Individual*> &topsort)
{
  for(Vector<Individual*>::iterator cur=topsort.begin(), end=topsort.end() ;
      cur!=end ; ++cur) {
    Individual *ind=*cur;
    if(ind->isRoot()) continue;
    Array1D<Genotype> &genotypes=ind->getGenotypes();
    MaternalPaternal fromMother=ind->getInherit(MOTHER);
    MaternalPaternal fromFather=ind->getInherit(FATHER);
    Individual *mother=ind->getParent(MOTHER);
    Individual *father=ind->getParent(FATHER);
    Array1D<Genotype> &motherGenotypes=mother->getGenotypes();
    Array1D<Genotype> &fatherGenotypes=father->getGenotypes();
    if(!mother || !father) throw RootException("Missing parent");
    const int numVar=genotypes.size();
    for(int i=0 ; i<numVar ; ++i) {
      Genotype &genotype=genotypes[i], &motherGT=motherGenotypes[i],
	&fatherGT=fatherGenotypes[i];
      genotype[MAT]=motherGT[fromMother];
      genotype[PAT]=fatherGT[fromFather];
    }
  }
}



/****************************************************************
                         Individual methods
 ****************************************************************/
Individual::Individual(const String &id,Sex sex,const String &motherID,
		       const String &fatherID)
  : ID(id), sex(sex), parentID(2), parents(2), inherit(2)
{
  parentID[MOTHER]=motherID;
  parentID[FATHER]=fatherID;
  parents.setAllTo(NULL);
  inherit.setAllTo(MATPAT_UNKNOWN);
}



void Individual::resizeGenotypes(int s)
{
  genotypes.resize(s);

  // Intialize to -1/-1 so the vector is the right size
  for(int i=0 ; i<s ; ++i) {
    Genotype &g=genotypes[i];
    g.addAllele(-1);
    g.addAllele(-1);
  }
}



Array1D<Genotype> &Individual::getGenotypes()
{
  return genotypes;
}



MaternalPaternal &Individual::getInherit(MotherFather parent)
{
  return inherit[parent];
}



Vector<Individual*> &Individual::getChildren()
{
  return children;
}



const String &Individual::getID() const
{
  return ID;
}



Individual *Individual::getParent(MotherFather i) const
{
  return parents[i];
}



void Individual::setParent(MotherFather which,Individual *parent)
{
  parents[which]=parent;
}



const String &Individual::getParentID(MotherFather i) const
{
  return parentID[i];
}



Sex Individual::getSex() const
{
  return sex;
}



bool Individual::isRoot() const
{
  return false;
}



bool Individual::isLeaf() const
{
  return children.size()==0;
}



void Individual::addChild(Individual *ind)
{
  children.push_back(ind);
}



/****************************************************************
                         Pedigree methods
 ****************************************************************/
Pedigree::Pedigree()
{
  // ctor
}



Pedigree::~Pedigree()
{
  for(Vector<Individual*>::iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur)
    delete *cur;
}



void Pedigree::getRoots(Vector<Root*> &into)
{
  for(Vector<Individual*>::iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur) {
    Individual *ind=*cur;
    if(ind->isRoot()) into.push_back(dynamic_cast<Root*>(ind));
  }
}



void Pedigree::printOn(ostream &os) const
{
  os<<"ID\tSex\tMother\tFather"<<endl;
  for(Vector<Individual*>::const_iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur) {
    Individual *ind=*cur;
    Individual *mother=ind->getParent(MOTHER);
    Individual *father=ind->getParent(FATHER);
    os<<ind->getID()<<"\t"<<ind->getSex()
      <<"\t"<<(mother ? mother->getID() : ".")
      <<"\t"<<(father ? father->getID() : ".")
      <<endl;
  }
}



ostream &operator<<(ostream &os,const Pedigree &ped)
{
  ped.printOn(os);
  return os;
}



Pedigree *Pedigree::loadFromTextFile(const String &filename)
{
  File file(filename);
  String header=file.getline();
  Vector<String> fields;
  Pedigree *pedigree=new Pedigree();
  while(!file.eof()) {
    String line=file.getline();
    line.getFields(fields);
    if(fields.size()<4) continue;
    const String ID=fields[0];
    const Sex sex=stringToSex(fields[1]);
    const String motherID=fields[2];
    const String fatherID=fields[3];
    const bool isRoot=motherID=="." && fatherID==".";
    Individual *ind=
      isRoot ? new Root(ID,sex) : new Individual(ID,sex,motherID,fatherID);
    pedigree->addIndividual(ind);
  }
  pedigree->installPointers();
  return pedigree;
}



void Pedigree::addIndividual(Individual *ind)
{
  individuals.push_back(ind);
}



int Pedigree::size() const
{
  return individuals.size();
}



Individual *Pedigree::operator[](int i)
{
  return individuals[i];
}



Individual *Pedigree::findIndiv(const String ID) const
{
  // Linear search: inefficient, but only used a few times
  for(Vector<Individual*>::const_iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur) {
    Individual *ind=*cur;
    if(ind->getID()==ID) return ind;
  }
  return NULL;
}



void Pedigree::installPointers()
{
  for(Vector<Individual*>::iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur) {
    Individual *ind=*cur;
    Individual *mother=findIndiv(ind->getParentID(MOTHER));
    Individual *father=findIndiv(ind->getParentID(FATHER));
    if(mother) { ind->setParent(MOTHER,mother); mother->addChild(ind); }
    if(father) { ind->setParent(FATHER,father); father->addChild(ind); }
  }
}



void Pedigree::topologicalSort(Vector<Individual*> &into)
{
  Set<Individual*> seen;
  Vector<Individual*> stack;
  for(Vector<Individual*>::iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur)
    dfs(*cur,seen,stack);
  while(!stack.empty()) { into.push_back(stack.back()); stack.pop_back(); }
}



void Pedigree::dfs(Individual *ind,Set<Individual*> &seen,
		   Vector<Individual*> &stack)
{
  if(seen.isMember(ind)) return;
  seen+=ind;
  Vector<Individual*> &children=ind->getChildren();
  for(Vector<Individual*>::iterator cur=children.begin(), end=children.end() ;
      cur!=end ; ++cur)
    dfs(*cur,seen,stack);
  stack.push_back(ind);
}



void Pedigree::simInherit()
{
  // This method stochastically chooses which haplotype each individual
  // inherits from each parent: the parent's maternal haplotype, or the
  // parent's paternal haplotype
  for(Vector<Individual*>::iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur) {
    Individual *ind=*cur;
    ind->getInherit(MOTHER)=Random0to1()<0.5 ? MAT : PAT;
    ind->getInherit(FATHER)=Random0to1()<0.5 ? MAT : PAT;
  }
}



void Pedigree::resizeGenotypes(int s) {
  for(Vector<Individual*>::iterator cur=individuals.begin(),
	end=individuals.end() ; cur!=end ; ++cur)
    (*cur)->resizeGenotypes(s);
}


/****************************************************************
                          Root methods
 ****************************************************************/
Root::Root(const String &id,Sex s)
  : Individual(id,s,"","")
{
  //ctor
}



void Root::setVcfID(const String &id)
{
  vcfID=id;
}



const String &Root::getVcfID() const
{
  return vcfID;
}



bool Root::isRoot() const
{
  return true;
}



/****************************************************************
                         utility functions
 ****************************************************************/
Sex stringToSex(const String &s)
{
  String sex=s; sex.toupper();
  if(sex=="FEMALE" || sex=="F") return FEMALE;
  else if(sex=="MALE" || sex=="M") return MALE;
  return SEX_UNKNOWN;
}



ostream &operator<<(ostream &os,Sex s)
{
  switch(s) {
  case FEMALE:
    os<<"female"; break;
  case MALE:
    os<<"male"; break;
  case SEX_UNKNOWN:
    os<<"."; break;
  }
  return os;
}

/****************************************************************
 sim1.C
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
private:
  String ID;
  Sex sex;
  Array1D<String> parentID;
  Array1D<Individual*> parents;
  Vector<Individual*> children;
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

  // Load pedigree
  Pedigree *pedigree=Pedigree::loadFromTextFile(PED_FILE);
  Vector<Individual*> topsort;
  pedigree->topologicalSort(topsort);
  cout<<"TOPOLOGICAL SORT:"<<endl;
  for(Vector<Individual*>::iterator cur=topsort.begin(), end=topsort.end() ;
      cur!=end ; ++cur)
    cout<<(*cur)->getID()<<endl;
  
  // Parse VCF identifier list
  ID_LIST.getFields(vcfIDs,",");
  
  // Process VCF file
  ofstream phasedFile(phasedFilename), unphasedFile(unphasedFilename);
  VcfReader reader(VCF_FILE);
  reader.hashSampleIDs();
  //motherIndex=reader.getSampleIndex(MOTHER_ID);
  //fatherIndex=reader.getSampleIndex(FATHER_ID);
  for(int geneNum=0 ; geneNum<NUM_GENES ; ++geneNum) {
    cout<<"Simulating gene "<<(geneNum+1)<<endl;
    
    //Vector<VariantAndGenotypes> variants;
    for(int varNum=0 ; varNum<VARIANTS_PER_GENE ; ++varNum) {
      //variants.push_back(simNext(reader));
      
      }
  }

  return 0;
}



/****************************************************************
                         Individual methods
 ****************************************************************/
Individual::Individual(const String &id,Sex sex,const String &motherID,
		       const String &fatherID)
  : ID(id), sex(sex), parentID(2), parents(2)
{
  parentID[MOTHER]=motherID;
  parentID[FATHER]=fatherID;
  parents.setAllTo(NULL);
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

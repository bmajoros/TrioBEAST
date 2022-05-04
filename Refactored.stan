//====================================================================
// TrioBEAST.stan
// (C)2022 W.H. Majoros (bmajoros@alumni.duke.edu)
// This is OPEN SOURCE software, released under the GPL 3.0 license.
//
// This version has been refactored to construct the likelihood
// terms automatically for each requested mode of inheritance.
//
// Indexing of arrays:
//   Individuals: 1=mother, 2=father, 3=child
//   Haplotypes:  1=maternal, 2=paternal
//====================================================================

functions {

real binom_lpmf(int x,int n,int het,real p) 
{
   return het ? binomial_lpmf(x|n,p) : 0.0;
}

real computeElem(int[,,] count,int[,] het,int[] isPhased,int site,
   int MI,int FI,int CI,real MP,real FP,real CP) 
{
   int MN=count[site,1,1]+count[site,1,2]; // Mother
   int FN=count[site,2,1]+count[site,2,2]; // Father
   int CN=count[site,3,1]+count[site,3,2]; // Child
   int MC=count[site,1,MI];
   int FC=count[site,2,FI];
   int CC=count[site,3,CI];
   if(isPhased[site]) 
      return 
           binom_lpmf(MC | MN,het[site,1],MP) // Mother
         + binom_lpmf(FC | FN,het[site,2],FP) // Father
         + binom_lpmf(CC | CN,het[site,3],CP);// Child
   else {
      real phase1=log(0.5) +
           binom_lpmf(MC | MN,het[site,1],MP) // Mother
         + binom_lpmf(FC | FN,het[site,2],FP) // Father
         + binom_lpmf(CC | CN,het[site,3],CP);// Child
      real phase2=log(0.5) +
           binom_lpmf(MC | MN,het[site,1],1-MP) // Mother
         + binom_lpmf(FC | FN,het[site,2],1-FP) // Father
         + binom_lpmf(CC | CN,het[site,3],1-CP);// Child
      return log_sum_exp(phase1,phase2);
   }       
}

int countDenovos(int parent,int[,] V) {
   // Count as de novo only if parent has no affected copies and child does
   if(V[parent,1]==0 && V[parent,2]==0 && V[3,parent]==1) return 1;
   return 0;
}

int countRecomb(int parent,int[,] V) {
   // If the child's inherited allele differs from the parent's transmitted
   // allele, there has been a recombination in the parent (as long as
   // the parent has an affected allele -- otherwise it may be a de novo)
   if(V[parent,1]==0 && V[parent,2]==1 && V[3,parent]==1 ||
      V[parent,1]==1 && V[parent,2]==0 && V[3,parent]==0) return 1;
   return 0;
}

int affectedCopy(int indiv,int[,] V) {
   if(V[indiv,1]==1) return 1;
   if(V[indiv,2]==1) return 2;
   return 1;
}

real getP(int indiv,int[,] V,real p) {
   if(V[indiv,1]==V[indiv,2]) return 0.5;
   return p;
}

real modeLik(int site,int[,] mode,int[,,] count,int[,] het,real logAffected,
   real logUnaffected,real logDenovo,real logNoDenovo,real logRecomb,
   real logNoRecomb,real p,int N_SITES,int[] isPhased) 
{
   // Count some things and initialize variables
   real s=0; 
   int numDenovos=countDenovos(1,mode)+countDenovos(2,mode);
   int numRecomb=countRecomb(1,mode)+countRecomb(2,mode);
   int mCopy=affectedCopy(1,mode); int fCopy=affectedCopy(2,mode);
   int cCopy=affectedCopy(3,mode); real mP=getP(1,mode,p);
   real fP=getP(2,mode,p); real cP=getP(3,mode,p);

   // Priors for parent affected status
   int numParentAffected=mode[1,1]+mode[1,2]+mode[2,1]+mode[2,2];
   int numParentUnaffected=4-numParentAffected;
   s+=numParentAffected*logAffected+numParentUnaffected*logUnaffected;
   
   // Compute binomial terms
   s+=computeElem(count,het,isPhased,site,mCopy,fCopy,cCopy,mP,fP,cP);

   // Handle recombinations (parents) and de novos (child)
   s+=numRecomb*logRecomb + (2-numRecomb)*logNoRecomb;
   s+=numDenovos*logDenovo + (2-numDenovos)*logNoDenovo;

   return s;
}

real likelihoods(int N_MODES,int[,,] modes,int[,,] count,int[,] het,
   real logAffected,real logUnaffected,real logDenovo,real logNoDenovo,
   real logRecomb,real logNoRecomb,real p,int N_SITES,int[] isPhased) 
{
   real array[N_MODES]; for(i in 1:N_MODES) array[i]=0.0;
   for(site in 1:N_SITES) {
      for(i in 1:N_MODES)
         array[i]+=modeLik(site,modes[i],count,het,logAffected,logUnaffected,
            logDenovo,logNoDenovo,logRecomb,logNoRecomb,p,N_SITES,isPhased);
   }
   return log_sum_exp(array);
}
}

data {
   int N_SITES; // number of sites (variants)
   int N_MODES; // number of modes of inheritance
   int modes[N_MODES,3,2]; // [mode,individual,haplotype]
   int<lower=0,upper=1> het[N_SITES,3]; // [site,indvidual]
   int<lower=0> count[N_SITES,3,2]; // [site,individual,haplotype]
   int<lower=0,upper=1> isPhased[N_SITES]; // triple hets are unphased
   real<lower=0,upper=1> probAffected; // prior prob of 1 parent copy affected
}

transformed data {
   real logAffected=log(probAffected);
}

parameters 
{
   real<lower=0.000001,upper=1> theta; // amount of ASE
   real<lower=0,upper=1> probRecomb;
   real<lower=0,upper=1> probDenovo;
}

transformed parameters 
{
   real p = theta / (1 + theta);
   real logUnaffected=log(1-probAffected);
   real logDenovo=log(probDenovo);
   real logNoDenovo=log(1-probDenovo);
   real logRecomb=log(probRecomb);
   real logNoRecomb=log(1-probRecomb);
}

model 
{
   // Priors:
   log2(theta) ~ normal(0, 1);
   target += -log(theta * log(2)); // Jacobian
   probRecomb ~ beta(1,99);
   probDenovo ~ beta(1,999);

   // Likelihoods:
   target+=likelihoods(N_MODES,modes,count,het,logAffected,logUnaffected,
      logDenovo,logNoDenovo,logRecomb,logNoRecomb,p,N_SITES,isPhased);
}

generated quantities 
{
   real numerator[N_MODES]; real denominator;
   for(i in 1:N_MODES) numerator[i]=0.0;
   for(site in 1:N_SITES)
      for(i in 1:N_MODES)
         numerator[i]+=modeLik(site,modes[i],count,het,logAffected,
		logUnaffected,logDenovo,logNoDenovo,logRecomb,logNoRecomb,
		p,N_SITES,isPhased);
   denominator=log_sum_exp(numerator);
}

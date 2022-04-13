//====================================================================
// TrioBEAST.stan
// (C)2022 W.H. Majoros (bmajoros@alumni.duke.edu)
// This is OPEN SOURCE software, released under the GPL 3.0 license.
//
// This version of the model uses latent variables for the prior
// probabilities, rather than expecting them as fixed inputs.
//
// Indexing of arrays:
//   Individuals: 0=mother, 1=father, 3=child
//   Haplotypes:  0=maternal, 1=paternal
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


real likelihoods(int[,,] count,int[,] het,real logAffected,
   real logUnaffected,real logDenovo,real logNoDenovo,real logRecomb,
   real logNoRecomb,real p,int N_SITES,int[] isPhased) 
{
   // Sum over all possible assignments of affected status:
   real array[27];
   for(i in 1:27) array[i]=0.0;

   for(i in 1:N_SITES) {
      // MM FF CC (Mother Father Child)
      // 00 00 00 = all unaffected
      array[1]+=computeElem(count,het,isPhased,i, 1,1,1, 0.5,0.5,0.5)
         + 4*logUnaffected + 2*logNoDenovo;

      // 00 00 10 = child has a de novo in the causal variant
      array[2]+=computeElem(count,het,isPhased,i, 1,1,1, 0.5,0.5,p)
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 00 00 01 = child has a de novo in the causal variant
      array[3]+=computeElem(count,het,isPhased,i, 1,1,2, 0.5,0.5,p)
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 01 00 00 = mother affected, child doesn't inherit
      array[4]+=computeElem(count,het,isPhased,i, 2,1,1, p,0.5,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 01 00 10 = mother affected and recombines, child inherits
      array[5]+=computeElem(count,het,isPhased,i, 2,1,1, p,0.5,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 01 00 = father affected, child doesn't inherit
      array[6]+=computeElem(count,het,isPhased,i, 1,2,1, 0.5,p,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 01 00 = father affected and recombines, child inherits
      array[7]+=computeElem(count,het,isPhased,i, 1,2,2, 0.5,p,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 10 00 10 = mother affected, child inherits
      array[8]+=computeElem(count,het,isPhased,i, 1,1,1, p,0.5,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 10 00 00 = mother affected and recombines, child doesn't inherit
      array[9]+=computeElem(count,het,isPhased,i, 1,1,1, p,0.5,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 10 01 = father affected, child inherits
      array[10]+=computeElem(count,het,isPhased,i, 1,1,2, 0.5,p,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 10 00 = father affected and recombines, child doesn't inherit
      array[11]+=computeElem(count,het,isPhased,i, 1,1,1, 0.5,p,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      //--------------------- Multiple affected copies --------------------

      // 01 01 00 = both parents affected, child doesn't inherit
      array[12]+=computeElem(count,het,isPhased,i, 2,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 01 01 10 = both parents affected, one recombines, child inherits
      array[13]+=computeElem(count,het,isPhased,i, 2,2,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 01 01 = both parents affected, one recombines, child inherits
      array[14]+=computeElem(count,het,isPhased,i, 2,2,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 01 11 = both parents affected, both recombine, child inherits 2
      array[15]+=computeElem(count,het,isPhased,i, 2,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logRecomb;

      // 10 10 11 = both parents affected, child inherits 2 bad copies
      array[16]+=computeElem(count,het,isPhased,i, 1,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 10 10 10 = both parents affected, one recombines, child inherits 1
      array[17]+=computeElem(count,het,isPhased,i, 1,1,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 10 01 = both parents affected, one recombines, child inherits 1
      array[18]+=computeElem(count,het,isPhased,i, 1,1,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 10 00 = both parents affected, both recombine, child inherits 0
      array[19]+=computeElem(count,het,isPhased,i, 1,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logRecomb;

      // 10 01 10 = both parents affected, child inherits 1 copy
      array[20]+=computeElem(count,het,isPhased,i, 1,2,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 10 01 11 = both parents affected, 1 recombines, child inherits 2
      array[21]+=computeElem(count,het,isPhased,i, 1,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 01 00 = both parents affected, 1 recombines, child inherits 0
      array[22]+=computeElem(count,het,isPhased,i, 1,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 01 01 = both parents affected, 2 recombine, child inherits 1
      array[23]+=computeElem(count,het,isPhased,i, 1,2,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logRecomb;

      // 01 10 01 = both parents affected, child inherits 1 copy
      array[24]+=computeElem(count,het,isPhased,i, 2,1,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 01 10 11 = both parents affected, 1 recombines child inherits 2
      array[25]+=computeElem(count,het,isPhased,i, 2,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 10 00 = both parents affected, 1 recombines child inherits 0
      array[26]+=computeElem(count,het,isPhased,i, 2,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 10 10 = both parents affected, child inherits 1 copy
      array[27]+=computeElem(count,het,isPhased,i, 2,1,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;
   }

   return log_sum_exp(array);
   }
}
data {
   int N_SITES;
   int<lower=0,upper=1> het[N_SITES,3]; // [site,indvidual]
   int<lower=0> count[N_SITES,3,2]; // [site,individual,haplotype]
   int<lower=0,upper=1> isPhased[N_SITES]; // triple hets are unphased
   real<lower=0,upper=1> probAffected; // prior prob of 1 parent copy affected
   //real<lower=0,upper=1> probDenovo; // de novo mutation rate, per copy
   //real<lower=0,upper=1> probRecomb; // recombination rate
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
   target+=likelihoods(count,het,logAffected,logUnaffected,logDenovo,
      logNoDenovo,logRecomb,logNoRecomb,p,N_SITES,isPhased);
}
generated quantities 
{
   real numerator[27];
   real denominator;

   for(i in 1:11) numerator[i]=0.0;

   for(i in 1:N_SITES) {
      //if(!isPhased[i]) continue; // ### Need to relax this (later)

      // MM FF CC (Mother Father Child)
      // 00 00 00 = all unaffected
      numerator[1]+=computeElem(count,het,isPhased,i, 1,1,1, 0.5,0.5,0.5)
         + 4*logUnaffected + 2*logNoDenovo;

      // 00 00 10 = child has a de novo in the causal variant
      numerator[2]+=computeElem(count,het,isPhased,i, 1,1,1, 0.5,0.5,p)
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 00 00 01 = child has a de novo in the causal variant
      numerator[3]+=computeElem(count,het,isPhased,i, 1,1,2, 0.5,0.5,p)
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 01 00 00 = mother affected, child doesn't inherit
      numerator[4]+=computeElem(count,het,isPhased,i, 2,1,1, p,0.5,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 01 00 10 = mother affected and recombines, child inherits
      numerator[5]+=computeElem(count,het,isPhased,i, 2,1,1, p,0.5,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 01 00 = father affected, child doesn't inherit
      numerator[6]+=computeElem(count,het,isPhased,i, 1,2,1, 0.5,p,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 01 00 = father affected and recombines, child inherits
      numerator[7]+=computeElem(count,het,isPhased,i, 1,2,2, 0.5,p,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 10 00 10 = mother affected, child inherits
      numerator[8]+=computeElem(count,het,isPhased,i, 1,1,1, p,0.5,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 10 00 00 = mother affected and recombines, child doesn't inherit
      numerator[9]+=computeElem(count,het,isPhased,i, 1,1,1, p,0.5,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 10 01 = father affected, child inherits
      numerator[10]+=computeElem(count,het,isPhased,i, 1,1,2, 0.5,p,p)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 10 00 = father affected and recombines, child doesn't inherit
      numerator[11]+=computeElem(count,het,isPhased,i, 1,1,1, 0.5,p,0.5)
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;
      
      //--------------------- Multiple affected copies --------------------

      // 01 01 00 = both parents affected, child doesn't inherit
      numerator[12]+=computeElem(count,het,isPhased,i, 2,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 01 01 10 = both parents affected, one recombines, child inherits
      numerator[13]+=computeElem(count,het,isPhased,i, 2,2,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 01 01 = both parents affected, one recombines, child inherits
      numerator[14]+=computeElem(count,het,isPhased,i, 2,2,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 01 11 = both parents affected, both recombine, child inherits 2
      numerator[15]+=computeElem(count,het,isPhased,i, 2,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logRecomb;

      // 10 10 11 = both parents affected, child inherits 2 bad copies
      numerator[16]+=computeElem(count,het,isPhased,i, 1,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 10 10 10 = both parents affected, one recombines, child inherits 1
      numerator[17]+=computeElem(count,het,isPhased,i, 1,1,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 10 01 = both parents affected, one recombines, child inherits 1
      numerator[18]+=computeElem(count,het,isPhased,i, 1,1,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 10 00 = both parents affected, both recombine, child inherits 0
      numerator[19]+=computeElem(count,het,isPhased,i, 1,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logRecomb;

      // 10 01 10 = both parents affected, child inherits 1 copy
      numerator[20]+=computeElem(count,het,isPhased,i, 1,2,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 10 01 11 = both parents affected, 1 recombines, child inherits 2
      numerator[21]+=computeElem(count,het,isPhased,i, 1,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 01 00 = both parents affected, 1 recombines, child inherits 0
      numerator[22]+=computeElem(count,het,isPhased,i, 1,2,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 10 01 01 = both parents affected, 2 recombine, child inherits 1
      numerator[23]+=computeElem(count,het,isPhased,i, 1,2,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logRecomb;

      // 01 10 01 = both parents affected, child inherits 1 copy
      numerator[24]+=computeElem(count,het,isPhased,i, 2,1,2, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;

      // 01 10 11 = both parents affected, 1 recombines child inherits 2
      numerator[25]+=computeElem(count,het,isPhased,i, 2,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 10 00 = both parents affected, 1 recombines child inherits 0
      numerator[26]+=computeElem(count,het,isPhased,i, 2,1,1, p,p,0.5)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + logNoRecomb
      + logRecomb;

      // 01 10 10 = both parents affected, child inherits 1 copy
      numerator[27]+=computeElem(count,het,isPhased,i, 2,1,1, p,p,p)
      + 2*logUnaffected + 2*logAffected + 2*logNoDenovo + 2*logNoRecomb;
   }
   denominator=log_sum_exp(numerator);
}

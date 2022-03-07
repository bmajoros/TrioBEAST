//====================================================================
// TrioBEAST.stan
// (C)2022 W.H. Majoros (bmajoros@alumni.duke.edu)
// This is OPEN SOURCE software, released under the GPL 3.0 license.
//
// Indexing of arrays:
//   Individuals: 0=mother, 1=father, 3=child
//   Haplotypes:  0=maternal, 1=paternal
//====================================================================

functions {

real binom_lpmf(int x,int n,int het,real p) {
   return het ? binomial_lpmf(x|n,p) : 0.0;
}

real likelihoods(int[,,] count,int[,] het,real logAffected,
   real logUnaffected,real logDenovo,real logNoDenovo,real logRecomb,
   real logNoRecomb,real p,int N_SITES,int[] isPhased) {

   // Sum over all possible assignments of affected status:
   real array[11];
      
   for(i in 1:11) array[i]=0.0;

   for(i in 1:N_SITES) {
      if(!isPhased[i]) continue; // ### Need to relax this (later)

      // MM FF CC (Mother Father Child)
      // 00 00 00 = all unaffected
      array[1]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5)//M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5)//F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5)//C
      + 4*logUnaffected + 2*logNoDenovo;

      // 00 00 10 = child has a de novo in the causal variant
      array[2]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5)//M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5)//F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],p)//C
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 00 00 01 = child has a de novo in the causal variant
      array[3]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,2]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 01 00 00 = mother affected, child doesn't inherit
      array[4]+=
        binom_lpmf(count[i,1,2]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 01 00 10 = mother affected and recombines, child inherits
      array[5]+=
        binom_lpmf(count[i,1,2]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 01 00 = father affected, child doesn't inherit
      array[6]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,2]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 01 00 = father affected and recombines, child inherits
      array[7]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,2]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,2]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 10 00 10 = mother affected, child inherits
      array[8]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 10 00 00 = mother affected and recombines, child doesn't inherit
      array[9]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 10 01 = father affected, child inherits
      array[10]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,2]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 10 00 = father affected and recombines, child doesn't inherit
      array[11]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;
   }

   //for(i in 1:11) print(array[i]);
 
   return log_sum_exp(array);
   }
}
data {
   int N_SITES;
   int<lower=0,upper=1> het[N_SITES,3]; // [site,indvidual]
   int<lower=0> count[N_SITES,3,2]; // [site,individual,haplotype]
   int<lower=0,upper=1> isPhased[N_SITES]; // triple hets are unphased
   real<lower=0,upper=1> probDenovo; // de novo mutation rate, per copy
   real<lower=0,upper=1> probRecomb; // recombination rate
   real<lower=0,upper=1> probAffected; // prior prob of 1 parent copy affected
}
transformed data {
   real logAffected=log(probAffected);
   real logUnaffected=log(1-probAffected);
   real logDenovo=log(probDenovo);
   real logNoDenovo=log(1-probDenovo);
   real logRecomb=log(probRecomb);
   real logNoRecomb=log(1-probRecomb);
}
parameters {
   real<lower=0> theta; // amount of ASE
}
transformed parameters {
   real p = theta / (1 + theta);
}
model {
   // Priors:
   log2(theta) ~ normal(0, 1);
   target += -log(theta * log(2)); // Jacobian

   // Likelihoods:
   target+=likelihoods(count,het,logAffected,logUnaffected,logDenovo,
      logNoDenovo,logRecomb,logNoRecomb,p,N_SITES,isPhased);
}
generated quantities {
   real numerator[11];
   real denominator;

   for(i in 1:11) numerator[i]=0.0;

   for(i in 1:N_SITES) {
      if(!isPhased[i]) continue; // ### Need to relax this (later)

      // MM FF CC (Mother Father Child)
      // 00 00 00 = all unaffected
      numerator[1]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5)//M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5)//F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5)//C
      + 4*logUnaffected + 2*logNoDenovo;

      // 00 00 10 = child has a de novo in the causal variant
      numerator[2]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5)//M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5)//F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],p)//C
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 00 00 01 = child has a de novo in the causal variant
      numerator[3]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,2]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 4*logUnaffected + logDenovo + logNoDenovo;

      // 01 00 00 = mother affected, child doesn't inherit
      numerator[4]+=
        binom_lpmf(count[i,1,2]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 01 00 10 = mother affected and recombines, child inherits
      numerator[5]+=
        binom_lpmf(count[i,1,2]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 01 00 = father affected, child doesn't inherit
      numerator[6]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,2]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 01 00 = father affected and recombines, child inherits
      numerator[7]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,2]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,2]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 10 00 10 = mother affected, child inherits
      numerator[8]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 10 00 00 = mother affected and recombines, child doesn't inherit
      numerator[9]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],p) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],0.5) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 10 01 = father affected, child inherits
      numerator[10]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,2]|count[i,3,1]+count[i,3,2],het[i,3],p) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 10 00 = father affected and recombines, child doesn't inherit
      numerator[11]+=
        binom_lpmf(count[i,1,1]|count[i,1,1]+count[i,1,2],het[i,1],0.5) //M
      + binom_lpmf(count[i,2,1]|count[i,2,1]+count[i,2,2],het[i,2],p) //F
      + binom_lpmf(count[i,3,1]|count[i,3,1]+count[i,3,2],het[i,3],0.5) //C
      + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;
   }
   denominator=log_sum_exp(numerator);
}

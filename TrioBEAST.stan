data {
   int N_SITES;
   int<lower=0> count[N_SITES,3,2]; // site,individual,haplotype
   int<lower=0,upper=1> genotype[N_SITES,3,2]; // site,individual,haplotype
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
   for(i in 1:N_SITES) {
      if(!isPhased[i]) continue; // ### Need to relax this (later)
      
      // Sum over all possible assignments of affected status:
      real array[11];
      
      // 00 00 00 = all unaffected
      array[1]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],0.5) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],0.5) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],0.5) //C
         + 4*logUnaffected + 2*logNoDenovo;

      // 00 00 10 = child has a de novo in the causal variant
      array[2]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],0.5) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],0.5) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],p) //C
         + 4*logUnaffected + logDenovo + logNoDenovo;

      // 00 00 01 = child has a de novo in the causal variant
      array[3]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],0.5) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],0.5) //F
         + binomial_lpmf(count[i][3][2]|count[i][3][1]+count[i][3][2],p) //C
         + 4*logUnaffected + logDenovo + logNoDenovo;

      // 01 00 00 = mother affected, child doesn't inherit
      array[4]=
           binomial_lpmf(count[i][1][2]|count[i][1][1]+count[i][1][2],p) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],0.5) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],0.5) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 01 00 10 = mother affected and recombines, child inherits
      array[5]=
           binomial_lpmf(count[i][1][2]|count[i][1][1]+count[i][1][2],p) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],0.5) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],p) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 01 00 = father affected, child doesn't inherit
      array[6]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],0.5) //M
         + binomial_lpmf(count[i][2][2]|count[i][2][1]+count[i][2][2],p) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],0.5) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 01 00 = father affected and recombines, child inherits
      array[7]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],0.5) //M
         + binomial_lpmf(count[i][2][2]|count[i][2][1]+count[i][2][2],p) //F
         + binomial_lpmf(count[i][3][2]|count[i][3][1]+count[i][3][2],p) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 10 00 10 = mother affected, child inherits
      array[8]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],p) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],0.5) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],p) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 10 00 00 = mother affected and recombines, child doesn't inherit
      array[9]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],p) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],0.5) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],0.5) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;

      // 00 10 01 = father affected, child inherits
      array[10]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],0.5) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],p) //F
         + binomial_lpmf(count[i][3][2]|count[i][3][1]+count[i][3][2],p) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logNoRecomb;

      // 00 10 00 = father affected and recombines, child doesn't inherit
      array[11]=
           binomial_lpmf(count[i][1][1]|count[i][1][1]+count[i][1][2],0.5) //M
         + binomial_lpmf(count[i][2][1]|count[i][2][1]+count[i][2][2],p) //F
         + binomial_lpmf(count[i][3][1]|count[i][3][1]+count[i][3][2],0.5) //C
         + 3*logUnaffected + logAffected + 2*logNoDenovo + logRecomb;
 
      target+=log_sum_exp(array);
   }
}


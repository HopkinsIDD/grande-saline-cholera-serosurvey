data {
  int<lower = 0> N;
  int y[N]; // an array with N elements
  int<lower = 0> K;
  vector[K] intervalBreaks; // a vector with K elements
  int<lower = 0> nPriorPar;
  vector[nPriorPar] prior;
}

parameters {
  ordered[2] mu; // this must be ordered (an ordered vector) so the chains do not constantly switch between the two distributions
  real<lower=0> sigma[2];
  real<lower=0, upper=1> theta;
}

transformed parameters {
  
  //do the interval censoring here (assign probabilities to tau)

  vector[K+1] tau; //the probabilities associated to the categories 1:K+1 with breaks intervalBreaks associated to them. 

  //using Phi (with scaled and translated x) here instead of normal_cdf because normal_cdf returns 1 for (y-mu)/sigma > 8.25 (which is a very low limit). Also has some advantages to use the actual probabilities instead of logs here (because log_diff_exponential returns nan if a==b). See Stan Manual p523 (version 2.17.0-1).
  //if this does not work I could also use Phi_approx, which has even a larger support.


  for (k in 1:K+1){
    if (k == 1)
          tau[k] =      theta * Phi_approx((intervalBreaks[1]-mu[1])/sigma[1])
                      + (1-theta) * Phi_approx((intervalBreaks[1]-mu[2])/sigma[2]);
    else if (k > 1 && k <= K)
          tau[k] =      theta * Phi_approx((intervalBreaks[k]-mu[1])/sigma[1])
                      + (1-theta) * Phi_approx((intervalBreaks[k]-mu[2])/sigma[2])
                      - theta * Phi_approx((intervalBreaks[k-1]-mu[1])/sigma[1])
                      - (1-theta) * Phi_approx((intervalBreaks[k-1]-mu[2])/sigma[2]);
    else if (k == (K+1))
          tau[k] =      1
                      - theta * Phi_approx((intervalBreaks[K]-mu[1])/sigma[1])
                      - (1-theta) * Phi_approx((intervalBreaks[K]-mu[2])/sigma[2]);
  }
  
  tau = tau/sum(tau); //make sure tau sums to 1

  //print(tau)

}

model {
  mu[1] ~ normal(prior[1], prior[2]);
  mu[2] ~ normal(prior[3], prior[4]);
  sigma ~ normal(prior[5], prior[6]);
  theta ~ beta(prior[7], prior[8]);

  y ~ categorical(tau);
  
}

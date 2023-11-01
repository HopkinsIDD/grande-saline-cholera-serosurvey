
data {
  
  // actual data 
  int<lower=1> L; // number of individuals
  int<lower=1> N; // number of observations
  row_vector[N]  t; //time of observation
  real y[N]; // value of observation
  int cens[N]; // value of censorship
  int<lower=1,upper=L> id_new[N]; // individual id
  
  real delay ; //delay in response

 
  // Measurement model type
  int <lower=1, upper =3> measure ;// 1= Luminex, 2= Vibriocidal, 3= ELISA
 
}
parameters {
  
  // population average baseline and boost
  real mu_omega;
  real<lower=0> mu_lambda;
  real<lower=log(delay),upper=log(365)> mu_logD ; //delay in response
  
  // individual baseline and boost
  vector[L] omega_ind; // individual specific baseline and boost
  vector<lower=0>[L] lambda_ind; // individual specific baseline and boost
  vector<lower=log(delay),upper=log(365)>[L] logD_ind; // individual specific baseline and boost

  // log halflife for everyone
  real<lower=0> log_halflife;
  real<lower=0> SIGMA; //random error
  
  cov_matrix[3] params_sigma; //covariance matrix for baseline and boost


}

transformed parameters {
    
    real halflife = exp(log_halflife);
    real delta = log(2)/halflife;
    
    real MU[N];
    vector[L] D_ind = exp(logD_ind); // individual specific baseline and boost
    
    real foldchange;
    real log_foldchange;


// calculate expected value for each observation
  for (n in 1:N){
    MU[n] = omega_ind[id_new[n]]; // baseline level
      if ((t[n] >=delay) && (t[n]<D_ind[id_new[n]])){ // linear increase increase
         MU[n] =MU[n] + (t[n]-delay)*lambda_ind[id_new[n]]/(D_ind[id_new[n]]-delay) ;
      }
      if (t[n]>=D_ind[id_new[n]]){ // boost if after D days
         MU[n] =MU[n] + lambda_ind[id_new[n]] * exp(-delta * (t[n]-D_ind[id_new[n]]));
      }
  }
  
  // calculate fold change based on measurement
  if (measure ==1) foldchange= 10^(mu_lambda);
  if (measure ==2) foldchange= 2^(mu_lambda);
  if (measure ==3) foldchange= 10^(mu_lambda);
  
  log_foldchange = log(foldchange);
  
}


model {
        
  target += normal_lpdf(mu_logD | log(14), log(2)) - 1 * normal_lccdf(0 | log(14),log(2));
  target += normal_lpdf(log_halflife | log(30), log(5)) ;//- 1 * normal_lccdf(0 | log(30),log(5));
  target += normal_lpdf(log_foldchange | log(8), 1) - 1 * normal_lccdf(0 | log(8),1);
        
    for (l in 1:L){
        target += multi_normal_lpdf([omega_ind[l],lambda_ind[l],logD_ind[l]]|[mu_omega,mu_lambda,mu_logD], params_sigma);

}
    
  
  // Measurement model to deal with censorship for RAU
  if (measure ==1){
    for (n in 1:N){
      if(cens[n]==-1) target += normal_lcdf(y[n]| MU[n], SIGMA);
      if(cens[n] == 0) target += normal_lpdf(y[n]| MU[n], SIGMA);
      if(cens[n]==1) target += normal_lccdf(y[n]| MU[n], SIGMA);
    }
  }

  // Measurement model to deal with censorship for Vibriocidal
  if (measure ==2){
  for (n in 1:N){
      if(cens[n]==0) target += normal_lcdf(y[n]| MU[n], SIGMA);
      if ((cens[n] > 0) && (cens[n] < 11)){
        target += log(normal_cdf((y[n]+1), MU[n], SIGMA) -normal_cdf(y[n], MU[n], SIGMA));
      }
      if(cens[n]==11) target += normal_lccdf(y[n]| MU[n], SIGMA);
    }
  }
  
  // Measurment model to deal with censorship for ELISA
  if (measure ==3) target += normal_lpdf(y| MU,SIGMA);
}

generated quantities{

  vector[N] log_lik;
  
  // calculate log likelihood for each measurement for loo
  for (n in 1:N) {
      log_lik[n] = normal_lpdf(y[n] | MU[n], SIGMA);
  }
  
}


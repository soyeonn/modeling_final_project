data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  real r[N, T];
  int a[N, T];
  int s[N, T];
  int Avoid[N, T];
}

transformed data {
  vector[7] initV;
  //real<lower=0, upper=1> alpha0;
  initV = rep_vector(0.0, 7);
  //alpha0 = normal_rng(0, 1);

}

parameters {
  vector[4] mu_p;
  vector<lower=0>[4] sigma;
  
  vector[N] kappa_pr;
  vector[N] eta_pr;
  vector[N] beta_pr;
  vector[N] bias_pr;
  
}

transformed parameters {
  vector<lower=0>[N] kappa;
  vector<lower=0>[N] eta;
  vector<lower=0>[N] beta;
  vector<lower=0>[N] bias;
  
  for (i in 1:N) {
    kappa[i] = Phi_approx(mu_p[1] + sigma[1] * kappa_pr[i]);
    eta[i] = Phi_approx(mu_p[2] + sigma[2] * eta_pr[i]);
    beta[i]    = Phi_approx(mu_p[3] + sigma[3] * beta_pr[i]) * 10;
    bias[i]    = Phi_approx(mu_p[4] + sigma[4] * bias_pr[i]);
  }
}

model {
   mu_p  ~ normal(0, 1);
  sigma  ~ normal(0, 0.2);

  // individual parameters
  kappa_pr    ~ normal(0, 1);
  eta_pr      ~ normal(0, 1);
  beta_pr     ~ normal(0, 1);
  bias_pr     ~ normal(0, 1);
  
 for (i in 1:N) { 
   vector[7] ev;
     real delta;
     real PE;
     real alpha;
     
     ev = initV;
     alpha = 0.5;
     delta = 0;
     
   for (t in 1:Tsubj[i]) {
     delta = - beta[i] * (ev[s[i, t]] + 0.2 + bias[i]);
     Avoid[i, t] ~ bernoulli_logit(delta);
     
     if (a[i, t] == 0) {
       PE = (r[i, t] - ev[s[i, t]]);
       
       ev[s[i, t]] = ev[s[i, t]] + kappa[i] * alpha * PE;
     } else {
       PE = 0;
     }
     alpha = eta[i] * (fabs(PE)) + (1-eta[i]) * alpha;
    }
  }
}

generated quantities {
  real<lower=0, upper=1> mu_kappa;
  real<lower=0, upper=1> mu_eta;
  real<lower=0> mu_beta;
  real<lower=0, upper=1> mu_bias;
  
  real log_lik[N];


    mu_kappa = Phi_approx( mu_p[1] );
    mu_eta = Phi_approx( mu_p[2] );
    mu_beta = Phi_approx( mu_p[3] ) * 10;
    mu_bias = Phi_approx( mu_p[4] );

{for (i in 1:N) { 
  vector[7] ev;
     real delta;
     real PE;
     real alpha;
     
     ev = initV;
     alpha = 0.5;
     delta = 0;
     
     log_lik[i] = 0.0;
     
   for (t in 1:Tsubj[i]) {
     delta = - beta[i] * (ev[s[i, t]] + 0.2 + bias[i]);
     log_lik[i] = log_lik[i] + bernoulli_logit_lpmf(Avoid[i, t] | delta);
     
     if (a[i, t] == 0) {
       PE = (r[i, t] - ev[s[i, t]]);
       
       ev[s[i, t]] = ev[s[i, t]] + kappa[i] * alpha * PE;
     } else {
       PE = 0;
     }
     alpha = eta[i] * (fabs(PE)) + (1-eta[i]) * alpha;
     }
   }
  }
} 


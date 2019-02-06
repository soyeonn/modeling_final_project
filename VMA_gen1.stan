data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  real<lower=-1, upper=0> r[N, T];
  int<lower=1, upper=7> s[N, T];
  int a[N, T];
  int Avoid[N, T];
}

transformed data {
  vector[7] initV;
  vector<lower=0, upper=1>[7] rhos;
  //real<lower=0, upper=1> alpha0;
  initV = rep_vector(0.0, 7);
  rhos = [0.1850, 0.2500, 0.3150, 0.6200, 0.6850, 0.7500, 0.8150]';
  //alpha0 = normal_rng(0, 1);

}

parameters {
  vector[5] mu_p;
  vector<lower=0>[5] sigma;
  
  vector[N] kappa_pr;
  vector[N] eta_pr;
  vector[N] theta_pr;
  vector[N] beta_pr;
  vector[N] bias_pr;
  
}

transformed parameters {
  vector<lower=0>[N] kappa;
  vector<lower=0>[N] eta;
  vector<lower=0>[N] theta;
  vector<lower=0>[N] beta;
  vector<lower=0>[N] bias;
  
  for (i in 1:N) {
    kappa[i]   = Phi_approx(mu_p[1] + sigma[1] * kappa_pr[i]);
    eta[i]     = Phi_approx(mu_p[2] + sigma[2] * eta_pr[i]);
    theta[i]   = Phi_approx(mu_p[3] + sigma[3] * theta_pr[i]);
    beta[i]    = Phi_approx(mu_p[4] + sigma[4] * beta_pr[i]) * 10;
    bias[i]    = Phi_approx(mu_p[5] + sigma[5] * bias_pr[i]);
  }
}

model {
   mu_p  ~ normal(0, 1);
  sigma  ~ normal(0, 0.2);

  // individual parameters
  kappa_pr    ~ normal(0, 1);
  eta_pr      ~ normal(0, 1);
  theta_pr    ~ normal(0, 1);
  beta_pr     ~ normal(0, 1);
  bias_pr     ~ normal(0, 1);
  
 for (i in 1:N) { 
   vector[7] ev;
   real delta;
   real PE;
   real alpha;
   int currShape;
   real currShapeRho;
   real otherShapeRho;
   real G;
     
   ev = initV;
   alpha = 0.5;
   delta = 0;
     
   for (t in 1:Tsubj[i]) {
     delta = -(ev[s[i, t]] + 0.2 + bias[i]);
     Avoid[i, t] ~ bernoulli_logit(beta[i] * delta);
 
     if (a[i, t] == 0) {
       currShape = s[i, t];
       currShapeRho = rhos[currShape];
       PE = r[i, t] - ev[s[i, t]];
     
       for (j in 1:7) {
         otherShapeRho = rhos[j];
         G = 1/exp((otherShapeRho - currShapeRho)^2 / (2*theta[i]^2));
         ev[j] = ev[j] + kappa[i] * alpha * PE * G;
       }
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
  real<lower=0, upper=1> mu_theta;
  real<lower=0> mu_beta;
  real<lower=0, upper=1> mu_bias;
  
  real log_lik[N];

    mu_kappa = Phi_approx( mu_p[1] );
    mu_eta = Phi_approx( mu_p[2] );
    mu_theta = Phi_approx( mu_p[3] );
    mu_beta = Phi_approx( mu_p[4] ) * 10;
    mu_bias = Phi_approx( mu_p[5] );
    
     {
     for (i in 1:N) { 
   vector[7] ev;
   real delta;
   real PE;
   real alpha;
   int currShape;
   real currShapeRho;
   real otherShapeRho;
   real G;
     
   ev = initV;
   alpha = 0.5;
   delta = 0;
   
   log_lik[i] = 0.0;
     
   for (t in 1:Tsubj[i]) {
     delta = -(ev[s[i, t]] + 0.2 + bias[i]);
     log_lik[i] = log_lik[i] + bernoulli_logit_lpmf(Avoid[i, t] | beta[i] * delta);
 
     if (a[i, t] == 0) {
       currShape = s[i, t];
       currShapeRho = rhos[currShape];
       PE = r[i, t] - ev[s[i, t]];
     
       for (j in 1:7) {
         otherShapeRho = rhos[j];
         G = 1/exp((otherShapeRho - currShapeRho)^2 / (2*theta[i]^2));
         ev[j] = ev[j] + kappa[i] * alpha * PE * G;
       }
     } else {
       PE = 0;
     }
     
     alpha = eta[i] * (fabs(PE)) + (1-eta[i]) * alpha;
      }
    }
  }
} 

// input data
data {
  int<lower=1> N; // number of observations
  int<lower=1> S; //number of species
  int<lower=1> max; //max number of draws
  int <lower=0,upper=max> presence [N]; //obs of presence absence
  int <lower=0,upper=max> draw [N]; 
  int <lower=1,upper=S> species [N];
  // real prior_t_fsm;
  // real prior_t_hsm; // vector[S] prior_s;
  vector[N] fsm; // predictor : frost safety margins
  vector[N] hsm; // predictor : hydraulic safety margins
}

// the parameters to be estimated from the data
parameters {
  real <lower=0,upper=1> K_int;
  vector <lower=0,upper=1> [S] K_sp; // plateau of the second segment, or threshold value
  real <lower=0> s_k;
  real <lower=0.1> r_fsm;
  real<lower=0.1> r_hsm;
  real <lower=min(fsm),upper=max(fsm)> t_fsm; //point where regression changes
  real <lower=min(hsm),upper=max(hsm)> t_hsm; //point where regression changes
}

transformed parameters {
  vector <lower=0,upper=1> [N] proba;
  proba = K_sp[species] ./
                (
                  (1 + exp(-r_fsm * (fsm - t_fsm))).*
                  (1 + exp(-r_hsm * (hsm - t_hsm)))
                  );
}
// The model to be estimated.
// presence are distributed according to Bernoulli law
// logit used as link-function
model {
  //priors
  // K_int~beta(2.5,3.5);
  // K_sp~beta(K_int,s_k); 
  // s_k ~ gamma(1,0.5);
  // K_int~gamma(3,0.5);
  // K_sp~beta(K_int,s_k);
  // s_k ~ gamma(10,0.6);
  K_int~normal(0.1,0.05);
  K_sp~normal(K_int,s_k); 
  s_k ~ gamma(0.3,0.5);
  t_fsm~normal(0,0.3);
  t_hsm~normal(0,0.3);

  // How the data are expected to have been generated:
  presence ~ binomial(draw,proba);
}

// // predictions
// generated quantities {
//   int presence_pred [N];
//   for (i in 1:N){
//     presence_pred[i] = bernoulli_rng(K_sp[species[i]] /(
//                           (1 + exp(-r_fsm * (fsm[i] - t_fsm)))*
//                           (1+exp(-r_hsm * (hsm[i] - t_hsm)))
//                           ));
//   }
// }

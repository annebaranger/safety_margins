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
  real <lower=1.1,upper=9> b_2;
  real <lower=0,upper=b_2/2> b_1;
  vector <lower=0,upper=1> [S] K_sp; // plateau of the second segment, or threshold value
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
  b_1 ~ gamma(1.5,0.5);
  b_2 ~ gamma(6,0.5);

  
  K_sp~beta(b_1,b_2);
  
  t_fsm~normal(0,0.3);
  t_hsm~normal(0,0.3);

  // How the data are expected to have been generated:
  presence ~ binomial(draw,proba);
}

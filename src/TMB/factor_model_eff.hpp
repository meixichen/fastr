#ifndef factor_model_eff_hpp
#define factor_model_eff_hpp

#include "fastr/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
  
template<class Type>
Type factor_model_eff(objective_function<Type>* obj){
  // data inputs
  DATA_INTEGER(n_factor); // number of factors 
  DATA_SCALAR(dt); // length of the time bin 
  DATA_ARRAY(Y); //  q x n x r array of 0s and 1s.
  DATA_SCALAR(lam); // regularization parameter
  DATA_SCALAR(nu); // sigmoid function scale parameter
  DATA_INTEGER(held_out_cell); // index of the neuron to be held out. If zero, all are included.

  // parameters
  PARAMETER_VECTOR(log_k); // log thresholds (k>0)
  PARAMETER_VECTOR(log_a); // log drift (0<drift<1)
  PARAMETER_VECTOR(Lt); // dq-d(d-1)/2 vector of the lower triangular elements of L (by column)
  PARAMETER_ARRAY(x); // q x n x r neuron paths
 
  using namespace density;
  using namespace fastr;
  // transformed data
  vector<int> Y_dim = Y.dim; // dimension of data Y (a 3d array)
  int n_cell = Y_dim(0); // number of neurons
  int n_bin = Y_dim(1); // number of time bins
  int n_trial = Y_dim(2); // number of trials
  Type nu_x; // nu*x
  int held_out_ind = held_out_cell-1; // convert to cpp index that starts from 0

  // transformed parameters
  vector<Type> k = exp(log_k);
  vector<Type> alpha = exp(log_a);
  // negative log-likelihood computation
  MVN_FA<Type> latent_nll(Lt, n_cell, n_factor, dt); // efficiently compute density of MVN using the Woodbury formula
  vector<Type> reg = Lt*Lt; // l2 (ridge) regularization 
  Type nll = lam*(reg.sum()); // negative log-likelihood penalization term
  vector<Type> mu(n_cell); // drift vector
  vector<Type> Nt(n_cell); // count the number of spikes up to and including time t for each neuron

  for (int u=0;u<n_trial;u++){
    mu = x.col(u).col(0) - alpha * dt;
    nll += latent_nll(mu);
    for(int w=1;w<n_bin;w++){
      mu = x.col(u).col(w) - x.col(u).col(w-1) - alpha * dt;
      nll += latent_nll(mu);
    }
    Nt.fill(1);

    auto start = std::chrono::high_resolution_clock::now();
    for(int j=0;j<n_bin;j++){
      for (int i=0;i<n_cell;i++){
        if (i != held_out_ind){
          nu_x = nu * (x(i,j,u) - Nt[i] * k[i]);
          nll -= Y(i,j,u) * nu_x - log(1 + exp(nu_x)); // log Bernoulli dens with sigmoid prob
          Nt(i) += Y(i,j,u);
        }
      }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Data likelihood evaluation took (in ms): " << std::endl;
    std::cout << duration.count() << std::endl;
  }

  matrix<Type> Sig = latent_nll.cov();
  ADREPORT(Sig);
  
  SIMULATE{
    matrix<Type> x_exceed(n_bin, n_trial);
    matrix<Type> y_pred(n_bin, n_trial);
    Type Nt_pred;
    if (held_out_ind >= Type(0)){
      for (int u=0;u<n_trial;u++){
        Nt_pred = 1;
        for (int j=0;j<n_bin;j++){
          x_exceed(j,u) = x(held_out_ind, j, u) - Nt_pred * k[held_out_ind];
          if (x_exceed(j,u) >= 1e-9){
            y_pred(j,u) = Type(1);
            Nt_pred += 1;
          } else{
            y_pred(j,u) = Type(0);
          }
        }
      }
    }
    REPORT(x_exceed);
    REPORT(y_pred);
  }
  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

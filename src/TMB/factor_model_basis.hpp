#ifndef factor_model_basis_hpp
#define factor_model_basis_hpp

#include "fastr/utils.hpp"
#include <chrono>

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
  
template<class Type>
Type factor_model_basis(objective_function<Type>* obj){
  // data inputs
  DATA_INTEGER(n_factor); // number of factors 
  DATA_SCALAR(dt); // length of the time bin 
  DATA_ARRAY(Y); //  q x n x r array of 0s and 1s.
  DATA_SCALAR(lam); // regularization parameter
  DATA_SCALAR(nu); // sigmoid function scale parameter
  DATA_ARRAY(Phi); // m x n basis functions
  DATA_SCALAR(T); // length of the time window 

  // parameters
  PARAMETER_VECTOR(log_k); // log thresholds (k>0)
  PARAMETER_VECTOR(log_a); // log drift (0<drift<1)
  PARAMETER_VECTOR(Lt); // dq-d(d-1)/2 vector of the lower triangular elements of L (by column)
  PARAMETER_ARRAY(Xi); // m x q x r orthonormal basis function coefficients 
 
  using namespace density;
  using namespace fastr;
  // transformed data
  vector<int> Y_dim = Y.dim; // dimension of data Y (a 3d array)
  vector<int> Phi_dim = Phi.dim; // dimension of basis function array (3d)
  int n_cell = Y_dim(0); // number of neurons
  int n_bin = Y_dim(1); // number of time bins
  int n_trial = Y_dim(2); // number of trials
  int n_basis = Phi_dim(0); //
  Type nu_x; // nu*x

  // transformed parameters
  vector<Type> k = exp(log_k);
  vector<Type> alpha = exp(log_a);
  // negative log-likelihood computation
  MVN_FA<Type> latent_nll(Lt, n_cell, n_factor, 1.); // efficiently compute density of MVN using the Woodbury formula
  vector<Type> reg = Lt*Lt; // l2 (ridge) regularization 
  Type nll = lam*(reg.sum()); // negative log-likelihood penalization term
  vector<Type> mu(n_cell); // drift vector
  vector<Type> Nt(n_cell); // count the number of spikes up to and including time t for each neuron
  
  for (int u=0;u<n_trial;u++){
    for(int w=0;w<n_basis;w++){
      mu = Xi.col(u).transpose().col(w);   
      nll += latent_nll(mu);
    }
    Nt.fill(1); 

    auto start = std::chrono::high_resolution_clock::now();
    for(int j=0;j<n_bin;j++){
      for (int i=0;i<n_cell;i++){
        nu_x = nu * ((Xi.col(u).col(i) * Phi.col(j)).sum() * sqrt(T) + 
			alpha[i]*dt*(j+1) - Nt[i] * k[i]);  // *****
        nll -= Y(i,j,u) * nu_x - log(1 + exp(nu_x)); // log Bernoulli dens with sigmoid prob
        Nt(i) += Y(i,j,u);
      }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Data likelihood evaluation took (in ms): " << std::endl;
    std::cout << duration.count() << std::endl;
  }
  
  matrix<Type> Sig = latent_nll.cov();
  ADREPORT(Sig);
  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

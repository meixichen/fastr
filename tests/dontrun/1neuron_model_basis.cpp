#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
  // data inputs                                                                                       
  DATA_SCALAR(dt); // length of the time bin
  DATA_ARRAY(Y); //  n x r array of 0s and 1s    
  DATA_SCALAR(nu); // sigmoid function scale parameter  
  DATA_ARRAY(Phi); // m x n basis function

  // parameters                                                                                        
  PARAMETER(log_k); // log thresholds (k>0)                                                   
  PARAMETER(log_a); // log drift (0<drift<1)                                                  
  PARAMETER_ARRAY(xi); // m x r coefficients                                                     
  
  using namespace density;         
  // transformed data                                                                                  
  vector<int> Y_dim = Y.dim; // dimension of data Y     
  vector<int> Phi_dim = Phi.dim;  // dimension of basis function                                   
  int n_bin = Y_dim(0); // number of time bins                                                   
  int n_trial = Y_dim(1); // number of trials                                                       
  int n_basis = Phi_dim(0); // number of basis functions

  // transformed parameters                                                                         
  Type k = exp(log_k);    
  Type alpha = exp(log_a);
  Type nu_x; // nu*x
 
  // negative log-likelihood computation  
  Type nll = 0;
  Type Nt; // count the number of spikes up to and including time t   
  for(int w=1;w<n_basis;w++){ 
    mu = xi(w);               
    nll -= dnorm(mu, Type(0), sqrt(dt), true);
  }
  for(int u=0;u<n_trial;u++){
    Nt = 1;   
    for(int j=0;j<n_bin;j++){
      nu_x = nu * ((Phi.col(j) * xi.col(r)).sum() * sqrt(T) + alpha*dt*(j+1) - Nt*k);            
      nll -= Y(j,u) * nu_x - log(1 + exp(nu_x)); 
      Nt += Y(j,u);
    }
  }
  return nll;
}


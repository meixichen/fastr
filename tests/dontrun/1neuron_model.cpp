#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
  // data inputs                                                                                       
  DATA_SCALAR(dt); // length of the time bin                                                           
  DATA_ARRAY(Y); //  n x r array of 0s and 1s                                                   
  
  // parameters                                                                                        
  PARAMETER(log_k); // log thresholds (k>0)                                                   
  PARAMETER(log_a); // log drift (0<drift<1)                                                  
  PARAMETER_ARRAY(x); // n x r neuron paths                                                      
  
  using namespace density;         
  // transformed data                                                                                  
  vector<int> Y_dim = Y.dim; // dimension of data Y (a 3d array)                                       
  int n_bin = Y_dim(0); // number of time bins                                                         
  int n_trial = Y_dim(1); // number of trials                                                          
  
  // transformed parameters                                                                            
  Type k = exp(log_k);    
  Type alpha = exp(log_a);

  // negative log-likelihood computation  
  Type nll = 0;
  Type mu;
  Type Nt; // count the number of spikes up to and including time t   
  for(int u=0;u<n_trial;u++){
    mu = x(0,u) - alpha * dt;       
    nll -= dnorm(mu, Type(0), sqrt(dt), true);
    for(int w=1;w<n_bin;w++){ 
      mu = x(w,u) - x(w-1,u) - alpha * dt;               
      nll -= dnorm(mu, Type(0), sqrt(dt), true);
    }
    Nt = 1;   
    for(int j=0;j<n_bin;j++){
      nll -= dbinom_robust(Y(j,u), Type(1), x(j,u) - Nt * k, true);            
      Nt += Y(j,u);
    }
  }
  return nll;
}


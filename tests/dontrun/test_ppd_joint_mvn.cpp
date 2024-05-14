#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
  // data inputs
  DATA_VECTOR(y); 
  DATA_SCALAR(log_a);

  // parameters
  PARAMETER(mu);
  PARAMETER(theta);
  PARAMETER(log_b);
  PARAMETER(log_c);
  
  Type nll = -sum(dnorm(y, mu, exp(log_c), true));
  nll -= dnorm(mu, theta, exp(log_b), true);
  nll -= dnorm(theta, Type(0), exp(log_a), true);
  
  return nll;

}


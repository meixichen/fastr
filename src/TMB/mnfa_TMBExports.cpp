// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_mnfa_TMBExports
#include <TMB.hpp>
#include "factor_model_big.hpp"
#include "factor_model.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "factor_model_big") {
    return factor_model_big(this);
  } else if(model == "factor_model") {
    return factor_model(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}

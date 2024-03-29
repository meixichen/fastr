#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
  // data inputs                                                                                     
  DATA_INTEGER(n_factor); // number of factors                                                       
  DATA_SCALAR(dt); // length of the time bin                                                         
  DATA_ARRAY(Y); //  q x n x r array of 0s and 1s.                                                   
  DATA_SCALAR(lam); // regularization parameter                                                      
  DATA_SCALAR(nu); 

  // parameters                                                                                      
  PARAMETER_VECTOR(log_k); // log thresholds (k>0)                                                   
  PARAMETER_VECTOR(log_a); // log drift (0<drift<1)                                                  
  PARAMETER_VECTOR(Lt); // dq vector of the lower triangular elements of L (by column)      
  PARAMETER_ARRAY(x); // q x n x r neuron paths                                                      
  
  using namespace density;         
  // transformed data                                                                                
  vector<int> Y_dim = Y.dim; // dimension of data Y (a 3d array)                                     
  int n_cell = Y_dim(0); // number of neurons                                                        
  int n_bin = Y_dim(1); // number of time bins                                                       
  int n_trial = Y_dim(2); // number of trials                                                        
  Type nu_x; 

  // transformed parameters                                                                          
  vector<Type> k = exp(log_k);                                                                       
  vector<Type> alpha = exp(log_a);                                                                   
  
  //----------Modularize below -------------                                                        
  matrix<Type> L(n_cell, n_factor); // loading matrix                                               
  matrix<Type> Sig(n_cell, n_cell); // latent path increment cov matrix                             
  matrix<Type> Psi(n_cell, n_cell); // error term cov matrix                                        
  int indx=0;
  L.setIdentity();
  Psi.setIdentity();
  for (int j=0;j<n_factor;j++){ // column                                                           
    for (int i=0;i<n_cell;i++){ // row                                                              
      L(i,j) = Lt(indx++);
    }
  }
  
  for (int i=0;i<n_cell;i++){
    Type norm2 = Psi(i,i);
    for (int j=0;j<n_factor;j++){
      norm2 += L(i,j) * L(i,j);
    }
    Psi(i,i) = Psi(i,i)/norm2;
    for (int j=0;j<n_factor;j++){
      L(i,j) /= sqrt(norm2);
    }
  }
  
  Sig = L * L.transpose() + Psi;
  //--------------------------------       
  // negative log-likelihood computation  
  MVNORM_t<Type> latent_nll(Sig*dt); // MVN distn for the increments of latent path x           
  vector<Type> reg = Lt*Lt; // l2 (ridge) regularization                                             
  parallel_accumulator<Type> nll(this); // initialize parallel accumulation of negative log-lik      
  nll += lam*(reg.sum()); // negative log-likelihood penalization term
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
    for(int j=0;j<n_bin;j++){                                                                        
      for (int i=0;i<n_cell;i++){                                                                    
        //nll -= dbinom_robust(Y(i,j,u), Type(1), x(i,j,u) - Nt[i] * k[i], true);                      
        nu_x = nu * (x(i,j,u) - Nt[i] * k[i]);
        nll -= Y(i,j,u) * nu_x - log(1 + exp(nu_x));
	Nt(i) += Y(i,j,u);
      } 
    } 
  } 
  return nll; 
}

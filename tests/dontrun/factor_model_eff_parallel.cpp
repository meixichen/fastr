#include <TMB.hpp>

/// @typedef
  /// @brief Standard typedefs for arguments to Eigen functions.
  template <class Type>
  using Matrix_t = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;
  template <class Type>
  using Vector_t = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
  template <class Type>
  using RefMatrix_t = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
  template <class Type>
  using cRefMatrix_t = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >  ;
  template <class Type>
  using RefVector_t = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, 1> >;
  template <class Type>
  using cRefVector_t = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, 1> >;
  template <class Type>
  using RefRowVector_t = Eigen::Ref <Eigen::Matrix<Type, 1, Eigen::Dynamic> >;
  template <class Type>
  using cRefRowVector_t = const Eigen::Ref <const Eigen::Matrix<Type, 1, Eigen::Dynamic> >;


/// Function to efficiently compute the negative log density of a MVN distribution.
  template <class Type>
  class MVN_FA {
  private:
    int n_cell_; ///> Number of cells.
    int n_factor_; ///> Number of neurons.
    Type scale_; ///> Scale parameter multiplied by the MVN variance matrix, e.g., time scale.
    Matrix_t<Type> L_; ///> Normalized loading matrix.
    Vector_t<Type> psi_; ///> Diagonal vector of the uniqueness matrix.
    Vector_t<Type> ipsi_; ///> Inverse of psi.
    Matrix_t<Type> LP_; ///> tLambda * iPsi
    Eigen::LDLT<Matrix_t<Type> > chol_; ///> Cholesky decomp
    Type logdetS_; ///> Log determinant of MVN Sigma.
  public:
    MVN_FA(cRefVector_t<Type>& Lt, int n_cell, int n_factor, Type scale);
    Matrix_t<Type> cov();
    Type operator()(Vector_t<Type> x);
    Type xQx(Vector_t<Type> x);
    void precomp();
    void normalize(cRefVector_t<Type>& Lt);
  };

  template <class Type>
  inline MVN_FA<Type>::MVN_FA(cRefVector_t<Type>& Lt, int n_cell, int n_factor, Type scale){
    n_cell_ = n_cell;
    n_factor_ = n_factor;
    scale_ = scale;
    // memory allocation.
    L_ = Matrix_t<Type>::Zero(n_cell_, n_factor_);
    psi_ = Vector_t<Type>::Zero(n_cell_);
    ipsi_ = Vector_t<Type>::Zero(n_cell_);
    LP_ = Matrix_t<Type>::Zero(n_factor_, n_cell_);
    // compute quantities needed for density evaluation.
    normalize(Lt);
    precomp();
  }

  /// Output the covariance matrix
  template <class Type>
  inline Matrix_t<Type> MVN_FA<Type>::cov(){
    Matrix_t<Type> Sig(n_cell_, n_cell_);
    Sig.noalias() = L_ * L_.transpose();
    Sig += psi_.asDiagonal();
    return Sig;
  }
 /// Create normalized loading matrix and diagonal vector of the uniqueness matrix.
  /// @param[in] Lt Vector of lower triangular elements of the computational basis of the loading ma  trix.
  template <class Type>
  inline void MVN_FA<Type>::normalize(cRefVector_t<Type>& Lt){
    int indx=0;
    int i,j;
    Type norm2;
    // populate loading matrix
    for (j=0;j<n_factor_;j++){ // column
      for (i=j;i<n_cell_;i++){ // row
        L_(i,j) = Lt(indx++);
      }
    }
    // normalize
    for (i=0;i<n_cell_;i++){
      norm2 = Type(1.0);
      for (j=0;j<n_factor_;j++){
        norm2 += L_(i,j) * L_(i,j);
      }
      psi_(i) = Type(1.0)/norm2;
      ipsi_(i) = norm2;
      for (j=0;j<n_factor_;j++){
        L_(i,j) /= sqrt(norm2);
      }
    }
  }
  /// Initialize components for calculating the negative log density.
  template <class Type>
  inline void MVN_FA<Type>::precomp(){
    // Question: do I need to move the following declarations to under Private?
    Matrix_t<Type> Omega(n_factor_, n_factor_); // Omega = I + L' iPsi L
    Matrix_t<Type> iPsi(n_cell_, n_cell_); // inverse of Psi
    Matrix_t<Type> I_k(n_factor_, n_factor_);

    I_k.setIdentity();
    iPsi = ipsi_.asDiagonal();
    Omega.noalias() = L_.transpose() * iPsi * L_;
    Omega += I_k;
    LP_.noalias() = L_.transpose() * iPsi;

    // matrix decomposition
    chol_ = Omega.ldlt();
    Vector_t<Type> Omega_D = chol_.vectorD();
    logdetS_ = Omega_D.array().log().sum() + psi_.array().log().sum(); // log|Omega| + log|Psi|
  }

  /// Calculate the quadratic component of the MVN density.
  /// @param[in] x Vector of length `n_cell` to be evaluated at.
  template <class Type>
  inline Type MVN_FA<Type>::xQx(Vector_t<Type> x){
    Vector_t<Type> z(n_factor_);
    Vector_t<Type> iOmega_z(n_factor_);
    z.noalias() = LP_ * x;
    iOmega_z = chol_.solve(z);
    Type res1 = iOmega_z.dot(z);
    Vector_t<Type> x2 = (x.array() * x.array()).matrix();
    Type res2 = x2.dot(ipsi_);
    return res2 - res1;
  }

  /// Evaluate the negative log density of a zero-mean MVN with cov matrix Sigma=LL'+Psi.
  /// @param[in] x Vector of length `n_cell` to be evaluated at.
  template <class Type>
  inline Type MVN_FA<Type>::operator()(Vector_t<Type> x){
    return Type(0.5)*logdetS_ + Type(0.5)* x.size() *log(scale_) + Type(0.5)/scale_ *xQx(x) + x.size    ()*Type(log(sqrt(2.0*M_PI))); // scaling by time scale is applied here
  }

////////////////////// Start TMB template ///////////////////////
template<class Type>
Type objective_function<Type>::operator() (){
  // data inputs                                                                                       
  DATA_INTEGER(n_factor); // number of factors                                                         
  DATA_SCALAR(dt); // length of the time bin                                                           
  DATA_ARRAY(Y); //  q x n x r array of 0s and 1s.                                                     
  DATA_SCALAR(lam); // regularization parameter                                                        

  // parameters                                                                                        
  PARAMETER_VECTOR(log_k); // log thresholds (k>0)                                                     
  PARAMETER_VECTOR(log_a); // log drift (0<drift<1)                                                    
  PARAMETER_VECTOR(Lt); // dq-d(d-1)/2 vector of the lower triangular elements of L (by column)        
  PARAMETER_ARRAY(x); // q x n x r neuron paths                                                        

  using namespace density;
  // transformed data                                                                                  
  vector<int> Y_dim = Y.dim; // dimension of data Y (a 3d array)                                       
  int n_cell = Y_dim(0); // number of neurons                                                          
  int n_bin = Y_dim(1); // number of time bins                                                         
  int n_trial = Y_dim(2); // number of trials                                                          

  // transformed parameters                                                                            
  vector<Type> k = exp(log_k);                                                                       
  vector<Type> alpha = exp(log_a);                                                                   
  // negative log-likelihood computation  
  MVN_FA<Type> latent_nll(Lt, n_cell, n_factor, dt); // efficiently compute density of MVN 
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
        nll -= dbinom_robust(Y(i,j,u), Type(1), x(i,j,u) - Nt[i] * k[i], true);                      
        Nt(i) += Y(i,j,u);
      }
    }
  }
  return nll;
}



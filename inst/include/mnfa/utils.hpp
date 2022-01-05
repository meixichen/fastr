/// @file utils.hpp
///
/// @brief Utilities for `mnfa`.

#ifndef MNFA_UTILS_HPP
#define MNFA_UTILS_HPP

namespace mnfa {

  /// @typedef
  /// @brief Standard typedefs for arguments to Eigen functions.
  template <class Type>
  using Matrix_t = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;
  template <class Type>
  using Vector_t = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
  template <class Type>
  using RefMatrix_t = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
  template <class Type>
  using cRefMatrix_t = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
  template <class Type>
  using RefVector_t = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, 1> >;
  template <class Type>
  using cRefVector_t = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, 1> >;
  template <class Type>
  using RefRowVector_t = Eigen::Ref <Eigen::Matrix<Type, 1, Eigen::Dynamic> >;
  template <class Type>
  using cRefRowVector_t = const Eigen::Ref <const Eigen::Matrix<Type, 1, Eigen::Dynamic> >;
  
  /// Function object for computing the variance matrix of a latent path increment.
  template <class Type>
  class SigmaFA {
  private:
    int n_cell_; ///> Number of cells.
    int n_factor_; ///> Number of neurons.
    Matrix_t<Type> L_; ///> Loading matrix.
    Vector_t<Type> psi_; ///> Diagonal vector.
  public:
    /// Class constructor.
    SigmaFA(int n_cell, int n_factor);
    void create(RefMatrix_t<Type> Sig, cRefVector_t<Type>& Lt);
  };

  /// @param[in] n_cell Number of cells.
  /// @param[in] n_factor Number of neurons.
  template <class Type>
  inline SigmaFA<Type>::SigmaFA(int n_cell, int n_factor) {
    n_cell_ = n_cell;
    n_factor_ = n_factor;
    // memory allocation
    L_ = Matrix_t<Type>::Zero(n_cell_, n_factor_);
    psi_ = Vector_t<Type>::Zero(n_cell_);
  }

  template <class Type>
  inline void SigmaFA<Type>::create(RefMatrix_t<Type> Sig,
				    cRefVector_t<Type>& Lt) {
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
      for (j=0;j<n_factor_;j++){
        L_(i,j) /= sqrt(norm2);
      }
    }
    // evaluate Sig
    Sig.noalias() = L_ * L_.transpose();
    Sig += psi_.asDiagonal();
    return;
  }

  ////////////////* NEW: MVN_FA will be used to replace SigmaFA
  //
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
    Matrix_t<Type> Q_; ///> Inverse of MVN cov mat Sigma.
    Type logdetS_; ///> Log determinant of MVN Sigma.
  public:
    MVN_FA(cRefVector_t<Type>& Lt, int n_cell, int n_factor, 
	Type scale, bool use_atomic=true);
    Matrix_t<Type> cov();
    Type operator()(Vector_t<Type> x);
    Type xQx(Vector_t<Type> x);
    void setQ(bool use_atomic=true);
    void normalize(cRefVector_t<Type>& Lt);  
  };

  template <class Type>
  inline MVN_FA<Type>::MVN_FA(cRefVector_t<Type>& Lt, int n_cell, int n_factor, 
      Type scale, bool use_atomic){
    n_cell_ = n_cell;
    n_factor_ = n_factor;
    scale_ = scale;
    // memory allocation.
    L_ = Matrix_t<Type>::Zero(n_cell_, n_factor_);
    psi_ = Vector_t<Type>::Zero(n_cell_);
    ipsi_ = Vector_t<Type>::Zero(n_cell_); 
    Q_ = Matrix_t<Type>::Zero(n_cell_, n_cell_);
    // compute quantities needed for density evaluation.
    normalize(Lt);
    setQ(use_atomic);
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
  /// @param[in] Lt Vector of lower triangular elements of the computational basis of the loading matrix.
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
  inline void MVN_FA<Type>::setQ(bool use_atomic){
    // Question: do I need to move the following declarations to under Private?
    Matrix_t<Type> Omega(n_factor_, n_factor_); // Omega = I + L' iPsi L
    Matrix_t<Type> iOmega(n_factor_, n_factor_); // inverse of Omega
    Matrix_t<Type> iPsi(n_cell_, n_cell_); // inverse of Psi
    Matrix_t<Type> PLOLP(n_cell_, n_cell_); // PLOLP = iPsi*L*iOmega*Lt*iPsi. Q=iPsi-PLOLP.
    Matrix_t<Type> I_k(n_factor_, n_factor_);
    
    I_k.setIdentity();
    iPsi = ipsi_.asDiagonal();
    Omega.noalias() = L_.transpose() * iPsi * L_; 
    Omega += I_k;
    
    // matrix decomposition
    Eigen::LDLT<Matrix_t<Type> > ldlt(Omega);
    iOmega = ldlt.solve(I_k);
    Vector_t<Type> Omega_D = ldlt.vectorD();
    logdetS_ = Omega_D.array().log().sum() + psi_.array().log().sum(); // log|Omega| + log|Psi|
    PLOLP = iPsi * L_ * iOmega * L_.transpose() * iPsi; 
    Q_ = iPsi - PLOLP;
  }

  /// Calculate the quadratic component of the MVN density.
  /// @param[in] x Vector of length `n_cell` to be evaluated at.
  template <class Type>
  inline Type MVN_FA<Type>::xQx(Vector_t<Type> x){
    Vector_t<Type> Qx = Q_ * x;
    return x.dot(Qx);
  }
  
  /// Evaluate the negative log density of a zero-mean MVN with cov matrix Sigma=LL'+Psi.
  /// @param[in] x Vector of length `n_cell` to be evaluated at.
  template <class Type>
  inline Type MVN_FA<Type>::operator()(Vector_t<Type> x){
    return Type(0.5)*logdetS_ + Type(0.5)* x.size() *log(scale_) + Type(0.5)/scale_ *xQx(x) + x.size()*Type(log(sqrt(2.0*M_PI))); // scaling by time scale is applied here
  }

  /// Get the lower triangular elements (by column) of a correlation matrix as a vector.
  ///
  /// @param[out] rho Vector into which to store the output.
  /// @param[in] Sig A correlation matrix.
  /// @param[in] n_cell Number of cells.
  template <class Type>
  inline void Sig2rho(RefVector_t<Type> rho, cRefMatrix_t<Type>& Sig, 
		      int n_cell){
    int indx = 0;
    int i,j;
    for (j=0;j<n_cell;j++){ // column
      for (i=j+1;i<n_cell;i++){ // row
        rho(indx++) = Sig(i,j);
      }
    }
    return;
  }

} // end namespace mnfa

#endif

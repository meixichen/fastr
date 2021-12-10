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

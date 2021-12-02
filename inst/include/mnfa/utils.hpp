/// @file utils.hpp
///
/// @brief Utilities for `mnfa`.

#ifndef MNFA_UTILS_HPP
#define MNFA_UTILS_HPP

namespace mnfa {

  /// @typedef
  /// @brief Standard typedefs for arguments to Eigen functions.
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
  
  
  /// Compute the variance matrix for the latent path increment.
  ///
  /// @param[out] Sig Matrix into which to store the output. Sig = LL' + uniqueness so that diag(Sig)=1.
  /// @param[in] Lt dq-d(d-1)/2 vector of the lower triangular elements of L (by column).
  /// @param[in] n_cell Number of cells.
  /// @param[in] n_factor Number of neurons.
  template <class Type>
  inline void create_Sig(RefMatrix_t<Type> Sig, cRefVector_t<Type>& Lt,
  	              Type n_cell, Type n_factor) {
    matrix<Type> L(n_cell, n_factor); // loading matrix
    matrix<Type> Psi(n_cell, n_cell); // error term cov matrix
    int indx=0;
    int i,j;
    L.setIdentity();
    Psi.setIdentity();
    for (j=0;j<n_factor;j++){ // column
      for (i=j;i<n_cell;i++){ // row
        L(i,j) = Lt(indx++);
      }
    }
    
    for (i=0;i<n_cell;i++){
      Type norm2 = Psi(i,i);
      for (j=0;j<n_factor;j++){
        norm2 += L(i,j) * L(i,j);
      }
      Psi(i,i) = Psi(i,i)/norm2;
      for (j=0;j<n_factor;j++){
        L(i,j) /= sqrt(norm2);
      }
    }
    
    Sig = L * L.transpose() + Psi;
    return;
  }
  
  /// Get the lower triangular elements (by column) of a correlation matrix as a vector.
  ///
  /// @param[out] rho Vector into which to store the output.
  /// @param[in] Sig A correlation matrix                                                         
  template <class Type>
  inline void Sig2rho(RefVector_t<Type> rho, cRefMatrix_t<Type>& Sig, 
               Type n_cell){
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

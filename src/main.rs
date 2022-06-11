use num_complex::Complex64;
use qrusty::SimplePauli::* ;
use ndarray::Array2 ;
use sprs::* ;

fn main() {
    print!("{:?}", I.to_matrix()) ;
    let sp_mat = kronecker_product(I.to_matrix().view(),X.to_matrix().view()) ;
    let mut dmat = sp_mat.to_dense() ;
    print!("{}", dmat) ;
    
}

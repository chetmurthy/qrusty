use num_complex::Complex64;
use qrusty::SimplePauli::* ;
use ndarray::Array2 ;


fn main() {
    print!("{:?}", I.to_matrix()) ;
    let sp_mat : sprs::CsMat<Complex64> = I.to_matrix() ;
    let mut dmat = sp_mat.to_dense() ;
    let mut dmat2 = sp_mat.to_dense() ;
    assert!(dmat == dmat2) ;
    print!("{}", dmat) ;
    
}

use num_complex::Complex64;
use qrusty::* ;
use qrusty::SimplePauli::* ;
use ndarray::Array2 ;
use sprs::* ;

fn main() {
    let p = Pauli::new("I").unwrap() ;

    let sp_mat = p.to_matrix() ;
    
    println!("I={:?}\n", sp_mat) ;

    let sp_mat_accel = p.to_matrix_accel() ;

    println!("accel I={:?}\n", sp_mat_accel) ;

    let (data, indices, indptr) = p.to_triplets() ;
    println!("data={:?}\nindices={:?}\nindptr={:?}\n", data, indices, indptr) ;

}

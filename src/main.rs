use num_complex::Complex64;
use qrusty::* ;
use qrusty::SimplePauli::* ;
use ndarray::Array2 ;
use sprs::* ;

use qrusty::util::* ;
use qrusty::util::list::* ;
use qrusty::util::BinaryTreeFold ;

fn main() {
    let mut bt = BinaryTreeFold::begin(atom(0), cons) ;
    (1..8).for_each(|n| bt.add(atom(n))) ;
    let rv = bt.end() ;
    println!("{}", (*rv).to_string()) ;
}

fn main2() {
    let p = Pauli::new("Y").unwrap() ;

    let sp_mat = p.to_matrix() ;
    
    println!("Y={:?}\n", sp_mat) ;

    let sp_mat_accel = p.to_matrix_accel() ;

    println!("accel Y={:?}\n", sp_mat_accel) ;

    let (data, indices, indptr) = p.to_triplets() ;
    println!("data={:?}\nindices={:?}\nindptr={:?}\n", data, indices, indptr) ;

}

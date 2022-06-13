use std::ops::Add;
use num_complex::Complex64;
use qrusty::* ;
use qrusty::SimplePauli::* ;
use ndarray::Array2 ;
use sprs::* ;

use qrusty::util::* ;
use qrusty::util::list::* ;
use qrusty::util::BinaryTreeFold ;

fn main() {
    let spop = SparsePauliOp::new(
        PauliList::from_labels(&vec![&"I".to_string(),&"X".to_string()]).unwrap(),
        vec![Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)]).unwrap() ;
    let sp_mat = spop.to_matrix() ;
    println!("{}",sp_mat.to_dense()) ;
}

fn main2() {
    let p = Pauli::new(&"Y".to_string()).unwrap() ;

    let sp_mat = p.to_matrix() ;
    
    println!("Y={:?}\n", sp_mat) ;

    let sp_mat_accel = p.to_matrix_accel() ;

    println!("accel Y={:?}\n", sp_mat_accel) ;

    let (data, indices, indptr) = p.to_triplets() ;
    println!("data={:?}\nindices={:?}\nindptr={:?}\n", data, indices, indptr) ;

}
fn main3() {
    let mut bt = BinaryTreeFold::begin(atom(0), cons) ;
    (1..8).for_each(|n| bt.add(atom(n))) ;
    let rv = bt.end() ;
    println!("{}", (*rv).to_string()) ;
}
fn main4() {
    use sprs::{CsMat, CsVec};
    let eye : CsMat<f64> = CsMat::eye(2);
    let b = &eye + &eye;
    println!("{}",b.to_dense()) ;
}

fn main5() {
    let p = Pauli::new(&"I".to_string()).unwrap() ;
    let sp_mat = p.to_matrix() ;
    let b = &sp_mat + &sp_mat ;
    println!("{}",b.to_dense()) ;
}

fn main6() {
    let p = Pauli::new(&"I".to_string()).unwrap() ;
    let q = Pauli::new(&"X".to_string()).unwrap() ;
    let p_mat = p.to_matrix() ;
    let q_mat = q.to_matrix() ;
    let a = Complex64::new(2.0, 0.0);
    let b = Complex64::new(4.0, 0.0);
    let b  = sprs::binop::csmat_binop(p_mat.view(),  q_mat.view(), |x,y| a * x + b * y) ;
    println!("{}",b.to_dense()) ;
}

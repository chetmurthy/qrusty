use num_complex::Complex64;
use qrusty::SimplePauli::* ;
use ndarray::Array2 ;
use sprs::* ;

fn main() {
    print!("{:?}", I.to_matrix()) ;
    let mut sp_mat = I.to_matrix() ;
    sp_mat.scale(Complex64::new(2.0, 0.0)) ;
    print!("{:?}", sp_mat) ;
    
let a = vec![1, 2, 3, 4];

// the sum of all of the elements of the array
    let sum = a[0..a.len()-1].iter().rev().fold(a[a.len()-1], |acc, x| { println!("{}", x); acc + x });

assert_eq!(sum, 10);

}

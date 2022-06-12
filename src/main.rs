use num_complex::Complex64;
use qrusty::SimplePauli::* ;
use ndarray::Array2 ;
use sprs::* ;

fn main() {
    print!("{:?}\n", I.to_matrix()) ;
    print!("{:?}\n", kronecker_product(I.to_matrix().view(), X.to_matrix().view())) ;
    
let a = vec![1, 2, 3, 4];

// the sum of all of the elements of the array
    let sum = a[0..a.len()-1].iter().rev().fold(a[a.len()-1], |acc, x| { println!("{}", x); acc + x });

assert_eq!(sum, 10);

}

#[macro_use]
extern crate lazy_static;

use num_complex::Complex64;
use regex::Regex;
use sprs::{CsMat, TriMat, TriMatI, CsMatI, kronecker_product};

#[derive(Debug, PartialEq)]
pub enum SimplePauli {
    I, X, Y, Z,
}

use crate::SimplePauli::* ;
impl SimplePauli {
    pub fn new(c : char)-> Result<SimplePauli, String> {
        match c {
            'I' => Ok(I),
            'X' => Ok(X),
            'Y' => Ok(Y),
            'Z' =>  Ok(Z),
            _ => Err(String::from("internal error: malformed pauli"))
        }
    }

/*
>>> Pauli("IXYZ").x
array([False,  True,  True, False])
>>> Pauli("IXYZ").z
array([ True,  True, False, False])
>>> 
*/
    pub fn x(&self) -> bool {
        match self {
            I|Z => false,
            X|Y => true
        }
    }
    pub fn z(&self) -> bool {
        match self {
            Y|Z => true,
            I|X => false
        }
    }

    pub fn to_matrix(&self) -> sprs::CsMatI<Complex64,  u64> {
        let mut a = TriMatI::<Complex64, u64>::new((2, 2));
        match self {
            I => {
                a.add_triplet(0, 0, Complex64::new(1.0, 0.0)) ;
                a.add_triplet(1, 1, Complex64::new(1.0, 0.0)) ;
            },
            Z => {
                a.add_triplet(0, 0, Complex64::new(1.0, 0.0)) ;
                a.add_triplet(1, 1, Complex64::new(-1.0, 0.0)) ;
            },
            X => {
                a.add_triplet(0, 1, Complex64::new(1.0, 0.0)) ;
                a.add_triplet(1, 0, Complex64::new(1.0, 0.0)) ;
            },
            Y => {
                a.add_triplet(0, 1, Complex64::new(0.0, -1.0)) ;
                a.add_triplet(1, 0, Complex64::new(0.0, 1.0)) ;
            },
        } ;
        a.to_csr()
    }
}

mod accel ;

#[derive(Debug)]
pub struct Pauli {
    phase : usize,
    paulis : Vec<SimplePauli>,
}
impl Pauli {

    pub fn num_qubits(&self) -> usize { self.paulis.len() }
    fn parse_label(s : &str) -> Result<(usize, Vec<SimplePauli>),  String> {
        lazy_static! {
            static ref RE : Regex = Regex::new(r"^([+-]?)1?([ij]?)([IXYZ]+)$").unwrap();
        }
        let caps = RE.captures(s).ok_or("error: malformed label")? ;

        let sign = caps.get(1).map_or("", |m| m.as_str());
        let imag = caps.get(2).map_or("", |m| m.as_str());
        let paulistr = caps.get(3).map_or("", |m| m.as_str());

        let sign = match sign { ""|"+" => Ok(true), "-" => Ok(false), _ => Err("internal error: malformed sign") } ?;
        let phase : usize = match imag {
            "" => Ok(0),
            "i"|"j" => Ok(1),
            _ => Err("internal error: malformed phase")
        }?;
        let phase : usize = if sign { phase } else { phase+2 } ;
        
        let mut paulis = Vec::new() ;
        for c in paulistr.chars().rev() {
            paulis.push(SimplePauli::new(c)?)
        }
        if 0 == paulis.len() { Err(String::from("internal error: no paulis")) }
        else { Ok((phase, paulis)) }
    }

    pub fn new(s : &str) -> Result<Pauli, String> {
        let (phase, paulis) = Pauli::parse_label(s)? ;
        Ok(Pauli{ phase, paulis })
    }

    pub fn phase(&self) -> usize { self.phase }
    pub fn coeff(&self) -> Complex64 {
        match self.phase {
            0 => Complex64::new(1.0, 0.0),
            1 => Complex64::new(0.0, 1.0),
            2 => Complex64::new(-1.0, 0.0),
            3 => Complex64::new(0.0, -1.0),
            _ => panic!("internal error: phase should never be outside [0..3]")
        }
    }
    pub fn xs(&self) -> Vec<bool> { self.paulis.iter().map(|p| p.x()).collect() }
    pub fn zs(&self) -> Vec<bool> { self.paulis.iter().map(|p| p.z()).collect() }

    pub fn to_matrix(&self) -> sprs::CsMatI<Complex64, u64> {
        let l = &(self.paulis) ;
        let mut sp_mat = l[1..]
            .iter()
            .fold(l[0].to_matrix(),
                  |acc, x| {
                      println!("{:?}", x);
                      kronecker_product(x.to_matrix().view(), acc.view())
                  });
        sp_mat.scale(self.coeff()) ;
        sp_mat
    }

    pub fn to_triplets(&self) -> (Vec<Complex64>, Vec<u64>, Vec<u64>) {
        let zs = self.zs() ;
        let xs = self.xs() ;
        let phase = self.phase() as i64 ;
        let coeff = Complex64::new(1.0, 0.0) ;
        let group_phase = false ;
        let (data, indices, indptr) = accel::rust_make_data(&zs, &xs, coeff, phase, group_phase, false).expect("accelerated matrix creation failed") ;
        (data, indices, indptr)
    }

    pub fn to_triplets_ffi(&self) -> (Vec<Complex64>, Vec<u64>, Vec<u64>) {
        let zs = self.zs() ;
        let xs = self.xs() ;
        let phase = self.phase() as i64 ;
        let coeff = Complex64::new(1.0, 0.0) ;
        let group_phase = false ;
        let (data, indices, indptr) = accel::rust_make_data(&zs, &xs, coeff, phase, group_phase, true).expect("accelerated matrix creation failed") ;
        (data, indices, indptr)
    }

    pub fn to_matrix_accel(&self) -> sprs::CsMatI<Complex64, u64> {
        let (data, indices, indptr) = self.to_triplets() ;
        let dim = 1 << self.num_qubits() ;
        let trimat = TriMatI::<Complex64, u64>::from_triplets((dim, dim),
                                                              indptr,
                                                              indices,
                                                              data) ;
        let sp_mat : CsMatI<Complex64, u64> = trimat.to_csr() ;
        sp_mat
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;
    use ndarray::array ;
    use ndarray::prelude::*;
    use sprs::*;

    use crate::Pauli ;
    use crate::SimplePauli::* ;

    #[test]
    fn parse_labels() {
        assert!(Pauli::parse_label("W").is_err());
        assert!(Pauli::parse_label("").is_err());
        assert!(Pauli::parse_label("+i").is_err());
        assert_eq!(Pauli::parse_label("I"), Ok( (0 as usize, vec![I]) ));
        assert_eq!(Pauli::parse_label("+I"), Ok( (0 as usize, vec![I]) ));
        assert_eq!(Pauli::parse_label("+iI"), Ok( (1 as usize, vec![I]) ));
        assert_eq!(Pauli::parse_label("+jI"), Ok( (1 as usize, vec![I]) ));
        assert_eq!(Pauli::parse_label("-1jI"), Ok( (3 as usize, vec![I]) ));
        assert_eq!(Pauli::parse_label("-1I"), Ok( (2 as usize, vec![I]) ));
        assert_eq!(Pauli::parse_label("-1jIX"), Ok( (3 as usize, vec![X, I]) ));
        assert_eq!(Pauli::parse_label("IXYZ"), Ok( (0 as usize, vec![Z, Y, X, I]) ));

/*
>>> (Pauli("I").z, Pauli("I").x)
(array([False]), array([False]))
>>> (Pauli("X").z, Pauli("X").x)
(array([False]), array([ True]))
>>> (Pauli("Y").z, Pauli("Y").x)
(array([ True]), array([ True]))
>>> (Pauli("Z").z, Pauli("Z").x)
(array([ True]), array([False]))
>>> 
*/
        assert_eq!(I.z(), false) ;
        assert_eq!(I.x(), false) ;
        assert_eq!(X.z(), false) ;
        assert_eq!(X.x(), true) ;
        assert_eq!(Y.z(), true) ;
        assert_eq!(Y.x(), true) ;
        assert_eq!(Z.z(), true) ;
        assert_eq!(Z.x(), false) ;


/*
>>> Pauli("IXYZ").x
array([False,  True,  True, False])
>>> Pauli("IXYZ").z
array([ True,  True, False, False])
>>> 
*/

        {
            let p = Pauli::new("IXYZ") ;
            assert!(p.is_ok()) ;
            let p = p.unwrap() ;
            assert_eq!(p.xs(),
                       vec![false,  true,  true, false]);
            assert_eq!(p.zs(),
                       vec![true,  true, false, false]);
        }

/*
Pauli::parse_label("-IIIIIIIIIIIIIIIIIIYXXY"), Ok(2 as usize, vec![ IIIIIIIIIIIIIIIIIIYXXY ]
z=[ True False False  True False False False False False False False False
 False False False False False False False False False False]
x=[ True  True  True  True False False False False False False False False
 False False False False False False False False False False]
phase=2

*/
        {
            let p = Pauli::new("-IIIIIIIIIIIIIIIIIIYXXY") ;
            assert!(p.is_ok()) ;
            let p = p.unwrap() ;
            assert_eq!(p.phase(), 2 as usize) ;
            assert_eq!(p.zs(),
                       vec![ true, false, false,  true,
                             false, false, false, false,
                             false, false, false, false,
                             false, false, false, false,
                             false, false, false, false,
                             false, false,]);
            assert_eq!(p.xs(),
                       vec![ true,  true,  true,  true,
                             false, false, false, false,
                             false, false, false, false,
                             false, false, false, false,
                             false, false, false, false,
                             false, false,]);
        }

    }
    fn simple_pauli_matrices() {
        let one = Complex64::new(1.0, 0.0) ;
        let minus_one = Complex64::new(-1.0, 0.0) ;
        let zero = Complex64::new(0.0, 0.0) ;
        let i = Complex64::new(0.0, 1.0) ;
        let minus_i = Complex64::new(0.0, -1.0) ;
        assert_eq!(I.to_matrix().to_dense(),
                   array![[one, zero],
                          [zero, one]]) ;
        assert_eq!(X.to_matrix().to_dense(),
                   array![[zero, one],
                          [one, zero]]) ;
        assert_eq!(Y.to_matrix().to_dense(),
                   array![[zero, minus_i],
                          [i, zero]]) ;
        assert_eq!(Z.to_matrix().to_dense(),
                   array![[one, zero],
                          [zero, minus_one]]) ;
        assert_eq!(kronecker_product(I.to_matrix().view(), X.to_matrix().view()).to_dense(),
                   array![[zero, one, zero, zero],
                          [one, zero, zero, zero],
                          [zero, zero, zero, one],
                          [zero, zero, one, zero]]) ;
    }

    fn pauli_matrices() {
        let one = Complex64::new(1.0, 0.0) ;
        let minus_one = Complex64::new(-1.0, 0.0) ;
        let zero = Complex64::new(0.0, 0.0) ;
        let i = Complex64::new(0.0, 1.0) ;
        let minus_i = Complex64::new(0.0, -1.0) ;

        let p = Pauli::new("IX") ;
        assert!(p.is_ok()) ;
        let p = p.unwrap() ;
        assert_eq!(p.to_matrix().view().to_dense(),
                   kronecker_product(I.to_matrix().view(), X.to_matrix().view()).to_dense()) ;

        let p = Pauli::new("XI") ;
        assert!(p.is_ok()) ;
        let p = p.unwrap() ;
        assert_eq!(p.to_matrix().view().to_dense(),
                   kronecker_product(X.to_matrix().view(), I.to_matrix().view()).to_dense()) ;

    }

    fn accel() {
        let p = Pauli::new("I").unwrap() ;

        assert_eq!(p.to_matrix().view().to_dense(),
                   p.to_matrix_accel().view().to_dense()) ;
    }

}

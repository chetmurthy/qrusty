#[macro_use]
extern crate lazy_static;

use std::cmp::min;

#[macro_use] extern crate impl_ops;
use std::ops;

use std::ops::Index;
use num_complex::Complex64;
use regex::Regex;
use rayon::prelude::*;
use sprs::{CompressedStorage, TriMatI, CsMatI, kronecker_product};

use pyo3::PyErr ;
use pyo3::exceptions::PyException;

mod accel ;
pub mod util ;
pub mod fixtures ;
pub mod rawio ;


#[derive(Debug, PartialEq, Eq, Clone)]
pub struct QrustyErr {
    s : &'static str
}

impl QrustyErr {
    pub fn new(s : &'static str) -> QrustyErr {
        QrustyErr { s }
    }
    pub fn message(&self) -> &'static str {
        self.s
    }
}

impl std::convert::From<QrustyErr> for pyo3::PyErr {
    fn from(err: QrustyErr) -> PyErr {
        PyException::new_err(err.message())
    }
}


#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum SimplePauli {
    I, X, Y, Z,
}

use crate::SimplePauli::* ;
impl SimplePauli {
    pub fn new(c : char)-> Result<SimplePauli, QrustyErr> {
        match c {
            'I' => Ok(I),
            'X' => Ok(X),
            'Y' => Ok(Y),
            'Z' =>  Ok(Z),
            _ => Err(QrustyErr::new("internal error: malformed pauli"))
        }
    }
    pub fn label(self)-> char {
        match self {
            I => 'I',
            X => 'X',
            Y => 'Y',
            Z => 'Z',
        }
    }

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


#[derive(Debug, Clone)]
pub struct Pauli {
    base_phase : usize,
    paulis : Vec<SimplePauli>,
}
impl Pauli {

    pub fn num_qubits(&self) -> usize { self.paulis.len() }
    fn parse_label(s : &str) -> Result<(usize, Vec<SimplePauli>),  QrustyErr> {
        lazy_static! {
            static ref RE : Regex = Regex::new(r"^([+-]?)1?([ij]?)([IXYZ]+)$").unwrap();
        }
        let caps = RE.captures(s).ok_or(QrustyErr::new("error: malformed label"))? ;

        let sign = caps.get(1).map_or("", |m| m.as_str());
        let imag = caps.get(2).map_or("", |m| m.as_str());
        let paulistr = caps.get(3).map_or("", |m| m.as_str());

        let sign = match sign { ""|"+" => Ok(true), "-" => Ok(false), _ => Err(QrustyErr::new("internal error: malformed sign")) } ?;
        let phase : usize = match imag {
            "" => Ok(0),
            "i"|"j" => Ok(1),
            _ => Err(QrustyErr::new("internal error: malformed phase"))
        }?;
        let phase : usize = if sign { phase } else { phase+2 } ;
        
        let mut paulis = Vec::new() ;
        for c in paulistr.chars().rev() {
            paulis.push(SimplePauli::new(c)?)
        }
        if 0 == paulis.len() { Err(QrustyErr::new("internal error: no paulis")) }
        else { Ok((phase, paulis)) }
    }
    pub fn label(&self) -> String {
        self.paulis.iter().rev().into_iter()
            .map(|p| p.label())
            .collect::<String>()
    }

    pub fn new(s : &str) -> Result<Pauli, QrustyErr> {
        let (base_phase, paulis) = Pauli::parse_label(s)? ;
        Ok(Pauli{ base_phase, paulis })
    }

    pub fn x_indices(&self) -> u64 {
        let x_indices : u64 =
            self.xs().iter().enumerate()
            .filter(|(_,x)| **x)
            .map(|(i,_)| 1<<i)
            .sum() ;
        x_indices
    }
    pub fn z_indices(&self) -> u64 {
        let z_indices : u64 =
            self.zs().iter().enumerate()
            .filter(|(_,z)| **z)
            .map(|(i,_)| 1<<i)
            .sum() ;
        z_indices
    }
    
    pub fn base_phase(&self) -> usize { self.base_phase }
    pub fn phase(&self) -> usize {
        (self.base_phase() + self.paulis.iter().fold(0, |acc, p| acc + match p { Y => 1, _ => 0 })) % 4
    }
    pub fn base_coeff(&self) -> Complex64 {
        match self.base_phase() {
            0 => Complex64::new(1.0, 0.0),
            1 => Complex64::new(0.0, 1.0),
            2 => Complex64::new(-1.0, 0.0),
            3 => Complex64::new(0.0, -1.0),
            _ => panic!("internal error: base_phase should never be outside [0..3]")
        }
    }
    pub fn coeff(&self) -> Complex64 {
        match self.phase() {
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
                      kronecker_product(x.to_matrix().view(), acc.view())
                  });
        sp_mat.scale(self.base_coeff()) ;
        sp_mat
    }

    pub fn to_unsafe_vectors(&self) -> accel::UnsafeVectors {
        let zs = self.zs() ;
        let xs = self.xs() ;
        let phase = self.phase() as i64 ;
        let coeff = Complex64::new(1.0, 0.0) ;
        let group_phase = false ;
        let usv = accel::make_unsafe_vectors(&zs, &xs, coeff, phase, group_phase, false).expect("accelerated matrix creation failed") ;
        usv
    }

    pub fn to_unsafe_vectors_ffi(&self) -> accel::UnsafeVectors {
        let zs = self.zs() ;
        let xs = self.xs() ;
        let phase = self.phase() as i64 ;
        let coeff = Complex64::new(1.0, 0.0) ;
        let group_phase = false ;
        let usv = accel::make_unsafe_vectors(&zs, &xs, coeff, phase, group_phase, true).expect("accelerated matrix creation failed") ;
        usv
    }

    pub fn to_matrix_accel(&self) -> sprs::CsMatI<Complex64, u64> {
        let usv = self.to_unsafe_vectors_ffi() ;
        let dim = 1 << self.num_qubits() ;

        unsafe {
            CsMatI::<Complex64, u64, u64>::new_unchecked(
                CompressedStorage::CSR,
                (dim, dim),
                usv.indptr,
                usv.indices,
                usv.data,
            )
        }
    }
    pub fn to_matrix_rowwise(&self) -> sprs::CsMatI<Complex64, u64> {
        let v = vec![(self.clone(), Complex64::new(1.0, 0.0))] ;
        let members = &v[..] ;
        let trimat = accel::rowwise::make_trimat(&members) ;
        trimat.to_csr()
    }

}
pub type PauliSummand = (Pauli, Complex64) ;
pub type PauliSum = Vec<PauliSummand> ;

impl<'a> Index<isize> for SparsePauliOp {
    type Output = PauliSummand ;
    fn index(self : &SparsePauliOp, index: isize) -> &PauliSummand {
        let self_len = self.members.len() as isize;
        let idx = (((index % self_len) + self_len) % self_len) as usize;
        &self.members[idx]
    }
}

#[derive(Debug, Clone)]
pub struct SparsePauliOp {
    members : PauliSum
}
impl SparsePauliOp {
    pub fn members(&self) -> &PauliSum {
        &self.members
    }
    pub fn num_qubits(&self) -> usize { self.members[0].0.num_qubits() }
    pub fn new(paulis : Vec<Pauli>, coeffs : &[Complex64]) -> Result<SparsePauliOp, QrustyErr> {
        if paulis.len() != coeffs.len() {
            Err(QrustyErr::new("SparsePauliOp::new: paulis and coeffs must have same length"))
        }
        else if paulis.len() == 0 {
            Err(QrustyErr::new("SparsePauliOp::new: at least one pauli must be supplied"))
        }
        else {
            let num_qubits = paulis[0].num_qubits() ;
            if paulis.iter().any(|p| p.num_qubits() != num_qubits) {
            Err(QrustyErr::new("SparsePauliOp::new: all supplied paulis must have the same #qubits"))
            }
            else {
            let mut members : PauliSum = Vec::new() ;
            paulis.iter()
                .zip(coeffs.iter())
                .for_each(|(p,c)| {
                    members.push((p.clone(),*c)) ;
                }) ;
            Ok(SparsePauliOp { members })
            }
        }
    }
    pub fn from_slice(summands : &[PauliSummand]) -> SparsePauliOp {
        let members : Vec<PauliSummand> =
            summands
            .iter()
            .map(|s| s.clone())
            .collect() ;
        SparsePauliOp { members }
    }

    pub fn from_labels(l : &[&str], coeffs : &[Complex64]) -> Result<SparsePauliOp, QrustyErr> {
        let mut v = Vec::new() ;
        for s in l.iter() {
            let p = Pauli::new(s)? ;
            v.push(p) ;
        }
        SparsePauliOp::new(v, coeffs)
    }
    pub fn from_labels_simple(l : &[&str]) -> Result<SparsePauliOp, QrustyErr> {
        let coeffs : Vec<Complex64> = l.iter()
            .map(|_| Complex64::new(1.0, 0.0))
            .collect() ;
        SparsePauliOp::from_labels(l, &coeffs[..])
    }

    pub fn to_matrix(&self) -> sprs::CsMatI<Complex64, u64> {
        let zeroth = &self.members[0] ;
        let mut sum = zeroth.0.to_matrix().clone() ;
        sum.scale(zeroth.1) ;
        for i in 1..self.members.len() {
            let ith = &self.members[i] ;
            let m = ith.0.to_matrix() ;
            let coeff = ith.1 ;
            sum = sprs::binop::csmat_binop(sum.view(), m.view(), |x,y| x + coeff * y) ;
        }
        sum
    }

    pub fn to_matrix_binary(&self) -> sprs::CsMatI<Complex64, u64> {
        fn addmul(l: (Complex64, sprs::CsMatI<Complex64, u64>),
                  r: (Complex64, sprs::CsMatI<Complex64, u64>))
                  -> (Complex64, sprs::CsMatI<Complex64, u64>) {
            (Complex64::new(1.0,  0.0),
             sprs::binop::csmat_binop(l.1.view(), r.1.view(), |x,y| l.0 * x + r.0 * y))
        }

        let zeroth = &self.members[0] ;
        let p0 = zeroth.0.to_matrix().clone() ;
        let coeff0 = zeroth.1 ;
        let mut bt = util::BinaryTreeFold::begin((coeff0, p0), |x,y| addmul(x,y)) ;

        for i in 1..self.members.len() {
            let ith = &self.members[i] ;
            let p = ith.0.to_matrix() ;
            let coeff = ith.1 ;
            bt.add((coeff, p)) ;
        }
        let mut p = bt.end() ;
        p.1.scale(p.0) ;
        p.1
    }

    pub fn to_matrix_accel(&self) -> sprs::CsMatI<Complex64, u64> {
        fn addmul(l: (Complex64, sprs::CsMatI<Complex64, u64>),
                  r: (Complex64, sprs::CsMatI<Complex64, u64>))
                  -> (Complex64, sprs::CsMatI<Complex64, u64>) {
            (Complex64::new(1.0,  0.0),
             sprs::binop::csmat_binop(l.1.view(), r.1.view(), |x,y| l.0 * x + r.0 * y))
        }

        let zeroth = &self.members[0] ;
        let p0 = zeroth.0.to_matrix_accel().clone() ;
        let coeff0 = zeroth.1 ;
        let mut bt = util::BinaryTreeFold::begin((coeff0, p0), |x,y| addmul(x,y)) ;

        for i in 1..self.members.len() {
            let ith = &self.members[i] ;
            let p = ith.0.to_matrix_accel() ;
            let coeff = ith.1 ;
            bt.add((coeff, p)) ;
        }
        let mut p = bt.end() ;
        p.1.scale(p.0) ;
        p.1
    }

    pub fn to_matrix_reduce(&self) -> sprs::CsMatI<Complex64, u64> {

        self.members.iter()
            .map(|(p,c)| (p.to_matrix(),*c))
            .reduce(|(a,acoeff),(b,bcoeff)| {
                (sprs::binop::csmat_binop(a.view(), b.view(), |x,y| acoeff * x + bcoeff * y),
                 Complex64::new(1.0, 0.0))
            })
            .unwrap().0
    }

    pub fn to_matrix_rayon(&self) -> sprs::CsMatI<Complex64, u64> {

        self.members
            .par_iter()
            .map(|(p,c)| {
                let mut m = p.to_matrix() ;
                m.scale(*c) ;
                Some(m)
            })
            .reduce(|| None,
                    |l,r| {
                        match (l,r) {
                            (None, None) => None,
                            (None,  Some(r)) => Some(r),
                            (Some(l), None) => Some(l),
                            (Some(l), Some(r)) => Some(&l + &r)
                        }
                    })
            .unwrap()
    }

    pub fn to_matrix_rayon_chunked(&self, step : usize) -> sprs::CsMatI<Complex64, u64> {

        let chunks : Vec<SparsePauliOp> =
            (0..self.members.len())
            .into_iter()
            .step_by(step)
            .map(|n| SparsePauliOp::from_slice(&self.members[n .. min(n + step, self.members.len())]))
            .collect();

        if chunks.len() == 1 {
            self.to_matrix_accel()
        }
        else {
            chunks
                .par_iter()
                .map(|p| {
                    let m = p.to_matrix_accel() ;
                    Some(m)
                })
                .reduce(|| None,
                        |l,r| {
                            match (l,r) {
                                (None, None) => None,
                                (None,  Some(r)) => Some(r),
                                (Some(l), None) => Some(l),
                                (Some(l), Some(r)) => Some(&l + &r)
                            }
                        })
                .unwrap()
        }
    }
    pub fn to_matrix_rowwise(&self) -> sprs::CsMatI<Complex64, u64> {
        let members = &self.members[..] ;
        let trimat = accel::rowwise::make_trimat(&members) ;
        trimat.to_csr()
    }
}

impl_op_ex!(+ |a: &SparsePauliOp, b: &SparsePauliOp| -> SparsePauliOp {
    let mut l = Vec::new();
    a.members.iter().for_each(|m| l.push(m.clone())) ;
    b.members.iter().for_each(|m| l.push(m.clone())) ;
    SparsePauliOp { members : l }
});

#[cfg(test)]
mod tests {
    use num_complex::Complex64;
    use ndarray::array ;
    use ndarray::Array2;
    use sprs::*;
    use ::approx::{ abs_diff_eq, assert_abs_diff_eq, abs_diff_ne, assert_abs_diff_ne, } ;

    use crate::* ;

    #[test]
    fn parse_labels() {
        assert!(Pauli::parse_label(&"W".to_string()).is_err());
        assert!(Pauli::parse_label(&"".to_string()).is_err());
        assert!(Pauli::parse_label(&"+i".to_string()).is_err());
        assert_eq!(Pauli::parse_label(&"I".to_string()), Ok( (0_usize, vec![I]) ));
        assert_eq!(Pauli::parse_label(&"+I".to_string()), Ok( (0_usize, vec![I]) ));
        assert_eq!(Pauli::parse_label(&"+iI".to_string()), Ok( (1_usize, vec![I]) ));
        assert_eq!(Pauli::parse_label(&"+jI".to_string()), Ok( (1_usize, vec![I]) ));
        assert_eq!(Pauli::parse_label(&"-1jI".to_string()), Ok( (3_usize, vec![I]) ));
        assert_eq!(Pauli::parse_label(&"-1I".to_string()), Ok( (2_usize, vec![I]) ));
        assert_eq!(Pauli::parse_label(&"-1jIX".to_string()), Ok( (3_usize, vec![X, I]) ));
        assert_eq!(Pauli::parse_label(&"IXYZ".to_string()), Ok( (0_usize, vec![Z, Y, X, I]) ));
        assert_eq!(Pauli::new("IXYZ").unwrap().label(), "IXYZ");

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
            let p = Pauli::new(&"IXYZ".to_string()) ;
            assert!(p.is_ok()) ;
            let p = p.unwrap() ;
            assert_eq!(p.xs(),
                       vec![false,  true,  true, false]);
            assert_eq!(p.zs(),
                       vec![true,  true, false, false]);
        }

/*
Pauli::parse_label("-IIIIIIIIIIIIIIIIIIYXXY"), Ok(2_usize, vec![ IIIIIIIIIIIIIIIIIIYXXY ]
z=[ True False False  True False False False False False False False False
 False False False False False False False False False False]
x=[ True  True  True  True False False False False False False False False
 False False False False False False False False False False]
phase=2

*/
        {
            let p = Pauli::new(&"-IIIIIIIIIIIIIIIIIIYXXY".to_string()) ;
            assert!(p.is_ok()) ;
            let p = p.unwrap() ;
            assert_eq!(p.base_phase(), 2_usize) ;
            assert_eq!(p.phase(), 0_usize) ;
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

    #[test]
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

    #[test]
    fn pauli_matrices() {
        let p = Pauli::new(&"IX".to_string()) ;
        assert!(p.is_ok()) ;
        let p = p.unwrap() ;
        assert_eq!(p.to_matrix().view().to_dense(),
                   kronecker_product(I.to_matrix().view(), X.to_matrix().view()).to_dense()) ;
        assert_eq!(p.to_matrix_rowwise().to_dense(), p.to_matrix().to_dense()) ;

        let p = Pauli::new(&"XI".to_string()) ;
        assert!(p.is_ok()) ;
        let p = p.unwrap() ;
        assert_eq!(p.to_matrix().view().to_dense(),
                   kronecker_product(X.to_matrix().view(), I.to_matrix().view()).to_dense()) ;
        assert_eq!(p.to_matrix_rowwise().to_dense(), p.to_matrix().to_dense()) ;

    }

    #[test]
    fn accel() {
        let p = Pauli::new(&"I".to_string()).unwrap() ;

        assert_eq!(p.to_matrix().view().to_dense(),
                   p.to_matrix_accel().view().to_dense()) ;
        assert_eq!(p.to_matrix_rowwise().to_dense(), p.to_matrix().to_dense()) ;

        let p = Pauli::new(&"Y".to_string()).unwrap() ;

        assert_eq!(p.to_matrix().view().to_dense(),
                   p.to_matrix_accel().view().to_dense()) ;
        assert_eq!(p.to_matrix_rowwise().to_dense(), p.to_matrix().to_dense()) ;

        let p = Pauli::new(&"YY".to_string()).unwrap() ;

        assert_eq!(p.to_matrix().view().to_dense(),
                   p.to_matrix_accel().view().to_dense()) ;
        assert_eq!(p.to_matrix_rowwise().to_dense(), p.to_matrix().to_dense()) ;
    }

    #[test]
    fn pauli_list() {
    }

    #[test]
    fn sparse_pauli_op() {
        let one = Complex64::new(1.0, 0.0) ;
        let two = 2.0 * one ;
        let _minus_one = Complex64::new(-1.0, 0.0) ;
        let _zero = Complex64::new(0.0, 0.0) ;
        let _i = Complex64::new(0.0, 1.0) ;
        let _minus_i = Complex64::new(0.0, -1.0) ;

        assert!(SparsePauliOp::from_labels_simple(&["I"][..]).is_ok()) ;
        assert!(SparsePauliOp::from_labels_simple(&["I", "II"][..]).is_err()) ;

        let spop = SparsePauliOp::from_labels(
            &["I","X"][..],
            &[Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)][..]) ;
        assert!(spop.is_ok()) ;
        let spop = spop.unwrap() ;
        assert_eq!(spop.to_matrix().to_dense(),
                   array![[one, two],
                          [two, one]]) ;
        assert_eq!(spop.to_matrix_binary().to_dense(),
                   array![[one, two],
                          [two, one]]) ;
        assert_eq!(spop.to_matrix_accel().to_dense(),
                   array![[one, two],
                          [two, one]]) ;
        assert_eq!(spop.to_matrix_reduce().to_dense(),
                   array![[one, two],
                          [two, one]]) ;
        assert_eq!(spop.to_matrix_rayon().to_dense(),
                   array![[one, two],
                          [two, one]]) ;
        assert_eq!(spop.to_matrix_rowwise().to_dense(),
                   spop.to_matrix().to_dense(),
        ) ;
    }

    #[test]
    fn h2() {
        let tc = &crate::fixtures::H2 ;
	let ll = &tc.labels[..] ;
	let cl = &tc.coeffs[..] ;
	let spop = SparsePauliOp::from_labels(ll, cl).unwrap() ;
        assert_abs_diff_eq!(spop.to_matrix_rayon(),
                     spop.to_matrix(),
                     epsilon = 1e-7
        ) ;
    }

    #[test]
    fn h4() {
        let tc = &crate::fixtures::H4 ;
	let ll = &tc.labels[..] ;
	let cl = &tc.coeffs[..] ;
	let spop = SparsePauliOp::from_labels(ll, cl).unwrap() ;
        assert_abs_diff_eq!(spop.to_matrix_rayon(),
                     spop.to_matrix(),
                     epsilon = 1e-7
        ) ;
    }

    #[test]
    fn h6() {
        let tc = &crate::fixtures::H6 ;
	let ll = &tc.labels[..] ;
	let cl = &tc.coeffs[..] ;
	let spop = SparsePauliOp::from_labels(ll, cl).unwrap() ;
        assert_abs_diff_eq!(spop.to_matrix_rayon(),
                     spop.to_matrix(),
                     epsilon = 1e-7
        ) ;
    }
}

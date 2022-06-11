#[macro_use]
extern crate lazy_static;

use num_complex::Complex64;
use regex::Regex;

pub fn greet(s : &str) {
    println!("Hello, {}!",  s) ;
}

#[derive(Debug, PartialEq)]
pub enum PauliKind {
    I, X, Y, Z,
}

use crate::PauliKind::* ;
impl PauliKind {
    pub fn new(c : char)-> Result<PauliKind, String> {
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
}

#[derive(Debug)]
pub struct Pauli {
    phase : usize,
    paulis : Vec<PauliKind>,
}
impl Pauli {

    fn parse_label(s : &str) -> Result<(usize, Vec<PauliKind>),  String> {
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
            paulis.push(PauliKind::new(c)?)
        }
        if 0 == paulis.len() { Err(String::from("internal error: no paulis")) }
        else { Ok((phase, paulis)) }
    }

    pub fn new(s : &str) -> Result<Pauli, String> {
        let (phase, paulis) = Pauli::parse_label(s)? ;
        Ok(Pauli{ phase, paulis })
    }

    pub fn phase(&self) -> usize { self.phase }
    pub fn xs(&self) -> Vec<bool> { self.paulis.iter().map(|p| p.x()).collect() }
    pub fn zs(&self) -> Vec<bool> { self.paulis.iter().map(|p| p.z()).collect() }

}

#[cfg(test)]
mod tests {
    use crate::Pauli ;
    use crate::PauliKind::* ;
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

}


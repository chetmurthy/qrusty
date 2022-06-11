#[macro_use]
extern crate lazy_static;

use num_complex::Complex64;
use regex::Regex;

pub fn greet(s : &str) {
    println!("Hello, {}!",  s) ;
}

pub enum PauliKind {
    I, X, Y, Z,
}

impl PauliKind {
    pub fn new(c : char)-> Result<PauliKind, String> {
        match c {
            'I' => Ok(PauliKind::I),
            'X' => Ok(PauliKind::X),
            'Y' => Ok(PauliKind::Y),
            'Z' =>  Ok(PauliKind::Z),
            _ => Err(String::from("internal error: malformed pauli"))
        }
    }
}

pub struct Pauli {
    coeff : Complex64,
    paulis : Vec<PauliKind>,
}
impl Pauli {    

    fn parse_label(s : &str) -> Result<(Complex64, Vec<PauliKind>),  String> {
        lazy_static! {
            static ref RE : Regex = Regex::new(r"^([+-]?)1?([ij]?)([IXYZ]+)$").unwrap();
        }
        let caps = RE.captures(s).unwrap();

        let sign = caps.get(1).map_or("", |m| m.as_str());
        let imag = caps.get(2).map_or("", |m| m.as_str());
        let paulistr = caps.get(3).map_or("", |m| m.as_str());

        let sign = match sign { ""|"+" => Ok(true), "-" => Ok(false), _ => Err("internal error: malformed sign") } ?;
        let coeff = match imag {
            "" => Ok(Complex64::new(1.0,  0.0)),
            "i"|"j" => Ok(Complex64::new(0.0, 1.0)),
            _ => Err("internal error: malformed imag")
        }?;
        let coeff = if sign { coeff } else { - coeff } ;
        
        let mut paulis = Vec::new() ;
        for c in paulistr.chars() {
            paulis.push(PauliKind::new(c)?)
        }
        if 0 == paulis.len() { Err(String::from("internal error: no paulis")) }
        else { Ok((coeff, paulis)) }
    }

    pub fn new(s : &str) -> Result<Pauli, String> {
        let (coeff, paulis) = Pauli::parse_label(s)? ;
        Ok(Pauli{ coeff, paulis })
    }

}

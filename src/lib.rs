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
    phase : usize,
    paulis : Vec<PauliKind>,
}
impl Pauli {    

    fn parse_label(s : &str) -> Result<(usize, Vec<PauliKind>),  String> {
        lazy_static! {
            static ref RE : Regex = Regex::new(r"^([+-]?)1?([ij]?)([IXYZ]+)$").unwrap();
        }
        let caps = RE.captures(s).unwrap();

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
        for c in paulistr.chars() {
            paulis.push(PauliKind::new(c)?)
        }
        if 0 == paulis.len() { Err(String::from("internal error: no paulis")) }
        else { Ok((phase, paulis)) }
    }

    pub fn new(s : &str) -> Result<Pauli, String> {
        let (phase, paulis) = Pauli::parse_label(s)? ;
        Ok(Pauli{ phase, paulis })
    }

}

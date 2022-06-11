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
            _ => Err(String::from("internal error: malformed paulis"))
        }
    }
}

pub struct PauliData {
    coeff : Complex64,
    paulis : Vec<PauliKind>,
}
impl PauliData {    

pub fn parse_label(s : &str) -> PauliData {
    lazy_static! {
        static ref RE : Regex = Regex::new(r"^([+-]?)1?([ij]?)([IXYZ]+)$").unwrap();
    }
    let caps = RE.captures(s).unwrap();

    let sign = caps.get(1).map_or("", |m| m.as_str());
    let imag = caps.get(2).map_or("", |m| m.as_str());
    let paulis = caps.get(3).map_or("", |m| m.as_str());

    let sign = match sign { ""|"+" => true, "-" => false, _ => panic!("internal error: malformed sign") } ;
    let coeff = match imag {
        "" => Complex64::new(1.0,  0.0),
        "i"|"j" => Complex64::new(0.0, 1.0),
        _ => panic!("internal error: malformed imag")
    };
    let coeff = if sign { coeff } else { - coeff } ;
            
    let mut paulis = Vec::new() ;
    for c in s.chars() {
        paulis.push(match PauliKind::new(c) {
            Ok(v) => v,
            _ => panic!("internal error: malformed paulis")
        })
    }
    PauliData{ coeff, paulis }
}
}

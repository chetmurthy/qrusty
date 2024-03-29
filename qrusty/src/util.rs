// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

#![allow(non_upper_case_globals)]
use regex::Regex;
use num_complex::*;
use std::fs::File;
use std::io;
use std::io::Write;
use sprs::SpIndex;

static DIGIT : &str = "[0-9]" ;
static  NZDIGIT : &str = "[1-9]" ;

fn unamp(l : Vec<&str>) -> Vec<String> {
    l.iter().map(|s| String::from(*s)).collect()
}

fn conc(l : Vec<&str>) -> String {
    let rv = unamp(l).join("") ;
    format!("(?:{})", rv)
}

fn disj(l : Vec<&str>) -> String {
    let rv = unamp(l).join("|") ;
    format!("(?:{})", rv)
}

fn top(s : &str) -> String {
    format!("^{}$", s)
}
fn capt(s : &str) -> String {
    format!("({})", s)
}
fn star(s : &str) -> String {
    format!("(?:{})*", s)
}
fn plus(s : &str) -> String {
    format!("(?:{})+", s)
}
fn opt(s : &str) -> String {
    format!("(?:{})?", s)
}

lazy_static! {
    /*
    let int = [%sedlex.regexp? '0' | ( ('1'..'9') , (Star digit) )]
    let frac = [%sedlex.regexp? '.' , (Star digit)]
    let ne_frac = [%sedlex.regexp? '.' , (Plus digit)]
    let exp = [%sedlex.regexp? ('e' | 'E') , (Opt ('-' | '+')) , (Plus digit)]
    let decimal_float_number = [%sedlex.regexp? (Opt '-') , ((int , (Opt frac) , (Opt exp)) | (ne_frac, Opt exp))]
    let json_number = [%sedlex.regexp? (Opt '-') , int, Opt ne_frac, Opt exp]
     */

    static ref UnsignedFloat : String = {
        let digit = DIGIT ;
        let nzdigit = NZDIGIT ;
        let int = disj(vec!["0", &conc(vec![&nzdigit, &star(&digit)])]) ;
        let frac = conc(vec!["\\.", &star(&digit)]) ;
        let ne_frac = conc(vec!["\\.", &plus(&digit)]) ;
        let exp = conc(vec!["[eE]", &opt("[-+]"), &plus(&digit)]) ;
        let unsigned_float = disj(vec![&conc(vec![&int , &opt(&frac) , &opt(&exp)]),
                                       &conc(vec![&ne_frac, &opt(&exp)])]) ;
        unsigned_float
    } ;

    static ref SignedFloat : String = {
        let spaces = star(" ") ;
        let signed_float = conc(vec![ &capt(&opt("[-+]")), &spaces, &capt(&UnsignedFloat) ]) ;
        signed_float
    } ;



    static ref Complex : String = {
        let spaces = star(" ") ;
        /* {SFLOAT} (({SFLOAT}J)? | J)? */
        let complex = conc(vec![&spaces,
                                &opt(&SignedFloat),
                                &spaces,
                                &opt(&disj(vec![
                                    &conc(vec![&SignedFloat, "j"]),
                                    "j"
                                ]))]) ;
        complex
    } ;


    static ref FloatRE : Regex = {
        Regex::new(&top(&SignedFloat)).unwrap()
    } ;
    static ref ComplexRE : Regex = {
        Regex::new(&top(&Complex)).unwrap()
    } ;
}

pub fn float64_from_string(s : &str) -> Result<f64,String> {
    let caps = FloatRE.captures(s).ok_or(format!("error: malformed float: {}", s))? ;
    
    let signop : &str = caps.get(1).map_or("+", |m| m.as_str());
    let s : &str = caps.get(2).map_or("", |m| m.as_str());

    let n : f64 = s.parse().map_err(|_| format!("error: malformed float: {}", s)) ? ;
    Ok(if signop != "-" { n } else { - n })
}

pub fn parse_complex64<'a>(s : &'a str) -> Result<(bool, &'a str, bool, &'a str), String> {
    let caps = ComplexRE.captures(s).ok_or(format!("error: malformed complex: {}", s))? ;
    
    let signop1 : &str = caps.get(1).map_or("+", |m| m.as_str());
    let sreal = caps.get(2).map_or("", |m| m.as_str());
    let signop2 : &str = caps.get(3).map_or("+", |m| m.as_str());
    let simag = caps.get(4).map_or("0", |m| m.as_str());
    Ok(((signop1 != "-"), sreal, (signop2 != "-"), simag))
}

pub fn complex64_from_string(s : &str) -> Result<Complex64,String> {
    let (resign, res,  imsign, ims) = parse_complex64(s) ? ;
    let real : f64 = float64_from_string(res) ? ;
    let imag : f64 = float64_from_string(ims) ? ;
    Ok(Complex64::new(if resign { real } else { - real },
		      if imsign { imag } else { - imag }))
}

pub fn complex64_list_from_string_list(l : &Vec<&'static str>) -> Result<Vec<Complex64>, String> {
    let l: Result<Vec<Complex64>, _> = l
        .iter()
        .map(|s| complex64_from_string(s))
        .collect::<Result<Vec<Complex64>, _>>() ;
    l
}

pub fn complex64_to_string(c : Complex64) -> String {
    format!("{}+{}j", c.re, c.im)
}

pub fn csmatrix_nz(it : &sprs::CsMatI<Complex64,  u64>, tolerance : f64) -> usize {
    it.view()
        .iter_rbr()
        .filter(|(c,_)| {
            let c : Complex64 = **c ;
            c.norm() <= tolerance
        })
        .count()
}

pub fn csmatrix_eliminate_zeroes(it : &sprs::CsMatI<Complex64,  u64>, tolerance : f64) -> sprs::CsMatI<Complex64,  u64> {
    let mut rows = Vec::new() ;
    let mut cols = Vec::new() ;
    let mut data = Vec::new() ;
    it.view()
        .iter_rbr()
        .filter(|(c,_)| {
            let c : Complex64 = **c ;
            c.norm() > tolerance
        })
        .for_each(|(c,(row,col))| {
            rows.push(row) ;
            cols.push(col) ;
            data.push(*c) ;
        }) ;
    let trimat = sprs::TriMatBase::<Vec<u64>, Vec<Complex64>>::from_triplets(it.shape(), rows, cols, data) ;
    trimat.to_csr()
}
pub fn write_matrix_market_for_scipy(
    path: &str,
    mat: &sprs::CsMatI<Complex64,  u64>,
) -> Result<(), io::Error>
{
    let (rows, cols, nnz) = (mat.rows(), mat.cols(), mat.nnz());
    let f = File::create(path)?;
    let mut writer = io::BufWriter::new(f);

    // header
    let data_type = "complex";

    writeln!(
        writer,
        "%%MatrixMarket matrix coordinate {} general",
        data_type
    )?;
    writeln!(writer, "% written by sprs")?;

    // dimensions and nnz
    writeln!(writer, "{} {} {}", rows, cols, nnz)?;

    // entries
    for (val, (row, col)) in mat {
        writeln!(writer, "{} {} {} {}", row.index() + 1, col.index() + 1, val.re, val.im)?;
    }
    Ok(())
}

pub struct BinaryTreeFold<T> {
    stk : Vec<(usize, T)>,
    op : fn(l : T, r: T) -> T
}

/* a BinaryTreeFold is a stack that keeps track of where we are in the
 * binary-tree-based folding process.
 *
 * The operations are:
 *
 * begin: given a value of type T and a binary operator over T, create
 * the BTF
 *
 * add: add a new T to the BTF, invoking the binary operator, perhaps
 * more than once
 *
 * end: extract the final value, again perhaps invoking the binary
 * operator more than once
 */

impl<T> BinaryTreeFold<T> {
    pub fn begin(init : T, op : fn(l: T,r: T) -> T) -> BinaryTreeFold<T> {
        let stk = vec![ (1, init) ] ;
        BinaryTreeFold { stk, op }
    }
/*
            def addp(self, p):
                if not self.stack:
                    self.stack.append(p)
                    return
                q = self.stack[-1]
                if q[0] == p[0]:
                    height = 1 + q[0]
                    newv = self.f(height, q[1], p[1])
                    self.stack.pop()
                    self.addp((height, newv))
                else:
                    self.stack.append(p)

            def add(self, v):
                self.addp((1,v))

            def flush(self):
                if not self.stack:
                    raise Exception("flush: stack was empty")
                v = self.stack.pop()[1]
                while self.stack:
                    q = self.stack.pop()
                    v = self.f(q[0], q[1], v)
                return v

*/
    /* invariant: the stack is always nonempty */
    fn addp(&mut self, v : (usize, T)) {
        if self.stk.len() == 0 {
            self.stk.push(v) ;
        }
        else {
            let q = self.stk.last().unwrap() ;
            if q.0 == v.0 {
                let q = self.stk.pop().unwrap() ;
                let newv = (self.op)(q.1, v.1) ;
                self.addp((q.0 + 1, newv)) ;
            } else {
                self.stk.push(v) ;
            }
        }
    }

    pub fn add(&mut self, v : T) {
        self.addp((1,v))
    }

    pub fn end(&mut self) -> T {
        let mut v = self.stk.pop().unwrap().1 ;
        while self.stk.len() > 0 {
            let w = self.stk.pop().unwrap() ;
            v = (self.op)(v, w.1) ;
        }
        drop(self) ;
        v
    }
}

pub mod list {
    use std::fmt;

    #[derive(Debug)]
    pub enum ListNode<T> {
        Cons(List<T>, List<T>),
        Atom(T)
    }
    type List<T> = Box<ListNode<T>> ;


    pub fn cons<T>(l : List<T>, r : List<T>) -> List<T> {
        Box::new(ListNode::Cons(l,r))
    }

    pub fn atom<T>(a : T) -> List<T> {
        Box::new(ListNode::Atom(a))
    }

    impl<T: fmt::Display> fmt::Display for ListNode<T> {
        // This trait requires `fmt` with this exact signature.
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            match self {
                ListNode::Atom(t) => write!(f, "{}", t),
                ListNode::Cons(l,r) =>
                    write!(f, "({} {})", l,  r)
            }
        }
    }
}

pub mod fileio {
    use flate2::bufread::GzDecoder;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io;
    use std::fs::File;
    use std::path::Path;

    pub fn with_input_file<R, E>(p : &Path,
                           f : &dyn Fn(&mut dyn io::BufRead) -> Result<R, E>,
    ) -> Result<R, E>
    where E : From<io::Error>
    {
        match p.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => {
                    let fp = File::open(p)?;
                    let buf = io::BufReader::new(fp);
                    let reader = GzDecoder::new(buf);
                    let mut buf = io::BufReader::new(reader);
                    let rv = f(&mut buf) ;
                    return rv ;
            },
            _ => {
                    let fp = File::open(p)?;
                    let mut reader = io::BufReader::new(fp);
                    let rv = f(&mut reader) ;
                    return rv ;
            }
        }
    }

    pub fn with_output_file<R, E>(p : &Path,
                                  f : &dyn Fn(&mut dyn io::Write) -> Result<R, E>,
    ) -> Result<R, E>
    where E : From<io::Error>
    {
        match p.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => {
                let fp = File::create(p)?;
                let e = GzEncoder::new(fp, Compression::default());
                let mut writer = io::BufWriter::new(e);
                let rv = f(&mut writer) ;
                return rv ;
            },
            _ => {
                let fp = File::create(p)?;
                let mut writer = io::BufWriter::new(fp);
                let rv = f(&mut writer) ;
                return rv ;
            }
        }
    }

}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;
    use crate::util::BinaryTreeFold ;
    use crate::util::csmatrix_nz;
    use crate::util::csmatrix_eliminate_zeroes ;
    use crate::util::list::* ;

    #[test]
    fn plus() {
        let mut bt = BinaryTreeFold::begin(1, |x,y| x+y) ;
        (2..10).for_each(|n| bt.add(n)) ;

        assert_eq!((1..10).fold(0,|a,b| a+b), bt.end()) ;
    }

    #[test]
    fn test_cons() {
        let mut bt = BinaryTreeFold::begin(atom(0), cons) ;
        (1..8).for_each(|n| bt.add(atom(n))) ;
        let rv = bt.end() ;
        assert_eq!("(((0 1) (2 3)) ((4 5) (6 7)))",
                   (*rv).to_string()) ;
    }

    #[test]
    fn test_parse_f64() {
        assert_eq!(crate::util::float64_from_string("0"), Ok(0.0)) ;
        assert_eq!(crate::util::float64_from_string("-0"), Ok(-0.0)) ;
        assert_eq!(crate::util::float64_from_string("+0"), Ok(0.0)) ;
        assert_eq!(crate::util::float64_from_string("-1.0"), Ok(-1.0)) ;
        assert_eq!(crate::util::float64_from_string("+1.0"), Ok(1.0)) ;
        assert_eq!(crate::util::float64_from_string("1.0"), Ok(1.0)) ;
        assert!(crate::util::float64_from_string("foo").is_err()) ;
    }

    #[test]
    fn test_parse_complex64() {
        assert_eq!(crate::util::complex64_from_string("0"), Ok(Complex64::new(0.0,  0.0))) ;
        assert_eq!(crate::util::parse_complex64("1.0"), Ok((true, "1.0", true,  "0"))) ;
        assert_eq!(crate::util::complex64_from_string("1.0"), Ok(Complex64::new(1.0,  0.0))) ;
        assert!(crate::util::complex64_from_string("foo").is_err()) ;
        assert_eq!(crate::util::parse_complex64("1+2j"), Ok((true, "1", true,  "2"))) ;
        assert_eq!(crate::util::complex64_from_string("1+2j"), Ok(Complex64::new(1.0,  2.0))) ;
        assert_eq!(crate::util::complex64_from_string("1-2j"), Ok(Complex64::new(1.0,  -2.0))) ;
        assert_eq!(crate::util::complex64_from_string("-1.56404847e01 + 0.0j"), Ok(Complex64::new(-1.56404847e01,  0.0))) ;
    }

    #[test]
    fn test_print_complex64() {
        fn roundtrip(c : Complex64) -> Complex64 {
            crate::util::complex64_from_string(&crate::util::complex64_to_string(c)).unwrap()
        }
        assert_eq!(roundtrip(Complex64::new(0.0,  0.0)), Complex64::new(0.0,  0.0)) ;
    }

    #[test]
    fn test_csmatrix_nz() {
        let eye_100 = sprs::CsMatI::<Complex64, u64>::eye(100) ;
        assert_eq!(csmatrix_nz(&eye_100, 1e-7), 0) ;
        let mut scaled = sprs::CsMatI::<Complex64, u64>::eye(100) ;
        scaled.scale(Complex64::new(1e-8, 0.0)) ;
        assert_eq!(csmatrix_nz(&scaled, 1e-7), 100) ;
        let no_zeros = csmatrix_eliminate_zeroes(&scaled, 1e-7) ;
        assert_eq!(no_zeros.nnz(), 0) ;
    }
}

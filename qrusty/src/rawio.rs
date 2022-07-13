// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use std::io ;
use std::io::{ Read, Write } ;
use num_complex::Complex64;
use cpu_endian::{Endian, working};
use num_traits::Zero;
use sprs::{CompressedStorage, CsMatI};

pub trait Swab
where
    Self: Sized,
{
    fn swab(&mut self) ;
}

impl Swab for u64 {
    fn swab(&mut self) {
        *self = self.swap_bytes() ;
    }
}

impl Swab for f64 {
    fn swab(&mut self) {
        *self = f64::from_bits(self.to_bits().swap_bytes()) ;
    }
}

impl Swab for Complex64 {
    fn swab(&mut self) {
        self.re.swab() ;
        self.im.swab() ;
    }
}

impl<T : Swab> Swab for Vec<T> {
    fn swab(&mut self) {
        for i in 0..self.len() {
            self[i].swab() ;
        }
    }
}

fn read_endianmark<R : io::Read>(r : &mut R) -> Result<u16, io::Error> {
    let mut buffer : [u8; 2] = [0u8; 2];
    r.read_exact(&mut buffer)?;
    let readmark = u16::from_ne_bytes(buffer) ;
    Ok(readmark)
}

fn mark() -> u16 {
    let buffer = ['M' as u8, 'I' as u8] ;
    u16::from_ne_bytes(buffer)
}

fn write_endianmark<W : io::Write>(w : &mut W) -> Result<(), io::Error> {
    let mark = mark() ;
    let buffer = mark.to_ne_bytes() ;
    w.write_all(&buffer[..])
}

fn read_u64<R : io::Read>(r : &mut R, need_swab : bool) -> Result<u64, io::Error> {
    let mut buffer = [0u8; std::mem::size_of::<u64>()];
    r.read_exact(&mut buffer)?;
    let mut n = match working() {
        Endian::Little => u64::from_le_bytes(buffer),
        Endian::Big => u64::from_be_bytes(buffer),
        _ => panic!("rawio::read_u64: Minor endian detected, completely unsupported")
    } ;
    if need_swab { n.swab() ; }
    Ok(n)
}

fn write_u64<W : io::Write>(w : &mut W, v : u64) -> Result<(), io::Error> {
    let buffer = v.to_ne_bytes() ;
    w.write_all(&buffer[..])
}

fn fread<R: Read, T>(
    reader: &mut R,
    length: usize, // estimated via a seek or other mechanism
    need_swab: bool,
) -> Result<Vec<T>,  io::Error>
    where T : Sized + Zero + Clone + Swab
{
    // it is undefined behavior to call read_exact on un-initialized, https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read
    // see also https://github.com/MaikKlein/ash/issues/354#issue-781730580
    let mut buffer = vec![T::zero(); length];

    unsafe {
        // transmute u64 to bytes.
        let slice = std::slice::from_raw_parts_mut(
            buffer.as_mut_ptr() as *mut u8,
            length * std::mem::size_of::<T>(),
        );
        reader.read_exact(slice)?;
    }

    if need_swab { buffer.swab() ; }
    Ok(buffer)
}

fn fwrite<W: Write, T>(
    w: &mut W,
    v : &[T],
) -> Result<(),  io::Error>
    where T : Sized + Zero + Clone + Swab,
{
    unsafe {
        // transmute u64 to bytes.
        let slice = std::slice::from_raw_parts(
            v.as_ptr() as *const u8,
            v.len() * std::mem::size_of::<T>(),
        );
        w.write_all(slice)?;
    }
    Ok(())
}

pub fn write<W : io::Write>(w : &mut W, m : &sprs::CsMatI<Complex64,  u64>) -> Result<(), io::Error> {
    write_endianmark(w)?;
    let (rows, cols) = m.shape() ;
    match m.storage() {
        CompressedStorage::CSR => write_u64(w, 0x0),
        CompressedStorage::CSC => write_u64(w, 0x1),
    }? ;
    write_u64(w, rows as u64)? ;
    write_u64(w, cols as u64)? ;
    let indptr_view = m.indptr() ;
    let indptr = indptr_view.raw_storage() ;
    write_u64(w, indptr.len() as u64)? ;
    fwrite(w, indptr)? ;
    let indices = m.indices() ;
    write_u64(w, indices.len() as u64)? ;
    fwrite(w, indices)? ;
    let data = m.data() ;
    write_u64(w, data.len() as u64)? ;
    fwrite(w, data)
}

pub fn read<R : io::Read>(r : &mut R) -> Result<sprs::CsMatI<Complex64,  u64>, io::Error> {
    let read_mark = read_endianmark(r)? ;
    let need_swab = mark() != read_mark ;
    let n = read_u64(r, need_swab)?;
    let storage = match n {
        0x0 => CompressedStorage::CSR,
        0x1 => CompressedStorage::CSC,
        _ => panic!("rawio::read: bad storage found")
    } ;
    let rows = read_u64(r, need_swab)?;
    let cols = read_u64(r, need_swab)?;

    let indptr_len = read_u64(r, need_swab)?;
    let indptr = fread(r, indptr_len as usize, need_swab)?;

    let indices_len = read_u64(r, need_swab)?;
    let indices = fread(r, indices_len as usize, need_swab)?;

    let data_len = read_u64(r, need_swab)?;
    let data = fread(r, data_len as usize, need_swab)?;

    unsafe {
        Ok(CsMatI::<Complex64, u64, u64>::new_unchecked(
            storage,
            (rows as usize, cols as usize),
            indptr,
            indices,
            data,
        ))
    }
}

#[cfg(test)]
mod test {
    use std::io::Cursor ;
    use num_complex::Complex64;
    use crate::rawio::{ Swab, write, read };
    use sprs::{ TriMatI , CsMatI };
    const U64_CONSTANT : u64 = 0x0123456789abcdef ;
    const SWABBED_U64_CONSTANT : u64 = 0xefcdab8967452301 ;
    const F64_CONSTANT : f64 = 3.14159 ;

    fn u64_hex(n : u64) -> String {
        format!("{:X}",n)
    }

    #[test]
    fn prim_u64() {
        let mut n : u64 = U64_CONSTANT ;
        n.swab() ;
        assert_eq!(u64_hex(SWABBED_U64_CONSTANT), u64_hex(n)) ;
    }

    #[test]
    fn prim_f64() {
        let mut n : f64 = F64_CONSTANT ;
        n.swab() ;
        n.swab() ;
        assert_eq!(F64_CONSTANT, n) ;
    }

    #[test]
    fn vec_u64() {
        let mut v : Vec<u64> = vec![1,2,3,4] ;

        v.swab() ;
        assert_eq!(vec![1<<56,2<<56,3<<56,4<<56].into_iter().map(u64_hex).collect::<Vec<String>>(),
                   v.into_iter().map(u64_hex).collect::<Vec<String>>()) ;
    }

    #[test]
    fn write_matrix() {
        let mut trimat : TriMatI<Complex64, u64> = TriMatI::new((2,2)) ;
        trimat.add_triplet(0,0,Complex64::new(1.0,0.0)) ;
        trimat.add_triplet(1,1,Complex64::new(1.0,0.0)) ;
        let csr : CsMatI<Complex64, u64> = trimat.to_csr() ;
        let mut wfile = Vec::<u8>::new() ;
        assert!(write(&mut wfile, &csr).is_ok()) ;

        
        let mut rfile = Cursor::new(wfile);
        let csr2 = read(&mut rfile).unwrap() ;
        assert_eq!(csr.to_dense(), csr2.to_dense()) ;
    }

}

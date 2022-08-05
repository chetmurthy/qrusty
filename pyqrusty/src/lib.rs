// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

#![allow(unused_imports)]
use std::fs::File;
use std::io;
use std::path::Path;
use pyo3::prelude::*;
use ndarray::{ArrayViewD, ArrayD, ArrayView, Dim, Array, Array1, par_azip};
use numpy::{IntoPyArray, PyReadonlyArray1, PyArrayDyn, PyReadonlyArrayDyn, PyArray1};
use num_complex::Complex64;
use pyo3::wrap_pyfunction;
use pyo3::Python;
use pyo3::PyErr ;
use pyo3::types::{ PySlice } ;
use pyo3::exceptions::PyOSError;
use pyo3::exceptions::PyException;

use sprs::{CompressedStorage, CsMatI, TriMatI};
use sprs::io::{ read_matrix_market_from_bufread, write_matrix_market_to_bufwrite, SymmetryMode::* } ;


use qrusty::util::fileio::{ with_output_file, with_input_file } ;
use qrusty::AccelMode;

#[pyclass]
#[repr(transparent)]
pub struct SpMat {
    pub it: Box< Option< sprs::CsMatI<Complex64,  u64> > >,
}

impl SpMat {
    fn binop(&self, other: &SpMat, name: &str, f : fn (&sprs::CsMatI<Complex64,  u64>, &sprs::CsMatI<Complex64,  u64>) -> sprs::CsMatI<Complex64,  u64>) -> PyResult<SpMat> {

        match (&*(self.it), &*(other.it)) {
            (Some(lhs), Some(rhs)) => {
                (lhs.shape() == rhs.shape()).then(|| ())
                    .ok_or(PyException::new_err(format!("cannot {} matrices with different shapes", name)))? ;
                let mat : sprs::CsMatI<Complex64,  u64> = f(lhs, rhs) ;
                Ok(SpMat::new_from_csmatrix(mat))
            },
            _ => Err(PyException::new_err(format!("cannot {} already-exported sparse matrices", name)))
        }
    }

    #[allow(dead_code)]
    fn unop(&self, name: &str, f : fn (&sprs::CsMatI<Complex64,  u64>) -> sprs::CsMatI<Complex64,  u64>) -> PyResult<SpMat> {

        match &*(self.it) {
            Some(lhs) => {
                let mat : sprs::CsMatI<Complex64,  u64> = f(lhs) ;
                Ok(SpMat { it : Box::new(Option::Some(mat)) })
            },
            _ => Err(PyException::new_err(format!("cannot {} already-exported sparse matrix", name)))
        }
    }

    fn new_from_csmatrix(sp_mat : sprs::CsMatI<Complex64,  u64>) -> SpMat {
        SpMat { it : Box::new(Option::Some(sp_mat)) }
    }

    fn map_immut<R,OKF, ERRF>(&self, errf : ERRF, mut okf : OKF) -> R
    where OKF: FnMut (&CsMatI::<Complex64, u64, u64>) -> R,
    ERRF : Fn () -> R
    {
        match &*(self.it) {
            None => errf(),
            Some(spmat) => okf(spmat)
        }
    }
    fn iter_mut<E, OKF, ERRF>(&mut self, errf : ERRF, okf : OKF) -> Result<(), E>
    where OKF : Fn (&mut CsMatI::<Complex64, u64, u64>) -> (),
    ERRF : Fn () -> E
    {
        match &mut *(self.it) {
            None => Err(errf()),
            Some(spmat) => Ok(okf(spmat))
        }
    }

    fn parse_symmetry_mode(mode : &str) -> PyResult<sprs::io::SymmetryMode> {
        match mode {
            ""|"general" => Ok(General),
            "hermitian" => Ok(Hermitian),
            "symmetric" => Ok(Symmetric),
            "skew-symmetric" => Ok(SkewSymmetric),
            _ => Err(PyException::new_err(format!("matrixmarket_write_symmetric: unrecognized mode {}", mode)))
        }
    }

}

#[pymethods]
impl SpMat {

    #[staticmethod]
    pub fn new_unchecked(shape : (usize, usize), data : Vec<Complex64>, indices: Vec<u64>, indptr: Vec<u64>) -> PyResult<SpMat> {
        let sp_mat = 
            unsafe {
                CsMatI::<Complex64, u64, u64>::new_unchecked(
                    CompressedStorage::CSR,
                    shape,
                    indptr,
                    indices,
                    data,
                )
            } ;
        Ok(SpMat { it : Box::new(Option::Some(sp_mat)) })
    }

    fn diagonal<'py>(
	&self,
        py: Python<'py>,
    ) -> PyResult<PyObject> {
	let y = self.map_immut(|| Err(PyException::new_err("SpMat.diagonal(): already-exported sparse matrix")),
			|spmat| Ok(spmat.diag().to_dense())) ? ;
	Ok(y.into_pyarray(py).into())
    }

    fn shape(&self) -> PyResult<(usize, usize)> {
        self.map_immut(|| Err(PyException::new_err("SpMat.shape(): matrix is already dropped")),
                        |spmat| { Ok(spmat.shape()) })
    }

    fn __copy__(&self) -> SpMat {
        SpMat { it : Box::new((&*self.it).clone()) }
    }

    fn __repr__(&self) -> String {
        self.map_immut(|| "<already-dropped sparse matrix of type Complex64>".to_string(),
                        |spmat| {
                            let (rows, cols) = spmat.shape() ;
                            let nnz = spmat.nnz() ;
                            format!("<{}x{} sparse matrix of type Complex64\n\twith {} stored elements in Compressed Sparse Row format>",
                                    rows,  cols, nnz)
                        })
    }
    
    fn __str__(&self) -> String {
        self.__repr__()
    }
    
    fn __add__(&self, other: &SpMat) -> PyResult<SpMat> {
        self.binop(other, "add",  |x,y| x+y)
    }
    
    fn __sub__(&self, other: &SpMat) -> PyResult<SpMat> {
        self.binop(other, "subtract",  |x,y| x-y)
    }

    fn scale(&mut self, factor: Complex64) -> PyResult<()> {
        self.iter_mut(|| PyException::new_err(format!("cannot scale already-exported sparse matrix")),
                      |lhs| {
                          lhs.scale(factor) ;
                      })
    }

/*    
    fn __neg__(&self) -> PyResult<SpMat> {
        self.binop("negate",  |y: &CsMatI::<Complex64, u64, u64>| -y)
    }
  */  
    #[args(tolerance = "1e-7")]
    fn count_zeros(&self,  tolerance: f64) -> PyResult<usize> {
        self.map_immut(|| Err(PyException::new_err("cannot count zeroes of an exported sparse matrix")),
                       |spmat| Ok(qrusty::util::csmatrix_nz(&spmat, tolerance)))
    }
    
    #[args(tolerance = "1e-7")]
    fn eliminate_zeros(&self,  tolerance: f64) -> PyResult<SpMat> {
        self.map_immut(|| Err(PyException::new_err("cannot eliminate zeroes of an exported sparse matrix")),
                       |spmat| {
                           let mat = qrusty::util::csmatrix_eliminate_zeroes(&spmat, tolerance) ;
                           Ok(SpMat { it : Box::new(Option::Some(mat)) })
                       })
    }

    fn nnz(&self) -> PyResult<usize> {
        self.map_immut(|| Err(PyException::new_err("cannot get NNZ of an exported sparse matrix")),
                       |spmat| Ok(spmat.nnz()))
    }

    pub fn export(
        &mut self,
        py: Python,
    ) -> PyResult<((usize,  usize), PyObject, PyObject, PyObject)> {

        match &*(self.it) {
            None =>
                Err(PyException::new_err("cannot export from an already-exported sparse matrix")),
            Some(_) => {
                match std::mem::take(&mut *(self.it)) {
                    None => panic!("internal error in SpMat::export"),
                    Some(spmat) => {
                        let shape = spmat.shape() ;
                        let (indptr, indices, data) = spmat.into_raw_storage();
                        Ok((
                            shape,
	                    data.into_pyarray(py).into(),
	                    indices.into_pyarray(py).into(),
	                    indptr.into_pyarray(py).into(),
                        ))
                    }
                }
            }
        }
    }

    pub fn matrixmarket_write(
        &mut self,
        save_path : &str,
    ) -> PyResult<()> {
        let path = Path::new(save_path) ;

        with_output_file(&path,
                         &|writer|{
                             self.map_immut(|| Err(PyException::new_err("cannot write an already-destroyed sparse matrix")),
                                            |spmat|
                                            sprs::io::write_matrix_market_to_bufwrite(writer, spmat)
                                            .map_err(|e| PyOSError::new_err(e)))
                         })
    }

    #[staticmethod]
    pub fn matrixmarket_read(
        read_path : &str,
    ) -> PyResult<SpMat> {
        let path = Path::new(read_path) ;

        let trimat : TriMatI<Complex64, u64> =
            with_input_file(path,
                            &|reader : &mut dyn io::BufRead|{
                                sprs::io::read_matrix_market_from_bufread(reader)
                                    .map_err(|e| PyException::new_err(format!("matrixmarket_read: {}", e)))
                            })?;

        let spmat = trimat.to_csr() ;
        Ok(SpMat::new_from_csmatrix(spmat))
    }

    #[args(mode = "\"general\"")]
    pub fn matrixmarket_write_symmetric(
        &mut self,
        save_path : &str,
        mode : &str,
    ) -> PyResult<()> {
        let mode = SpMat::parse_symmetry_mode(mode) ? ;
        self.map_immut(|| Err(PyException::new_err("cannot write an already-destroyed sparse matrix")),
                       |spmat|
                       sprs::io::write_matrix_market_sym(&save_path, spmat, mode)
                       .map_err(|e| PyOSError::new_err(e)))
    }
}

#[pyclass]
#[repr(transparent)]
#[derive(Clone)]
pub struct Pauli {
    pub it: qrusty::Pauli,
}

impl From<qrusty::Pauli> for Pauli {
    fn from(it: qrusty::Pauli) -> Self {
        Pauli { it : it.clone() }
    }
}

#[pymethods]
impl Pauli {

    #[new]
    pub fn new(label: String) -> PyResult<Pauli> {
        let it = qrusty::Pauli::new(&label) ? ;
        Ok(Pauli::from(it))
    }

    pub fn num_qubits(&self) -> usize {
        self.it.num_qubits()
    }

    pub fn label(&self) -> String {
        self.it.label()
    }

   fn __repr__(&self) -> String {
       format!("Pauli('{}')", self.it.label())
   }

   fn __str__(&self) -> String {
       self.it.label()
   }

    fn to_matrix(&self) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix())
    }
}

#[pyclass]
#[repr(transparent)]
#[derive(Clone)]
pub struct SparsePauliOp {
    pub it: qrusty::SparsePauliOp,
}

impl From<qrusty::SparsePauliOp> for SparsePauliOp {
    fn from(it: qrusty::SparsePauliOp) -> Self {
        SparsePauliOp { it : it.clone() }
    }
}

#[derive(FromPyObject)]
enum SliceOrInt<'a> {
    Slice(&'a PySlice),
    Int(isize),
}


enum SliceResult<T> {
    It(T),
    Slice(Vec<T>)
}
impl<T : IntoPy<PyObject>> IntoPy<PyObject> for SliceResult<T> {
    fn into_py(self, py: Python<'_>) -> PyObject {
        match self {
            SliceResult::It(it) => it.into_py(py),
            SliceResult::Slice(v) => {
                v.into_py(py)
            }
        }
    }
}

#[pymethods]
impl SparsePauliOp {

    #[new]
    pub fn new(paulis : Vec<Pauli>, coeffs : Vec<Complex64>) -> PyResult<SparsePauliOp> {
        let paulis = paulis.iter().map(|p| p.it.clone()).collect() ;
        let coeffs : &[Complex64] = &coeffs[..] ;
        let it = qrusty::SparsePauliOp::new(paulis, coeffs) ? ;
        Ok(SparsePauliOp { it })
    }

    pub fn __len__(&self) -> PyResult<usize> {
        Ok(self.it.members().len())
    }

    fn __getitem__(&self, idx: SliceOrInt) -> PyResult<SliceResult<(Pauli, Complex64)>> {
        match idx {
            SliceOrInt::Slice(slice) => {
                let psi = slice.indices(self.it.members().len() as i64)? ;
                let (start, stop, step) = (psi.start, psi.stop, psi.step) ;
                let m : Vec<(Pauli, Complex64)> =
                    self.it.members()[start as usize..stop as usize].iter()
                    .step_by(step as usize)
                    .map(|p| (Pauli::from(p.0.clone()), p.1))
                    .collect() ;
                let m = SliceResult::Slice(m) ;
                Ok(m)
            }
            SliceOrInt::Int(idx) => {
                (0 <= idx && idx < self.it.members().len() as isize).then(|| ())
                    .ok_or(PyException::new_err(format!("__getitem__ called on invalid index {}", idx))) ? ;
                let m = &self.it.members()[idx as usize] ;
                let m = SliceResult::It((Pauli::from(m.0.clone()), m.1)) ;
                Ok(m)
            }
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.it.num_qubits()
    }

    pub fn to_matrix(&self) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix())
    }
/*
    pub fn to_matrix_binary(&self) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix_binary())
    }

    pub fn to_matrix_accel(&self) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix_accel())
    }

    pub fn to_matrix_reduce(&self) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix_reduce())
    }

    pub fn to_matrix_rayon(&self) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix_rayon())
    }

    pub fn to_matrix_rayon_chunked(&self, step : usize) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix_rayon_chunked(step))
    }
*/
    #[args(mode = "\"\"")]
    fn to_matrix_mode(&self, mode: &str) -> PyResult<SpMat> {
        let mode = if mode == "" { None }
        else {
            let mode = AccelMode::try_from(mode)
                .map_err(|_| PyException::new_err(format!("to_matrix_mode: unrecognized mode {}", mode))) ? ;
            Some(mode)
        } ;
        Ok(SpMat::new_from_csmatrix(self.it.to_matrix_mode(&mode)))
    }

    fn __repr__(&self) -> String {
        let labels : Vec<String> = self.it.members()
            .iter()
            .map(|(p,_)| p.label())
            .collect() ;
        let coeffs : Vec<String> = self.it.members()
            .iter()
            .map(|(_,c)| qrusty::util::complex64_to_string(*c))
            .collect() ;
        format!("SparsePauliOp('{}', [{}])", labels.join("','"), coeffs.join(", "))
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }

}


fn reg(x: &ArrayView<Complex64, Dim<[usize; 1]>>, tol: f64) -> Array<Complex64, Dim<[usize; 1]>> {
    x.iter()
	.map(|x| if x.norm() < tol { Complex64::new(tol, 0.0) } else { *x })
	.collect()
}

fn precond(
    spmat : &sprs::CsMatI<Complex64,  u64>,
    dx : &ArrayView<Complex64, Dim<[usize; 1]>>,
    e : Complex64,
    tol : f64) -> Array<Complex64, Dim<[usize; 1]>> {
    let x : Array<Complex64, Dim<[usize; 1]>> =
	spmat.diag()
	.to_dense()
	.iter()
	.map(|x| x - e)
	.collect() ;
    let y = dx / reg(&x.view(), tol) ;
    y
}

fn precond2(
    diag : &ArrayView<Complex64, Dim<[usize; 1]>>,
    dx : &ArrayView<Complex64, Dim<[usize; 1]>>,
    e : Complex64,
    tol : f64) -> Array<Complex64, Dim<[usize; 1]>> {
    let x : Array<Complex64, Dim<[usize; 1]>> =
	diag.iter()
	.map(|x| x - e)
	.collect() ;
    let y = dx / reg(&x.view(), tol) ;
    y
}


/// A Python module implemented in Rust.
#[pymodule]
fn pyqrusty(_py: Python, m: &PyModule) -> PyResult<()> {

    // example using complex numbers
    fn spmat_dot_densevec(spmat: &SpMat, x: ArrayView<Complex64, Dim<[usize; 1]>>) -> PyResult<Array<Complex64, Dim<[usize; 1]>>> {
        spmat.map_immut(|| Err(PyException::new_err("cannot multiply with an exported sparse matrix")),
			|spmat| Ok(qrusty::accel::rowwise::spmat_dot_densevec(&spmat, &x)))
    }

    // wrapper of `conj`
    #[pyfn(m)]
    #[pyo3(name = "spmat_dot_densevec")]
    fn spmat_dot_densevec_py<'py>(
        py: Python<'py>,
	spmat: &SpMat,
        x: PyReadonlyArray1<'_, Complex64>,
    ) -> PyResult<PyObject> {
        let y = spmat_dot_densevec(spmat, x.as_array()) ? ;
	Ok(y.into_pyarray(py).into())
    }

    #[pyfn(m)]
    #[pyo3(name = "axpy")]
    fn axpy_py<'py>(
        py: Python<'py>,
        a: Complex64,
        x: PyReadonlyArray1<'_, Complex64>,
        y: PyReadonlyArray1<'_, Complex64>,
    ) -> PyResult<PyObject> {
        let x = x.as_array();
        let y = y.as_array();
        let z = qrusty::accel::axpy(a, &x, &y);
	Ok(z.into_pyarray(py).into())
    }

    #[pyfn(m)]
    #[pyo3(name = "ax")]
    fn ax_py<'py>(
        py: Python<'py>,
        a: Complex64,
        x: PyReadonlyArray1<'_, Complex64>,
    ) -> PyResult<PyObject> {
        let x = x.as_array();
        let z = qrusty::accel::ax(a, &x);
	Ok(z.into_pyarray(py).into())
    }


    #[pyfn(m)]
    #[pyo3(name = "precond")]
    fn precond_py<'py>(
        py: Python<'py>,
	m : &SpMat,
        dx: PyReadonlyArray1<'_, Complex64>,
        e: Complex64,
	tol: f64,
    ) -> PyResult<PyObject> {
        let dx = dx.as_array();
        let z = m.map_immut(|| Err(PyException::new_err("cannot call precond with an exported sparse matrix")),
			    |spmat| Ok(precond(spmat, &dx, e, tol))) ? ;
	Ok(z.into_pyarray(py).into())
    }

    #[pyfn(m)]
    #[pyo3(name = "precond2")]
    fn precond2_py<'py>(
        py: Python<'py>,
        diag: PyReadonlyArray1<'_, Complex64>,
        dx: PyReadonlyArray1<'_, Complex64>,
        e: Complex64,
	tol: f64,
    ) -> PyResult<PyObject> {
        let diag = diag.as_array();
        let dx = dx.as_array();
        let z = precond2(&diag, &dx, e, tol) ;
	Ok(z.into_pyarray(py).into())
    }

    m.add_class::<Pauli>()?;
    m.add_class::<SparsePauliOp>()?;
    m.add_class::<SpMat>()?;
    Ok(())
}

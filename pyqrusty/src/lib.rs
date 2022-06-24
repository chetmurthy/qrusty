#![allow(unused_imports)]
use pyo3::prelude::*;
use numpy::{IntoPyArray, PyReadonlyArray1};
use num_complex::Complex64;
use pyo3::wrap_pyfunction;
use pyo3::Python;
use sprs::{CompressedStorage, CsMatI};
use pyo3::PyErr ;
use pyo3::exceptions::PyOSError;
use pyo3::exceptions::PyException;

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

    fn __repr__(&self) -> String {
        match &*(self.it) {
            None => "<already-dropped sparse matrix of type Complex64>".to_string(),
            Some(spmat) => {
                let (rows, cols) = spmat.shape() ;
                let nnz = spmat.nnz() ;
                format!("<{}x{} sparse matrix of type Complex64\n\twith {} stored elements in Compressed Sparse Row format>",
                        rows,  cols, nnz)
            }
        }
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
      match &mut *(self.it) {
          Some(lhs) => {
              lhs.scale(factor) ;
              Ok(())
          },
          _ => Err(PyException::new_err(format!("cannot scale already-exported sparse matrix")))
      }
    }

/*    
    fn __neg__(&self) -> PyResult<SpMat> {
        self.binop("negate",  |y: &CsMatI::<Complex64, u64, u64>| -y)
    }
  */  
    #[args(tolerance = "1e-7")]
    fn count_zeros(&self,  tolerance: f64) -> PyResult<usize> {
        match &*(self.it) {
            None => Err(PyException::new_err("cannot count zeroes of an exported sparse matrix")),
            Some(spmat) => {
                Ok(qrusty::util::csmatrix_nz(&spmat, tolerance))
            }
        }
    }
    
    #[args(tolerance = "1e-7")]
    fn eliminate_zeros(&self,  tolerance: f64) -> PyResult<SpMat> {
        match &*(self.it) {
            None => Err(PyException::new_err("cannot eliminate zeroes of an exported sparse matrix")),
            Some(spmat) => {
                let mat = qrusty::util::csmatrix_eliminate_zeroes(&spmat, tolerance) ;
                Ok(SpMat { it : Box::new(Option::Some(mat)) })
            }
        }
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

    pub fn write_to_file_for_scipy(
        &mut self,
        save_path : &str,
    ) -> PyResult<()> {

        match &*(self.it) {
            None =>
                Err(PyException::new_err("cannot write an already-destroyed sparse matrix")),
            Some(spmat) => {
                qrusty::util::write_matrix_market_for_scipy(&save_path, &spmat)
                    .map_err(|e| PyOSError::new_err(e))
            }
        }
    }

    pub fn write_to_file(
        &mut self,
        save_path : &str,
    ) -> PyResult<()> {

        match &*(self.it) {
            None =>
                Err(PyException::new_err("cannot write an already-destroyed sparse matrix")),
            Some(spmat) => {
                sprs::io::write_matrix_market(&save_path, spmat)
                    .map_err(|e| PyOSError::new_err(e))
            }
        }
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
        Ok(Pauli { it })
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

 #[pymethods]
impl SparsePauliOp {

    #[new]
    pub fn new(paulis : Vec<Pauli>, coeffs : Vec<Complex64>) -> PyResult<SparsePauliOp> {
        let paulis = paulis.iter().map(|p| p.it.clone()).collect() ;
        let coeffs : &[Complex64] = &coeffs[..] ;
        let it = qrusty::SparsePauliOp::new(paulis, coeffs) ? ;
        Ok(SparsePauliOp { it })
    }

    pub fn num_qubits(&self) -> usize {
        self.it.num_qubits()
    }

    pub fn to_matrix(&self) -> SpMat {
        SpMat::new_from_csmatrix(self.it.to_matrix())
    }

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

/// A Python module implemented in Rust.
#[pymodule]
fn pyqrusty(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Pauli>()?;
    m.add_class::<SparsePauliOp>()?;
    m.add_class::<SpMat>()?;
    Ok(())
}

#![allow(unused_imports)]
use pyo3::prelude::*;
use numpy::{IntoPyArray, PyReadonlyArray1};
use num_complex::Complex64;
use pyo3::wrap_pyfunction;
use pyo3::Python;

#[pyclass]
#[repr(transparent)]
#[derive(Clone)]
pub struct PyPauli {
    pub it: qrusty::Pauli,
}

impl From<qrusty::Pauli> for PyPauli {
    fn from(it: qrusty::Pauli) -> Self {
        PyPauli { it : it.clone() }
    }
}

#[pymethods]
impl PyPauli {

    #[new]
    pub fn new(label: String) -> PyResult<PyPauli> {
        let it = qrusty::Pauli::new(&label) ? ;
        Ok(PyPauli { it })
    }

    pub fn num_qubits(&self) -> usize {
        self.it.num_qubits()
    }

    pub fn label(&self) -> String {
        self.it.label()
    }

}

#[pyclass]
#[repr(transparent)]
#[derive(Clone)]
pub struct PySparsePauliOp {
    pub it: qrusty::SparsePauliOp,
}

impl From<qrusty::SparsePauliOp> for PySparsePauliOp {
    fn from(it: qrusty::SparsePauliOp) -> Self {
        PySparsePauliOp { it : it.clone() }
    }
}

#[pymethods]
impl PySparsePauliOp {

    #[new]
    pub fn new(paulis : Vec<PyPauli>, coeffs : Vec<Complex64>) -> PyResult<PySparsePauliOp> {
        let paulis = paulis.iter().map(|p| p.it.clone()).collect() ;
        let coeffs : &[Complex64] = &coeffs[..] ;
        let it = qrusty::SparsePauliOp::new(paulis, coeffs) ? ;
        Ok(PySparsePauliOp { it })
    }

    pub fn num_qubits(&self) -> usize {
        self.it.num_qubits()
    }

}

/// A Python module implemented in Rust.
#[pymodule]
fn pyqrusty(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPauli>()?;
    m.add_class::<PySparsePauliOp>()?;
    Ok(())
}

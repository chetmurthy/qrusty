#![allow(unused_imports)]
use pyo3::prelude::*;
use numpy::{IntoPyArray, PyReadonlyArray1};
use num_complex::Complex64;
use pyo3::wrap_pyfunction;
use pyo3::Python;

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
    Ok(())
}

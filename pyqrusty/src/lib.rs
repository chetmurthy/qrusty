use pyo3::prelude::*;

#[pyclass]
pub struct Pauli {
    pub it: qrusty::Pauli,
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

}


/// A Python module implemented in Rust.
#[pymodule]
fn pyqrusty(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Pauli>()?;
    Ok(())
}

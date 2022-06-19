use pyo3::prelude::*;

#[pyclass]
pub struct RustStruct {
    #[pyo3(get, set)]
    pub data: String,
    #[pyo3(get, set)]
    pub vector: Vec<u8>,
}

#[pymethods]
impl RustStruct {
    #[new]
    pub fn new(data: String, vector: Vec<u8>) -> RustStruct {
        RustStruct { data, vector }
    }
    pub fn printer(&self) {
        println!("{}", self.data);
        for i in &self.vector {
            println!("{}", i);
        }
    }
    pub fn extend_vector(&mut self, extension: Vec<u8>) {
        println!("{}", self.data);
        for i in extension {
            self.vector.push(i);
        }
    }
}

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
}


/// A Python module implemented in Rust.
#[pymodule]
fn pyqrusty(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<RustStruct>()?;
    m.add_class::<Pauli>()?;
    Ok(())
}

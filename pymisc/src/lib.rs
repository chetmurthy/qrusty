
use pyo3::prelude::*;
use num_complex::Complex64;
use ndarray::{ArrayViewD, ArrayD, ArrayView, Dim, Array, Array1, par_azip};
use numpy::{IntoPyArray, PyReadonlyArray1, PyArrayDyn, PyReadonlyArrayDyn, PyArray1};
use pyo3::wrap_pyfunction;
use pyo3::Python;
use pyo3::PyErr ;
use pyo3::types::{ PySlice } ;
use pyo3::exceptions::PyOSError;
use pyo3::exceptions::PyException;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}



/// A Python module implemented in Rust.
#[pymodule]
fn pymisc(_py: Python, m: &PyModule) -> PyResult<()> {

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
        let z : Array<Complex64, Dim<[usize; 1]>> = x.iter().zip(y.iter()).map(|(x,y)| a * x + y).collect() ;
	Ok(z.into_pyarray(py).into())
    }

    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    Ok(())
}

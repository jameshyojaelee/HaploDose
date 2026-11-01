//! Python bindings for the Haplodose data plane. The GIL is released while rayon
//! threads run, so Python can read the next chunk concurrently.

mod aggregate;

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn slice<'a, T: numpy::Element>(a: &'a PyReadonlyArray1<'_, T>, name: &str) -> PyResult<&'a [T]> {
    a.as_slice()
        .map_err(|_| PyValueError::new_err(format!("{name} must be a contiguous 1-D array")))
}

/// Sum CSR rows into genomic bins. Returns (rows x bins float32, row totals float64).
#[pyfunction]
fn bin_sums<'py>(
    py: Python<'py>,
    indptr: PyReadonlyArray1<'py, i64>,
    indices: PyReadonlyArray1<'py, i32>,
    data: PyReadonlyArray1<'py, f32>,
    gene_bin: PyReadonlyArray1<'py, i32>,
    n_bins: usize,
) -> PyResult<(Bound<'py, PyArray2<f32>>, Bound<'py, PyArray1<f64>>)> {
    let (p, i, d, b) = (
        slice(&indptr, "indptr")?,
        slice(&indices, "indices")?,
        slice(&data, "data")?,
        slice(&gene_bin, "gene_bin")?,
    );
    let n_rows = p.len().saturating_sub(1);
    let (out, totals) = py
        .detach(|| aggregate::bin_sums(p, i, d, b, n_bins))
        .map_err(PyValueError::new_err)?;
    let out = Array2::from_shape_vec((n_rows, n_bins), out)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok((out.into_pyarray(py), totals.into_pyarray(py)))
}

/// Column sums per row group. Returns groups x columns float64.
#[pyfunction]
fn grouped_column_sums<'py>(
    py: Python<'py>,
    indptr: PyReadonlyArray1<'py, i64>,
    indices: PyReadonlyArray1<'py, i32>,
    data: PyReadonlyArray1<'py, f32>,
    row_group: PyReadonlyArray1<'py, i32>,
    n_groups: usize,
    n_cols: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let (p, i, d, g) = (
        slice(&indptr, "indptr")?,
        slice(&indices, "indices")?,
        slice(&data, "data")?,
        slice(&row_group, "row_group")?,
    );
    let out = py
        .detach(|| aggregate::grouped_column_sums(p, i, d, g, n_groups, n_cols))
        .map_err(PyValueError::new_err)?;
    let out = Array2::from_shape_vec((n_groups, n_cols), out)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(out.into_pyarray(py))
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bin_sums, m)?)?;
    m.add_function(wrap_pyfunction!(grouped_column_sums, m)?)?;
    Ok(())
}

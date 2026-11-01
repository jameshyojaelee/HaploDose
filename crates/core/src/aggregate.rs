//! Row-parallel aggregation of cell-major CSR count chunks.
//!
//! Both kernels read a bounded chunk handed over from Python and never densify
//! the cells x genes matrix. Outputs are cells x bins or groups x genes.

use rayon::prelude::*;

/// Errors are plain strings; lib.rs maps them to Python ValueError.
pub type Result<T> = std::result::Result<T, String>;

pub fn check_csr(indptr: &[i64], indices: &[i32], data: &[f32], n_cols: usize) -> Result<usize> {
    if indptr.is_empty() || indptr[0] != 0 {
        return Err("indptr must start at 0".into());
    }
    let nnz = *indptr.last().unwrap();
    if nnz < 0 || nnz as usize != indices.len() || indices.len() != data.len() {
        return Err("indptr, indices and data lengths disagree".into());
    }
    if indptr.windows(2).any(|w| w[1] < w[0]) {
        return Err("indptr must be nondecreasing".into());
    }
    if indices.iter().any(|&j| j < 0 || j as usize >= n_cols) {
        return Err("column index out of range".into());
    }
    if data.iter().any(|&x| !x.is_finite() || x < 0.0) {
        return Err("counts must be finite and nonnegative".into());
    }
    Ok(indptr.len() - 1)
}

/// out[r, b] = sum of row r over genes with gene_bin == b; genes with -1 are skipped.
/// Also returns each row's total over all genes (the library size).
pub fn bin_sums(
    indptr: &[i64],
    indices: &[i32],
    data: &[f32],
    gene_bin: &[i32],
    n_bins: usize,
) -> Result<(Vec<f32>, Vec<f64>)> {
    let n_rows = check_csr(indptr, indices, data, gene_bin.len())?;
    if gene_bin.iter().any(|&b| b < -1 || b >= n_bins as i32) {
        return Err("gene_bin must be -1 or a valid bin index".into());
    }
    let mut out = vec![0f32; n_rows * n_bins];
    let mut totals = vec![0f64; n_rows];
    if n_bins == 0 {
        // Library sizes only.
        totals.par_iter_mut().enumerate().for_each(|(r, total)| {
            let (lo, hi) = (indptr[r] as usize, indptr[r + 1] as usize);
            *total = data[lo..hi].iter().map(|&x| x as f64).sum();
        });
        return Ok((out, totals));
    }
    out.par_chunks_mut(n_bins)
        .zip(totals.par_iter_mut())
        .enumerate()
        .for_each(|(r, (row, total))| {
            let (lo, hi) = (indptr[r] as usize, indptr[r + 1] as usize);
            let mut sum = 0f64;
            for k in lo..hi {
                let x = data[k];
                sum += x as f64;
                let b = gene_bin[indices[k] as usize];
                if b >= 0 {
                    row[b as usize] += x;
                }
            }
            *total = sum;
        });
    Ok((out, totals))
}

/// out[g, j] = sum over rows with row_group == g of X[r, j]; rows with -1 are skipped.
/// Thread-local accumulators bound memory to threads x groups x columns.
pub fn grouped_column_sums(
    indptr: &[i64],
    indices: &[i32],
    data: &[f32],
    row_group: &[i32],
    n_groups: usize,
    n_cols: usize,
) -> Result<Vec<f64>> {
    let n_rows = check_csr(indptr, indices, data, n_cols)?;
    if row_group.len() != n_rows {
        return Err("row_group length must equal the number of rows".into());
    }
    if row_group.iter().any(|&g| g < -1 || g >= n_groups as i32) {
        return Err("row_group must be -1 or a valid group index".into());
    }
    let width = n_groups * n_cols;
    let out = (0..n_rows)
        .into_par_iter()
        .fold(
            || vec![0f64; width],
            |mut acc, r| {
                let g = row_group[r];
                if g >= 0 {
                    let base = g as usize * n_cols;
                    for k in indptr[r] as usize..indptr[r + 1] as usize {
                        acc[base + indices[k] as usize] += data[k] as f64;
                    }
                }
                acc
            },
        )
        .reduce(
            || vec![0f64; width],
            |mut a, b| {
                a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                a
            },
        );
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Rows: [1 0 2], [0 3 0], [4 0 0]
    fn toy() -> (Vec<i64>, Vec<i32>, Vec<f32>) {
        (vec![0, 2, 3, 4], vec![0, 2, 1, 0], vec![1., 2., 3., 4.])
    }

    #[test]
    fn bins_and_totals() {
        let (p, i, d) = toy();
        let (out, tot) = bin_sums(&p, &i, &d, &[0, 1, -1], 2).unwrap();
        assert_eq!(out, vec![1., 0., 0., 3., 4., 0.]);
        assert_eq!(tot, vec![3., 3., 4.]);
    }

    #[test]
    fn grouped_sums_skip_unassigned_rows() {
        let (p, i, d) = toy();
        let out = grouped_column_sums(&p, &i, &d, &[0, -1, 0], 1, 3).unwrap();
        assert_eq!(out, vec![5., 0., 2.]);
    }

    #[test]
    fn rejects_bad_input() {
        let (p, i, d) = toy();
        assert!(bin_sums(&p, &i, &d, &[0, 5, 0], 2).is_err());
        assert!(bin_sums(&p, &[0, 9, 1, 0], &d, &[0, 1, 0], 2).is_err());
        assert!(bin_sums(&p, &i, &[1., -2., 3., 4.], &[0, 1, 0], 2).is_err());
    }
}

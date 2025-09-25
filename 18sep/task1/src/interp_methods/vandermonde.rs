use nalgebra::{DMatrix, DVector};

pub fn vandermonde_interpolation(
    n: u32,
    x_list_input: &[f64],
    f_list_input: &[f64],
    x_list_with_midpoints: &[f64],
) -> Vec<f64> {
    let a_list = solve_matrix(n, x_list_input, f_list_input);
    x_list_with_midpoints
        .iter()
        .map(|x| {
            a_list
                .iter()
                .enumerate()
                .map(|(i, a_i)| x.powi(i as i32) * a_i)
                .sum()
        })
        .collect()
}

fn solve_matrix(n: u32, x_list: &[f64], f_list: &[f64]) -> Vec<f64> {
    let x_matrix: Vec<f64> = (0..=n)
        .into_iter()
        .map(|i| x_list.iter().map(move |x| x.powi(i as i32)))
        .flatten()
        .collect();

    let x_matrix = DMatrix::from_vec((n + 1) as usize, (n + 1) as usize, x_matrix);
    let f_vector = DVector::from_row_slice(f_list);

    x_matrix
        .lu()
        .solve(&f_vector)
        .expect("Failed to solve")
        .into_iter()
        .map(|a| *a)
        .collect()
}

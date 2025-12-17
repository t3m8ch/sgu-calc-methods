use nalgebra::{DMatrix, DVector};

pub fn vandermonde_interpolation(x_nodes: &[f64], f_nodes: &[f64], x_targets: &[f64]) -> Vec<f64> {
    // TODO: Validate len for x_nodes and f_nodes
    let n = x_nodes.len() - 1;
    let a_list = solve_matrix(n, x_nodes, f_nodes);
    x_targets
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

fn solve_matrix(n: usize, x_list: &[f64], f_list: &[f64]) -> Vec<f64> {
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
        .expect("Failed to solve") // TODO: Change to Result
        .into_iter()
        .map(|a| *a)
        .collect()
}

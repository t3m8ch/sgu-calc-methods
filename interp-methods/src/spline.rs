use nalgebra::{DMatrix, DVector};
use shared::cli::print_matrix;

pub fn cubic_spline_interpolation(x_nodes: &[f64], f_nodes: &[f64], x_targets: &[f64]) -> Vec<f64> {
    let (matrix_data, rhs_data) = build_spline_system(x_nodes, f_nodes);
    let n = x_nodes.len() - 1;
    let size = 4 * n;

    println!("Матрица системы ({}x{}):", size, size);
    print_matrix(&matrix_data, size, size);

    println!("\nВектор правой части:");
    println!("{:?}\n", rhs_data);

    let coefficients = solve_spline_coefficients(x_nodes, f_nodes);

    println!("Коэффициенты сплайнов [a_i, b_i, c_i, d_i] для каждого интервала:");
    for (i, coef) in coefficients.iter().enumerate() {
        println!(
            "  Интервал {}: a={:.6}, b={:.6}, c={:.6}, d={:.6}",
            i, coef[0], coef[1], coef[2], coef[3]
        );
    }
    println!();

    x_targets
        .iter()
        .map(|&x| evaluate_spline(x, x_nodes, &coefficients))
        .collect()
}

fn solve_spline_coefficients(x_nodes: &[f64], f_nodes: &[f64]) -> Vec<[f64; 4]> {
    let n = x_nodes.len() - 1;
    let size = 4 * n;
    let (matrix, rhs) = build_spline_system(x_nodes, f_nodes);

    let matrix = DMatrix::from_row_slice(size, size, &matrix);
    let rhs = DVector::from_row_slice(&rhs);

    matrix
        .lu()
        .solve(&rhs)
        .expect("Failed to solve spline system")
        .as_slice()
        .chunks(4)
        .map(|chunk| [chunk[0], chunk[1], chunk[2], chunk[3]])
        .collect()
}

fn build_spline_system(x_nodes: &[f64], f_nodes: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = x_nodes.len() - 1;
    let size = 4 * n;
    let mut matrix = vec![0.0; size * size];
    let mut rhs = vec![0.0; size];

    let h: Vec<f64> = (0..n).map(|i| x_nodes[i + 1] - x_nodes[i]).collect();

    (0..n).for_each(|i| {
        matrix[i * size + 4 * i] = 1.0;
        rhs[i] = f_nodes[i];
    });

    (0..n).for_each(|i| {
        let row = n + i;
        matrix[row * size + 4 * i] = 1.0;
        matrix[row * size + 4 * i + 1] = h[i];
        matrix[row * size + 4 * i + 2] = h[i].powi(2);
        matrix[row * size + 4 * i + 3] = h[i].powi(3);
        rhs[row] = f_nodes[i + 1];
    });

    (0..n - 1).for_each(|i| {
        let row = 2 * n + i;
        matrix[row * size + 4 * i + 1] = 1.0;
        matrix[row * size + 4 * i + 2] = 2.0 * h[i];
        matrix[row * size + 4 * i + 3] = 3.0 * h[i].powi(2);
        matrix[row * size + 4 * (i + 1) + 1] = -1.0;
    });

    (0..n - 1).for_each(|i| {
        let row = 3 * n - 1 + i;
        matrix[row * size + 4 * i + 2] = 2.0;
        matrix[row * size + 4 * i + 3] = 6.0 * h[i];
        matrix[row * size + 4 * (i + 1) + 2] = -2.0;
    });

    let row1 = 4 * n - 2;
    matrix[row1 * size + 2] = 2.0;

    let row2 = 4 * n - 1;
    matrix[row2 * size + 4 * (n - 1) + 2] = 2.0;
    matrix[row2 * size + 4 * (n - 1) + 3] = 6.0 * h[n - 1];

    (matrix, rhs)
}

fn evaluate_spline(x: f64, x_nodes: &[f64], coefficients: &[[f64; 4]]) -> f64 {
    let interval = x_nodes
        .windows(2)
        .position(|w| x >= w[0] && x <= w[1])
        .unwrap_or(coefficients.len() - 1);

    let [a, b, c, d] = coefficients[interval];
    let dx = x - x_nodes[interval];

    a + b * dx + c * dx.powi(2) + d * dx.powi(3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_spline() {
        let x_nodes = vec![0.0, 1.0, 2.0, 3.0];
        let f_nodes = vec![0.0, 1.0, 4.0, 9.0];
        let x_targets = vec![0.5, 1.5, 2.5];

        let result = cubic_spline_interpolation(&x_nodes, &f_nodes, &x_targets);
        println!("Spline values: {:?}", result);
    }
}

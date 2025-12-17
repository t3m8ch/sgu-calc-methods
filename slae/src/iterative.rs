use shared::cli::print_matrix;

use crate::gauss::calculate_determinant;

pub fn solve_iterative(matrix: &[f64], b: &[f64], epsilon: f64, max_iterations: usize) -> Vec<f64> {
    let n = b.len();

    println!("Матрица A:");
    print_matrix(matrix, n, n);

    let det = calculate_determinant(matrix, n);
    println!("\nОпределитель матрицы A: {:.6}", det);

    println!("\nСтолбец b: {:.2?}", b);

    let alpha = compute_alpha(matrix, n);
    println!("\nМатрица alpha:");
    print_matrix(&alpha, n, n);

    let beta = compute_beta(matrix, b, n);
    println!("\nСтолбец beta: {:.10?}", beta);

    println!("\nСчитаем до точности epsilon = {:.0e}", epsilon);

    let mut xk = vec![0.0; n];
    println!(
        "x^(0) = [ {} ]",
        xk.iter()
            .map(|_| "0.0".to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );

    for iteration in 0..max_iterations {
        let xkp1 = compute_next_iteration(&alpha, &xk, &beta, n);

        print!("x^({}) = [ ", iteration + 1);
        xkp1.iter().for_each(|&val| print!("{:.10} ", val));
        println!("]");

        if norm_stop(&xk, &xkp1, epsilon) {
            println!("\nСходимость достигнута за {} итераций", iteration + 1);
            return xkp1;
        }

        xk = xkp1;
    }

    println!(
        "\nДостигнуто максимальное количество итераций: {}",
        max_iterations
    );
    xk
}

fn compute_alpha(matrix: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .flat_map(|i| {
            (0..n).map(move |j| {
                if i == j {
                    0.0
                } else {
                    -matrix[i * n + j] / matrix[i * n + i]
                }
            })
        })
        .collect()
}

fn compute_beta(matrix: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    (0..n).map(|i| b[i] / matrix[i * n + i]).collect()
}

fn compute_next_iteration(alpha: &[f64], xk: &[f64], beta: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let sum: f64 = (0..n).map(|j| alpha[i * n + j] * xk[j]).sum();
            sum + beta[i]
        })
        .collect()
}

fn norm_stop(xk: &[f64], xkp1: &[f64], epsilon: f64) -> bool {
    xk.iter()
        .zip(xkp1.iter())
        .map(|(x, xp)| (x - xp).abs())
        .fold(0.0, f64::max)
        < epsilon
}

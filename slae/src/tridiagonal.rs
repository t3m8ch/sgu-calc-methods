use shared::cli::print_matrix;

pub fn solve_tridiagonal(matrix: &[f64], b: &[f64]) -> Vec<f64> {
    let n = b.len();

    println!("Матрица A:");
    print_matrix(matrix, n, n);
    println!("\nСтолбец b: {:?}", b);

    let (p_coeffs, q_coeffs) = forward_sweep(matrix, b, n);

    println!("\nПрямая прогонка");
    println!("Список Pi и Qi");
    print!("P: [");
    p_coeffs.iter().for_each(|&val| print!(" {:.10}", val));
    println!(" ]");

    print!("Q: [");
    q_coeffs.iter().for_each(|&val| print!(" {:.10}", val));
    println!(" ]");

    let solution = backward_sweep(&p_coeffs, &q_coeffs, n);

    println!("\nОбратная прогонка");
    println!("x {} = {:.10}", n, solution[n - 1]);
    (0..n - 1).rev().for_each(|i| {
        println!("x {} = {:.10}", i + 1, solution[i]);
    });

    solution
}

fn forward_sweep(matrix: &[f64], b: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut p = vec![0.0; n - 1];
    let mut q = vec![0.0; n];

    let a0 = matrix[0 * n + 0];
    let c0 = matrix[0 * n + 1];
    p[0] = -c0 / a0;
    q[0] = b[0] / a0;

    for i in 1..n - 1 {
        let a_i = matrix[i * n + i - 1];
        let b_i = matrix[i * n + i];
        let c_i = matrix[i * n + i + 1];

        p[i] = c_i / (-b_i - a_i * p[i - 1]);
    }

    for i in 1..n {
        let a_i = matrix[i * n + i - 1];
        let b_i = matrix[i * n + i];

        q[i] = (a_i * q[i - 1] - b[i]) / (-b_i - a_i * p[i - 1]);
    }

    (p, q)
}

fn backward_sweep(p: &[f64], q: &[f64], n: usize) -> Vec<f64> {
    let mut x = vec![0.0; n];

    x[n - 1] = q[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    x
}

use comfy_table::Table;

pub fn print_task_11(v: f64) {
    let (x, y_calc, y_exact, error) = fredholm_solver(v, 3);

    let mut table = Table::new();
    table.set_header(vec!["x", "y_мет", "y_точн", "погрешность"]);

    for i in 0..x.len() {
        table.add_row(vec![
            format!("{:7.3}", x[i]),
            format!("{:7.3}", y_calc[i]),
            format!("{:7.3}", y_exact[i]),
            format!("{:7.3}", error[i]),
        ]);
    }

    println!("{table}");
}

fn fredholm_solver(v: f64, rank: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let alpha = build_alpha(rank);
    let gamma = build_gamma(rank, v);

    let mut system_matrix = vec![0.0; rank * rank];
    for i in 0..rank {
        for j in 0..rank {
            system_matrix[i * rank + j] = if i == j { 1.0 } else { 0.0 } + alpha[i * rank + j];
        }
    }

    let coeffs = solve_slae_silent(&system_matrix, &gamma);

    let step = 0.1;
    let x_vals: Vec<f64> = (0..=10).map(|i| i as f64 * step).collect();
    let mut y_num = vec![0.0; x_vals.len()];

    for (i, &x) in x_vals.iter().enumerate() {
        let mut y_val = rhs_func(x, v);
        for j in 0..rank {
            y_val -= coeffs[j] * basis_func(j, x);
        }
        y_num[i] = y_val;
    }

    let y_exact: Vec<f64> = x_vals.iter().map(|&x| v * x).collect();
    let error: Vec<f64> = y_num
        .iter()
        .zip(y_exact.iter())
        .map(|(yn, ye)| (yn - ye).abs())
        .collect();

    (x_vals, y_num, y_exact, error)
}

fn rhs_func(x: f64, v: f64) -> f64 {
    v * (4.0 / 3.0 * x + 0.25 * x.powi(2) + 0.2 * x.powi(3))
}

fn basis_func(index: usize, x: f64) -> f64 {
    match index {
        0 => x,
        1 => x.powi(2),
        2 => x.powi(3),
        _ => panic!("Неверный индекс базисной функции"),
    }
}

fn integrate_trapezoidal<F>(f: F, a: f64, b: f64, n: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));

    for i in 1..n {
        sum += f(a + i as f64 * h);
    }

    sum * h
}

fn build_alpha(size: usize) -> Vec<f64> {
    let mut alpha = vec![0.0; size * size];

    for i in 0..size {
        for j in 0..size {
            alpha[i * size + j] =
                integrate_trapezoidal(|t| basis_func(i, t) * basis_func(j, t), 0.0, 1.0, 1000);
        }
    }

    alpha
}

fn build_gamma(size: usize, v: f64) -> Vec<f64> {
    let mut gamma = vec![0.0; size];

    for i in 0..size {
        gamma[i] = integrate_trapezoidal(|t| rhs_func(t, v) * basis_func(i, t), 0.0, 1.0, 1000);
    }

    gamma
}

pub fn solve_slae_silent(matrix: &[f64], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut a = matrix.to_vec();
    let mut b_work = b.to_vec();

    forward_elimination(&mut a, &mut b_work, n);
    backward_substitution(&a, &b_work, n)
}

fn forward_elimination(a: &mut [f64], b: &mut [f64], n: usize) {
    for k in 0..n - 1 {
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];

            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            a[i * n + k] = 0.0;
            b[i] -= factor * b[k];
        }
    }
}

fn backward_substitution(u: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    (0..n).rev().fold(vec![0.0; n], |mut x, i| {
        let sum: f64 = (i + 1..n).map(|j| u[i * n + j] * x[j]).sum();

        x[i] = (b[i] - sum) / u[i * n + i];
        x
    })
}

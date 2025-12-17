use shared::cli::print_matrix;

pub fn solve_slae_gauss(matrix: &[f64], b: &[f64]) -> Vec<f64> {
    let n = b.len();

    let mut a = matrix.to_vec();
    let mut b_work = b.to_vec();

    println!("Матрица A:");
    print_matrix(&a, n, n);

    let det = calculate_determinant(&a, n);
    println!("\nОпределитель матрицы A = {:.6}", det);

    println!("\nКолонка b: {:?}", b);

    forward_elimination(&mut a, &mut b_work, n);

    println!("\nМатрица после прямого прохода:");
    print_matrix(&a, n, n);

    println!("\nВектор b после прямого прохода: {:.6?}", b_work);

    let solution = backward_substitution(&a, &b_work, n);

    solution
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

pub fn calculate_determinant(matrix: &[f64], n: usize) -> f64 {
    let mut a = matrix.to_vec();
    let mut det = 1.0;

    for k in 0..n - 1 {
        if a[k * n + k].abs() < 1e-10 {
            return 0.0;
        }

        det *= a[k * n + k];

        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
        }
    }

    det * a[(n - 1) * n + (n - 1)]
}

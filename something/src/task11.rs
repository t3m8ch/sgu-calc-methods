use comfy_table::Table;

pub fn print_task_11(v: f64, num_points: usize) {
    let (a, b) = (0.0, 1.0);
    let n = 3;
    let lambda = 1.0;

    let y_exact = |x: f64| v * x;
    let f = |x: f64| v * (4.0 / 3.0 * x + 1.0 / 4.0 * x.powi(2) + 1.0 / 5.0 * x.powi(3));

    let a_funcs: Vec<Box<dyn Fn(f64) -> f64>> = vec![
        Box::new(|x| x),
        Box::new(|x| x.powi(2)),
        Box::new(|x| x.powi(3)),
    ];

    let mut alpha = vec![vec![0.0; n]; n];
    for i in 0..n {
        for k in 0..n {
            alpha[i][k] = 1.0 / ((i + k + 3) as f64);
        }
    }

    let mut gamma = vec![0.0; n];
    for i in 0..n {
        gamma[i] = integrate_fredholm(|x| f(x) * x.powi((i + 1) as i32), a, b, 1000);
    }

    let mut mat_a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            mat_a[i][j] = if i == j { 1.0 } else { 0.0 } + lambda * alpha[j][i];
        }
    }

    let q = solve_linear_system(&mat_a, &gamma);

    let y_numerical = |x: f64| -> f64 {
        let mut result = f(x);
        for i in 0..n {
            result -= lambda * q[i] * a_funcs[i](x);
        }
        result
    };

    let x_test: Vec<f64> = (0..num_points)
        .map(|i| a + (b - a) * i as f64 / (num_points - 1) as f64)
        .collect();

    let mut table = Table::new();
    table.set_header(vec!["x", "y_метода", "y_точн", "eps"]);

    for &x in &x_test {
        let y_num = y_numerical(x);
        let y_ex = y_exact(x);
        let error = (y_num - y_ex).abs();

        table.add_row(vec![
            format!("{:.7}", x),
            format!("{:.7}", y_num),
            format!("{:.7}", y_ex),
            format!("{:.2e}", error),
        ]);
    }

    println!("{table}");
}

fn integrate_fredholm<F>(f: F, a: f64, b: f64, n: usize) -> f64
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

fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut mat = a.iter().map(|row| row.clone()).collect::<Vec<_>>();
    let mut vec = b.to_vec();

    for i in 0..n {
        let pivot = mat[i][i];
        for j in i..n {
            mat[i][j] /= pivot;
        }
        vec[i] /= pivot;

        for j in (i + 1)..n {
            let factor = mat[j][i];
            for k in i..n {
                mat[j][k] -= factor * mat[i][k];
            }
            vec[j] -= factor * vec[i];
        }
    }

    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = vec[i];
        for j in (i + 1)..n {
            x[i] -= mat[i][j] * x[j];
        }
    }

    x
}

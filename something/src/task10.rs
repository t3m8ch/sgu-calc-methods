use comfy_table::Table;

pub fn print_task_10(v: f64, h: f64) {
    let n = (v / h) as usize;

    let y_exact = |x: f64| v * x.powi(2) * (x - v);
    let f =
        |x: f64| 4.0 * v * x.powi(4) - 3.0 * v.powi(2) * x.powi(3) + 6.0 * v * x - 2.0 * v.powi(2);
    let p = |x: f64| x.powi(2);
    let q = |x: f64| x;

    let phi_k = |x: f64, k: i32| x.powi(k) * (x - v);
    let dphi_k = |x: f64, k: i32| (k + 1) as f64 * x.powi(k) - v * k as f64 * x.powi(k - 1);
    let ddphi_k = |x: f64, k: i32| {
        k as f64 * (k + 1) as f64 * x.powi(k - 1) - v * k as f64 * (k - 1) as f64 * x.powi(k - 2)
    };

    let xk: Vec<f64> = (0..=n).map(|i| i as f64 * h).collect();

    println!("Проверка точного решения:");
    let mut check_table = Table::new();
    check_table.set_header(vec!["x", "y_toch"]);

    for &x in &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0] {
        if x <= v {
            check_table.add_row(vec![format!("{}", x), format!("{:.2}", y_exact(x))]);
        }
    }
    println!("{check_table}");

    let mut a = vec![vec![0.0; n]; n];
    let mut b = vec![0.0; n];

    for i in 1..=n {
        b[i - 1] = f(xk[i]);
        for k in 1..=n {
            a[i - 1][k - 1] = ddphi_k(xk[i], k as i32)
                + p(xk[i]) * dphi_k(xk[i], k as i32)
                + q(xk[i]) * phi_k(xk[i], k as i32);
        }
    }

    println!("\nРазмерность: A = {}x{}, b = {}", n, n, n);

    let c = gauss_with_pivot(&mut a, &mut b);

    println!("\nКоэффициенты a_k:");
    let mut coef_table = Table::new();
    coef_table.set_header(vec!["k", "a_k"]);

    for i in 0..c.len().min(10) {
        coef_table.add_row(vec![format!("{}", i + 1), format!("{:.6e}", c[i])]);
    }
    println!("{coef_table}");

    let y_approx = |x: f64, coeff: &[f64]| -> f64 {
        coeff
            .iter()
            .enumerate()
            .map(|(i, &ci)| ci * phi_k(x, (i + 1) as i32))
            .sum()
    };

    println!("\nСравнение решений:");
    let mut result_table = Table::new();
    result_table.set_header(vec!["x", "Точное y", "Приближенное y", "Относ. погрешн."]);

    for i in 0..=(v as usize) {
        let x = i as f64;
        let y_ex = y_exact(x);
        let y_ap = y_approx(x, &c);
        let rel_error = if y_ex.abs() > 1e-12 {
            ((y_ap - y_ex) / y_ex).abs()
        } else {
            y_ap.abs()
        };

        result_table.add_row(vec![
            format!("{}", x),
            format!("{:.4}", y_ex),
            format!("{:.4}", y_ap),
            format!("{:.2e}", rel_error),
        ]);
    }
    println!("{result_table}");
}

fn gauss_with_pivot(a: &mut [Vec<f64>], b: &mut [f64]) -> Vec<f64> {
    let n = b.len();

    for i in 0..n {
        let mut max_row = i;
        let mut max_val = a[i][i].abs();

        for j in (i + 1)..n {
            if a[j][i].abs() > max_val {
                max_val = a[j][i].abs();
                max_row = j;
            }
        }

        if max_row != i {
            a.swap(i, max_row);
            b.swap(i, max_row);
        }

        if a[i][i].abs() < 1e-12 {
            a[i][i] = 1e-12;
        }

        let pivot = a[i][i];
        for j in i..n {
            a[i][j] /= pivot;
        }
        b[i] /= pivot;

        for j in (i + 1)..n {
            let factor = a[j][i];
            for k in i..n {
                a[j][k] -= factor * a[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    let mut c = vec![0.0; n];
    for i in (0..n).rev() {
        c[i] = b[i];
        for j in (i + 1)..n {
            c[i] -= a[i][j] * c[j];
        }
    }

    c
}

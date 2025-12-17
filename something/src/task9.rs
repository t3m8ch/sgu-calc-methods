use comfy_table::Table;

pub fn print_task_9(v: f64) {
    let n = 10;
    let x0 = 0.0;
    let h = v / n as f64;

    let derivative = |x: f64| {
        -(4.0 * v * x.powi(4) - 3.0 * v.powi(2) * x.powi(3) + 6.0 * v * x - 2.0 * v.powi(2))
    };
    let y_exact = |x: f64| v * x.powi(2) * (x - v);
    let p = |x: f64| -x.powi(2);
    let q = |x: f64| -x;

    let x: Vec<f64> = (0..=n).map(|i| x0 + i as f64 * h).collect();
    let exact: Vec<f64> = x.iter().map(|&xi| y_exact(xi)).collect();

    let mut f = vec![0.0; n + 1];
    let mut s = vec![0.0; n + 1];
    let mut t = vec![0.0; n + 1];
    let mut r = vec![0.0; n + 1];
    let mut f1 = vec![0.0; n + 1];
    let mut s1 = vec![0.0; n + 1];
    let mut y = vec![0.0; n + 1];

    for i in 1..n {
        f[i] = 0.5 * (1.0 + 0.5 * h * p(x[i]));
        s[i] = 0.5 * (1.0 - 0.5 * h * p(x[i]));
        t[i] = 1.0 + 0.5 * h.powi(2) * q(x[i]);
        r[i] = 0.5 * h.powi(2) * derivative(x[i]);
    }

    f1[1] = 0.0;
    s1[1] = 0.0;

    for j in 1..n {
        let denom = t[j] - f[j] * f1[j];
        f1[j + 1] = s[j] / denom;
        s1[j + 1] = (r[j] + f[j] * s1[j]) / denom;
    }

    y[n] = 0.0;
    for j in (1..n).rev() {
        y[j] = f1[j + 1] * y[j + 1] + s1[j + 1];
    }

    let errors: Vec<f64> = y
        .iter()
        .zip(exact.iter())
        .map(|(yi, ei)| (yi - ei).abs())
        .collect();

    let (max_e, max_e_index) = errors
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();

    let mut table = Table::new();
    table.set_header(vec!["x", "y", "exact", "e"]);

    for i in 0..=n {
        table.add_row(vec![
            format!("{:.2}", x[i]),
            format!("{:.6}", y[i]),
            format!("{:.6}", exact[i]),
            format!("{:.8}", errors[i]),
        ]);
    }

    println!("{table}");
    println!("Максимальный e: {:.8}, Номер: {}", max_e, max_e_index);
}

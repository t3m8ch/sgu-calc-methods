fn f(x: f64, y: f64, v: f64) -> f64 {
    2.0 * v * x + v * x.powi(2) - y
}

fn exact_solution(x: f64, v: f64) -> f64 {
    v * x.powi(2)
}

fn euler_method(x0: f64, y0: f64, h: f64, n: usize, v: f64) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n + 1];
    let mut y = vec![0.0; n + 1];

    x[0] = x0;
    y[0] = y0;

    for i in 0..n {
        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + h * f(x[i], y[i], v);
    }

    (x, y)
}

pub fn improved_euler_method(x0: f64, y0: f64, h: f64, n: usize, v: f64) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n + 1];
    let mut y = vec![0.0; n + 1];

    x[0] = x0;
    y[0] = y0;

    for i in 0..n {
        x[i + 1] = x[i] + h;

        let x_half = x[i] + h / 2.0;
        let y_half = y[i] + (h / 2.0) * f(x[i], y[i], v);

        y[i + 1] = y[i] + h * f(x_half, y_half, v);
    }

    (x, y)
}

pub fn print_task_8(x0: f64, y0: f64, h: f64, n: usize, v: f64) {
    let (x_euler, y_euler) = euler_method(x0, y0, h, n, v);
    let (x_improved, y_improved) = improved_euler_method(x0, y0, h, n, v);

    let y_exact: Vec<f64> = x_euler.iter().map(|&x| exact_solution(x, v)).collect();

    let error_euler: Vec<f64> = y_euler
        .iter()
        .zip(y_exact.iter())
        .map(|(ym, yt)| (ym - yt).abs())
        .collect();

    let error_improved: Vec<f64> = y_improved
        .iter()
        .zip(y_exact.iter())
        .map(|(ym, yt)| (ym - yt).abs())
        .collect();

    println!("\nМетод Эйлера:");
    println!("{}", "-".repeat(130));

    print!("x:       ");
    x_euler.iter().for_each(|&x| print!(" {:>10.7}", x));
    println!();

    print!("y_M:     ");
    y_euler.iter().for_each(|&y| print!(" {:>10.7}", y));
    println!();

    print!("y_T:     ");
    y_exact.iter().for_each(|&y| print!(" {:>10.7}", y));
    println!();

    print!("Погрешн: ");
    error_euler.iter().for_each(|&e| print!(" {:>10.7}", e));
    println!();

    println!("{}", "-".repeat(130));

    println!("\nУсовершенствованный метод Эйлера:");
    println!("{}", "-".repeat(130));

    print!("x:       ");
    x_improved.iter().for_each(|&x| print!(" {:>10.7}", x));
    println!();

    print!("y_M:     ");
    y_improved.iter().for_each(|&y| print!(" {:>10.7}", y));
    println!();

    print!("y_T:     ");
    y_exact.iter().for_each(|&y| print!(" {:>10.7}", y));
    println!();

    print!("Погрешн: ");
    error_improved.iter().for_each(|&e| print!(" {:>10.7}", e));
    println!();

    println!("{}", "-".repeat(130));
}

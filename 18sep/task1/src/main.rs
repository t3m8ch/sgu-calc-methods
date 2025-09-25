use std::io::{self, Write};

use comfy_table::Table;
use nalgebra::{DMatrix, DVector};

fn main() {
    let mut n = String::new();

    print!("Введите число n: ");
    io::stdout().flush().expect("Failed to flush");
    io::stdin().read_line(&mut n).expect("Failed to read line");

    let n: u32 = n.trim().parse().expect("Failed to parse n");

    let mut xv_interpolated = Vec::new();
    let mut fv_interpolated = Vec::new();

    for i in 0..=n {
        let mut x = String::new();

        print!("Введите x{i}: ");
        io::stdout().flush().expect("Failed to flush");
        io::stdin().read_line(&mut x).expect("Failed to read line");
        let x: f64 = x.trim().parse().expect("Failed to parse x");

        xv_interpolated.push(x);
    }

    for i in 0..=n {
        let mut f = String::new();

        print!("Введите f{i}: ");
        io::stdout().flush().expect("Failed to flush");
        io::stdin().read_line(&mut f).expect("Failed to read line");
        let f: f64 = f.trim().parse().expect("Failed to parse f");

        fv_interpolated.push(f);
    }

    println!("\nВходная таблица:");
    print_table(&xv_interpolated, &fv_interpolated);

    let x_interim_vector: Vec<f64> = [xv_interpolated[0]]
        .into_iter()
        .chain(
            xv_interpolated
                .iter()
                .skip(1)
                .zip(xv_interpolated.iter())
                .map(|(x0, x1)| [(x0 + x1) / 2.0, *x0])
                .flatten(),
        )
        .collect();

    println!("x_interim_vector: {:#?}", x_interim_vector);

    let f_vector =
        interpolation_polynomial(n, &xv_interpolated, &fv_interpolated, &x_interim_vector);

    println!("\nВыходная таблица");
    print_table(&x_interim_vector, &f_vector);
}

fn interpolation_polynomial(
    n: u32,
    xv_interpolated: &[f64],
    fv_interpolated: &[f64],
    xv: &[f64],
) -> Vec<f64> {
    let x_matrix: Vec<f64> = (0..=n)
        .into_iter()
        .map(|i| xv_interpolated.iter().map(move |x| x.powi(i as i32)))
        .flatten()
        .collect();
    let x_matrix = DMatrix::from_vec((n + 1) as usize, (n + 1) as usize, x_matrix);
    let f_vector = DVector::from_row_slice(fv_interpolated);

    let Some(a_vector) = x_matrix.lu().solve(&f_vector) else {
        // TODO: Remove panic
        panic!("Система не имеет решений!");
    };

    xv.iter()
        .map(|x| {
            (0..=n)
                .map(|i| x.powi(i as i32) * a_vector[(i as usize, 0)])
                .sum()
        })
        .collect()
}

fn print_table(xs: &[f64], fs: &[f64]) {
    let mut output_table = Table::new();
    output_table.set_header(
        ["x".to_string()]
            .into_iter()
            .chain(xs.iter().map(|x| x.to_string()))
            .collect::<Vec<String>>(),
    );
    output_table.add_row(
        ["f".to_string()]
            .into_iter()
            .chain(fs.iter().map(|f| f.to_string()))
            .collect::<Vec<String>>(),
    );

    println!("{output_table}");
}

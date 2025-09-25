use std::io::{self, Write};

use comfy_table::Table;
use itertools::Itertools;
use nalgebra::{DMatrix, DVector};

fn main() {
    let n = input_n();
    let x_list_input = input_vector(n, "x");
    let f_list_input = input_vector(n, "f");

    println!("\nВходная таблица:");
    print_table(&x_list_input, &f_list_input);

    let x_list_with_midpoints = x_list_with_midpoints(&x_list_input);

    println!("\nИнтерполяционный многолчен P_{n}");
    print_table(
        &x_list_with_midpoints,
        &vandermonde_interpolation(n, &x_list_input, &f_list_input, &x_list_with_midpoints),
    );

    println!("\nИнтерполяционный многолчен в форме Лагранжа l_{n}");
    print_table(
        &x_list_with_midpoints,
        &lagrange_interpolation(n, &x_list_input, &f_list_input, &x_list_with_midpoints),
    );
}

fn x_list_with_midpoints(x_list: &[f64]) -> Vec<f64> {
    x_list
        .iter()
        .tuple_windows()
        .flat_map(|(x0, x1)| [*x0, (x0 + x1) / 2.0])
        .chain(x_list.last().copied())
        .collect()
}

fn vandermonde_interpolation(
    n: u32,
    x_list_input: &[f64],
    f_list_input: &[f64],
    x_list_with_midpoints: &[f64],
) -> Vec<f64> {
    let x_matrix: Vec<f64> = (0..=n)
        .into_iter()
        .map(|i| x_list_input.iter().map(move |x| x.powi(i as i32)))
        .flatten()
        .collect();

    let x_matrix = DMatrix::from_vec((n + 1) as usize, (n + 1) as usize, x_matrix);
    let f_vector = DVector::from_row_slice(f_list_input);
    let a_vector: Vec<f64> = x_matrix
        .lu()
        .solve(&f_vector)
        .expect("Failed to solve")
        .into_iter()
        .map(|a| *a)
        .collect();

    x_list_with_midpoints
        .iter()
        .map(|x| {
            a_vector
                .iter()
                .enumerate()
                .map(|(i, a_i)| x.powi(i as i32) * a_i)
                .sum()
        })
        .collect()
}

fn lagrange_interpolation(
    n: u32,
    xv_interpolated: &[f64],
    fv_interpolated: &[f64],
    xv: &[f64],
) -> Vec<f64> {
    xv.iter()
        .map(|x| {
            fv_interpolated
                .into_iter()
                .zip(xv_interpolated)
                .enumerate()
                .map(|(f_idx, (fk, xk))| {
                    let xv_interpolated: Vec<f64> = xv_interpolated
                        .iter()
                        .enumerate()
                        .filter(|(x_idx, _)| *x_idx != f_idx)
                        .map(|(_, xi)| *xi)
                        .collect();

                    fk * xv_interpolated
                        .iter()
                        .map(|xi| x - xi)
                        .fold(1.0, |acc, m| acc * m)
                        / xv_interpolated
                            .iter()
                            .map(|xi| xk - xi)
                            .fold(1.0, |acc, m| acc * m)
                })
                .sum()
        })
        .collect()
}

fn input_n() -> u32 {
    let mut n = String::new();

    print!("Введите число n: ");
    io::stdout().flush().expect("Failed to flush");
    io::stdin().read_line(&mut n).expect("Failed to read line");

    n.trim().parse().expect("Failed to parse n")
}

fn input_vector(n: u32, letter: &str) -> Vec<f64> {
    let mut result = Vec::new();

    for i in 0..=n {
        let mut x = String::new();

        print!("Введите {letter}_{i}: ");
        io::stdout().flush().expect("Failed to flush");
        io::stdin().read_line(&mut x).expect("Failed to read line");
        let x: f64 = x.trim().parse().expect("Failed to parse x");

        result.push(x);
    }

    result
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

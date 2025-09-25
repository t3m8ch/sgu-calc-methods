use std::io::{self, Write};

use comfy_table::Table;
use nalgebra::{DMatrix, DVector};

fn main() {
    let mut n = String::new();

    print!("Введите число n: ");
    io::stdout().flush().expect("Failed to flush");
    io::stdin().read_line(&mut n).expect("Failed to read line");

    let n: u32 = n.trim().parse().expect("Failed to parse n");

    let mut x_vector = Vec::new();
    let mut f_vector = Vec::new();

    for i in 0..=n {
        let mut x = String::new();

        print!("Введите x{i}: ");
        io::stdout().flush().expect("Failed to flush");
        io::stdin().read_line(&mut x).expect("Failed to read line");
        let x: f64 = x.trim().parse().expect("Failed to parse x");

        x_vector.push(x);
    }

    for i in 0..=n {
        let mut f = String::new();

        print!("Введите f{i}: ");
        io::stdout().flush().expect("Failed to flush");
        io::stdin().read_line(&mut f).expect("Failed to read line");
        let f: f64 = f.trim().parse().expect("Failed to parse f");

        f_vector.push(f);
    }

    println!("\nВходная таблица:");
    print_table(&x_vector, &f_vector);

    let x_matrix: Vec<f64> = (0..=n)
        .into_iter()
        .map(|i| x_vector.iter().map(move |x| x.powi(i as i32)))
        .flatten()
        .collect();

    let x_matrix = DMatrix::from_vec((n + 1) as usize, (n + 1) as usize, x_matrix);
    let f_vector = DVector::from_vec(f_vector);

    let Some(a_vector) = x_matrix.lu().solve(&f_vector) else {
        eprintln!("Система не имеет решений!");
        return;
    };

    let x_vector: Vec<f64> = [x_vector[0]]
        .into_iter()
        .chain(
            x_vector
                .iter()
                .skip(1)
                .zip(x_vector.iter())
                .map(|(x0, x1)| [(x0 + x1) / 2.0, *x0])
                .flatten(),
        )
        .collect();

    let f_vector: Vec<f64> = x_vector
        .iter()
        .map(|x| {
            (0..=n)
                .map(|i| x.powi(i as i32) * a_vector[(i as usize, 0)])
                .sum()
        })
        .collect();

    println!("\nВыходная таблица");
    print_table(&x_vector, &f_vector);
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

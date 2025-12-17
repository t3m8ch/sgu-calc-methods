use std::io::{self, Write};

use comfy_table::Table;

pub fn input_n() -> u32 {
    let mut n = String::new();

    print!("Введите число n: ");
    io::stdout().flush().expect("Failed to flush");
    io::stdin().read_line(&mut n).expect("Failed to read line");

    n.trim().parse().expect("Failed to parse n")
}

pub fn input_vector(n: u32, letter: &str) -> Vec<f64> {
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

pub fn print_table(xs: &[f64], fs: &[f64]) {
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

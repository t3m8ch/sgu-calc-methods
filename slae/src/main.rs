use shared::cli::{input_matrix, input_n, input_vector};

use crate::gauss::solve_slae_gauss;

mod gauss;

fn main() {
    let n = input_n();
    let matrix = input_matrix(n - 1);
    let b = input_vector(n - 1, "b");

    println!("\n=== Метод Гаусса ===");
    println!(
        "Решение методом Гаусса: {:.6?}",
        solve_slae_gauss(&matrix, &b)
    );
}

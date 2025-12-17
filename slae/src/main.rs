use shared::cli::{input_matrix, input_n, input_vector};

use crate::{gauss::solve_slae_gauss, tridiagonal::solve_tridiagonal};

mod gauss;
mod tridiagonal;

fn main() {
    println!("=== Метод Гаусса ===");
    let n = input_n();
    let matrix = input_matrix(n - 1);
    let b = input_vector(n - 1, "b");

    println!("Решение: {:.6?}", solve_slae_gauss(&matrix, &b));

    println!("\n\n=== Метод прогонки решения СЛАУ (трехдиагональных) ===");
    let n = input_n();
    let matrix = input_matrix(n - 1);
    let b = input_vector(n - 1, "b");

    println!("\nРешение: {:.10?}", solve_tridiagonal(&matrix, &b));
}

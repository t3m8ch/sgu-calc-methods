pub fn lagrange_interpolation(
    n: u32,
    x_list_input: &[f64],
    f_list_input: &[f64],
    x_list_with_midpoints: &[f64],
) -> Vec<f64> {
    x_list_with_midpoints
        .iter()
        .map(|x| {
            f_list_input
                .into_iter()
                .zip(x_list_input)
                .enumerate()
                .map(|(f_idx, (fk, xk))| {
                    let x_list_without_k: Vec<f64> = x_list_input
                        .iter()
                        .enumerate()
                        .filter(|(x_idx, _)| *x_idx != f_idx)
                        .map(|(_, xi)| *xi)
                        .collect();

                    let numerator: f64 = x_list_without_k.iter().map(|xi| x - xi).product();
                    let denominator: f64 = x_list_without_k.iter().map(|xi| xk - xi).product();

                    fk * numerator / denominator
                })
                .sum()
        })
        .collect()
}

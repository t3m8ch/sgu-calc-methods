pub fn lagrange_interpolation(x_nodes: &[f64], f_nodes: &[f64], x_targets: &[f64]) -> Vec<f64> {
    x_targets
        .iter()
        .map(|x| {
            f_nodes
                .into_iter()
                .zip(x_nodes)
                .enumerate()
                .map(|(f_idx, (fk, xk))| {
                    let x_list_without_k: Vec<f64> = x_nodes
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

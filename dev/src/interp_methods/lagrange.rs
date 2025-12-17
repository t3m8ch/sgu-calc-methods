pub fn lagrange_interpolation(x_nodes: &[f64], f_nodes: &[f64], x_targets: &[f64]) -> Vec<f64> {
    x_targets
        .iter()
        .map(|x| {
            f_nodes
                .into_iter()
                .zip(x_nodes)
                .enumerate()
                .map(|(k, (fk, xk))| {
                    let x_list_without_k = x_nodes[..k].iter().chain(&x_nodes[k + 1..]);

                    let numerator: f64 = x_list_without_k.clone().map(|xi| x - xi).product();
                    let denominator: f64 = x_list_without_k.map(|xi| xk - xi).product();

                    fk * numerator / denominator
                })
                .sum()
        })
        .collect()
}

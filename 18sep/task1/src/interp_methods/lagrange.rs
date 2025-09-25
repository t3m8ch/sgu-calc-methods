pub fn lagrange_interpolation(
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

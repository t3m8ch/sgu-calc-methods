pub fn neuton_interpolation(x_nodes: &[f64], f_nodes: &[f64], x_targets: &[f64]) -> Vec<f64> {
    todo!()
}

fn separated_differences(x_nodes: &[f64], f_nodes: &[f64]) -> Vec<Vec<f64>> {
    let n = x_nodes.len() - 1;
    let mut diffs = Vec::with_capacity(n + 1);
    diffs.push(f_nodes.to_vec());

    for l in 1..=n {
        diffs.push(
            (0..=n - l)
                .map(|k| (diffs[l - 1][k + 1] - diffs[l - 1][k]) / (x_nodes[k + l] - x_nodes[k]))
                .collect(),
        );
    }

    diffs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_separated_differences() {
        let x_nodes = vec![5.0, 6.0, 9.0, 11.0];
        let f_nodes = vec![12.0, 13.0, 14.0, 16.0];

        let result = separated_differences(&x_nodes, &f_nodes);

        // Порядок 0
        assert_eq!(result[0], vec![12.0, 13.0, 14.0, 16.0]);

        // Порядок 1
        assert!((result[1][0] - 1.0).abs() < 1e-10);
        assert!((result[1][1] - 1.0 / 3.0).abs() < 1e-10);
        assert!((result[1][2] - 1.0).abs() < 1e-10);

        // Порядок 2
        assert!((result[2][0] - (-1.0 / 6.0)).abs() < 1e-10);
        assert!((result[2][1] - (2.0 / 15.0)).abs() < 1e-10);

        // Порядок 3
        assert!((result[3][0] - 0.05).abs() < 1e-10);

        println!("separated differences {:#?}!", result);
    }
}

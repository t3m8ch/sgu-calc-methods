use itertools::Itertools;

pub fn x_list_with_midpoints(x_list: &[f64]) -> Vec<f64> {
    x_list
        .iter()
        .tuple_windows()
        .flat_map(|(x0, x1)| [*x0, (x0 + x1) / 2.0])
        .chain(x_list.last().copied())
        .collect()
}

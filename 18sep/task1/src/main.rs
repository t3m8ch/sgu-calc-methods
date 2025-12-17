use crate::{
    cli::{input_n, input_vector, print_table},
    interp_methods::{
        lagrange::lagrange_interpolation, newton::newton_interpolation,
        spline::cubic_spline_interpolation, vandermonde::vandermonde_interpolation,
    },
    midpoints::x_list_with_midpoints,
};

mod cli;
mod interp_methods;
mod midpoints;

fn main() {
    let n = input_n();
    let x_list_input = input_vector(n, "x");
    let f_list_input = input_vector(n, "f");

    println!("\nВходная таблица:");
    print_table(&x_list_input, &f_list_input);

    let x_list_with_midpoints = x_list_with_midpoints(&x_list_input);

    println!("\nИнтерполяционный многолчен P_{n}");
    print_table(
        &x_list_with_midpoints,
        &vandermonde_interpolation(&x_list_input, &f_list_input, &x_list_with_midpoints),
    );

    println!("\nИнтерполяционный многолчен в форме Лагранжа l_{n}");
    print_table(
        &x_list_with_midpoints,
        &lagrange_interpolation(&x_list_input, &f_list_input, &x_list_with_midpoints),
    );

    println!("\nИнтерполяционный многолчен в форме Ньютона N_{n}");
    print_table(
        &x_list_with_midpoints,
        &newton_interpolation(&x_list_input, &f_list_input, &x_list_with_midpoints),
    );

    println!("\nИнтерполяция кубическими сплайнами");
    print_table(
        &x_list_with_midpoints,
        &cubic_spline_interpolation(&x_list_input, &f_list_input, &x_list_with_midpoints),
    );
}

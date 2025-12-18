#import "conf.typ" : conf
#show: conf.with(
  title: [Отчёт\ по практической подготовке],
  type: "pract",
  info: (
      author: (
        name: [Кудяков Артём Александрович],
        faculty: [компьютерных наук и информационных технологий],
        group: "351",
        sex: "male"
      ),
      inspector: (
        degree: "",
        name: ""
      )
  ),
  settings: (
    title_page: (
      enabled: true
    ),
    contents_page: (
      enabled: true
    )
  )
)

= Построение интерполяционного многочлена в общем виде

*Условие*

Необходимо найти интерполяционный многочлен в общем виде.

#table(rows: 2, columns: (1fr, 1fr, 1fr, 1fr, 1fr))[*$x$*][
  0][1][2][3][*$f(x)$*][1][2][9][28]

*Код*

```rust
use nalgebra::{DMatrix, DVector};

pub fn vandermonde_interpolation(x_nodes: &[f64],
                                 f_nodes: &[f64],
                                 x_targets: &[f64]) -> Vec<f64> {
    let n = x_nodes.len() - 1;
    let a_list = solve_matrix(n, x_nodes, f_nodes);
    x_targets
        .iter()
        .map(|x| {
            a_list
                .iter()
                .enumerate()
                .map(|(i, a_i)| x.powi(i as i32) * a_i)
                .sum()
        })
        .collect()
}

fn solve_matrix(n: usize, x_list: &[f64], f_list: &[f64]) -> Vec<f64> {
    let x_matrix: Vec<f64> = (0..=n)
        .into_iter()
        .map(|i| x_list.iter().map(move |x| x.powi(i as i32)))
        .flatten()
        .collect();

    let x_matrix = DMatrix::from_vec((n + 1) as usize, (n + 1) as usize, x_matrix);
    let f_vector = DVector::from_row_slice(f_list);

    x_matrix
        .lu()
        .solve(&f_vector)
        .expect("Failed to solve")
        .into_iter()
        .map(|a| *a)
        .collect()
}
```

*Результат*

```
Интерполяционный многолчен P_3
+---+---+-------+---+-------+---+--------+----+
| x | 0 | 0.5   | 1 | 1.5   | 2 | 2.5    | 3  |
+=============================================+
| f | 1 | 1.125 | 2 | 4.375 | 9 | 16.625 | 28 |
+---+---+-------+---+-------+---+--------+----+
```


= Интерполяционный многочлен в форме Лагранжа

*Условие*

По данным интерполяции из предыдущего задания построить интерполяционный многочлен
в форме Лагранжа.

#table(rows: 2, columns: (1fr, 1fr, 1fr, 1fr, 1fr))[*$x$*][
  0][1][2][3][*$f(x)$*][1][2][9][28]

*Код*

```rust
pub fn lagrange_interpolation(x_nodes: &[f64],
                              f_nodes: &[f64],
                              x_targets: &[f64]) -> Vec<f64> {
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
```

*Результат*

```
Интерполяционный многолчен в форме Лагранжа l_3
+---+---+-------+---+-------+---+--------+----+
| x | 0 | 0.5   | 1 | 1.5   | 2 | 2.5    | 3  |
+=============================================+
| f | 1 | 1.125 | 2 | 4.375 | 9 | 16.625 | 28 |
+---+---+-------+---+-------+---+--------+----+
```



= Интерполяционный многочлен в форме Ньютона

*Условие*

По данным интерполяции из предыдущего задания построить интерполяционный многочлен
в форме Ньютона.

#table(rows: 2, columns: (1fr, 1fr, 1fr, 1fr, 1fr))[*$x$*][
  0][1][2][3][*$f(x)$*][1][2][9][28]

*Код*

```rust
pub fn newton_interpolation(x_nodes: &[f64], f_nodes: &[f64], x_targets: &[f64]) -> Vec<f64> {
    let diffs = separated_differences(x_nodes, f_nodes);
    x_targets
        .iter()
        .map(|x| {
            f_nodes[0]
                + (1..x_nodes.len())
                    .map(|i| diffs[i][0] * (0..=i - 1).map(|j| (x - x_nodes[j])).product::<f64>())
                    .sum::<f64>()
        })
        .collect()
}

fn separated_differences(x_nodes: &[f64], f_nodes: &[f64]) -> Vec<Vec<f64>> {
    let n = x_nodes.len() - 1;
    (1..=n).fold(vec![f_nodes.to_vec()], |mut acc, l| {
        acc.push(
            (0..=n - l)
                .map(|k| (acc[l - 1][k + 1] - acc[l - 1][k]) / (x_nodes[k + l] - x_nodes[k]))
                .collect(),
        );
        acc
    })
}
```

*Результат*

```
Интерполяционный многолчен в форме Ньютона N_3
+---+---+-------+---+-------+---+--------+----+
| x | 0 | 0.5   | 1 | 1.5   | 2 | 2.5    | 3  |
+=============================================+
| f | 1 | 1.125 | 2 | 4.375 | 9 | 16.625 | 28 |
+---+---+-------+---+-------+---+--------+----+
```


= Интерполяция кубическими сплайнами

*Условие*

Необходимо построить интерполяционный многочлен с помощью кубических сплайнов
(алгебраических многочленов третьей степени, где сплайн --- фрагмент, отрезок чего-либо).

#table(rows: 2, columns: (1fr, 1fr, 1fr, 1fr, 1fr))[*$x$*][
  0][1][2][3][*$f(x)$*][1][2][9][28]

*Код*

```rust
use nalgebra::{DMatrix, DVector};

pub fn cubic_spline_interpolation(x_nodes: &[f64], f_nodes: &[f64], x_targets: &[f64]) -> Vec<f64> {
    let (matrix_data, rhs_data) = build_spline_system(x_nodes, f_nodes);
    let n = x_nodes.len() - 1;
    let size = 4 * n;

    println!("Матрица системы ({}x{}):", size, size);
    print_matrix(&matrix_data, size, size);

    println!("\nВектор правой части:");
    println!("{:?}\n", rhs_data);

    let coefficients = solve_spline_coefficients(x_nodes, f_nodes);

    println!("Коэффициенты сплайнов [a_i, b_i, c_i, d_i] для каждого интервала:");
    for (i, coef) in coefficients.iter().enumerate() {
        println!(
            "  Интервал {}: a={:.6}, b={:.6}, c={:.6}, d={:.6}",
            i, coef[0], coef[1], coef[2], coef[3]
        );
    }
    println!();

    x_targets
        .iter()
        .map(|&x| evaluate_spline(x, x_nodes, &coefficients))
        .collect()
}

fn solve_spline_coefficients(x_nodes: &[f64], f_nodes: &[f64]) -> Vec<[f64; 4]> {
    let n = x_nodes.len() - 1;
    let size = 4 * n;
    let (matrix, rhs) = build_spline_system(x_nodes, f_nodes);

    let matrix = DMatrix::from_row_slice(size, size, &matrix);
    let rhs = DVector::from_row_slice(&rhs);

    matrix
        .lu()
        .solve(&rhs)
        .expect("Failed to solve spline system")
        .as_slice()
        .chunks(4)
        .map(|chunk| [chunk[0], chunk[1], chunk[2], chunk[3]])
        .collect()
}

fn build_spline_system(x_nodes: &[f64], f_nodes: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = x_nodes.len() - 1;
    let size = 4 * n;
    let mut matrix = vec![0.0; size * size];
    let mut rhs = vec![0.0; size];

    let h: Vec<f64> = (0..n).map(|i| x_nodes[i + 1] - x_nodes[i]).collect();

    (0..n).for_each(|i| {
        matrix[i * size + 4 * i] = 1.0;
        rhs[i] = f_nodes[i];
    });

    (0..n).for_each(|i| {
        let row = n + i;
        matrix[row * size + 4 * i] = 1.0;
        matrix[row * size + 4 * i + 1] = h[i];
        matrix[row * size + 4 * i + 2] = h[i].powi(2);
        matrix[row * size + 4 * i + 3] = h[i].powi(3);
        rhs[row] = f_nodes[i + 1];
    });

    (0..n - 1).for_each(|i| {
        let row = 2 * n + i;
        matrix[row * size + 4 * i + 1] = 1.0;
        matrix[row * size + 4 * i + 2] = 2.0 * h[i];
        matrix[row * size + 4 * i + 3] = 3.0 * h[i].powi(2);
        matrix[row * size + 4 * (i + 1) + 1] = -1.0;
    });

    (0..n - 1).for_each(|i| {
        let row = 3 * n - 1 + i;
        matrix[row * size + 4 * i + 2] = 2.0;
        matrix[row * size + 4 * i + 3] = 6.0 * h[i];
        matrix[row * size + 4 * (i + 1) + 2] = -2.0;
    });

    let row1 = 4 * n - 2;
    matrix[row1 * size + 2] = 2.0;

    let row2 = 4 * n - 1;
    matrix[row2 * size + 4 * (n - 1) + 2] = 2.0;
    matrix[row2 * size + 4 * (n - 1) + 3] = 6.0 * h[n - 1];

    (matrix, rhs)
}

fn evaluate_spline(x: f64, x_nodes: &[f64], coefficients: &[[f64; 4]]) -> f64 {
    let interval = x_nodes
        .windows(2)
        .position(|w| x >= w[0] && x <= w[1])
        .unwrap_or(coefficients.len() - 1);

    let [a, b, c, d] = coefficients[interval];
    let dx = x - x_nodes[interval];

    a + b * dx + c * dx.powi(2) + d * dx.powi(3)
}
```

*Результат*

```
Интерполяция кубическими сплайнами
Матрица системы (12x12):
  [   1.00    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0 ]
  [    0.0    0.0    0.0    0.0   1.00    0.0    0.0    0.0    0.0    0.0    0.0    0.0 ]
  [    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0   1.00    0.0    0.0    0.0 ]
  [   1.00   1.00   1.00   1.00    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0 ]
  [    0.0    0.0    0.0    0.0   1.00   1.00   1.00   1.00    0.0    0.0    0.0    0.0 ]
  [    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0   1.00   1.00   1.00   1.00 ]
  [    0.0   1.00   2.00   3.00    0.0  -1.00    0.0    0.0    0.0    0.0    0.0    0.0 ]
  [    0.0    0.0    0.0    0.0    0.0   1.00   2.00   3.00    0.0  -1.00    0.0    0.0 ]
  [    0.0    0.0   2.00   6.00    0.0    0.0  -2.00    0.0    0.0    0.0    0.0    0.0 ]
  [    0.0    0.0    0.0    0.0    0.0    0.0   2.00   6.00    0.0    0.0  -2.00    0.0 ]
  [    0.0    0.0   2.00    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0 ]
  [    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0   2.00   6.00 ]

Вектор правой части:
[1.0, 2.0, 9.0, 2.0, 9.0, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

Коэффициенты сплайнов [a_i, b_i, c_i, d_i] для каждого интервала:
  Интервал 0: a=1.000000, b=0.200000, c=0.000000, d=0.800000
  Интервал 1: a=2.000000, b=2.600000, c=2.400000, d=2.000000
  Интервал 2: a=9.000000, b=13.400000, c=8.400000, d=-2.800000

+---+---+--------------------+---+--------------------+---+-------+----+
| x | 0 | 0.5                | 1 | 1.5                | 2 | 2.5   | 3  |
+======================================================================+
| f | 1 | 1.1999999999999997 | 2 | 4.1499999999999995 | 9 | 17.45 | 28 |
+---+---+--------------------+---+--------------------+---+-------+----+
```


= Метод Гаусса решения СЛАУ

*Условие*

Решить следующую СЛАУ методом Гаусса
Метод Гаусса должен решать уравнения вида $A x = b$, где $A$ - матрица.
Для упрощения тестирования матрица $А$ примет вид:

$
  A = mat(10, 0.10, 0.10, 0.10, 0.10;
      0.11, 11, 0.11, 0.11, 0.11;
      0.12, 0.12, 12, 0.12, 0.12;
      0.13, 0.13, 0.13, 13, 0.13;
      0.14, 0.14, 0.14, 0.14, 14) quad
  b = mat(10;11;12;13;14).
$

*Код*

```rust
use shared::cli::print_matrix;

pub fn solve_slae_gauss(matrix: &[f64], b: &[f64]) -> Vec<f64> {
    let n = b.len();

    let mut a = matrix.to_vec();
    let mut b_work = b.to_vec();

    println!("1) Матрица A:");
    print_matrix(&a, n, n);

    let det = calculate_determinant(&a, n);
    println!("\nОпределитель матрицы A = {:.6}", det);

    println!("\nКолонка b: {:?}", b);

    forward_elimination(&mut a, &mut b_work, n);

    println!("\nМатрица после прямого прохода:");
    print_matrix(&a, n, n);

    println!("\nВектор b после прямого прохода: {:.6?}", b_work);

    let solution = backward_substitution(&a, &b_work, n);

    solution
}

fn forward_elimination(a: &mut [f64], b: &mut [f64], n: usize) {
    for k in 0..n - 1 {
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];

            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            a[i * n + k] = 0.0;

            b[i] -= factor * b[k];
        }
    }
}

fn backward_substitution(u: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    (0..n).rev().fold(vec![0.0; n], |mut x, i| {
        let sum: f64 = (i + 1..n).map(|j| u[i * n + j] * x[j]).sum();

        x[i] = (b[i] - sum) / u[i * n + i];
        x
    })
}

fn calculate_determinant(matrix: &[f64], n: usize) -> f64 {
    let mut a = matrix.to_vec();
    let mut det = 1.0;

    for k in 0..n - 1 {
        if a[k * n + k].abs() < 1e-10 {
            return 0.0;
        }

        det *= a[k * n + k];

        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
        }
    }

    det * a[(n - 1) * n + (n - 1)]
}
```

*Результат*

```
Матрица A:
  [  10.00   0.10   0.10   0.10   0.10 ]
  [   0.11  11.00   0.11   0.11   0.11 ]
  [   0.12   0.12  12.00   0.12   0.12 ]
  [   0.13   0.13   0.13  13.00   0.13 ]
  [   0.14   0.14   0.14   0.14  14.00 ]

Определитель матрицы A = 240004.528860

Колонка b: [10.0, 11.0, 12.0, 13.0, 14.0]

Матрица после прямого прохода:
  [  10.00   0.10   0.10   0.10   0.10 ]
  [    0.0  11.00   0.11   0.11   0.11 ]
  [    0.0    0.0  12.00   0.12   0.12 ]
  [    0.0    0.0    0.0  13.00   0.13 ]
  [    0.0    0.0    0.0    0.0  13.99 ]

Вектор b после прямого прохода: [10.000000, 10.890000, 11.762376, 12.617647, 13.456311]
Решение: [0.961538, 0.961538, 0.961538, 0.961538, 0.961538]
```


=	Метод прогонки решения СЛАУ (трехдиагональных)

*Условие*
В данном случае решается система линейных уравнений вида $A x = b$, где $A$ --- матрица вида:

$
    mat(
      -10, 0.10, 0, 0, 0;
      0.10, -11, 0.11, 0, 0;
      0, 0.11, -12, 0.12, 0;
      0, 0, 0.12, -13, 0.13;
      0, 0, 0, 0.13, -14) x = mat(10; 11; 12; 13; 14).
$

*Код*

```rust
use shared::cli::print_matrix;

pub fn solve_tridiagonal(matrix: &[f64], b: &[f64]) -> Vec<f64> {
    let n = b.len();

    println!("Матрица A:");
    print_matrix(matrix, n, n);
    println!("\nСтолбец b: {:?}", b);

    let (p_coeffs, q_coeffs) = forward_sweep(matrix, b, n);

    println!("\nПрямая прогонка");
    println!("Список Pi и Qi");
    print!("P: [");
    p_coeffs.iter().for_each(|&val| print!(" {:.10}", val));
    println!(" ]");

    print!("Q: [");
    q_coeffs.iter().for_each(|&val| print!(" {:.10}", val));
    println!(" ]");

    let solution = backward_sweep(&p_coeffs, &q_coeffs, n);

    println!("\nОбратная прогонка");
    println!("x {} = {:.10}", n, solution[n - 1]);
    (0..n - 1).rev().for_each(|i| {
        println!("x {} = {:.10}", i + 1, solution[i]);
    });

    solution
}

fn forward_sweep(matrix: &[f64], b: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut p = vec![0.0; n - 1];
    let mut q = vec![0.0; n];

    let a0 = matrix[0 * n + 0];
    let c0 = matrix[0 * n + 1];
    p[0] = -c0 / a0;
    q[0] = b[0] / a0;

    for i in 1..n - 1 {
        let a_i = matrix[i * n + i - 1];
        let b_i = matrix[i * n + i];
        let c_i = matrix[i * n + i + 1];

        p[i] = c_i / (-b_i - a_i * p[i - 1]);
    }

    for i in 1..n {
        let a_i = matrix[i * n + i - 1];
        let b_i = matrix[i * n + i];

        q[i] = (a_i * q[i - 1] - b[i]) / (-b_i - a_i * p[i - 1]);
    }

    (p, q)
}

fn backward_sweep(p: &[f64], q: &[f64], n: usize) -> Vec<f64> {
    let mut x = vec![0.0; n];

    x[n - 1] = q[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    x
}
```

*Результат*

```
Матрица A:
  [ -10.00   0.10    0.0    0.0    0.0 ]
  [   0.10 -11.00   0.11    0.0    0.0 ]
  [    0.0   0.11 -12.00   0.12    0.0 ]
  [    0.0    0.0   0.12 -13.00   0.13 ]
  [    0.0    0.0    0.0   0.13 -14.00 ]

Столбец b: [10.0, 11.0, 12.0, 13.0, 14.0]

Прямая прогонка
Список Pi и Qi
P: [ 0.0100000000 0.0100009092 0.0100009168 0.0100009232 ]
Q: [ -1.0000000000 -1.0091826530 -1.0093433725 -1.0094102006 -1.0094668396 ]

Обратная прогонка
x 5 = -1.0094668396
x 4 = -1.0195058010
x 3 = -1.0195393653
x 2 = -1.0193789736
x 1 = -1.0101937897

Решение: [-1.0101937897, -1.0193789736, -1.0195393653, -1.0195058010, -1.0094668396]
```

= Метод простой итерации

*Условие*

При решении СЛАУ вида $A x = b$, где $A$ --- квадратная матрица, мы можем преобразовать ее к эквивалентному виду:

$
  mat(0, - a_12/a_11, ..., -a_(1 n) / a_11;
      -a_21/a_22, 0, ..., -a_(2 n)/a_22;
      dots.v, dots.v, dots.down, dots.v;
    -a_(n 1)/a_(n n), -a_(n 2)/a_(n n), ..., 0
    )
  x = mat(b_1 / a_11; b_2 / a_22; dots.v; b_n / a_(n n)).
$

Таким образом исходная система допускает представление в виде:

$
 alpha x + beta = x,
$

а критерий остановки вычислений:

$
  ||x^(k) - x^(k-1)|| < e.
$

*Код*

```rust
use shared::cli::print_matrix;

use crate::gauss::calculate_determinant;

pub fn solve_iterative(matrix: &[f64], b: &[f64], epsilon: f64, max_iterations: usize) -> Vec<f64> {
    let n = b.len();

    println!("Матрица A:");
    print_matrix(matrix, n, n);

    let det = calculate_determinant(matrix, n);
    println!("\nОпределитель матрицы A: {:.6}", det);

    println!("\nСтолбец b: {:.2?}", b);

    let alpha = compute_alpha(matrix, n);
    println!("\nМатрица alpha:");
    print_matrix(&alpha, n, n);

    let beta = compute_beta(matrix, b, n);
    println!("\nСтолбец beta: {:.10?}", beta);

    println!("\nСчитаем до точности epsilon = {:.0e}", epsilon);

    let mut xk = vec![0.0; n];
    println!(
        "x^(0) = [ {} ]",
        xk.iter()
            .map(|_| "0.0".to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );

    for iteration in 0..max_iterations {
        let xkp1 = compute_next_iteration(&alpha, &xk, &beta, n);

        print!("x^({}) = [ ", iteration + 1);
        xkp1.iter().for_each(|&val| print!("{:.10} ", val));
        println!("]");

        if norm_stop(&xk, &xkp1, epsilon) {
            println!("\nСходимость достигнута за {} итераций", iteration + 1);
            return xkp1;
        }

        xk = xkp1;
    }

    println!(
        "\nДостигнуто максимальное количество итераций: {}",
        max_iterations
    );
    xk
}

fn compute_alpha(matrix: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .flat_map(|i| {
            (0..n).map(move |j| {
                if i == j {
                    0.0
                } else {
                    -matrix[i * n + j] / matrix[i * n + i]
                }
            })
        })
        .collect()
}

fn compute_beta(matrix: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    (0..n).map(|i| b[i] / matrix[i * n + i]).collect()
}

fn compute_next_iteration(alpha: &[f64], xk: &[f64], beta: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let sum: f64 = (0..n).map(|j| alpha[i * n + j] * xk[j]).sum();
            sum + beta[i]
        })
        .collect()
}

fn norm_stop(xk: &[f64], xkp1: &[f64], epsilon: f64) -> bool {
    xk.iter()
        .zip(xkp1.iter())
        .map(|(x, xp)| (x - xp).abs())
        .fold(0.0, f64::max)
        < epsilon
}
```

*Результат*

```
Матрица A:
  [  10.00   0.10   0.10   0.10   0.10 ]
  [   0.11  11.00   0.11   0.11   0.11 ]
  [   0.12   0.12  12.00   0.12   0.12 ]
  [   0.13   0.13   0.13  13.00   0.13 ]
  [   0.14   0.14   0.14   0.14  14.00 ]

Определитель матрицы A: 240004.528860

Столбец b: [10.00, 11.00, 12.00, 13.00, 14.00]

Матрица alpha:
  [    0.0  -0.01  -0.01  -0.01  -0.01 ]
  [  -0.01    0.0  -0.01  -0.01  -0.01 ]
  [  -0.01  -0.01    0.0  -0.01  -0.01 ]
  [  -0.01  -0.01  -0.01    0.0  -0.01 ]
  [  -0.01  -0.01  -0.01  -0.01    0.0 ]

Столбец beta: [1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000, 1.0000000000]

Считаем до точности epsilon = 1e-8
x^(0) = [ 0.0, 0.0, 0.0, 0.0, 0.0 ]
x^(1) = [ 1.0000000000 1.0000000000 1.0000000000 1.0000000000 1.0000000000 ]
x^(2) = [ 0.9600000000 0.9600000000 0.9600000000 0.9600000000 0.9600000000 ]
x^(3) = [ 0.9616000000 0.9616000000 0.9616000000 0.9616000000 0.9616000000 ]
x^(4) = [ 0.9615360000 0.9615360000 0.9615360000 0.9615360000 0.9615360000 ]
x^(5) = [ 0.9615385600 0.9615385600 0.9615385600 0.9615385600 0.9615385600 ]
x^(6) = [ 0.9615384576 0.9615384576 0.9615384576 0.9615384576 0.9615384576 ]
x^(7) = [ 0.9615384617 0.9615384617 0.9615384617 0.9615384617 0.9615384617 ]

Сходимость достигнута за 7 итераций

Решение: [0.9615384617, 0.9615384617, 0.9615384617, 0.9615384617, 0.9615384617]
```

= Задача Коши методами Эйлера
Решить задачу Коши:
+ методом Эйлера;
+ усовершенствованным методом Эйлера: \
  $𝑦′ = 2 ⋅ 𝑉 ⋅ 𝑥 + 𝑉 ⋅ 𝑥^2 − 𝑦, quad 𝑦(x_0) = 𝑉 ⋅ 𝑥^2$.

*Код*

```rust
fn f(x: f64, y: f64, v: f64) -> f64 {
    2.0 * v * x + v * x.powi(2) - y
}

fn exact_solution(x: f64, v: f64) -> f64 {
    v * x.powi(2)
}

fn euler_method(x0: f64, y0: f64, h: f64, n: usize, v: f64) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n + 1];
    let mut y = vec![0.0; n + 1];

    x[0] = x0;
    y[0] = y0;

    for i in 0..n {
        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + h * f(x[i], y[i], v);
    }

    (x, y)
}

pub fn improved_euler_method(x0: f64,
                             y0: f64,
                             h: f64,
                             n: usize,
                             v: f64) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n + 1];
    let mut y = vec![0.0; n + 1];

    x[0] = x0;
    y[0] = y0;

    for i in 0..n {
        x[i + 1] = x[i] + h;

        let x_half = x[i] + h / 2.0;
        let y_half = y[i] + (h / 2.0) * f(x[i], y[i], v);

        y[i + 1] = y[i] + h * f(x_half, y_half, v);
    }

    (x, y)
}
```

*Результат*

#image("task8_result.png")

= Решить краевую задачу разностным методом.
#align(center)[
  $
  cases(
    y'' + x^2y' + x y = 4 V x^4 - 3 V T x^3 + 6 V x -  2 V T, \
    y'(0) = y(T) = 0, \
    V = T = 5
  )
  $
]

*Код*
```rust
use comfy_table::Table;

pub fn print_task_9(v: f64) {
    let n = 10;
    let x0 = 0.0;
    let h = v / n as f64;

    let derivative = |x: f64| {
        -(4.0 * v * x.powi(4) - 3.0 * v.powi(2) * x.powi(3) + 6.0 * v * x - 2.0 * v.powi(2))
    };
    let y_exact = |x: f64| v * x.powi(2) * (x - v);
    let p = |x: f64| -x.powi(2);
    let q = |x: f64| -x;

    let x: Vec<f64> = (0..=n).map(|i| x0 + i as f64 * h).collect();
    let exact: Vec<f64> = x.iter().map(|&xi| y_exact(xi)).collect();

    let mut f = vec![0.0; n + 1];
    let mut s = vec![0.0; n + 1];
    let mut t = vec![0.0; n + 1];
    let mut r = vec![0.0; n + 1];
    let mut f1 = vec![0.0; n + 1];
    let mut s1 = vec![0.0; n + 1];
    let mut y = vec![0.0; n + 1];

    for i in 1..n {
        f[i] = 0.5 * (1.0 + 0.5 * h * p(x[i]));
        s[i] = 0.5 * (1.0 - 0.5 * h * p(x[i]));
        t[i] = 1.0 + 0.5 * h.powi(2) * q(x[i]);
        r[i] = 0.5 * h.powi(2) * derivative(x[i]);
    }

    f1[1] = 0.0;
    s1[1] = 0.0;

    for j in 1..n {
        let denom = t[j] - f[j] * f1[j];
        f1[j + 1] = s[j] / denom;
        s1[j + 1] = (r[j] + f[j] * s1[j]) / denom;
    }

    y[n] = 0.0;
    for j in (1..n).rev() {
        y[j] = f1[j + 1] * y[j + 1] + s1[j + 1];
    }

    let errors: Vec<f64> = y
        .iter()
        .zip(exact.iter())
        .map(|(yi, ei)| (yi - ei).abs())
        .collect();

    let (max_e, max_e_index) = errors
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();

    let mut table = Table::new();
    table.set_header(vec!["x", "y", "exact", "e"]);

    for i in 0..=n {
        table.add_row(vec![
            format!("{:.2}", x[i]),
            format!("{:.6}", y[i]),
            format!("{:.6}", exact[i]),
            format!("{:.8}", errors[i]),
        ]);
    }

    println!("{table}");
    println!("Максимальный e: {:.8}, Номер: {}", max_e, max_e_index);
}
```

*Результат*

```
+-------+--------------+--------------+--------------+
| x     | y            | exact        | e            |
+====================================================+
| 0.00  | 0.000000     | -0.000000    | 0.00000000   |
|-------+--------------+--------------+--------------|
| 1.00  | 144.844890   | -90.000000   | 234.84488975 |
|-------+--------------+--------------+--------------|
| 2.00  | -170.103407  | -320.000000  | 149.89659317 |
|-------+--------------+--------------+--------------|
| 3.00  | -565.051703  | -630.000000  | 64.94829658  |
|-------+--------------+--------------+--------------|
| 4.00  | -892.783676  | -960.000000  | 67.21632355  |
|-------+--------------+--------------+--------------|
| 5.00  | -1232.199397 | -1250.000000 | 17.80060322  |
|-------+--------------+--------------+--------------|
| 6.00  | -1405.215858 | -1440.000000 | 34.78414156  |
|-------+--------------+--------------+--------------|
| 7.00  | -1480.343490 | -1470.000000 | 10.34349008  |
|-------+--------------+--------------+--------------|
| 8.00  | -1265.131577 | -1280.000000 | 14.86842263  |
|-------+--------------+--------------+--------------|
| 9.00  | -841.813901  | -810.000000  | 31.81390086  |
|-------+--------------+--------------+--------------|
| 10.00 | 0.000000     | 0.000000     | 0.00000000   |
+-------+--------------+--------------+--------------+
Максимальный e: 1, Номер: 234.84488975179295
```

= Краевая задача методом неопределенных коэффициентов
#align(center)[
  $
  cases(
    y'' + x^2y' + x y = 4 V x^4 - 3 V T x^3 + 6 V x -  2 V T, \
    y'(0) = y(T) = 0, \
    V = T = 5
  )
  $
]

*Код*

```rust
use comfy_table::Table;

pub fn print_task_10(v: f64, h: f64) {
    let n = (v / h) as usize;

    let y_exact = |x: f64| v * x.powi(2) * (x - v);
    let f =
        |x: f64| 4.0 * v * x.powi(4) - 3.0 * v.powi(2) * x.powi(3) + 6.0 * v * x - 2.0 * v.powi(2);
    let p = |x: f64| x.powi(2);
    let q = |x: f64| x;

    let phi_k = |x: f64, k: i32| x.powi(k) * (x - v);
    let dphi_k = |x: f64, k: i32| (k + 1) as f64 * x.powi(k) - v * k as f64 * x.powi(k - 1);
    let ddphi_k = |x: f64, k: i32| {
        k as f64 * (k + 1) as f64 * x.powi(k - 1) - v * k as f64 * (k - 1) as f64 * x.powi(k - 2)
    };

    let xk: Vec<f64> = (0..=n).map(|i| i as f64 * h).collect();

    println!("Проверка точного решения:");
    let mut check_table = Table::new();
    check_table.set_header(vec!["x", "y_toch"]);

    for &x in &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0] {
        if x <= v {
            check_table.add_row(vec![format!("{}", x), format!("{:.2}", y_exact(x))]);
        }
    }
    println!("{check_table}");

    let mut a = vec![vec![0.0; n]; n];
    let mut b = vec![0.0; n];

    for i in 1..=n {
        b[i - 1] = f(xk[i]);
        for k in 1..=n {
            a[i - 1][k - 1] = ddphi_k(xk[i], k as i32)
                + p(xk[i]) * dphi_k(xk[i], k as i32)
                + q(xk[i]) * phi_k(xk[i], k as i32);
        }
    }

    println!("\nРазмерность: A = {}x{}, b = {}", n, n, n);

    let c = gauss_with_pivot(&mut a, &mut b);

    println!("\nКоэффициенты a_k:");
    let mut coef_table = Table::new();
    coef_table.set_header(vec!["k", "a_k"]);

    for i in 0..c.len().min(10) {
        coef_table.add_row(vec![format!("{}", i + 1), format!("{:.6e}", c[i])]);
    }
    println!("{coef_table}");

    let y_approx = |x: f64, coeff: &[f64]| -> f64 {
        coeff
            .iter()
            .enumerate()
            .map(|(i, &ci)| ci * phi_k(x, (i + 1) as i32))
            .sum()
    };

    println!("\nСравнение решений:");
    let mut result_table = Table::new();
    result_table.set_header(vec!["x", "Точное y", "Приближенное y", "Относ. погрешн."]);

    for i in 0..=(v as usize) {
        let x = i as f64;
        let y_ex = y_exact(x);
        let y_ap = y_approx(x, &c);
        let rel_error = if y_ex.abs() > 1e-12 {
            ((y_ap - y_ex) / y_ex).abs()
        } else {
            y_ap.abs()
        };

        result_table.add_row(vec![
            format!("{}", x),
            format!("{:.4}", y_ex),
            format!("{:.4}", y_ap),
            format!("{:.2e}", rel_error),
        ]);
    }
    println!("{result_table}");
}

fn gauss_with_pivot(a: &mut [Vec<f64>], b: &mut [f64]) -> Vec<f64> {
    let n = b.len();

    for i in 0..n {
        let mut max_row = i;
        let mut max_val = a[i][i].abs();

        for j in (i + 1)..n {
            if a[j][i].abs() > max_val {
                max_val = a[j][i].abs();
                max_row = j;
            }
        }

        if max_row != i {
            a.swap(i, max_row);
            b.swap(i, max_row);
        }

        if a[i][i].abs() < 1e-12 {
            a[i][i] = 1e-12;
        }

        let pivot = a[i][i];
        for j in i..n {
            a[i][j] /= pivot;
        }
        b[i] /= pivot;

        for j in (i + 1)..n {
            let factor = a[j][i];
            for k in i..n {
                a[j][k] -= factor * a[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    let mut c = vec![0.0; n];
    for i in (0..n).rev() {
        c[i] = b[i];
        for j in (i + 1)..n {
            c[i] -= a[i][j] * c[j];
        }
    }

    c
}
```

*Результат*

```
Проверка точного решения:
+---+----------+
| x | y_toch   |
+==============+
| 0 | -0.00    |
|---+----------|
| 1 | -90.00   |
|---+----------|
| 2 | -320.00  |
|---+----------|
| 3 | -630.00  |
|---+----------|
| 4 | -960.00  |
|---+----------|
| 5 | -1250.00 |
+---+----------+

Размерность: A = 100x100, b = 100

Коэффициенты a_k:
+----+---------------+
| k  | a_k           |
+====================+
| 1  | -6.036771e-11 |
|----+---------------|
| 2  | 1.000000e1    |
|----+---------------|
| 3  | -1.889591e-10 |
|----+---------------|
| 4  | 9.864554e-10  |
|----+---------------|
| 5  | -2.997600e-9  |
|----+---------------|
| 6  | 5.913662e-9   |
|----+---------------|
| 7  | -8.000079e-9  |
|----+---------------|
| 8  | 7.749350e-9   |
|----+---------------|
| 9  | -5.524570e-9  |
|----+---------------|
| 10 | 2.939595e-9   |
+----+---------------+

Сравнение решений:
+----+------------+----------------+-----------------+
| x  | Точное y   | Приближенное y | Относ. погрешн. |
+====================================================+
| 0  | -0.0000    | 0.0000         | 0.00e0          |
|----+------------+----------------+-----------------|
| 1  | -90.0000   | -90.0000       | 5.48e-12        |
|----+------------+----------------+-----------------|
| 2  | -320.0000  | -320.0000      | 1.37e-12        |
|----+------------+----------------+-----------------|
| 3  | -630.0000  | -630.0000      | 4.39e-13        |
|----+------------+----------------+-----------------|
| 4  | -960.0000  | -960.0000      | 2.15e-13        |
|----+------------+----------------+-----------------|
| 5  | -1250.0000 | -1250.0000     | 1.34e-13        |
|----+------------+----------------+-----------------|
| 6  | -1440.0000 | -1440.0000     | 9.68e-14        |
|----+------------+----------------+-----------------|
| 7  | -1470.0000 | -1470.0000     | 8.03e-14        |
|----+------------+----------------+-----------------|
| 8  | -1280.0000 | -1280.0000     | 8.10e-14        |
|----+------------+----------------+-----------------|
| 9  | -810.0000  | -810.0000      | 1.13e-13        |
|----+------------+----------------+-----------------|
| 10 | 0.0000     | 0.0000         | 0.00e0          |
+----+------------+----------------+-----------------+
```

= Решение интегральных уравнений при помощи метода неопределённых коэффициентов
$ y(x) + 1 * integral_0^1 (x t + x^2 t^2 + x^3 t^3) y(t) d t
  = V (4/3 x + 1/4 x^2 + 1/5 x^3) $

= Решение интегральных уравнений при помощи метода квадратур:
$ y(x) + 1 * integral_0^1 (x t + x^2 t^2 + x^3 t^3) y(t) d t
  = V (4/3 x + 1/4 x^2 + 1/5 x^3) $

*Код*

```rust
use comfy_table::Table;

pub fn print_task_11(v: f64) {
    let (x, y_calc, y_exact, error) = fredholm_solver(v, 3);

    let mut table = Table::new();
    table.set_header(vec!["x", "y_мет", "y_точн", "погрешность"]);

    for i in 0..x.len() {
        table.add_row(vec![
            format!("{:7.3}", x[i]),
            format!("{:7.3}", y_calc[i]),
            format!("{:7.3}", y_exact[i]),
            format!("{:7.3}", error[i]),
        ]);
    }

    println!("{table}");
}

fn fredholm_solver(v: f64, rank: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let alpha = build_alpha(rank);
    let gamma = build_gamma(rank, v);

    let mut system_matrix = vec![0.0; rank * rank];
    for i in 0..rank {
        for j in 0..rank {
            system_matrix[i * rank + j] = if i == j { 1.0 } else { 0.0 } + alpha[i * rank + j];
        }
    }

    let coeffs = solve_slae_silent(&system_matrix, &gamma);

    let step = 0.1;
    let x_vals: Vec<f64> = (0..=10).map(|i| i as f64 * step).collect();
    let mut y_num = vec![0.0; x_vals.len()];

    for (i, &x) in x_vals.iter().enumerate() {
        let mut y_val = rhs_func(x, v);
        for j in 0..rank {
            y_val -= coeffs[j] * basis_func(j, x);
        }
        y_num[i] = y_val;
    }

    let y_exact: Vec<f64> = x_vals.iter().map(|&x| v * x).collect();
    let error: Vec<f64> = y_num
        .iter()
        .zip(y_exact.iter())
        .map(|(yn, ye)| (yn - ye).abs())
        .collect();

    (x_vals, y_num, y_exact, error)
}

fn rhs_func(x: f64, v: f64) -> f64 {
    v * (4.0 / 3.0 * x + 0.25 * x.powi(2) + 0.2 * x.powi(3))
}

fn basis_func(index: usize, x: f64) -> f64 {
    match index {
        0 => x,
        1 => x.powi(2),
        2 => x.powi(3),
        _ => panic!("Неверный индекс базисной функции"),
    }
}

fn integrate_trapezoidal<F>(f: F, a: f64, b: f64, n: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));

    for i in 1..n {
        sum += f(a + i as f64 * h);
    }

    sum * h
}

fn build_alpha(size: usize) -> Vec<f64> {
    let mut alpha = vec![0.0; size * size];

    for i in 0..size {
        for j in 0..size {
            alpha[i * size + j] =
                integrate_trapezoidal(|t| basis_func(i, t) * basis_func(j, t), 0.0, 1.0, 1000);
        }
    }

    alpha
}

fn build_gamma(size: usize, v: f64) -> Vec<f64> {
    let mut gamma = vec![0.0; size];

    for i in 0..size {
        gamma[i] = integrate_trapezoidal(|t| rhs_func(t, v) * basis_func(i, t), 0.0, 1.0, 1000);
    }

    gamma
}

pub fn solve_slae_silent(matrix: &[f64], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut a = matrix.to_vec();
    let mut b_work = b.to_vec();

    forward_elimination(&mut a, &mut b_work, n);
    backward_substitution(&a, &b_work, n)
}

fn forward_elimination(a: &mut [f64], b: &mut [f64], n: usize) {
    for k in 0..n - 1 {
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];

            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            a[i * n + k] = 0.0;
            b[i] -= factor * b[k];
        }
    }
}

fn backward_substitution(u: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    (0..n).rev().fold(vec![0.0; n], |mut x, i| {
        let sum: f64 = (i + 1..n).map(|j| u[i * n + j] * x[j]).sum();

        x[i] = (b[i] - sum) / u[i * n + i];
        x
    })
}
```

*Результат*

```
+---------+---------+---------+-------------+
| x       | y_мет   | y_точн  | погрешность |
+===========================================+
|   0.000 |   0.000 |   0.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.100 |   1.000 |   1.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.200 |   2.000 |   2.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.300 |   3.000 |   3.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.400 |   4.000 |   4.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.500 |   5.000 |   5.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.600 |   6.000 |   6.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.700 |   7.000 |   7.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.800 |   8.000 |   8.000 |   0.000     |
|---------+---------+---------+-------------|
|   0.900 |   9.000 |   9.000 |   0.000     |
|---------+---------+---------+-------------|
|   1.000 |  10.000 |  10.000 |   0.000     |
+---------+---------+---------+-------------+
```

*Код*

```rust
use comfy_table::Table;

pub fn print_task_11(v: f64, num_points: usize) {
    let (a, b) = (0.0, 1.0);
    let n = 3;
    let lambda = 1.0;

    let y_exact = |x: f64| v * x;
    let f = |x: f64| v * (4.0 / 3.0 * x + 1.0 / 4.0 * x.powi(2) + 1.0 / 5.0 * x.powi(3));

    let a_funcs: Vec<Box<dyn Fn(f64) -> f64>> = vec![
        Box::new(|x| x),
        Box::new(|x| x.powi(2)),
        Box::new(|x| x.powi(3)),
    ];

    let mut alpha = vec![vec![0.0; n]; n];
    for i in 0..n {
        for k in 0..n {
            alpha[i][k] = 1.0 / ((i + k + 3) as f64);
        }
    }

    let mut gamma = vec![0.0; n];
    for i in 0..n {
        gamma[i] = integrate_fredholm(|x| f(x) * x.powi((i + 1) as i32), a, b, 1000);
    }

    let mut mat_a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            mat_a[i][j] = if i == j { 1.0 } else { 0.0 } + lambda * alpha[j][i];
        }
    }

    let q = solve_linear_system(&mat_a, &gamma);

    let y_numerical = |x: f64| -> f64 {
        let mut result = f(x);
        for i in 0..n {
            result -= lambda * q[i] * a_funcs[i](x);
        }
        result
    };

    let x_test: Vec<f64> = (0..num_points)
        .map(|i| a + (b - a) * i as f64 / (num_points - 1) as f64)
        .collect();

    let mut table = Table::new();
    table.set_header(vec!["x", "y_метода", "y_точн", "eps"]);

    for &x in &x_test {
        let y_num = y_numerical(x);
        let y_ex = y_exact(x);
        let error = (y_num - y_ex).abs();

        table.add_row(vec![
            format!("{:.7}", x),
            format!("{:.7}", y_num),
            format!("{:.7}", y_ex),
            format!("{:.2e}", error),
        ]);
    }

    println!("{table}");
}

fn integrate_fredholm<F>(f: F, a: f64, b: f64, n: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));

    for i in 1..n {
        sum += f(a + i as f64 * h);
    }

    sum * h
}

fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut mat = a.iter().map(|row| row.clone()).collect::<Vec<_>>();
    let mut vec = b.to_vec();

    for i in 0..n {
        let pivot = mat[i][i];
        for j in i..n {
            mat[i][j] /= pivot;
        }
        vec[i] /= pivot;

        for j in (i + 1)..n {
            let factor = mat[j][i];
            for k in i..n {
                mat[j][k] -= factor * mat[i][k];
            }
            vec[j] -= factor * vec[i];
        }
    }

    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = vec[i];
        for j in (i + 1)..n {
            x[i] -= mat[i][j] * x[j];
        }
    }

    x
}
```

*Результат*

```
+-----------+-----------+------------+---------+
| x         | y_метода  | y_точн     | eps     |
+==============================================+
| 0.0000000 | 0.0000000 | 0.0000000  | 0.00e0  |
|-----------+-----------+------------+---------|
| 0.1111111 | 1.1111109 | 1.1111111  | 1.89e-7 |
|-----------+-----------+------------+---------|
| 0.2222222 | 2.2222217 | 2.2222222  | 4.99e-7 |
|-----------+-----------+------------+---------|
| 0.3333333 | 3.3333324 | 3.3333333  | 9.70e-7 |
|-----------+-----------+------------+---------|
| 0.4444444 | 4.4444428 | 4.4444444  | 1.64e-6 |
|-----------+-----------+------------+---------|
| 0.5555556 | 5.5555530 | 5.5555556  | 2.56e-6 |
|-----------+-----------+------------+---------|
| 0.6666667 | 6.6666629 | 6.6666667  | 3.76e-6 |
|-----------+-----------+------------+---------|
| 0.7777778 | 7.7777725 | 7.7777778  | 5.28e-6 |
|-----------+-----------+------------+---------|
| 0.8888889 | 8.8888817 | 8.8888889  | 7.18e-6 |
|-----------+-----------+------------+---------|
| 1.0000000 | 9.9999905 | 10.0000000 | 9.48e-6 |
+-----------+-----------+------------+---------+
```

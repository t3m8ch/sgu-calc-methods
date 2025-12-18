use shared::cli::input_v;

use crate::{
    task8::print_task_8, task9::print_task_9, task10::print_task_10, task12::print_task_12,
};

mod task10;
mod task12;
mod task8;
mod task9;

fn main() {
    let v = input_v();

    println!("\n=== Задача Коши методами Эйлера ===");
    print_task_8(1.0, v, 0.001, 10, v);

    println!("\n=== Решить краевую задачу разностным методом ===");
    print_task_9(v);

    println!("\n=== Краевая задача методом неопределенных коэффициентов ===");
    print_task_10(v, 0.1);

    println!("\n=== Решение интегрального уравнения ===");
    print_task_12(v, 10);
}

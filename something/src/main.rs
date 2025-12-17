use shared::cli::input_v;

use crate::task8::print_task_8;

mod task8;

fn main() {
    let v = input_v();
    print_task_8(1.0, v, 0.001, 10, v);
}

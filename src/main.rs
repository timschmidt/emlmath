use emlmath::{EmlExpr, eval_real, x, y};

fn main() {
    let expr = EmlExpr::add(EmlExpr::sin(x()), EmlExpr::sqrt(y()));
    let value = eval_real(&expr, &[("x", 0.5), ("y", 9.0)]).expect("variables are defined");

    println!("expr  = {expr}");
    println!("nodes = {}", expr.node_count());
    println!("value = {value}");
}

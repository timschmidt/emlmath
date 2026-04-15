use emlmath::{
    ComplexBall, ScientificExpr, eval_ball_scientific, eval_real_scientific, sx, sy,
};

fn main() {
    let scientific = ScientificExpr::add(ScientificExpr::sin(sx()), ScientificExpr::sqrt(sy()));
    let expr = scientific.to_eml();
    let value =
        eval_real_scientific(&scientific, &[("x", 0.5), ("y", 9.0)]).expect("variables are defined");
    let bound = eval_ball_scientific(
        &scientific,
        &[
            ("x", ComplexBall::from_real(0.5, 1e-6)),
            ("y", ComplexBall::from_real(9.0, 1e-6)),
        ],
    )
    .expect("scientific ball inputs stay away from singularities and branch cuts");

    println!("scientific = {scientific:?}");
    println!("expr  = {expr}");
    println!("nodes = {}", expr.node_count());
    println!("value = {value}");
    println!("ball  = {bound}");
}

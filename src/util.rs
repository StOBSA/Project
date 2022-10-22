/// a module with utility functions.
use crate::{OPoint, Point};
use ordered_float::*;

/// turn a Point into an OPoint. *for example to hash it*.
pub fn to_graph(point: crate::Point) -> OPoint {
    (OrderedFloat(point.0), OrderedFloat(point.1))
}

/// turn an OPoint into a Point. *for example to perform calculations*
pub fn to_point(point: OPoint) -> Point {
    (*point.0, *point.1)
}

pub fn is_improvement_by_factor(current_value : f64, new_value : f64, factor : f64) -> bool {
    new_value < (current_value-current_value*factor)
}

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

pub fn is_improvement_by_factor(current_value : f32, new_value : f32, factor : f32) -> bool {
    new_value < (current_value-current_value*factor)
}

pub fn average_from_iterator<I:Iterator<Item=f32> + Clone>(values : I) -> f32 {
    let mut len = 0;
    let mut sum = 0.0;
    for number in values {
        sum += number;
        len += 1;
    }
    sum / len as f32
}
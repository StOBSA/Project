
use crate::{OPoint, Point};
use ordered_float::*;
pub fn to_graph(point: crate::Point) -> OPoint {
    (OrderedFloat(point.0), OrderedFloat(point.1))
}
pub fn to_point(point: OPoint) -> Point {
    (*point.0, *point.1)
}

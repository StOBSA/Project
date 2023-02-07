pub const RADIANS_120_DEGREE: f32 = 2.0 * std::f32::consts::PI / 3.0;

use std::f32::INFINITY;

use itertools::Itertools;

use crate::{Point, EPSILON};

#[derive(Debug, Clone)]
pub struct Bounds {
    pub min_x: f32,
    pub max_x: f32,
    pub min_y: f32,
    pub max_y: f32,
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            min_x: INFINITY,
            max_x: -INFINITY,
            min_y: INFINITY,
            max_y: -INFINITY,
        }
    }
}

pub fn euclidean_distance(a: Point, b: Point) -> f32 {
    ((a.0 - b.0).powf(2.0) + (a.1 - b.1).powf(2.0)).sqrt()
}

pub fn overlap(x1: f32, y1: f32, x2: f32, y2: f32, x3: f32, y3: f32, x4: f32, y4: f32) -> bool {
    !(x2 < x3 || x4 < x1 || y2 < y3 || y4 < y1)
}

pub fn segment_segment_intersection(
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    x3: f32,
    y3: f32,
    x4: f32,
    y4: f32,
    point_overlap: bool,
) -> Option<Point> {
    let denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if denom == 0.0 {
        return None;
    }
    let t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
    let u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / denom;

    let test = if point_overlap {
        |t, u| 0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0
    } else {
        |t, u| 0.0 < t && t < 1.0 && 0.0 < u && u < 1.0
    };
    if test(t, u) {
        let p = (x1 + t * (x2 - x1), y1 + t * (y2 - y1));
        if point_overlap {
            return Some(p);
        } else {
            if p.0 != x3 && p.1 != y3 && p.0 != x4 && p.1 != y4 {
                return Some(p);
            }
        }
    }
    None
}

pub fn segment_polygon_intersection(
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    polygon: &[Point],
    point_overlap: bool,
) -> Vec<Point> {
    let mut result = Vec::new();
    for i in -1..(polygon.len() - 1) as i32 {
        let (x3, y3) = polygon[if i == -1 {
            polygon.len() - 1
        } else {
            i as usize
        }];
        let (x4, y4) = polygon[(i + 1) as usize];
        let intersection =
            segment_segment_intersection(x1, y1, x2, y2, x3, y3, x4, y4, point_overlap);
        if intersection.is_some() {
            let new_intersection = intersection.unwrap();
            let mut add_new_intersection = true;
            for intersection in result.iter() {
                if euclidean_distance(new_intersection, *intersection) < EPSILON {
                    add_new_intersection = false;
                    break;
                }
            }
            if add_new_intersection {
                result.push(new_intersection);
            }
        }
    }
    let mut indices_to_delete = Vec::new();
    for index in 0..result.len() {
        if euclidean_distance(result[index], (x1,y1)) < EPSILON {
            indices_to_delete.push(index);
            continue;
        }
        if euclidean_distance(result[index], (x2,y2)) < EPSILON {
            indices_to_delete.push(index);
            continue;
        }
    }
    indices_to_delete.reverse();
    for index in indices_to_delete {
        result.remove(index);
    }
    // if result.contains(&(x1, y1)) {
    //     result.remove(result.iter().position(|i| i.0 == x1 && i.1 == y1).unwrap());
    // }
    // if result.contains(&(x2, y2)) {
    //     result.remove(result.iter().position(|i| i.0 == x2 && i.1 == y2).unwrap());
    // }
    result
        .iter()
        .sorted_by(|&&a, &&b| {
            euclidean_distance((x1, y1), a).total_cmp(&euclidean_distance((x1, y1), b))
        })
        .map(|i| i.clone())
        .collect()
}

pub fn _ray_segment_intersection(
    px1: f32,
    py1: f32,
    px2: f32,
    py2: f32,
    px3: f32,
    py3: f32,
    epsilon: f32,
) -> bool {
    let x1 = px1;
    let y1 = py1;
    let mut x2 = px2;
    let mut y2 = py2;
    let mut x3 = px3;
    let mut y3 = py3;

    if y2 > y3 {
        (y2, y3) = (y3, y2);
        (x2, x3) = (x3, x2);
    }
    let y1 = if y1 == y2 || y1 == y3 {
        y1 + epsilon
    } else {
        y1
    };
    if (y1 > y3 || y1 < y2) || (x1 > f32::max(x2, x3)) {
        return false;
    }
    if x1 < f32::min(x2, x3) {
        return true;
    } else {
        let m_red = if (x2 - x3).abs() > epsilon {
            (y3 - y2) / (x3 - x2)
        } else {
            f32::INFINITY
        };
        let m_blue = if (x2 - x1).abs() > epsilon {
            (y1 - y2) / (x1 - x2)
        } else {
            f32::INFINITY
        };
        return m_blue >= m_red;
    }
}

pub fn middle(x1: f32, y1: f32, x2: f32, y2: f32) -> Point {
    let dx = x2 - x1;
    let dy = y2 - y1;
    (x1 + dx / 2.0, y1 + dy / 2.0)
}

pub fn point_in_polygon(x1: f32, y1: f32, polygon: &[Point], bounds: &Bounds) -> bool {
    let intersections = segment_polygon_intersection(
        bounds.min_x - 1.0,
        y1,
        bounds.max_x + 1.0,
        y1,
        polygon,
        true,
    );
    let (mut left, mut right) = (0, 0);
    for cut in intersections {
        if cut.0 < x1 && significantly_different(cut.0, x1) {
            left += 1;
        }
        if cut.0 > x1 && significantly_different(cut.0, x1) {
            right += 1;
        }
    }
    return left % 2 == 1 && right % 2 == 1;
}

fn significantly_different(f1:f32, f2:f32) -> bool {
    (f1-f2).abs() > EPSILON
}

pub fn intersection_length(
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    polygon: &[Point],
    bounds: &Bounds,
) -> f32 {
    let mut cuts = segment_polygon_intersection(x1, y1, x2, y2, polygon, true);
    cuts.push((x2, y2));
    cuts.insert(0, (x1, y1));
    let mut distance = 0.0;
    for i in 0..cuts.len() - 1 {
        let (x3, y3) = (cuts[i].0, cuts[i].1);
        let (x4, y4) = (cuts[i + 1].0, cuts[i + 1].1);
        let (mx, my) = middle(x3, y3, x4, y4);
        if point_in_polygon(mx, my, polygon, bounds) {
            distance += euclidean_distance((x3, y3), (x4, y4));
        }
    }
    return distance;
}
pub fn fermat_point(a: Point, b: Point, c: Point, epsilon: f32) -> Point {
    use nalgebra::{Matrix2, Vector2};

    let va = Vector2::new(a.0, a.1);
    let vb = Vector2::new(b.0, b.1);
    let vc = Vector2::new(c.0, c.1);

    let ab = vb - va;
    let ac = vc - va;
    let ba = va - vb;
    let bc = vc - vb;

    let ang1 = ((ab.dot(&ac)) / (ab.norm() * ac.norm())).acos();
    let ang2 = ((ba.dot(&bc)) / (ba.norm() * bc.norm())).acos();
    let ang3 = std::f32::consts::PI - (ang1 + ang2);

    let deg_lim = RADIANS_120_DEGREE - epsilon;
    if ang1 >= deg_lim {
        return a;
    }
    if ang2 >= deg_lim {
        return b;
    }
    if ang3 >= deg_lim {
        return c;
    }
    if ab.norm() < epsilon {
        return a;
    }
    if bc.norm() < epsilon {
        return b;
    }
    if ac.norm() < epsilon {
        return c;
    }

    let theta = RADIANS_120_DEGREE / 2.0; // sixty degree
    let rot_a = Matrix2::from([[theta.cos(), -theta.sin()], [theta.sin(), theta.cos()]]);
    let theta = -RADIANS_120_DEGREE / 2.0;
    let rot_b = Matrix2::from([[theta.cos(), -theta.sin()], [theta.sin(), theta.cos()]]);

    let b_star1 = rot_a * ac + va;
    let c_star1 = rot_a * ab + va;
    let b_star2 = rot_b * ac + va;
    let c_star2 = rot_b * ab + va;

    let mut b_star = b_star1;
    if (vb - b_star1).norm() < (vb - b_star2).norm() {
        b_star = b_star2;
    }

    let mut c_star = c_star1;
    if (vc - c_star1).norm() < (vc - c_star2).norm() {
        c_star = c_star2;
    }

    let x1 = b_star.x;
    let y1 = b_star.y;
    let x2 = vb.x;
    let y2 = vb.y;
    let x3 = c_star.x;
    let y3 = c_star.y;
    let x4 = vc.x;
    let y4 = vc.y;

    match segment_segment_intersection(x1, y1, x2, y2, x3, y3, x4, y4, true) {
        Some(x) => x,
        None => (x1, y1),
    }
}

pub fn _centroid(a: Point, b: Point,c: Point) -> Point {
    ((a.0+b.0+c.0)/3.0,(a.1+b.1+c.1)/3.0)
}
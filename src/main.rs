use geometry::Bounds;
use geometry::overlap;
use itertools::Itertools;
use petgraph::data::FromElements;
use petgraph::visit::EdgeRef;

use rand::{distributions::Uniform, prelude::Distribution, Rng, SeedableRng};

mod geometry {
    pub const RADIANS_120_DEGREE: f64 = 2.0 * std::f64::consts::PI / 3.0;

    use std::f64::INFINITY;

    use itertools::Itertools;

    use crate::Point;

    #[derive(Debug, Clone)]
    pub struct Bounds {
        pub min_x: f64,
        pub max_x: f64,
        pub min_y: f64,
        pub max_y: f64,
    }

    impl Default for Bounds {
        fn default() -> Self {
            Self {
                min_x: INFINITY,
                max_x: 0.0,
                min_y: INFINITY,
                max_y: 0.0,
            }
        }
    }

    pub fn euclidean_distance(a: Point, b: Point) -> f64 {
        ((a.0 - b.0).powf(2.0) + (a.1 - b.1).powf(2.0)).sqrt()
    }

    pub fn overlap(x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64, x4: f64, y4: f64) -> bool {
        !(x2 < x3 || x4 < x1 || y2 < y3 || y4 < y1)
    }

    pub fn segment_segment_intersection(
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        x3: f64,
        y3: f64,
        x4: f64,
        y4: f64,
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
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
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
                // println!("line {:?} and line {:?} intersect at {:?}", (x1, y1, x2, y2), (x3, y3, x4, y4), intersection);
                if !result.contains(&intersection.unwrap()) {
                    result.push(intersection.unwrap())
                }
            }
        }
        if result.contains(&(x1, y1)) {
            result.remove(result.iter().position(|i| i.0 == x1 && i.1 == y1).unwrap());
        }
        if result.contains(&(x2, y2)) {
            result.remove(result.iter().position(|i| i.0 == x2 && i.1 == y2).unwrap());
        }
        result
            .iter()
            .sorted_by(|&&a, &&b| {
                euclidean_distance((x1, y1), a).total_cmp(&euclidean_distance((x1, y1), b))
            })
            .map(|i| i.clone())
            .collect()
    }

    pub fn _ray_segment_intersection(
        px1: f64,
        py1: f64,
        px2: f64,
        py2: f64,
        px3: f64,
        py3: f64,
        epsilon: f64,
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
        if (y1 > y3 || y1 < y2) || (x1 > f64::max(x2, x3)) {
            return false;
        }
        if x1 < f64::min(x2, x3) {
            return true;
        } else {
            let m_red = if (x2 - x3).abs() > epsilon {
                (y3 - y2) / (x3 - x2)
            } else {
                f64::INFINITY
            };
            let m_blue = if (x2 - x1).abs() > epsilon {
                (y1 - y2) / (x1 - x2)
            } else {
                f64::INFINITY
            };
            return m_blue >= m_red;
        }
    }

    pub fn middle(x1: f64, y1: f64, x2: f64, y2: f64) -> Point {
        let dx = x2 - x1;
        let dy = y2 - y1;
        (x1 + dx / 2.0, y1 + dy / 2.0)
    }

    pub fn point_in_polygon(x1: f64, y1: f64, polygon: &[Point], bounds: &Bounds) -> bool {
        // let mut intersections = 0;
        // for i in -1..(polygon.len() as i32 - 1) {
        //     let (x2, y2) = (polygon[if i == -1 {polygon.len() - 1} else {i as usize}].0,polygon[if i == -1 {polygon.len() - 1} else {i as usize}].1);
        //     let (x3, y3) = (polygon[(i + 1) as usize].0, polygon[(i + 1) as usize].1);
        //     if ray_segment_intersection(x1, y1, x2, y2, x3, y3, 1e-6) {
        //         intersections += 1
        //     }
        // }
        // return intersections % 2 == 1;

        let intersections = segment_polygon_intersection(
            bounds.min_x - 1.0,
            y1,
            bounds.max_x + 1.0,
            y1,
            polygon,
            true,
        );
        // println!("point at {:?}", (x1,y1));
        // println!("intersections at {:?}", intersections);
        let (mut left, mut right) = (0, 0);
        for cut in intersections {
            if cut.0 < x1 {
                left += 1;
            }
            if cut.0 > x1 {
                right += 1;
            }
        }
        return left % 2 == 1 && right % 2 == 1;

        // let (bx1,by1,bx2,by2) = bounds;
        // let intersections1 = segment_polygon_intersection(bx1-1.0, y1, bx2+1.0, y1, polygon, true);
        // let intersections2 = segment_polygon_intersection(x1, by1-1.0, x1, by2+1.0, polygon, true);
        // let (mut left, mut right, mut above, mut below) = (0,0,0,0);
        // for cut in intersections1.iter().chain(intersections2.iter()) {
        //     if cut.0 < x1 {
        //         left+=1;
        //     }
        //     if cut.0 > x1 {
        //         right+=1;
        //     }
        //     if cut.1 < y1 {
        //         below+=1;
        //     }
        //     if cut.1 > y1 {
        //         above+=1;
        //     }
        // }
        // return left%2==1 && right%2==1 && above%2==1 && below%2==1;
    }

    pub fn intersection_length(
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        polygon: &[Point],
        bounds: &Bounds,
    ) -> f64 {
        let mut cuts = segment_polygon_intersection(x1, y1, x2, y2, polygon, true);
        cuts.push((x2, y2));
        cuts.insert(0, (x1, y1));
        let mut distance = 0.0;
        for i in 0..cuts.len() - 1 {
            let (x2, y2) = (cuts[i].0, cuts[i].1);
            let (x3, y3) = (cuts[i + 1].0, cuts[i + 1].1);
            let (mx, my) = middle(x2, y2, x3, y3);
            if point_in_polygon(mx, my, polygon, bounds) {
                distance += euclidean_distance((x2, y2), (x3, y3));
            }
        }
        return distance;
    }
    pub fn fermat_point(a: Point, b: Point, c: Point, epsilon: f64) -> Point {
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
        let ang3 = std::f64::consts::PI - (ang1 + ang2);

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
}

use std::{collections::HashSet, rc::Rc};

type Point = (f64, f64);

const M_RANGE_MIN: f64 = 0.01;
const NUMBER_OFFSPRING: i32 = 500 / 3;
const P_FLIP_MOVE_MIN: f64 = 0.01;
const INF: f64 = 1e10;

struct SteinerProblem {
    terminals: Vec<Point>,
    obstacles: Vec<Obstacle>,
    obstacle_corners: Vec<Point>,
    centroids: Vec<Point>,
    bounds: Bounds,
    average_terminal_distance: f64,
}

impl SteinerProblem {
    fn new(terminals: Vec<Point>, obstacles: Vec<Obstacle>) -> Self {
        let mut obstacle_corners = Vec::new();
        for obstacle in &obstacles {
            for point in &obstacle.points {
                obstacle_corners.push(*point);
            }
        }
        let mut centroids = Vec::new();
        let mut v = terminals.clone();
        v.append(&mut obstacle_corners.clone());
        let v = v.iter().map(|p| [p.0, p.1]).collect::<Vec<[f64; 2]>>();
        use geos::Geom;
        let seq = geos::CoordSeq::new_from_vec(&v).expect("failed to create CoordSeq");
        let geo = geos::Geometry::delaunay_triangulation(
            &seq.create_line_string()
                .expect("could not create line string"),
            0.0001,
            false,
        )
        .expect("could not triangulate");

        for n in 0..geo.get_num_geometries().expect("failed to get n") {
            let geometry = geo.get_geometry_n(n);
            if geometry.is_err() {
                break;
            }
            let mut triangle = Vec::new();
            let geometry = geometry.unwrap();
            let ring = geometry.get_exterior_ring().unwrap();
            for i in 0..3 {
                let coord_seq = ring.get_coord_seq().unwrap();
                triangle.push([coord_seq.get_x(i).unwrap(), coord_seq.get_y(i).unwrap()]);
            }
            centroids.push(geometry::fermat_point(
                (triangle[0][0], triangle[0][1]),
                (triangle[1][0], triangle[1][1]),
                (triangle[2][0], triangle[2][1]),
                1e-6,
            ));
        }

        let mut bounds = Bounds {
            min_x: INF,
            max_x: 0.0,
            min_y: INF,
            max_y: 0.0,
        };
        for point in terminals.iter().chain(obstacle_corners.iter()) {
            if point.0 < bounds.min_x {
                bounds.min_x = point.0
            }
            if point.1 < bounds.min_y {
                bounds.min_y = point.1
            }
            if point.0 > bounds.max_x {
                bounds.max_x = point.0
            }
            if point.1 > bounds.max_y {
                bounds.max_y = point.1
            }
        }
        let mut average_terminal_distance = 0.0;
        let mut counter = 0;
        for i in 0..terminals.len() {
            for j in (i + 1)..terminals.len() {
                average_terminal_distance +=
                    (nalgebra::Vector2::new(terminals[i].0, terminals[i].1)
                        - nalgebra::Vector2::new(terminals[j].0, terminals[j].1))
                    .norm();
                counter += 1;
            }
        }
        average_terminal_distance /= counter as f64;

        SteinerProblem {
            terminals,
            obstacles,
            obstacle_corners,
            centroids,
            bounds,
            average_terminal_distance,
        }
    }

    fn coordinates_in_solid_obstacle(&self, coordinates: Point) -> bool {
        for obstacle in self.obstacles.iter() {
            if obstacle.weight == INF {
                if geometry::point_in_polygon(
                    coordinates.0,
                    coordinates.1,
                    &obstacle.points,
                    &obstacle.bounds,
                ) {
                    return true;
                }
            }
        }
        false
    }
}

#[derive(Clone)]
struct Chromosome {
    steiner_points: Vec<Point>,
    included_corners: HashSet<usize>,
}

impl std::fmt::Debug for Chromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let string = format!("{:?}", self.included_corners);
        let len = string.len();
        f.write_str(
            format!(
                "Chromosome(steinerPoints={:?}, includedObstacleCornersIndices=set([{}]))",
                self.steiner_points,
                string.chars().skip(1).take(len - 2).collect::<String>()
            )
            .as_str(),
        )
    }
}

#[derive(Clone)]
struct MinimumSpanningTree {
    // edge_weights: Vec<f64>,
    total_weight: f64,
    // mst: Vec<(usize, usize)>,
    graph: petgraph::graph::UnGraph<Point, f64, u32>,
}

#[derive(Clone)]
struct Instance {
    problem: Rc<SteinerProblem>,
    chromosome: Chromosome,
    minimum_spanning_tree: Option<MinimumSpanningTree>,
}

struct StOBGA {
    problem: Rc<SteinerProblem>,
    population: Vec<Instance>,
    random_generator: rand::rngs::StdRng,
    current_generation: usize,
    child_buffer: Vec<Instance>,
}

impl StOBGA {
    fn crossover(&mut self, parent_1_index: usize, parent_2_index: usize) {
        let min_x = self.problem.bounds.min_x;
        let max_x = self.problem.bounds.max_x;
        let random_x_value = self.random_generator.gen_range(min_x..max_x);

        let mut steiner_points_1 = Vec::new();
        let mut steiner_points_2 = Vec::new();

        let mut obstacle_corners_1 = HashSet::new();
        let mut obstacle_corners_2 = HashSet::new();

        for point in self.population[parent_1_index]
            .chromosome
            .steiner_points
            .iter()
        {
            if point.0 < random_x_value {
                steiner_points_1.push(point.clone());
            } else {
                steiner_points_2.push(point.clone());
            }
        }
        for point in self.population[parent_2_index]
            .chromosome
            .steiner_points
            .iter()
        {
            if point.0 > random_x_value {
                steiner_points_1.push(point.clone());
            } else {
                steiner_points_2.push(point.clone());
            }
        }

        for index in self.population[parent_1_index]
            .chromosome
            .included_corners
            .iter()
        {
            let point = self.population[parent_1_index].problem.obstacle_corners[*index];
            if point.0 < random_x_value {
                obstacle_corners_1.insert(*index);
            } else {
                obstacle_corners_2.insert(*index);
            }
        }

        for index in self.population[parent_2_index]
            .chromosome
            .included_corners
            .iter()
        {
            let point = self.population[parent_2_index].problem.obstacle_corners[*index];
            if point.0 > random_x_value {
                obstacle_corners_1.insert(*index);
            } else {
                obstacle_corners_2.insert(*index);
            }
        }

        self.child_buffer.push(Instance {
            chromosome: Chromosome {
                steiner_points: steiner_points_1,
                included_corners: obstacle_corners_1,
            },
            minimum_spanning_tree: None,
            problem: self.problem.clone(),
        });
        self.child_buffer.push(Instance {
            chromosome: Chromosome {
                steiner_points: steiner_points_2,
                included_corners: obstacle_corners_2,
            },
            minimum_spanning_tree: None,
            problem: self.problem.clone(),
        });
    }

    fn new(
        problem: Rc<SteinerProblem>,
        population_size: usize,
        t1: usize,
        t2: usize,
        t3: usize,
    ) -> Self {
        let mut rng = rand::rngs::StdRng::from_seed([0; 32]);
        let mut population = vec![];

        for _ in 0..t1 {
            population.push(Instance {
                problem: Rc::clone(&problem),
                chromosome: Chromosome {
                    steiner_points: problem.centroids.clone(),
                    included_corners: HashSet::new(),
                },
                minimum_spanning_tree: Option::None,
            });
        }

        let k = problem.obstacle_corners.len();
        let n = problem.terminals.len();
        let min_x = problem.bounds.min_x;
        let max_x = problem.bounds.max_x;
        let min_y = problem.bounds.min_y;
        let max_y = problem.bounds.max_y;
        let x_dist = Uniform::new(min_x, max_x);
        let y_dist = Uniform::new(min_y, max_y);
        for _ in 0..t2 {
            let mut steiner_points = Vec::new();
            let r = rng.gen_range(0..(n + k));
            for _ in 0..r {
                steiner_points.push((rng.sample(x_dist), rng.sample(y_dist)));
            }
            population.push(Instance {
                problem: Rc::clone(&problem),
                chromosome: Chromosome {
                    steiner_points: steiner_points,
                    included_corners: HashSet::new(),
                },
                minimum_spanning_tree: Option::None,
            });
        }

        for _ in 0..t3 {
            let distribution = Uniform::new(0, k + 1);
            let amount = rng.sample(distribution);
            let draws = rand::seq::index::sample(&mut rng, k, amount);
            let mut corners = HashSet::new();
            for elem in draws {
                corners.insert(elem);
            }

            population.push(Instance {
                problem: Rc::clone(&problem),
                chromosome: Chromosome {
                    steiner_points: Vec::new(),
                    included_corners: corners,
                },
                minimum_spanning_tree: Option::None,
            })
        }

        let mut stobga = StOBGA {
            problem,
            population,
            random_generator: rng,
            current_generation: 0,
            child_buffer: Vec::new(),
        };

        let mut parents = Vec::new();
        let mut children = Vec::new();
        for _ in 0..(population_size - (t1 + t2 + t3)) {
            parents.push(stobga.tournament_select(5, false));
            children.push(stobga.tournament_select(5, true));
        }
        let mut save: Vec<Instance> = children
            .iter()
            .map(|i| stobga.population[*i].clone())
            .collect();
        for i in 0..(parents.len() / 2) {
            stobga.crossover(parents[2 * i], parents[2 * i + 1]);
        }
        for (child_index, population_index) in children.iter().enumerate() {
            stobga.population[*population_index] = stobga.child_buffer[child_index].clone();
        }
        stobga.population.append(&mut save);
        stobga.child_buffer.clear();
        stobga
    }

    fn tournament_select(&mut self, size: usize, to_die: bool) -> usize {
        if to_die {
            return rand::seq::index::sample(
                &mut self.random_generator,
                self.population.len(),
                size,
            )
            .iter()
            .max_by(|i1, i2| {
                let w1 = self.population[*i1].get_mst().total_weight;
                let w2 = self.population[*i2].get_mst().total_weight;
                w1.total_cmp(&w2)
            })
            .unwrap();
        } else {
            return rand::seq::index::sample(
                &mut self.random_generator,
                self.population.len(),
                size,
            )
            .iter()
            .min_by(|i1, i2| {
                let w1 = self.population[*i1].get_mst().total_weight;
                let w2 = self.population[*i2].get_mst().total_weight;
                w1.total_cmp(&w2)
            })
            .unwrap();
        }
    }

    fn step(&mut self) {
        let mut parents = Vec::new();
        let mut children = Vec::new();
        for _ in 0..NUMBER_OFFSPRING {
            parents.push(self.tournament_select(5, false));
            children.push(self.tournament_select(5, true));
        }
        let p_flip_move = f64::max(
            1.0 - (self.current_generation as f64) / 1000.0,
            P_FLIP_MOVE_MIN,
        );
        self.child_buffer.clear();
        for pair in parents.iter().as_slice().windows(2) {
            self.crossover(pair[0], pair[1]);
        }
        for (child_index, population_index) in children.iter().enumerate() {
            self.child_buffer[child_index].mutate(
                &mut self.random_generator,
                p_flip_move,
                self.current_generation,
            );
            self.population[*population_index] = self.child_buffer[child_index].clone();
        }
        self.child_buffer.clear();
        for instance in self.population.iter_mut() {
            let _ = instance.get_mst();
        }
        self.population.sort_unstable_by(|i1, i2| {
            i1.minimum_spanning_tree
                .as_ref()
                .unwrap()
                .total_weight
                .total_cmp(&i2.minimum_spanning_tree.as_ref().unwrap().total_weight)
        });
        self.current_generation += 1;
    }
}

impl Instance {
    fn get_mst(&mut self) -> &MinimumSpanningTree {
        if self.minimum_spanning_tree.is_none() {
            let mut graph = petgraph::Graph::new_undirected();
            let vertices = self
                .problem
                .terminals
                .iter()
                .chain(self.chromosome.steiner_points.iter())
                .chain(
                    self.chromosome
                        .included_corners
                        .iter()
                        .map(|c| &self.problem.obstacle_corners[*c]),
                )
                .map(|&c| c)
                .collect::<Vec<_>>();

            for vertex in vertices.iter() {
                graph.add_node((*vertex).clone());
            }
            for i1 in 0..(vertices.len() - 1) {
                let t1 = vertices[i1];
                for i2 in (i1 + 1)..vertices.len() {
                    let t2 = vertices[i2];
                    let mut length = geometry::euclidean_distance(t1, t2);
                    let line_bounds = Bounds {
                        min_x: t1.0.min(t2.0),
                        min_y: t1.1.min(t2.1),
                        max_x: t1.0.max(t2.0),
                        max_y: t1.1.max(t2.1),
                    };
                    for obstacle in &self.problem.obstacles {
                        let bounds = &obstacle.bounds;
                        if overlap(
                            line_bounds.min_x,
                            line_bounds.min_y,
                            line_bounds.max_x,
                            line_bounds.max_y,
                            bounds.min_x,
                            bounds.min_y,
                            bounds.max_x,
                            bounds.max_y,
                        ) {
                            let intersection_len = geometry::intersection_length(
                                t1.0,
                                t1.1,
                                t2.0,
                                t2.1,
                                &obstacle.points,
                                &obstacle.bounds,
                            );
                            if intersection_len > 0.0 {
                                if obstacle.weight == INF {
                                    length = INF;
                                } else {
                                    length -= intersection_len;
                                    length += intersection_len * obstacle.weight;
                                }
                            }
                        }
                    }
                    graph.add_edge(
                        petgraph::graph::NodeIndex::new(i1),
                        petgraph::graph::NodeIndex::new(i2),
                        length,
                    );
                }
            }
            let mst = petgraph::graph::UnGraph::<_, _>::from_elements(
                petgraph::algo::min_spanning_tree(&graph),
            );
            let total_distance = mst.edge_weights().sum::<f64>();
            let mst = MinimumSpanningTree {
                total_weight: total_distance,
                graph: mst,
            };

            self.minimum_spanning_tree = Option::Some(mst);
        }
        return self.minimum_spanning_tree.as_ref().unwrap();
    }

    fn mutate(
        &mut self,
        rng: &mut rand::rngs::StdRng,
        probability_flip_move: f64,
        generation: usize,
    ) {
        if rng.gen_bool(probability_flip_move) {
            self.mutation_flip_move(rng, generation);
        }
        if rng.gen_bool(0.5) {
            self.mutation_add_steiner(rng);
        } else {
            self.mutation_remove_steiner(rng)
        }
    }

    fn mutation_remove_steiner(&mut self, rng: &mut rand::rngs::StdRng) {
        let _ = self.get_mst();
        let mut candidate_steiner_points = Vec::new();

        let graph = &self.minimum_spanning_tree.as_ref().unwrap().graph;
        for steiner_point in &self.chromosome.steiner_points {
            let id = graph
                .node_indices()
                .find(|id| graph[*id].0 == steiner_point.0 && graph[*id].1 == steiner_point.1)
                .unwrap();
            let edges = graph.edges(id);
            if edges.count() <= 2 {
                candidate_steiner_points.push(
                    self.chromosome
                        .steiner_points
                        .iter()
                        .position(|i| i.0 == steiner_point.0 && i.1 == steiner_point.1)
                        .unwrap(),
                );
            }
        }
        let mut candidate_corners = Vec::new();
        for index_corner in &self.chromosome.included_corners {
            let steiner_point = self.problem.obstacle_corners[*index_corner];
            let id = graph
                .node_indices()
                .find(|id| graph[*id].0 == steiner_point.0 && graph[*id].1 == steiner_point.1)
                .unwrap();
            let edges = graph.edges(id);
            if edges.count() <= 2 {
                candidate_corners.push(index_corner.clone());
            }
        }
        match (candidate_steiner_points.len(), candidate_corners.len()) {
            (0, 0) => {}
            (0, n) => {
                self.chromosome
                    .included_corners
                    .remove(&candidate_corners[if n > 1 { rng.gen_range(0..n) } else { 0 }]);
            }
            (n, 0) => {
                self.chromosome
                    .steiner_points
                    .remove(candidate_steiner_points[if n > 1 { rng.gen_range(0..n) } else { 0 }]);
            }
            (n, m) => {
                if rng.gen_bool((n as f64 / m as f64).clamp(0.0, 1.0)) {
                    self.chromosome.steiner_points.remove(
                        candidate_steiner_points[if n > 1 { rng.gen_range(0..n) } else { 0 }],
                    );
                } else {
                    self.chromosome
                        .included_corners
                        .remove(&candidate_corners[if m > 1 { rng.gen_range(0..m) } else { 0 }]);
                }
            }
        }
        self.minimum_spanning_tree = None;
    }

    fn mutation_add_steiner(&mut self, rng: &mut rand::rngs::StdRng) {
        self.get_mst(); // making sure the mst is present
        let mut candidates = Vec::new();
        let graph = &self.minimum_spanning_tree.as_ref().unwrap().graph;
        for i1 in graph.node_indices() {
            let connections = graph.edges(i1);
            let c1 = graph[i1];
            let v1 = nalgebra::Vector2::new(c1.0, c1.1);
            for edge in connections.combinations(2) {
                let i2 = edge[0].target();
                let i3 = edge[1].target();
                let c2 = graph[i2];
                let c3 = graph[i3];
                let v2 = nalgebra::Vector2::new(c2.0, c2.1);
                let v3 = nalgebra::Vector2::new(c3.0, c3.1);
                let v12 = v2 - v1;
                let v13 = v3 - v1;
                let dot = v12.dot(&v13);
                let den = v12.norm() * v13.norm();
                let angle = (dot / den).acos();
                if angle < geometry::RADIANS_120_DEGREE {
                    candidates.push((i1, i2, i3));
                }
            }
        }
        if candidates.len() == 0 {
            // add random steiner point
            let min_x = self.problem.bounds.min_x;
            let max_x = self.problem.bounds.max_x;
            let min_y = self.problem.bounds.min_y;
            let max_y = self.problem.bounds.max_y;
            let mut new_steiner = (rng.gen_range(min_x..max_x), rng.gen_range(min_y..max_y));
            while self.problem.coordinates_in_solid_obstacle(new_steiner) {
                new_steiner = (rng.gen_range(min_x..max_x), rng.gen_range(min_y..max_y));
            }
            self.chromosome.steiner_points.push(new_steiner);
        } else {
            let random_triple = candidates[if candidates.len() > 1 {
                rng.gen_range(0..candidates.len())
            } else {
                0
            }];
            let p1 = graph[random_triple.0];
            let p2 = graph[random_triple.1];
            let p3 = graph[random_triple.2];
            let p4 = geometry::fermat_point(p1, p2, p3, 1e-9);
            if !self.problem.coordinates_in_solid_obstacle(p4) {
                self.chromosome.steiner_points.push(p4);
            }
        }
        self.minimum_spanning_tree = None;
    }

    fn mutation_flip_move(&mut self, rng: &mut rand::rngs::StdRng, generation: usize) {
        let s = self.chromosome.steiner_points.len();
        let k = self.problem.obstacle_corners.len();
        let p_gene = if s + k == 0 {
            1.0
        } else {
            1.0 / ((s + k) as f64)
        };
        let m_range = self.problem.average_terminal_distance
            * f64::max(1.0 - (generation as f64) / 1000.0, M_RANGE_MIN);
        for i in 0..s {
            if rng.gen_bool(p_gene) {
                let x_sign = if rng.gen_bool(0.5) { 1.0 } else { -1.0 };
                let y_sign = if rng.gen_bool(0.5) { 1.0 } else { -1.0 };

                if m_range > M_RANGE_MIN {
                    let dist = Uniform::new(M_RANGE_MIN, m_range);
                    self.chromosome.steiner_points[i].0 = dist.sample(rng) * x_sign;
                    self.chromosome.steiner_points[i].1 = dist.sample(rng) * y_sign;
                } else {
                    self.chromosome.steiner_points[i].0 = M_RANGE_MIN * x_sign;
                    self.chromosome.steiner_points[i].1 = M_RANGE_MIN * y_sign;
                }
            }
        }
        for i in 0..k {
            if rng.gen_bool(p_gene) {
                if self.chromosome.included_corners.contains(&i) {
                    self.chromosome.included_corners.remove(&i);
                } else {
                    self.chromosome.included_corners.insert(i);
                }
            }
        }
        self.minimum_spanning_tree = None
    }
}

#[derive(Clone)]
struct Obstacle {
    weight: f64,
    bounds: Bounds,
    points: Vec<Point>,
}

impl std::fmt::Debug for Obstacle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Obstacle")
            .field("weight", &self.weight)
            .field("bounds", &self.bounds)
            .field("points", &self.points)
            .finish()
    }
}

impl Obstacle {
    fn new(weight: f64, points: Vec<Point>) -> Self {
        Self {
            weight,
            points,
            bounds: Bounds::default(),
        }
    }

    pub(crate) fn compute_bounds(mut self) -> Obstacle {
        let mut bounds = Bounds::default();
        for point in &self.points {
            if point.0 < bounds.min_x {
                bounds.min_x = point.0
            }
            if point.1 < bounds.min_y {
                bounds.min_y = point.1
            }
            if point.0 > bounds.max_x {
                bounds.max_x = point.0
            }
            if point.1 > bounds.max_y {
                bounds.max_y = point.1
            }
        }
        self.bounds = bounds;
        self
    }
}

fn main() {
    std::env::set_var("RUST_BACKTRACE", "full");
    let mut terminals = Vec::new();
    for line in std::fs::read_to_string(std::env::args().nth(1).expect("please specify terminal file"))
        .unwrap()
        .lines()
        .skip(1)
    {
        let coords = line
            .split(",")
            .map(|c| c.parse().unwrap())
            .collect::<Vec<_>>();
        terminals.push((coords[0], coords[1]));
    }

    let mut obstacles = Vec::new();
    {
        let mut current_obstacle = Obstacle::new(0.0, vec![]);
        for line in std::fs::read_to_string(std::env::args().nth(2).expect("please specify obstacle file"))
            .unwrap()
            .lines()
        {
            if line == "" || line == "," {
                obstacles.push(current_obstacle.compute_bounds());
                current_obstacle = Obstacle::new(0.0, vec![]);
            } else if line.starts_with("max") {
                current_obstacle.weight = INF
            } else {
                let fields = line.split(",").collect::<Vec<_>>();
                if fields.get(1) == Some(&"") || fields.len() < 2 {
                    current_obstacle.weight = fields[0].parse().unwrap();
                } else {
                    current_obstacle
                        .points
                        .push((fields[0].parse().unwrap(), fields[1].parse().unwrap()));
                }
            }
        }
        obstacles.push(current_obstacle.compute_bounds());
        // println!("{:?}", obstacles);
    }
    // return;
    let problem = SteinerProblem::new(terminals.clone(), obstacles.clone());
    let mut stobga = StOBGA::new(Rc::new(problem), 500, 166, 166, 166);

    let mut streak = (0, f64::INFINITY);
    println!("generation;average;best;chromosome");
    while stobga.current_generation < 1500 {
        stobga.step();
        if stobga.population[0].get_mst().total_weight < streak.1-1e-6 {
            streak = (0, stobga.population[0].get_mst().total_weight);

            println!(
                "{};{};{};{:?}",
                stobga.current_generation,
                {
                    let mut avg = 0.0;
                    for i in stobga.population.iter_mut() {
                        avg += i.get_mst().total_weight;
                    }
                    avg / 500.0
                },
                { stobga.population[0].get_mst().total_weight },
                stobga.population[0].chromosome
            );
        } else {
            streak.0 += 1
        }
        if streak.0 == 200 {
            break;
        }
    }
    println!(
        "{};{};{};{:?}",
        stobga.current_generation,
        {
            let mut avg = 0.0;
            for i in stobga.population.iter_mut() {
                avg += i.get_mst().total_weight;
            }
            avg / 500.0
        },
        { stobga.population[0].get_mst().total_weight },
        stobga.population[0].chromosome
    );
}

#[cfg(test)]
mod test {
    use crate::geometry::*;
    use petgraph::{data::FromElements, prelude::UnGraph};

    #[test]
    fn test_geometry() {
        assert_eq!(
            crate::geometry::point_in_polygon(
                0.0,
                0.0,
                &[(-1.0, -1.0), (1.0, 1.0), (0.0, 2.0)],
                &(-1.0, -1.0, 1.0, 2.0)
            ),
            false
        )
    }

    #[test]
    fn test_geometry2() {
        assert_eq!(
            crate::geometry::segment_polygon_intersection(
                0.0,
                0.0,
                2.0,
                0.0,
                &[(1.0, 0.0), (1.0, -1.0), (-1.0, -1.0)],
                true
            ),
            vec![(1.0, 0.0)]
        );
        assert_eq!(
            crate::geometry::intersection_length(
                0.0,
                0.0,
                2.0,
                0.0,
                &[(1.0, 0.0), (1.0, -1.0), (-1.0, -1.0)],
                &(-1.0, 0.0, 1.0, 1.0)
            ),
            0.0
        );
    }

    #[test]
    fn test_geometry3() {
        assert_eq!(
            crate::geometry::segment_polygon_intersection(
                0.0,
                0.0,
                1.0,
                1.0,
                &[(0.0, 0.0), (1.0, 1.0), (1.0, -1.0)],
                true
            ),
            Vec::new()
        )
    }

    #[test]
    fn test_geometry4() {
        assert_eq!(
            crate::geometry::intersection_length(
                3.0,
                1.0,
                4.0,
                5.0,
                &[(0.0, 0.0), (3.0, 1.0), (4.0, 5.0)],
                &(-1.0, -1.0, 4.0, 2.0)
            ),
            0.0
        )
    }

    #[test]
    fn test_geometry5() {
        assert_eq!(
            crate::geometry::segment_polygon_intersection(
                3.0,
                1.0,
                4.0,
                5.0,
                &[(0.0, 0.0), (3.0, 1.0), (4.0, 5.0)],
                true
            ),
            Vec::new()
        )
    }

    #[test]
    fn test_geometry6() {
        let middle = middle(3.0, 1.0, 4.0, 5.0);
        assert!(!point_in_polygon(
            middle.0,
            middle.1,
            &[(0.0, 0.0), (3.0, 1.0), (4.0, 5.0)],
            &(0.0, 0.0, 4.0, 5.0)
        ))
    }

    #[test]
    fn test_geometry7() {
        let middle = middle(0.0, 0.0, 4.0, 5.0);
        assert!(!point_in_polygon(
            middle.0,
            middle.1,
            &[(0.0, 0.0), (3.0, 1.0), (4.0, 5.0)],
            &(0.0, 0.0, 4.0, 5.0)
        ))
    }

    #[test]
    fn test_geometry8() {
        let middle = middle(0.0, 0.0, 3.0, 1.0);
        assert!(!point_in_polygon(
            middle.0,
            middle.1,
            &[(0.0, 0.0), (3.0, 1.0), (4.0, 5.0)],
            &(0.0, 0.0, 4.0, 5.0)
        ))
    }

    #[test]
    fn test_geometry9() {
        assert_eq!(
            crate::geometry::intersection_length(
                0.0,
                1.0,
                1.0,
                1.0,
                &[(0.0, 0.0), (1.0, 0.0), (0.5, -1.0)],
                &(0.0, 0.0, 1.0, 0.0)
            ),
            0.0
        )
    }

    #[test]
    fn test_geometry10() {
        assert!(
            crate::geometry::intersection_length(
                0.845641974,
                0.904959172,
                0.753467217,
                0.42431886,
                &[
                    (0.796, 0.898),
                    (0.804, 0.784),
                    (0.906, 0.792),
                    (0.908, 0.886),
                ],
                &(0.0, 0.0, 1.0, 0.0)
            ) > 0.0
        )
    }

    #[test]
    fn test_geometry11() {
        println!(
            "{}",
            crate::geometry::intersection_length(
                0.936640447,
                0.706594727,
                0.753467217,
                0.42431886,
                &[
                    (0.784, 0.522),
                    (0.798, 0.44799999999999995),
                    (0.906, 0.45199999999999996),
                    (0.9, 0.534),
                ],
                &(0.0, 0.0, 1.0, 0.0)
            )
        );
        assert!(
            crate::geometry::intersection_length(
                0.936640447,
                0.706594727,
                0.753467217,
                0.42431886,
                &[
                    (0.784, 0.522),
                    (0.798, 0.44799999999999995),
                    (0.906, 0.45199999999999996),
                    (0.9, 0.534),
                ],
                &(0.0, 0.0, 1.0, 0.0)
            ) > 0.0
        )
    }

    #[test]
    fn using_petgraph() {
        let mut graph = petgraph::Graph::new_undirected();
        let i1 = graph.add_node((1.0, 1.0));
        let i2 = graph.add_node((2.0, 2.0));
        graph.add_edge(i1, i2, 1.0);
        let g2 = UnGraph::<_, _>::from_elements(petgraph::algo::min_spanning_tree(&graph));
        assert!(g2.edge_weights().sum::<f64>() == 1.0)
    }
}

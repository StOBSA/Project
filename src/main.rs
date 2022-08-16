use core::panic;
use geos::Geom;
use pathfinding::prelude::kruskal;

use std::{
    collections::{HashMap, HashSet},
    f64::INFINITY,
    ops::Index,
    rc::Rc,
};

type Point = [f64; 2];

struct SteinerProblem<'t> {
    terminals: Vec<Point>,
    obstacles: Vec<Obstacle<'t>>,
    obstacle_corners: Vec<Point>,
}

impl<'t> SteinerProblem<'t> {
    fn new(terminals: Vec<Point>, obstacles: Vec<Obstacle<'t>>) -> Self {
        let mut obstacle_corners = Vec::new();
        for obstacle in &obstacles {
            for point in &obstacle.points {
                obstacle_corners.push(*point);
            }
        }
        SteinerProblem {
            terminals,
            obstacles,
            obstacle_corners,
        }
    }
}

struct Chromosome {
    steiner_points: Vec<Point>,
    included_corners: HashSet<usize>,
}

struct MinimumSpanningTree {
    edges: Vec<(usize, usize, u64)>,
    edge_weights: Vec<f64>,
    total_weight: f64,
}

struct Instance<'t> {
    problem: Rc<SteinerProblem<'t>>,
    chromosome: Chromosome,
    minimum_spanning_tree: Option<MinimumSpanningTree>,
}

impl<'t> Instance<'t> {
    fn get_mst(&mut self) -> &MinimumSpanningTree {
        if self.minimum_spanning_tree.is_none() {
            let mut edges = Vec::new();
            let mut distances = Vec::new();
            let mut counter = 0;

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
                .collect::<Vec<_>>();
            
            // let mut graph = Graph::new_undirected();
            
            for i in 0..vertices.len() {
                for j in (i + 1)..vertices.len() {
                    let t1 = vertices[i];
                    let t2 = vertices[j];
                    let line = geos::CoordSeq::new_from_vec(&[t1, t2])
                        .expect("could not create CoordSeq")
                        .create_line_string()
                        .expect("could not convert CoordSeq to LineString");
                    let mut length = line.length().unwrap();
                    for obstacle in &self.problem.obstacles {
                        match line.intersection(obstacle.polygon.as_ref().unwrap()) {
                            Ok(intersect) => {
                                if length > 0.0 {
                                    length -= intersect.length().unwrap();
                                    length += intersect.length().unwrap() * obstacle.weight;
                                }
                            }
                            Err(err) => println!("{}", err),
                        }
                    }
                    distances.push(length)
                }
            }
            let mut max = 0.0;
            for distance in &distances {
                if *distance > max {
                    max = *distance
                }
            }
            let scaled_distances = distances
                .iter()
                .map(|v| (u32::MAX as f64 * (*v) / max) as u64)
                .collect::<Vec<_>>();
            // println!("distance = {:?}", distances);
            // println!("scaled distances = {:?}", scaled_distances);

            for i in 0..vertices.len() {
                for j in (i + 1)..vertices.len() {
                    edges.push((i, j, scaled_distances[counter]));
                    counter += 1;
                }
            }

            let mst = kruskal(&edges);
            let mut total_distance = 0.0;
            for edge in mst {
                let mut index = 0;
                'outer: for i in 0..vertices.len() {
                    for j in (i + 1)..vertices.len() {
                        if i == *edge.0 && j == *edge.1 {
                            break 'outer;
                        }
                        index += 1;
                    }
                }
                // println!("{:?} to {:?}: {}", vertices[*edge.0], vertices[*edge.1], distances[index]);
                total_distance += distances[index];
                // total_distance += edge.2 as f64
            }
            let mst = MinimumSpanningTree {
                edge_weights: distances,
                edges: edges,
                total_weight: total_distance,
            };

            self.minimum_spanning_tree = Option::Some(mst);
        }
        return self.minimum_spanning_tree.as_ref().unwrap();
    }
}

#[derive(Clone)]
struct Obstacle<'t> {
    weight: f64,
    points: Vec<Point>,
    polygon: Option<geos::Geometry<'t>>,
}

impl<'t> Obstacle<'t> {
    fn new(weight: f64, points: Vec<[f64; 2]>) -> Self {
        Self {
            weight,
            points,
            polygon: Option::None,
        }
    }
}

fn main() {
    let mut reader = csv::ReaderBuilder::new()
        .from_path(
            "/home/lama/Dokumente/Uni/Amsterdam/sem\
        ester 22/Thesis/code/dataset/p618-rosenberg_suppl/\
        TestCases/SoftObstacles/terminals1.csv",
        )
        .expect("error");
    let terminals = reader
        .records()
        .map(|u| match u {
            Ok(a) => [
                a[0].parse::<f64>().expect("could not parse first field"),
                a[1].parse::<f64>().expect("could not parse second field"),
            ],
            Err(_) => panic!("could not read line in terminal file"),
        })
        .collect::<Vec<[f64; 2]>>();

    let mut obstacles = Vec::new();
    {
        let mut current_obstacle = Obstacle::new(0.0, vec![]);
        for line in csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(
                "/home/lama/Dokumente/Uni/Amsterdam/semester 22\
            /Thesis/code/dataset/p618-rosenberg_suppl/TestCases/SoftObs\
            tacles/obstacles1.csv",
            )
            .expect("could not open obstacle file")
            .records()
            .map(|u| u.unwrap())
        {
            if &line[0] == "" {
                obstacles.push(current_obstacle);
                current_obstacle = Obstacle::new(0.0, vec![]);
            } else if &line[1] == "" {
                if &line[0] == "max" {
                    current_obstacle.weight = INFINITY;
                } else {
                    current_obstacle.weight = (&line[0]).parse().expect("bad csv");
                }
            } else {
                current_obstacle.points.push([
                    (&line[0]).parse().expect("bad csv"),
                    (&line[1]).parse().expect("bad csv"),
                ]);
            }
        }
        obstacles.push(current_obstacle);
    }

    for obstacle in obstacles.iter_mut() {
        let mut points = obstacle.points.clone();
        points.push(points[0]);
        obstacle.polygon = Option::Some(
            geos::Geometry::create_polygon(
                geos::CoordSeq::new_from_vec(&points)
                    .unwrap()
                    .create_linear_ring()
                    .unwrap(),
                vec![],
            )
            .unwrap(),
        )
    }
    let prob = SteinerProblem::new(terminals.clone(), obstacles.clone());
    let mut ins = Instance {
        problem: Rc::from(prob),
        minimum_spanning_tree: Option::None,
        chromosome: Chromosome {
            steiner_points: vec![],
            included_corners: HashSet::new(),
        },
    };

    println!("{:?}", ins.get_mst().total_weight);

    // let line_a = geos::CoordSeq::new_from_vec(
    //     &[
    //       &[0.0,0.0],
    //       &[0.0,1.0],
    //       &[1.0,1.0],
    //       &[1.0,0.0],
    //       &[0.0,0.0]])
    //       .expect("failed")
    //       .create_linear_ring()
    //       .expect("failed");

    // let line_b = geos::Geometry::new_from_wkt("LINESTRING (-1 0.5, 2 0.5)")
    //     .expect("could not create point");

    // let polygon = geos::Geometry::create_polygon(line_a, vec![])
    //     .expect("err");

    // let intersection = polygon.intersection(&line_b).expect("failed");
    // println!("{:?}", intersection.get_length());
}

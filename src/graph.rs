use std::{hash::Hash, collections::{HashMap, HashSet}, cmp::Ordering};

use crate::{OPoint, Point};

#[derive(Debug, Clone, Copy)]
pub struct Edge {
    pub start : OPoint,
    pub end   : OPoint
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.start == other.start && self.end == other.end) || 
        (self.start == other.end && self.end == other.start)
    }
}

impl Eq for Edge {
    fn assert_receiver_is_total_eq(&self) {}
}

impl Hash for Edge {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        if self.start.1 < self.end.1 {
            self.start.hash(state);
            self.end.hash(state);
        } else if self.start.1 == self.end.1 {
            if self.start.0 < self.end.0 {
                self.start.hash(state);
                self.end.hash(state);    
            } else {
                self.end.hash(state);
                self.start.hash(state);    
            }
        } else {
            self.end.hash(state);
            self.start.hash(state);
        }
    }
}

#[derive(Debug)]
pub struct Graph {
    pub nodes : HashSet<OPoint>,
    pub edges : HashMap<Edge, f32>
}

impl Graph {
    pub fn new() -> Self {
        Graph {
            nodes : HashSet::new(),
            edges : HashMap::new(),
        }
    }
    pub fn add_node(&mut self, node : OPoint) {
        self.nodes.insert(node);
    }
    pub fn add_edge(&mut self, a: OPoint, b: OPoint, weight: f32) {
        let edge = Edge { start: a, end: b };
        self.add_node(a);
        self.add_node(b);
        self.edges.insert(edge, weight);
    }
    pub fn add_edge_from_points(&mut self, a: Point, b: Point, weight: f32) {
        let a = crate::util::to_graph(a);
        let b = crate::util::to_graph(b);
        self.add_edge(a, b, weight)
    }
    pub fn edges_connected_to(&self, node : OPoint) -> HashSet<Edge> {
        self.edges.iter()
        .filter(|(edge, &_)|edge.start==node||edge.end==node)
        .map(|(edge, _)| edge.clone())
        .collect()
    }
    pub fn edges_connected_to_point(&self, node : Point) -> HashSet<Edge> {
        let node = crate::util::to_graph(node);
        self.edges_connected_to(node)
    }
    // fn has_circle(&self, start_node : OPoint) -> bool {
    //     let mut seen = HashSet::new();
    //     seen.insert(start_node);
    //     let mut visible = HashSet::new();
    //     // 
    //     for self.edges_connected_to(node)
    // }
    pub fn minimum_spanning_tree(&self) -> Self {
        fn add_edges(accumulator: &mut Vec<Edge>, other: &HashSet<Edge>, graph: &Graph) {
            for node in other {
                accumulator.push(*node);
            }
            accumulator.sort_by(|e1, e2| if &graph.edges[e1] < &graph.edges[e2] {Ordering::Less} else {Ordering::Greater});
        }
        let first_node = self.nodes.iter().next().expect("graph has no nodes.");
        let mut visited = HashSet::new();
        visited.insert(*first_node);
        let mut available = Vec::new();
        add_edges(&mut available, &self.edges_connected_to(*first_node), self);
        let mut accepted_edges = HashMap::new();
        let target_len = self.nodes.len();
        let mut current_len = 1;
        while current_len < target_len {
            let edge = available.remove(0);
            match (visited.contains(&edge.start), visited.contains(&edge.end)) {
                (true, true) => {},
                (true, false) => {
                    visited.insert(edge.end);
                    add_edges(&mut available, &self.edges_connected_to(edge.end), self);
                    current_len+=1;
                    accepted_edges.insert(edge, self.edges[&edge]);
                },
                (false, true) => {
                    visited.insert(edge.start);
                    add_edges(&mut available, &self.edges_connected_to(edge.start), self);
                    current_len+=1;
                    accepted_edges.insert(edge, self.edges[&edge]);
                },
                (false, false) => panic!("got forrest"),
            }
        }
        Graph { nodes: visited, edges: accepted_edges }
    }
}
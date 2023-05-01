# StOBGA

## What is it?
StOBGA is a genetic algorithm introduced in 2021 by Rosenberg et al. 
It is a metaheuristic algorithm, searching for solutions
to instances of the Euclidean Steiner Tree Problem with Soft Obstacles.

## How to use it?
This is an independent replication of the algorithm described by Rosenberg et al..
The Steiner Tree Problem itself is an NP-Hard problem that plays a role in
a large number of fields, like transportation or communication network planning.

The replication was implemented in the [Rust](https://www.rust-lang.org/tools/install) programming language. To build
the executable, execute
    
    cargo build --release

Alternatively, the program can be run using

    cargo run --release -- <terminal-file.csv> <obstacle-file.csv> <random-seed : int>

## Variants
The code used for the the variants
M1, M2 ,M3, P1, P2, P3, SA1, SA2, and SA3
can be found on the corresponding branches of this git repository.

## Results
The results of ten runs of ten variants of StOBGA over 79 problem instances
can be found in the `results` folder.

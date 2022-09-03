# StOBGA

## What is it?
StOBGA is a genetic algorithm introduced in 2020 and described in a paper
authored by Rosenberg et al. It is a novel approach to find the Minimum 
Euclidean Steiner Tree, connecting a set of *Terminals* in the presence of 
polygonal *Obstacles* with arbitrary weights (anything goes, 0 to ∞).

## How does it work?
As a genetic algorithm usually does, the principles of evolution are employed
to shape the individuals over the course of multiple generations to form ever
better solutions. Read Rosenberg et al.'s paper if you are curious - or try
    
    $ git checkout graphics
    $ cargo run --release -- SoftObstacles/terminals2.csv SoftObstacles/obstacles2.csv

if you want to see it in action.

## Why?
This is an independent replication of Rosenberg et al.'s paper to verify the
authors' results.
The Steiner Tree Problem itself is an NP-Hard problem that plays a role in
a large number of fields, like transportation or communication network planning.
StOBGA shows how a theoretically hard problem can often be solved well enough
through the help of heuristic algorithms!

## Notes
 - geos.triangulation most likely introduces randomness into the program,
   making perfect replications (even through the use of the same seed) infeasible...
  - mh, no that wasn't it...
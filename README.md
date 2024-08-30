# CS 570 Traveling Salesman Problem Solver

## Description

This project uses three different algorithms to provide solutions for the Traveling Salesman Problem: Brute Force, Nearest Neighbor, and an Original Solution. The Brute Force algorithm calculates the distance of every possible tour for a graph using permutations, guaranteeing the most optimal tour. This is only practical for very small graphs. The Nearest Neighbor algorithm is a popular heuristic algorithm that quickly finds a "pretty good" tour of a graph, but this solution is rarely optimal. The Original Solution developed combines the Nearest Neighbor algorithm with some brute force optimizations an dother strategies to attempt to remove locality bias. This algorithm takes longer than the Nearest Neighbor algorithm, but it usually finds a more optimal solution. This project folder also contains a report that discusses each algorithm in depth. It also contains the solutions found by the Original algorithm for the four graphs used in the class competition.

## Getting Started

### Dependencies

* C++ compiler
* Standard C++ libraries

### Installing

* All of the files needed are in this folder
* The source code files are labelled 'TSP_bruteforce.cpp', 'TSP_nearestneighbor.cpp', and 'TSP_Original.cpp'
* A graph file is required. It is recommended to use one of the graph files used in the competition.

### Executing program

* Change the graph variables to match the graph file being used. For TSP_bruteforce.cpp, these are 'graphFile' in line 55 and numVertices in line 62. For TSP_nearestneighbor.cpp, these are 'graphFile' in line 66 and 'numVertices' in line 71. For TSP_Original, these are 'graphName' in line 32 and 'numberOfVertices' in line 33.
* Compile each program using a C++ compiler with the '-03' optimization flag. Example:
```
g++ TSP_Original.cpp -O3
```
* Run the executable file using './a.out' or './a.exe'

### Optional Parameters (for TSP_Original.cpp)

* There are five optional parameters that can be edited when running the Original algorithm. The goal of this feature was to make the program more flexible for quick changes during the competition. These parameters were all changed frequently during the competition.
* 'K' (line 30): This is the scalar value used for generating alternate tours. The default of this is 1.5, but it can be reduced to 1.0 for very large graphs or graphs with a lot of duplicate edge weights.
* 'runAnyStartingPoint' (line 34): This boolean value can be toggled to enable or disable running the Nearest Neighbor algorithm from each starting point. This can be useful for the first run of a graph, but it can be disabled for every subsequent run.
* 'minStartingPoint' (line 35): After running the program once, the previously mentioned boolean value can be changed to false and this integer value can be changed to the best starting point. This saves time and optimizes subsequent runs.
* 'runUntilFindSolutionLessThanX' (line 38): This boolean value can be enabled to continuously loop the Simulated Annealing section to find a tour shorter than a specified distance 'X'. Once a 'good' solution or a goal is found, it is recommended to enable this option so that the program can continuously loop the Simulated Annealing section to optimize the tour and find a more ideal tour.
* 'X' (line 39): This is the integer value that can be set to the current best known solution or a goal solution.

## Help

For more details, please refer to the report file 'Report.pdf'.

## Authors

Riley A. Durbin
radurbin@crimson.ua.edu
The University of Alabama
CS 570: Computer Algorithms

## Version History

* 1.0
    * Initial Release

## Acknowledgments

References:
* Carnegie Mellon Simulated Annealing Article: https://www.cs.cmu.edu/afs/cs.cmu.edu/project/learn-43/lib/photoz/.g/web/glossary/anneal.html
* Numerical Recipes in C, Second Edition

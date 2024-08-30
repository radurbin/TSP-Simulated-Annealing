#include <iostream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <ctime>

using namespace std;

// Original Solution by Riley Durbin for the TSP Competition on 2/16/2024

// compile with g++ TSP_Original.cpp -O3

// GENERAL PROCESS

// - run nearest neighbor from every point to find the best starting point (OPTIONAL)
// - store every alternate neighbor path within a K scalar from the actual nearest neighbor and run each of those alternate paths completely to find the shortest one
// - try to move the vertex at the longest edge to each position in the tour and move it if the total distance decreases
// - try to move every vertex to each position in the tour and move it if the total distance decreases (brute force)
// - run simulated annealing once on this optimal tour to try to remove locality bias
// - run two brute force optimizations again
// - optionally loop the last two points

// OPTIONALLY - loop the simulated annealing and brute force section until the program finds a solution lower than a distance X (for competition)

// scalar for next closest neighbors to check the path for
const double K = 1.5; // should be 1 for 25,0000+ or 10,000+ if there are lots of duplicate edge weights

const string graphName = "g250.graph";
const int numberOfVertices = 250;
const bool runAnyStartingPoint = false; // should be false for 10,000+ because it takes forever and ultimately does not optimize the score all that much
int minStartingPoint = 0; // after running the program once, change the above bool to false and change this value to the best starting point

// COMPETITIVE LOOP FUNCTION
const bool runUntilFindSolutionLessThanX = true;
const int X = 900;
// Best distances for each size and time to run algo once without competitive loop:
// 250 - 1772 (1.543889s) - [K = 1.5, RASP = true]
// 1K - 1612 (3.348727s) - [K = 1.5, RASP = true]
// 5K - 197079 (147.139350s) - [K = 1.5, RASP = true]
// 15K - 191625 (2591.250717s) - [K = 1.5, RASP = false]

// SIMULATED ANNEALING VARIABLES
const double initialTemperature = 1000.0;
const double coolingRate = 0.9;
const int numIterations = 10000000;
// the higher the number of iterations, the better the solution but the slower the process

// global variable to store the distance of the best tour
int nearestNeighborKPercentSolution = 0;

// optional function that runs Nearest Neighbor N times, once for each vertex, to find which starting point returns the most optimal nearest neighbor tour
int whichStartingPoint(vector<vector<int> > &edgeArray)
{
    int bestStartingPoint = 0;
    double bestMinDistance = numeric_limits<double>::infinity();
    bool breakCondition = false;

    int numberOfVertices = edgeArray.size();
    for (int startingPoint = 0; startingPoint < numberOfVertices; startingPoint++)
    {
        vector<int> unvisited; // to keep track of the unvisited vertices
        for (int i = 0; i < numberOfVertices; i++)
        {
            if (i != startingPoint)
            {
                unvisited.push_back(i);
            }
        }
        vector<int> tour; // to keep track of the tour
        tour.push_back(startingPoint);
        int total_distance = 0.0;
        int nearestNeighbor = -1;

        for (int i = 0; i < numberOfVertices - 1; i++)
        {
            int currentPoint = tour.back();
            double min_distance = numeric_limits<double>::infinity();
            nearestNeighbor = -1;
            int NNindex = -1;
            // for each unvisited point
            for (int k = 0; k < unvisited.size(); k++)
            {
                if (currentPoint != unvisited[k])
                {
                    if (edgeArray[currentPoint][unvisited[k]] < min_distance)
                    {
                        min_distance = edgeArray[currentPoint][unvisited[k]];
                        nearestNeighbor = unvisited[k];
                        NNindex = k;
                    }
                }
            }
            tour.push_back(nearestNeighbor);
            total_distance += min_distance;
            // small optimization to speed up this function if the current distance is already longer than the current most optimal tour
            if (total_distance >= bestMinDistance)
            {
                breakCondition = true;
                break;
            }
            unvisited.erase(unvisited.begin() + NNindex);
        }
        if (breakCondition)
        {
            continue;
        }
        total_distance += edgeArray[nearestNeighbor][startingPoint];
        tour.push_back(startingPoint);
        if (total_distance < bestMinDistance)
        {
            bestMinDistance = total_distance;
            bestStartingPoint = startingPoint;
        }
    }

    // return the index of the starting point
    return bestStartingPoint;
}

// runs the nearest neighbor algorithm from a given starting point, finding alternate paths to try as it goes
// this function returns the distance of the nearest neighbor tour
int nearestNeighborFindAlts(vector<vector<int> > &edgeArray, int startingPoint, vector<int> &tour, vector<vector<int> > &altTours, vector<vector<int> > &altUnvisiteds, vector<int> &altTotalDistances)
{
    vector<int> unvisited; // to keep track of the unvisited vertices
    for (int i = 0; i < numberOfVertices; i++)
    {
        if (i != startingPoint)
        {
            unvisited.push_back(i);
        }
    }
    tour.push_back(startingPoint); // tour will start and end with the starting point
    int total_distance = 0.0;
    int nearestNeighbor = -1;

    for (int i = 0; i < numberOfVertices - 1; i++)
    { // for each vertex besides the starting point
        int currentPoint = tour.back();
        double min_distance = numeric_limits<double>::infinity();
        nearestNeighbor = -1;
        int NNindex = -1;
        // for each unvisited point
        for (int k = 0; k < unvisited.size(); k++)
        {
            if (currentPoint != unvisited[k])
            {
                if (edgeArray[currentPoint][unvisited[k]] < min_distance)
                {
                    min_distance = edgeArray[currentPoint][unvisited[k]];
                    nearestNeighbor = unvisited[k];
                    NNindex = k;
                }
            }
        }

        // once we have found the nearest neighbor from the current point,
        // store all of the other neighbors within a K scalar of the nearest neighbor's edge distance
        for (int k = 0; k < unvisited.size(); k++)
        {
            if (currentPoint != unvisited[k] && nearestNeighbor != unvisited[k])
            {
                if (edgeArray[currentPoint][unvisited[k]] < (min_distance * K))
                {
                    vector<int> altTour = tour;
                    altTour.push_back(unvisited[k]);
                    altTours.push_back(altTour);
                    int altTotalDistance = total_distance + edgeArray[currentPoint][unvisited[k]];
                    altTotalDistances.push_back(altTotalDistance);
                    vector<int> altUnvisited = unvisited;
                    altUnvisited.erase(altUnvisited.begin() + k);
                    altUnvisiteds.push_back(altUnvisited);
                    // duplicate tour, unvisited, and total distance
                    // add this neighbor as the next stop in the tour
                    // save these to a vector to be referenced later
                }
            }
        }

        // add nearest neighbor to the tour
        tour.push_back(nearestNeighbor);
        total_distance += min_distance;
        unvisited.erase(unvisited.begin() + NNindex);
    }

    total_distance += edgeArray[nearestNeighbor][startingPoint];
    tour.push_back(startingPoint);

    cout << endl
         << "Any Starting Point Nearest Neighbor Total Distance: " << fixed << total_distance << endl << endl;
    return total_distance;
}

int checkAltTours(vector<vector<int> > &edgeArray, int startingPoint, vector<vector<int> > &altTours, vector<vector<int> > &altUnvisiteds, vector<int> &altTotalDistances, int nearestNeighborSolution)
{
    int whichAlt = -1;
    cout << "Alternate Paths found: " << altTours.size() << endl
         << endl;
    int numAltTours = altTours.size();

    for (int a = 0; a < numAltTours; a++)
    {
        bool breakCondition = false;

        int numberOfVerticesLeft = altUnvisiteds[a].size();
        int nearestNeighbor = -1;

        for (int i = 0; i < numberOfVerticesLeft; i++)
        {
            int currentPoint = altTours[a].back();
            double min_distance = numeric_limits<double>::infinity();
            nearestNeighbor = -1;
            int NNindex = -1;
            // for each unvisited point
            for (int k = 0; k < altUnvisiteds[a].size(); k++)
            {
                if (currentPoint != altUnvisiteds[a][k])
                {
                    if (edgeArray[currentPoint][altUnvisiteds[a][k]] < min_distance)
                    {
                        min_distance = edgeArray[currentPoint][altUnvisiteds[a][k]];
                        nearestNeighbor = altUnvisiteds[a][k];
                        NNindex = k;
                    }
                }
            }

            altTours[a].push_back(nearestNeighbor);
            altTotalDistances[a] += min_distance;
            if (altTotalDistances[a] >= nearestNeighborSolution)
            {
                breakCondition = true;
                break;
            }
            altUnvisiteds.at(a).erase(altUnvisiteds.at(a).begin() + NNindex);
        }
        if (breakCondition == true)
        {
            continue;
        }

        altTotalDistances[a] += edgeArray[nearestNeighbor][startingPoint];
        altTours[a].push_back(startingPoint);

        if (altTotalDistances[a] < nearestNeighborKPercentSolution)
        {
            nearestNeighborKPercentSolution = altTotalDistances[a];
            whichAlt = a;
        }
    }
    return whichAlt;
}

// small brute force function to try to remove the longest edges from the graph
void longestEdgeBruteForce(vector<vector<int> > &edgeArray, vector<int> &tour)
{
    double lastMax = numeric_limits<double>::infinity();
    for (int edgesTested = 0; edgesTested < tour.size() - 2; edgesTested++)
    {
        int maxEdge = -1;
        int from = -1;
        int to = -1;
        int before = -1;
        int after = -1;
        for (int i = 1; i < tour.size() - 1; i++)
        {
            if (edgeArray[tour[i]][tour[i + 1]] > maxEdge && edgeArray[tour[i]][tour[i + 1]] < lastMax)
            {
                maxEdge = edgeArray[tour[i]][tour[i + 1]];
                if (i != 1)
                {
                    before = i - 1;
                }
                else
                {
                    before = -1;
                }
                from = i;
                to = i + 1;
                if (i < tour.size() - 2)
                {
                    after = i + 2;
                }
                else
                {
                    after = -1;
                }
            }
        }
        int vertexToMoveIndex = -1;
        int vertexToMoveValue = -1;
        // once we have the longest edge, check the one after it and before it to figure out which to move
        // check the lengths of the edges before and after using before and after variables
        if (before != -1 && after != -1)
        {
            if (edgeArray[tour[before]][tour[from]] >= edgeArray[tour[to]][tour[after]])
            {
                vertexToMoveIndex = from;
                vertexToMoveValue = tour[from];
            }
            else
            {
                vertexToMoveIndex = to;
                vertexToMoveValue = tour[to];
            }
        }
        else if (before != -1 && after == -1)
        {
            vertexToMoveIndex = from;
            vertexToMoveValue = tour[from];
        }
        else if (before == -1 && after != -1)
        {
            vertexToMoveIndex = to;
            vertexToMoveValue = tour[to];
        }
        else
        {
            // Couldn't find an edge greater than maxEdge and less than lastMax
            break;
        }
        // do some math for the distance if we remove the current edge, then try adding that to each of the edges in between each vertex, subtracting their current edge and adding the two new edges
        // newMinDistance is the current best distance

        int distanceWithoutRemoveVertex = nearestNeighborKPercentSolution - edgeArray[tour[vertexToMoveIndex - 1]][tour[vertexToMoveIndex]] - edgeArray[tour[vertexToMoveIndex]][tour[vertexToMoveIndex + 1]] + edgeArray[tour[vertexToMoveIndex - 1]][tour[vertexToMoveIndex + 1]];
        if (distanceWithoutRemoveVertex > nearestNeighborKPercentSolution)
        {
            // trying to remove vertex results in a greater distance immediately
            // small time optimization
            lastMax = maxEdge;
            continue;
        }
        double newMinDistance = numeric_limits<double>::infinity();
        int insertAfter = -1;
        for (int i = 0; i < tour.size() - 1; i++)
        {
            if (i == vertexToMoveIndex - 1 || i == vertexToMoveIndex)
            {
                continue;
            }
            if (distanceWithoutRemoveVertex + edgeArray[tour[i]][vertexToMoveValue] + edgeArray[vertexToMoveValue][tour[i + 1]] - edgeArray[tour[i]][tour[i + 1]] < nearestNeighborKPercentSolution)
            {
                if (distanceWithoutRemoveVertex + edgeArray[tour[i]][vertexToMoveValue] + edgeArray[vertexToMoveValue][tour[i + 1]] - edgeArray[tour[i]][tour[i + 1]] < newMinDistance)
                {
                    newMinDistance = distanceWithoutRemoveVertex + edgeArray[tour[i]][vertexToMoveValue] + edgeArray[vertexToMoveValue][tour[i + 1]] - edgeArray[tour[i]][tour[i + 1]];
                    insertAfter = i;
                }
            }
        }
        if (insertAfter != -1)
        {
            tour.erase(tour.begin() + vertexToMoveIndex);
            if (insertAfter > vertexToMoveIndex)
            {
                tour.insert(tour.begin() + insertAfter, vertexToMoveValue);
            }
            else
            {
                tour.insert(tour.begin() + insertAfter + 1, vertexToMoveValue);
            }
            // cout << "NEW Total Distance: " << fixed << newMinDistance << endl;
            nearestNeighborKPercentSolution = newMinDistance;
            // cout << endl;
            edgesTested = -1;
            lastMax = numeric_limits<double>::infinity();
        }
        else
        {
            lastMax = maxEdge;
        }
    }
}

// another brute force function to try to move each vertex to each position to optimize the tour
void eachVertexBruteForce(vector<vector<int> > &edgeArray, vector<int> &tour)
{
    for (int edgesTested = 1; edgesTested < tour.size() - 2; edgesTested++)
    {
        int vertexToMoveIndex = edgesTested;
        int vertexToMoveValue = tour[vertexToMoveIndex];

        double distanceWithoutRemoveVertex = nearestNeighborKPercentSolution - edgeArray[tour[vertexToMoveIndex - 1]][tour[vertexToMoveIndex]] - edgeArray[tour[vertexToMoveIndex]][tour[vertexToMoveIndex + 1]] + edgeArray[tour[vertexToMoveIndex - 1]][tour[vertexToMoveIndex + 1]];
        if (distanceWithoutRemoveVertex > nearestNeighborKPercentSolution)
        {
            // trying to remove vertex results in a greater distance immediately
            // small time optimization
            continue;
        }
        double newMinDistance = numeric_limits<double>::infinity();
        int insertAfter = -1;
        for (int i = 0; i < tour.size() - 1; i++)
        {
            if (i == vertexToMoveIndex - 1 || i == vertexToMoveIndex)
            {
                continue;
            }
            if (distanceWithoutRemoveVertex + edgeArray[tour[i]][vertexToMoveValue] + edgeArray[vertexToMoveValue][tour[i + 1]] - edgeArray[tour[i]][tour[i + 1]] < nearestNeighborKPercentSolution)
            {
                if (distanceWithoutRemoveVertex + edgeArray[tour[i]][vertexToMoveValue] + edgeArray[vertexToMoveValue][tour[i + 1]] - edgeArray[tour[i]][tour[i + 1]] < newMinDistance)
                {
                    newMinDistance = distanceWithoutRemoveVertex + edgeArray[tour[i]][vertexToMoveValue] + edgeArray[vertexToMoveValue][tour[i + 1]] - edgeArray[tour[i]][tour[i + 1]];
                    insertAfter = i;
                }
            }
        }
        if (insertAfter != -1)
        {
            tour.erase(tour.begin() + vertexToMoveIndex);
            if (insertAfter > vertexToMoveIndex)
            {
                tour.insert(tour.begin() + insertAfter, vertexToMoveValue);
            }
            else
            {
                tour.insert(tour.begin() + insertAfter + 1, vertexToMoveValue);
            }
            // cout << "NEW Total Distance: " << fixed << newMinDistance << endl;
            nearestNeighborKPercentSolution = newMinDistance;
            // cout << endl;
            edgesTested = 0;
            continue;
        }
        else
        {
            continue;
        }
    }
}

// runs the simulated annealing process, then calls both brute force optimizations
void simulatedAnnealing(vector<vector<int> > &edgeArray, int currentDistance, vector<int> &tour)
{
    double currentTemperature = initialTemperature;
    srand(static_cast<unsigned>(time(nullptr)));

    for (int iteration = 0; iteration < numIterations; ++iteration)
    {
        // Randomly choose between swapping points or reversing a segment
        bool doSwap = rand() % 2 == 0;
        int point1;
        int point2;
        int start;
        int length;
        int newDistance;
        int end;

        // swap two points
        if (doSwap)
        {
            // INDEXES of the points to be swapped
            point1 = rand() % (tour.size() - 3) + 1; // Avoid the first and last elements
            point2 = rand() % (tour.size() - 3) + 1;

            // ensure the two points are distinct
            while (point1 == point2)
            {
                point2 = rand() % (tour.size() - 3) + 1;
            }

            // calculate new distance
            if (abs(point1-point2) == 1 && point1 < point2) { // if points are next to each other
                newDistance = currentDistance - edgeArray[tour[point1 - 1]][tour[point1]] - edgeArray[tour[point2]][tour[point2 + 1]] + edgeArray[tour[point1 - 1]][tour[point2]] + edgeArray[tour[point1]][tour[point2 + 1]];
            }
            else if (abs(point1-point2) == 1 && point1 > point2) {
                newDistance = currentDistance - edgeArray[tour[point2 - 1]][tour[point2]] - edgeArray[tour[point1]][tour[point1 + 1]] + edgeArray[tour[point2 - 1]][tour[point1]] + edgeArray[tour[point2]][tour[point1 + 1]];
            }
            else {
                newDistance = currentDistance - edgeArray[tour[point1 - 1]][tour[point1]] - edgeArray[tour[point1]][tour[point1 + 1]] - edgeArray[tour[point2 - 1]][tour[point2]] - edgeArray[tour[point2]][tour[point2 + 1]] + edgeArray[tour[point2 - 1]][tour[point1]] + edgeArray[tour[point1]][tour[point2 + 1]] + edgeArray[tour[point1 - 1]][tour[point2]] + edgeArray[tour[point2]][tour[point1 + 1]];
            }
        }
        // reverse a segment of the tour
        else
        {
            start = rand() % (tour.size() - 3) + 1;
            length = rand() % (tour.size() - start - 2) + 1; // ensure length is at least 1

            // ensure the segment doesn't wrap around
            if (start + length >= tour.size() - 1)
            {
                length = tour.size() - start - 2;
            }
            end = start + length - 1;

            // calculate new distance
            newDistance = currentDistance - edgeArray[tour[start - 1]][tour[start]] - edgeArray[tour[end]][tour[end + 1]] + edgeArray[tour[start - 1]][tour[end]] + edgeArray[tour[start]][tour[end + 1]];
        }

        // decide whether to accept the new solution
        if (newDistance < currentDistance || exp((currentDistance - newDistance) / currentTemperature) > (rand() % 1000) / 1000.0)
        {
            if (doSwap)
            {
                swap(tour[point1], tour[point2]);
            }
            else
            {
                reverse(tour.begin() + start, tour.begin() + start + length);
            }
            currentDistance = newDistance;
        }

        // update the temperature
        currentTemperature *= coolingRate;
    }
    nearestNeighborKPercentSolution = currentDistance;
    longestEdgeBruteForce(edgeArray, tour);
    eachVertexBruteForce(edgeArray, tour);
}

// main function that finds the optimal tour and distance
void findTour(vector<vector<int> > edgeArray, int startingPoint)
{
    int nearestNeighborSolution = 0;

    vector<int> tour; // to keep track of the tour

    vector<vector<int> > altTours;      // 2D vector holding the alternate tours
    vector<vector<int> > altUnvisiteds; // 2D vector holding the unvisited nodes at this point of the alternate tour
    vector<int> altTotalDistances;     // 2D vector holding the distance of each alternate tour at the point of the branch

    // run nearest neighbor algorithm again, finding and storing alternate paths as it goes
    nearestNeighborSolution = nearestNeighborFindAlts(edgeArray, startingPoint, tour, altTours, altUnvisiteds, altTotalDistances);
    nearestNeighborKPercentSolution = nearestNeighborSolution;

    // perform nearest neighbor with the alternate tours to remove small nearest neighbor bias to find a better alt tour
    int whichAlt = checkAltTours(edgeArray, startingPoint, altTours, altUnvisiteds, altTotalDistances, nearestNeighborSolution);

    if (whichAlt != -1)
    { // found an alternate tour with a better distance than the nearest neighbor solution
        // print the distance and tour of the best alt tour found
        cout << "Best ALT Tour Total Distance: " << fixed << nearestNeighborKPercentSolution << endl << endl;
        // find the longest edge in the graph, test it in each new position and find the smallest, continue this process
        longestEdgeBruteForce(edgeArray, altTours[whichAlt]);

        // try to brute force each vertex into every position, if it would result in a lower overall distance, move it there, then restart trying each vertex
        eachVertexBruteForce(edgeArray, altTours[whichAlt]);

        cout << "When K equals " << K << endl;
        cout << "Optimal Tour Before Annealing Total Distance: " << fixed << nearestNeighborKPercentSolution << endl;

        if (runUntilFindSolutionLessThanX)
        { // if we want to loop simulated annealing to try to beat a provided distance
            int solutionBeforeSimulatedAnnealing = nearestNeighborKPercentSolution;
            vector<int> tourBeforeSimulatedAnnealing = altTours[whichAlt];
            cout << endl
                 << "Looping simulated annealing section until I find a solution less than " << X << endl;
            cout << "Solution before annealing: " << solutionBeforeSimulatedAnnealing << endl;
            simulatedAnnealing(edgeArray, nearestNeighborKPercentSolution, altTours[whichAlt]);
            cout << "Solution after annealing once: " << nearestNeighborKPercentSolution << endl;
            while (nearestNeighborKPercentSolution >= X)
            {
                if (nearestNeighborKPercentSolution < solutionBeforeSimulatedAnnealing) { // if last loop found a shorter solution, run with that solution
                    cout << "This solution is better. Saving and running annealing with tour of distance " << nearestNeighborKPercentSolution << endl;
                    solutionBeforeSimulatedAnnealing = nearestNeighborKPercentSolution;
                    tourBeforeSimulatedAnnealing = altTours[whichAlt];

                    simulatedAnnealing(edgeArray, nearestNeighborKPercentSolution, altTours[whichAlt]);
                    cout << "Found tour of length " << nearestNeighborKPercentSolution << endl;
                }
                else { // if last loop did not find a better solution, run with tour from before that
                    cout << "This solution is worse. Running annealing with tour of distance " << solutionBeforeSimulatedAnnealing << endl;
                    altTours[whichAlt] = tourBeforeSimulatedAnnealing;
                    simulatedAnnealing(edgeArray, solutionBeforeSimulatedAnnealing, altTours[whichAlt]);
                    cout << "Found tour of length " << nearestNeighborKPercentSolution << endl;
                }
            }
        }
        else
        {
            // SIMLUATED ANNEALING SECTION ONLY PERFORMED ONCE
            int currentDistance = nearestNeighborKPercentSolution;
            simulatedAnnealing(edgeArray, currentDistance, altTours[whichAlt]);
        }
        ofstream outputFile2("s_radurbin.txt");
        outputFile2 << fixed << nearestNeighborKPercentSolution << endl;
        cout << endl << "Optimal Tour After Annealing Total Distance: " << fixed << nearestNeighborKPercentSolution << endl;
        for (int z = 0; z < altTours[whichAlt].size(); z++)
        {
            outputFile2 << altTours[whichAlt][z] << " ";
        }
        outputFile2.close();
        cout << endl;
    }
    else
    { // for when K = 1.0000 or there are no alternate tours
        // find the longest edge in the graph, test it in each new position and find the smallest, continue this process
        longestEdgeBruteForce(edgeArray, tour);

        // try to brute force each vertex into every position, if it would result in a lower overall distance, move it there, then restart trying each vertex
        eachVertexBruteForce(edgeArray, tour);

        cout << "When K equals " << K << endl;
        cout << "Optimal Tour Before Annealing Total Distance: " << fixed << nearestNeighborKPercentSolution << endl;

        if (runUntilFindSolutionLessThanX)
        {
            int solutionBeforeSimulatedAnnealing = nearestNeighborKPercentSolution;
            vector<int> tourBeforeSimulatedAnnealing = tour;
            cout << endl
                 << "Looping simulated annealing section until I find a solution less than " << X << endl;
            cout << "Solution before annealing: " << solutionBeforeSimulatedAnnealing << endl;
            simulatedAnnealing(edgeArray, nearestNeighborKPercentSolution, tour);
            cout << "Solution after annealing once: " << nearestNeighborKPercentSolution << endl;
            while (nearestNeighborKPercentSolution >= X)
            {
                if (nearestNeighborKPercentSolution < solutionBeforeSimulatedAnnealing) { // if last loop found a shorter solution, run with that solution
                    cout << "This solution is better. Saving and running annealing with tour of distance " << nearestNeighborKPercentSolution << endl;
                    solutionBeforeSimulatedAnnealing = nearestNeighborKPercentSolution;
                    tourBeforeSimulatedAnnealing = tour;
                    simulatedAnnealing(edgeArray, nearestNeighborKPercentSolution, tour);
                    cout << "Found tour of length " << nearestNeighborKPercentSolution << endl;
                }
                else { // if last loop did not find a better solution, run with tour from before that
                    cout << "This solution is worse. Running annealing with tour of distance " << solutionBeforeSimulatedAnnealing << endl;
                    tour = tourBeforeSimulatedAnnealing;
                    simulatedAnnealing(edgeArray, solutionBeforeSimulatedAnnealing, tour);
                    cout << "Found tour of length " << nearestNeighborKPercentSolution << endl;
                }
            }
        }
        else
        {
            // SIMLUATED ANNEALING SECTION ONLY PERFORMED ONCE
            int currentDistance = nearestNeighborKPercentSolution;
            simulatedAnnealing(edgeArray, currentDistance, tour);
        }
        ofstream outputFile2("s_radurbin.txt");
        outputFile2 << fixed << nearestNeighborKPercentSolution << endl;
        cout << endl << "Optimal Tour After Annealing Total Distance: " << fixed << nearestNeighborKPercentSolution << endl;
        for (int z = 0; z < tour.size(); z++)
        {
            outputFile2 << tour[z] << " ";
        }
        outputFile2.close();
        cout << endl;
    }
}

int main(int argc, char *argv[])
{

    using chrono::duration;
    using chrono::duration_cast;
    using chrono::high_resolution_clock;
    using chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    // open graph file and read into 2D vector o
    ifstream graphFile(graphName);
    if (!graphFile)
    {
        cerr << "Error opening graph file." << endl;
        return 1;
    }

    vector<vector<int> > edgeArray;
    for (int i = 0; i < numberOfVertices; i++)
    {
        vector<int> oneDvector;
        edgeArray.push_back(oneDvector);
    }

    // initializes the 2D vector to all 0s
    for (int i = 0; i < numberOfVertices; i++)
    {
        for (int j = 0; j < numberOfVertices; j++)
        {
            edgeArray[i].push_back(0.0);
        }
    }

    // reads edge weights from graph file into 2D vector
    for (int i = 0; i < numberOfVertices; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int distance;
            graphFile >> distance;
            edgeArray[i][j] = distance;
            edgeArray[j][i] = distance;
        }
    }
    graphFile.close();

    // if specified, run minStartingPoint to find the optimal nearest neighbor path to start with
    // should only be used with relatively small graphs
    if (runAnyStartingPoint == true)
    {
        minStartingPoint = whichStartingPoint(edgeArray);
    }

    // run algorithm
    findTour(edgeArray, minStartingPoint);
    auto t2 = high_resolution_clock::now();

    duration<double, milli> ms_double = t2 - t1;

    // print time taken by algorithm in seconds
    std::cout << (ms_double.count() / 1000.0) << "s\n";

    return 0;
}
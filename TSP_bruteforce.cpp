#include <iostream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

// Simple Brute Force solution for the TSP by Riley Durbin

// function to calculate the factorial of a number
double factorial(int x)
{
    double fact = (double)x;
    for (int i = x - 1; i > 0; i--)
    {
        fact *= i;
    }
    return fact;
}

// slightly optimized brute force algorithm with assumption that 0 is start and end
double bruteForce(int tourPerm[], vector<vector<double> > edgeArray, double min_distance)
{
    // total the distance of each possible tour
    int numVertices = edgeArray.size() - 1;
    double total_distance = 0.0;
    total_distance += edgeArray[0][tourPerm[0]];
    for (int i = 0; i < numVertices - 1; i++)
    {
        total_distance += edgeArray[tourPerm[i]][tourPerm[i + 1]];
        if (total_distance >= min_distance)
        {
            return numeric_limits<double>::infinity();
        }
    }
    total_distance += edgeArray[tourPerm[numVertices - 1]][0];

    return total_distance;
}

int main(int argc, char *argv[])
{

    using chrono::duration;
    using chrono::duration_cast;
    using chrono::high_resolution_clock;
    using chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    ifstream graphFile("g250.graph");
    if (!graphFile)
    {
        cerr << "Error opening graph file." << endl;
        return 1;
    }

    int numVertices = 13;
    // 10 - 0.446524s
    // 11 - 4.870739s
    // 12 - 59.107826s

    vector<vector<double> > edgeArray;
    for (int i = 0; i < numVertices; i++)
    {
        vector<double> oneDvector;
        edgeArray.push_back(oneDvector);
    }

    // initializes the 2d vector to all 0s
    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            edgeArray[i].push_back(0.0);
        }
    }

    for (int i = 0; i < numVertices; i++)
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

    int tourPerm[numVertices - 1];
    for (int i = 1; i < numVertices; i++)
    {
        tourPerm[i - 1] = i;
    }
    int optimalTour[numVertices - 1];
    for (int x = 0; x < numVertices - 1; x++)
    {
        optimalTour[x] = tourPerm[x];
    }
    double min_distance = bruteForce(tourPerm, edgeArray, numeric_limits<double>::infinity());
    cout << "Brute forcing. Running all " << fixed << factorial(numVertices - 1) << " possibilities." << endl;
    // loop calculate distance function for every permutation of the tour
    for (int i = 0; i < factorial(numVertices - 1) - 1; i++)
    {
        if (i % 10000000 == 0)
        { // periodically update the terminal with the smallest distance
            cout << i << ": Optimal distance so far is " << fixed << min_distance << endl;
        }
        next_permutation(tourPerm, tourPerm + numVertices - 1);
        double permDistance = bruteForce(tourPerm, edgeArray, min_distance);
        if (permDistance < min_distance)
        {
            min_distance = permDistance;
            // update the terminal when a new min distance is found
            cout << "New optimal tour at " << i << ", distance = " << fixed << min_distance << endl;
            cout << "0 ";
            for (int x = 0; x < numVertices - 1; x++)
            {
                optimalTour[x] = tourPerm[x];
                cout << optimalTour[x] << " ";
            }
            cout << "0" << endl;
        }
    }

    // print optimal tour and distance
    cout << "Optimal tour, distance = " << fixed << min_distance << endl;
    cout << "0 ";
    for (int x = 0; x < numVertices - 1; x++)
    {
        cout << optimalTour[x] << " ";
    }
    cout << "0" << endl;

    auto t2 = high_resolution_clock::now();

    duration<double, milli> ms_double = t2 - t1;

    // print time taken by algorithm in seconds
    std::cout << (ms_double.count() / 1000.0) << "s\n";

    return 0;
}
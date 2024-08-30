#include <iostream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

// nearest neighbor algorithm by Riley Durbin with the assumption that the starting point is the 0th vertex
void nearestNeighbor(vector<vector<double> > edgeArray) {
    int numVertices = edgeArray.size();
    vector<int> unvisited; // to keep track of the unvisited vertices
    for (int i = 0; i < numVertices - 1; i++) {
        unvisited.push_back(i + 1);
    }
    vector<int> tour; // to keep track of the tour
    tour.push_back(0);
    double total_distance = 0.0;
    int nearestNeighbor = -1;

    for (int i = 0; i < numVertices - 1; i++) {
        int currentPoint = tour.back();
        double min_distance = numeric_limits<double>::infinity();
        nearestNeighbor = -1;
        int NNindex = -1;
        // for each unvisited point
        for (int k = 0; k < unvisited.size(); k++) {
            if (currentPoint != unvisited[k]) {
                if (edgeArray[currentPoint][unvisited[k]] < min_distance) {
                    min_distance = edgeArray[currentPoint][unvisited[k]];
                    nearestNeighbor = unvisited[k];
                    NNindex = k;
                }
            }
        }
        // found nearest neighbor, add them to the tour and remove them from the unvisited list
        tour.push_back(nearestNeighbor);
        total_distance += min_distance;
        unvisited.erase(unvisited.begin() + NNindex);
    }

    // add 0 to the back of the tour
    total_distance += edgeArray[nearestNeighbor][0];
    tour.push_back(0);

    // print the total distance and tour
    cout << "Total Distance: " << fixed << total_distance << endl;
    cout << "Nearest Neighbor Tour: ";
    for (int i = 0; i < tour.size(); i++) {
        cout << tour[i] << " ";
    }
    cout << endl;
}

int main(int argc, char *argv[]) {

    using chrono::duration;
    using chrono::duration_cast;
    using chrono::high_resolution_clock;
    using chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    ifstream graphFile("g15000.graph");
    if (!graphFile) {
        cerr << "Error opening graph file." << endl;
        return 1;
    }
    int numVertices = 15000;
    // 250 - 0.009215s
    // 1K - 0.128493s
    // 5K - 4.443482s
    // 15K - 50.943934s

    vector<vector<double> > edgeArray;
    for (int i = 0; i < numVertices; i++) {
        vector<double> oneDvector;
        edgeArray.push_back(oneDvector);
    }

    // initializes the 2d vector to all 0s
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < numVertices; j++) {
            edgeArray[i].push_back(0.0);
        }
    }

    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j <= i; j++) {
            int distance;
            graphFile >> distance;
            edgeArray[i][j] = distance;
            edgeArray[j][i] = distance;
        }
    }
    graphFile.close();

    nearestNeighbor(edgeArray);

    auto t2 = high_resolution_clock::now();

    duration<double, milli> ms_double = t2 - t1;

    // print time taken by algorithm in seconds
    std::cout << (ms_double.count() / 1000.0) << "s\n";

    return 0;
}
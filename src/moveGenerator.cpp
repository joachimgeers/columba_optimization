#include "moveGenerator.h"

using namespace std;

//
// UNIFORM WEIGHT MOVE GENERATOR
//

UniformWeightMoveGenerator::UniformWeightMoveGenerator(vector<double> stepSizes, int numWeights = 0) {
    // moves that are dependent on step sizes need to be generated now we have the stepSizes
    for (auto step: stepSizes) {
        for (auto move: stepMoves) {
            moveList.push_back(Move(move, step));
        }

        for (auto move: positionMoves) {
            if (numWeights <= 0) {
                moveList.push_back(Move(move, step));
            } else {
                for (int i = -1; i < numWeights; i++) {
                    moveList.push_back(Move(move, step, i));
                }
            }
        }
    }
}

//
// UNIFORM SEED MOVE GENERATOR
//

UniformSeedMoveGenerator::UniformSeedMoveGenerator(vector<double> stepSizes, int numSeeds = 0) {
    // moves that are dependent on step sizes need to be generated now we have the stepSizes
    for (auto step: stepSizes) {
        for (auto move: stepMoves) {
            moveList.push_back(Move(move, step));
        }

        for (auto move: positionMoves) {
            if (numSeeds <= 0) {
                moveList.push_back(Move(move, step));
            } else {
                for (int i = -1; i < numSeeds; i++) {
                    moveList.push_back(Move(move, step, i));
                }
            }
        }
    }
}

//
// WEIGHTED MOVE GENERATOR
//

int WeightedMoveGenerator::binarySearch(vector<double> list, double element) {
    int lower = 0;
    int upper = list.size() - 1;
    int index;
    while (upper - lower > 1) {
        index = (lower + upper)/2;
        auto el = list[index];
        if (el < element) {
            lower = index;
        } else if (el > element) {
            upper = index;
        } else {
            return index;
        }
    }

    if (list[upper] < element) return upper;
    else return lower;
}

WeightedMoveGenerator::WeightedMoveGenerator(vector<double> weightStepSizes, vector<double> seedStepSizes, int numSeeds = 0, int numWeights = 0, double startWeight = 50) {
    vector<Move> seedMoveList;
    
    if (numSeeds == 0) {
        seedMoveList = vector<Move>(); 
        cout << "lol?" << endl;
    } else {
        seedMoveList = UniformSeedMoveGenerator(seedStepSizes, numSeeds).getMoveList();
    }
    
    vector<Move> weightMoveList = UniformWeightMoveGenerator(weightStepSizes, numWeights).getMoveList();
    total = 0;
    int index = 0;
    weights.resize(seedMoveList.size() + weightMoveList.size(), startWeight);

    // the move list contains the cumulative weight of moves
    // as to allow binary searching
    for (auto seedMove: seedMoveList) {
        moveList.push_back(seedMove);
        weights[index] = total; // cumulative weight at index
        total += startWeight; // add weight to the cumulative weight
        moveToIndex[seedMove.getId()] = index++;
    }

    for (auto weightMove: weightMoveList) {
        moveList.push_back(weightMove);
        weights[index] = total;
        total += startWeight;
        moveToIndex[weightMove.getId()] = index++;
    }
}

void WeightedMoveGenerator::setWeights(map<string, double> weightMap) {
    double tmpTotal = 0, tmp;

    for (size_t i = 0; i < weights.size(); i++) {
        auto move = moveList[i].getId();
        if (weightMap.count(move)) {
            // new weight is updated gradually
            tmp = weightMap[move] * 0.6 + getWeight(i) * 0.4;
        } else {
            tmp = getWeight(i);
        }
        weights[i] = tmpTotal; // cumulative weight at index
        tmpTotal += tmp; // add weight to the cumulative weight
    }

    total = tmpTotal; // new total weight
}
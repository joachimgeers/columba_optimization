#include "move.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cstdarg>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
Move::Move(enum LocalMove move, double stepSize, int position, bool validStepSize) : lower(0), upper(max(1., stepSize)), stepSize(stepSize), id(move), position(position), validStepSize(validStepSize){
    this->movePointer = this->moveMap[move];
}

Solution Move::apply(Solution solution) const {
    Solution next = (this->*movePointer)(solution);
    vector<double> seeds;
    for (auto seed: next.getSeeds()) {
        seeds.push_back(round(seed * 1000.0)/1000.0);
    }
    vector<int> weights;

    for (auto weight: next.getWeights()) {
        weights.push_back(min(weight, solution.getMaxWeight()));
    }
    next = next.setSeeds(seeds).setWeights(weights);
    auto v = next.isValid() ? next : solution;
    return v;
}

Solution Move::updateRandom(Solution solution, int amount) const {
    vector<int> weights = solution.getWeights();
    int r = rand() % weights.size();
    if (weights[r] <= -amount) weights = updateAll(solution, -(amount - weights[r] + 1)).getWeights();
    weights[r] += amount;
    return solution.setWeights(weights);
}

Solution Move::updateRandomRelative(Solution solution, double factor) const {
    vector<int> weights = solution.getWeights();
    int r;
    if (position < 0 || position - weights.size() >= 0) r = rand() % weights.size();
    else r = position;
    int amount = ceil(weights[r] * factor);
    if (amount == 0) amount = 1;
    
    // if the weight to be updated would end up being < 0, we increase all weights by the difference between the two + 1
    // so that our weight ends up as the lower bound (1)
    if (weights[r] <= -amount) weights = updateAllRelative(solution, -(amount - weights[r] + 1)).getWeights();
    weights[r] += amount;
    
    return solution.setWeights(weights);
}

Solution Move::updateAllRelative(Solution solution, double factor) const {
    vector<int> weights = solution.getWeights();
    int amount;
    for (size_t i = 0; i < weights.size(); i++) {
        amount = ceil(weights[i] * factor) + 1;
        if (weights[i] <= -amount) weights[i] = 1;
        else weights[i] += amount;
    }
    return solution.setWeights(weights);

}

Solution Move::updateAll(Solution solution, int amount) const {
    vector<int> weights = solution.getWeights();
    for (size_t i = 0; i < weights.size(); i++) {
        if ((int) weights[i] <= -amount) continue;
        else weights[i] += amount;
    }
    return solution.setWeights(weights);

}

Solution Move::raiseOuter(Solution solution) const {
    vector<int> weights = solution.getWeights();
    int amount =  ceil((weights[0] + weights[weights.size() - 1])/2 * stepSize) + 1;
    weights[0] += amount;
    weights[weights.size() - 1] += amount;
    return solution.setWeights(weights);
}

Solution Move::updateWeightsBySeeds(Solution solution) const {
    auto seeds = solution.getSeeds();
    auto weights = solution.getWeights();
    if (seeds.size() <= 0) return solution;
    
    seeds.insert(seeds.begin(), 0);
    seeds.push_back(1);
    vector<double> partSizes(seeds.size());
    auto weightSum = weights.front() + weights.back();
    
    partSizes.front() = seeds[1];
    for (size_t i = 1; i < seeds.size() - 1; i++) {
        weightSum += weights[i];
        partSizes[i] = (seeds[i+1] - seeds[i-1])/2;
    }
    partSizes.back() = (1 - seeds[seeds.size() - 2]);

    auto sortedPartSizes = partSizes;
    sort(sortedPartSizes.begin(), sortedPartSizes.end());

    auto median = sortedPartSizes[partSizes.size()/2];
    weightSum *= stepSize;
    
    for (size_t i = 0; i < partSizes.size(); i++) {
        auto p = partSizes[i]/median;
        if (p < 1) p = -1/p;
        weights[i] = max(weights[i] + p * weightSum, 1.);
    }

    return solution.setWeights(weights);;
}

Solution Move::lowerCenter(Solution solution) const {
    vector<int> weights = solution.getWeights();
    int amount;
    for (size_t i = 1; i < weights.size() - 1; i++) {
        amount = ceil(weights[i] * stepSize);
        if (amount == 0) amount = 1;
        // update all other weights so that our weight takes the lower bound, in case we were smaller than 0
        if (weights[i] < amount) weights = updateAll(solution, amount - weights[i] + 1).getWeights();
        weights[i] -= amount;
    }
    return solution.setWeights(weights);
}

Solution Move::normalize(Solution solution) const {
    vector<int> weights = solution.getWeights();
    int minimum = weights[std::distance(weights.begin(), min_element(weights.begin(), weights.end()))];
    for (size_t i = 0; i < weights.size(); i++) {
        weights[i] = ceil(weights[i]/minimum);
    }
    return solution.setWeights(weights);
}

Solution Move::multiplyAll(Solution solution) const {
    vector<int> weights = solution.getWeights();
    for (size_t i = 0; i < weights.size(); i++) {
        weights[i] *= stepSize;
    }
    return solution.setWeights(weights);
}

Solution Move::raiseRandom(Solution solution) const {
    return updateRandomRelative(solution, this->stepSize);
}

Solution Move::raiseAll(Solution solution) const {
    return updateAllRelative(solution, this->stepSize);
}

Solution Move::lowerRandom(Solution solution) const {
    return updateRandomRelative(solution, -this->stepSize);
}
Solution Move::lowerAll(Solution solution) const {
    return updateAllRelative(solution, -this->stepSize);
}

Solution Move::swapRandom(Solution solution) const {
    vector<int> weights = solution.getWeights();
    int r = rand() % weights.size();
    int rr = rand() % weights.size();

    if (rr == r) {
        rr = (rr + 1) % weights.size();
    }

    int temp = weights[r];
    weights[r] = weights[rr];
    weights[rr] = temp;
    return solution.setWeights(weights);
}

Solution Move::exchangeRandom(Solution solution) const {
    vector<int> weights = solution.getWeights();
    int r = rand() % weights.size();
    int rr = rand() % weights.size();
    
    if (rr == r) {
        rr = (rr + 1) % weights.size();
    }

    int amount = ceil((weights[r] + weights[rr])/2 * stepSize) + 1;

    if (amount > weights[r] || amount > weights[rr]) {
        weights = updateAll(solution, amount).getWeights();
    }
    weights[r] += amount;
    weights[rr] -= amount;
    return solution.setWeights(weights);
}

Solution Move::copyRandom(Solution solution) const {
    vector<int> weights = solution.getWeights();
    int r = rand() % weights.size();
    int rr = rand() % weights.size();
    weights[r] = weights[rr];
    return solution.setWeights(weights);
}

Solution Move::flip(Solution solution) const {
    vector<int> weights = solution.getWeights();
    reverse(weights.begin(), weights.end());

    return solution.setWeights(weights);
}

Solution Move::moveSeed(Solution solution, size_t index, double amount) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    bool isMove = false;
    // index needs to be valid
    if (index >= 0 && index < seeds.size()) {
        if (amount < 0) {
            double neighbour = index == 0 ? 0. : seeds[index-1];
            isMove = (seeds[index] + amount) > neighbour;
        } else {
            double neighbour = index == seeds.size() ? 1. : seeds[index+1];
            isMove = (seeds[index] + amount) < neighbour;
        }
    }

    if (isMove) seeds[index] += amount;
    return solution.setSeeds(seeds);
}

Solution Move::moveSeedLeft(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    size_t r;
    if (position < 0 || position - seeds.size() >= 0) r = rand() % seeds.size();
    else r = position;
    return moveSeed(solution, r, -stepSize);

}
Solution Move::moveSeedRight(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    size_t r;
    if (position < 0 || position - seeds.size() >= 0) r = rand() % seeds.size();
    else r = position;
    return moveSeed(solution, r, stepSize);
}

Solution Move::moveAllSeedsLeft(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    for (size_t i = 0; i < seeds.size(); i++) {
        solution = moveSeed(solution, i, -stepSize);
    }
    return solution;
}
Solution Move::moveAllSeedsRight(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    for (int i = seeds.size() - 1; i >= 0; i--) {
        solution = moveSeed(solution, i, stepSize);
    }
    return solution;
}

Solution Move::moveNeighboursCloser(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;

    if (position < 0 || position - seeds.size() >= 0) {
        size_t n = rand() % 2;
        size_t r = rand() % seeds.size();

        if (r == 0) n = 1;
        if (r == seeds.size() - 1) n = 0;

        r += n;

        solution = moveSeed(solution, r, -stepSize);
        solution = moveSeed(solution, r-1, stepSize);
    } else {
        solution = moveSeed(solution, position, stepSize);
        solution = moveSeed(solution, (position + 1) % seeds.size(), -stepSize);
    }

    return solution;

}

Solution Move::moveNeighboursFurther(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;

    if (position < 0 || position - seeds.size() >= 0) { 
        size_t n = rand() % 2; // randomly select left or right neighbour
        size_t r = rand() % seeds.size();

        // if we're on the side we only have on neighbour so select that one 
        if (r == 0) n = 1;
        if (r == seeds.size() - 1) n = 0;

        // if n == 0: increase our current position, decrease the position on our left
        // if n == 1: increase our right neighbour (r+n) and decrease our current position
        r += n;

        solution = moveSeed(solution, r, stepSize);
        solution = moveSeed(solution, r-1, -stepSize);
    } else {
        solution = moveSeed(solution, position, -stepSize);
        solution = moveSeed(solution, (position + 1) % seeds.size(), stepSize);
    }

    return solution;
}

Solution Move::moveAllInward(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    for (size_t i = 0; i < seeds.size(); i++) {
        // if less than center move right, otherwise left
        if (seeds[i] < 0.5) solution = moveSeed(solution, i, stepSize);
        else solution = moveSeed(solution, i, -stepSize);
    }
    return solution;
}

Solution Move::moveAllOutward(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    for (size_t i = 0; i < seeds.size(); i++) {
        // if less than center move left, otherwise right
        if (seeds[i] < 0.5) solution = moveSeed(solution, i, -stepSize);
        else solution = moveSeed(solution, i, stepSize);
    }
    return solution;
}

Solution Move::recenterSeed(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    size_t c = rand() % seeds.size();

    // if our seed is the first, our left side is 0
    // if our seed is the last, our left side is 1
    double l = c == 0 ? 0. : seeds[c - 1];
    double r = c == seeds.size() - 1 ? 1. : seeds[c + 1];

    // move seed to the new position (left + factor * diff) (0 <= factor <= 1)
    double value = l + (r - l) * 0.5;
    
    // round the value to 2 decimals
    seeds[c] = floor(value * 100.0) / 100.0;
    
    return solution.setSeeds(seeds);

}

Solution Move::invertSeeds(Solution solution) const {
    vector<double> seeds = solution.getSeeds();
    if (seeds.size() <= 0) return solution;
    size_t numSeeds = seeds.size() - 1;
    double tmp;
    for (size_t i = 0; i < floor(seeds.size()/2); i++) {
        tmp = seeds[i];
        seeds[i] = 1 - seeds[numSeeds - i];
        seeds[numSeeds - i] = 1 - tmp;
    }

    return solution.setSeeds(seeds);
}

Solution Move::none(Solution solution) const {
    return solution;
}


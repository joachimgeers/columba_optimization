#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>
#include <cstdlib>
#include <boost/functional/hash/hash.hpp>
#include "trainingSearchstrategy.h"
#include "inputProcessing.h"

using namespace std;

typedef struct TimeCost {
    size_t averageTime;
    size_t minimumTime;
    size_t cost;
    double sd; // standard deviation
} TimeCost;

//
// SOLUTION CLASS
//
// A solution represents a current set of weights and seeds
// We also store the fm index and search strategy so we can compute the cost
//

class Solution {
    private:
        vector<int> weights;
        vector<double> seeds;
        string readsFile;
        FMIndex *mapper;
        TrainingSearchStrategy* strategy;
        length_t ED;
        double minSeedDist;
        size_t chunkSize;
        int threads;
        int maxWeight;
        vector<pair<string, string>> reads;

    public:
        Solution(vector<int> weights, 
                 vector<double> seeds, 
                 string readsFile, 
                 FMIndex* mapper, TrainingSearchStrategy* strategy, 
                 length_t ED, double minSeedDist, size_t chunkSize, int threads, int maxWeight, vector<pair<string, string>> reads)
        : weights(weights), seeds(seeds), readsFile(readsFile), mapper(mapper), strategy(strategy), ED(ED), minSeedDist(minSeedDist), 
        chunkSize(chunkSize), threads(threads), maxWeight(maxWeight), reads(reads){};
        

        /**
         * Get the weights of this solution
         **/
        vector<int> getWeights() const {
            return weights;
        }
        
        /**
         * Get the seeds of this solution
         **/
        vector<double> getSeeds() const {
            return seeds;
        }

        /**
         * Get the maximum allowed weight
         **/
        int getMaxWeight() const {
            return maxWeight;
        }


        // SETTERS
        // note that under any circumstance we might want to ignore a new solution
        // which is why we always return a new solution instead of adjusting the old solution
        
        /**
         * Get a new solution with a new set of weights
         * @param weights the new set of weights
         **/
        Solution setWeights(vector<int> weights) {
            return Solution(weights, seeds, readsFile, mapper, strategy, ED, minSeedDist, chunkSize, threads, maxWeight, reads);
        }

        /**
         * Get a new solution with a new set of seeds
         * @param seeds the new set of seeds
         **/
        Solution setSeeds(vector<double> seeds) {
            return Solution(weights, seeds, readsFile, mapper, strategy, ED, minSeedDist, chunkSize, threads, maxWeight, reads);
        }


        /**
         * Compute and return cost of this solution
         * The cost is computed by matching the set of samples related to this solution
         **/
        size_t getCost() const;


        size_t getSequentialCost(vector<pair<string, string>> reads) const;



        /**
         * Compute and return the time cost of this solution
         * the time cost is represented 
         **/
        TimeCost getTimeCost(int iterations) const;


        /**
         * Check whether the solution is valid
         * A solution is valid when the seeds are appropriately positioned, and all weights are > 0
         **/
        bool isValid() {
            bool validSeeds = true;
            
            if (seeds.size() > 0) {
                validSeeds = (seeds.front() > minSeedDist * 1.5 && 1 - seeds.back() > minSeedDist * 1.5);
                // check if no seed violates the minimum seed distance
                for (size_t i = 0; i < seeds.size() - 1; i++) {
                    if ((seeds[i + 1] - seeds[i]) <= minSeedDist) validSeeds = false;
                }
            }


            // check if no weight is negative or 0
            bool validWeights = true;
            for (auto i: weights) if (i <= 0 || i > maxWeight) validWeights = false;

            return !weights.empty() &&
                   // the seeds on the side are not specified
                   // so check if the number of seeds and weights make sense
                   validSeeds && validWeights;
        }

        bool operator==(const Solution& other) const {
            return other.getWeights() == getWeights() && other.getSeeds() == this->getSeeds();
        }

        bool operator!=(const Solution& other) const{
            return other.getWeights() != getWeights() || other.getSeeds() != this->getSeeds();
        }
};

// SOLUTIONHASH CLASS
//
// This class is used to compute the hashvalue of a solution
// for this we just use teh boost library and combine the hash of the weight and seed vectors
// We need to be able to hash our solutions for the tabu list
//
class SolutionHash {
    public:
        size_t operator()(const Solution& solution) const {
            size_t seed = 0;
            auto weights = solution.getWeights();
            size_t weightHash = boost::hash_range(weights.begin(), weights.end());
            auto seeds = solution.getSeeds();
            size_t seedHash = boost::hash_range(seeds.begin(), seeds.end());
            boost::hash_combine(seed, weightHash);
            boost::hash_combine(seed, seedHash);
            return seed;
        }
};

#endif
#ifndef MOVE_GENERATOR_H
#define MOVE_GENERATOR_H

#include "move.h"

using namespace std;

/**
 * 
 *  A move generator selects a move out of a list of potential moves
 * 
**/
class MoveGenerator {
    protected:

    public:
        virtual ~MoveGenerator(){};
        virtual vector<Move> getMoveList() const = 0;
        virtual Move* getMove() = 0;
};

class StaticMoveGenerator : public MoveGenerator {
    private:
        Move move = Move(NONE, 0);
    public:
        /**
         * Generate a move
         **/
        virtual Move* getMove() {
            return &move;
        }

        /**
         * @return a copy of the list of moves
         */
        virtual vector<Move> getMoveList() {
            return {move};
        }
};

class UniformWeightMoveGenerator : public MoveGenerator {
    private:
        // define all the moves that aren't dependent on the step sizes
        vector<Move> moveList{
            Move(MULTIPLY_ALL, 2., false),
            Move(MULTIPLY_ALL, 3., false),
            Move(FLIP, NO_VALUE, false), 
            Move(COPY, NO_VALUE, false), 
            Move(NORMALIZE, NO_VALUE, false), 
            Move(SWAP_RANDOM, NO_VALUE, false)
        };
        
        // define all the moves that depend on stepsizes
        vector<LocalMove> stepMoves = {RAISE_ALL, LOWER_ALL, EXCHANGE_RANDOM, RAISE_OUTER, LOWER_CENTER, UPDATE_WEIGHT_BY_SEEDS};
        vector<LocalMove> positionMoves = {RAISE_RANDOM, LOWER_RANDOM};
    public:

        /**
         * A generator that uniformly selects between all weight moves
         * @param stepSizes all the step sizes that moves requiring step sizes must be initialized on (only use multiple if you are using static step sizes)
         * @param numWeights the number of weights in the solution that moves are generated for
         */
        UniformWeightMoveGenerator(vector<double> stepSizes, int numWeights);

        /**
         * @return a copy of the move list
         */
        virtual vector<Move> getMoveList() const {
            return moveList;
        }

        /**
         * randomly select a move from the move list
         * @return the selected move
         */
        virtual Move* getMove() {
            Move* move = &(*select_randomly(moveList.begin(), moveList.end()));
            return move;
        }
};


class UniformSeedMoveGenerator : public MoveGenerator {
    private:
        // define all the moves that aren't dependent on the step sizes
        vector<Move> moveList{
            Move(INVERT_ALL_SEEDS, NO_VALUE, false),
            Move(RECENTER_SEED, NO_VALUE, false),
        };
        
        // define all the moves that depend on stepsizes
        vector<LocalMove> stepMoves = {MOVE_ALL_SEEDS_LEFT, MOVE_ALL_SEEDS_RIGHT, MOVE_ALL_INWARD, MOVE_ALL_OUTWARD};
        vector<LocalMove> positionMoves = {MOVE_SEED_LEFT, MOVE_SEED_RIGHT, MOVE_NEIGHBOURS_CLOSER, MOVE_NEIGHBOURS_FURTHER};

    public:

        /**
         * A generator that uniformly selects between all seed moves
         * @param stepSizes all the step sizes that moves requiring step sizes must be initialized on (only use multiple if you are using static step sizes)
         * @param numWeights the number of seeds in the solution that moves are generated for
         */
        UniformSeedMoveGenerator(vector<double> stepSizes, int numSeeds);


        /**
         * @return a copy of the move list
         */
        virtual vector<Move> getMoveList() const {
            return moveList;
        }
    
        /**
         * randomly select a move from the move list
         * @return the selected move
         */
        virtual Move* getMove() {
            Move* move = &(*select_randomly(moveList.begin(), moveList.end()));
            return move;
        }


};


class UniformMoveGenerator : public MoveGenerator {
    private:
        // define all the moves that aren't dependent on the step sizes
        vector<Move> moveList;

    public:

        /**
         * A generator that uniformly selects between all defined moves
         * @param weightStepSizes the step sizes that the weight moves must be initialized with
         * @param seedStepSizes   the step sizes that the seed moves must be initialized with
         * @param numSeeds        the number of  seeds in the solution that moves must be generated for
         * @param numWeights      the number of  weights in the solution that moves must be generated for 
         */
        UniformMoveGenerator(vector<double> weightStepSizes, vector<double> seedStepSizes, int numSeeds = 0, int numWeights = 0) {
            moveList = UniformWeightMoveGenerator(weightStepSizes, numWeights).getMoveList();
            auto seedMoves = UniformSeedMoveGenerator(seedStepSizes, numSeeds).getMoveList();
            moveList.insert(moveList.end(), seedMoves.begin(), seedMoves.end());
        }


        /**
         * @return a copy of the move list
         */
        virtual vector<Move> getMoveList() const {
            return moveList;
        }


        /**
         * randomly select a move from the move list
         * @return the selected move
         */
        virtual Move* getMove() {
            Move* move = &(*select_randomly(moveList.begin(), moveList.end()));
            return move;
        }


};


class WeightedMoveGenerator : public MoveGenerator {
    private:
        // define all the moves that aren't dependent on the step sizes
        map<string, int> moveToIndex;
        // list of moves
        vector<Move> moveList;
        // list of cumulative weights (index corresponding with movelist)
        vector<double> weights;
        // the total cumulative weight
        double total = 0;
        // position in move list of the last generated move
        int lastIndex = -1;

        /**
         * simple binary search, preferring the lower bound, as to match the interval
         * @param list the list of elements to search in
         * @param element the element to be found
         * @return the index of the move matching the given cumulative weight
         */
        int binarySearch(vector<double> list, double element);

    public:

        /**
         * A generator that uniformly selects between all defined moves
         * @param weightStepSizes the step sizes that the weight moves must be initialized with
         * @param seedStepSizes   the step sizes that the seed moves must be initialized with
         * @param numSeeds        the number of  seeds in the solution that moves must be generated for
         * @param numWeights      the number of  weights in the solution that moves must be generated for 
         * @param maxWeight       the highest allowed weight for each move during initialization
         */
        WeightedMoveGenerator(vector<double> weightStepSizes, vector<double> seedStepSizes, int numSeeds, int numWeights, double maxWeight);


        /**
         * set the weights in accordance to the given weightMap (moves not included will retain the same weight if they are already in the generator)
         * @param weightMap a map of moveId::weight
         */
        void setWeights(map<string, double> weightMap);

        /**
         * set the weights in accordance to a vector (each index corresponding to the moveList)
         * @param newWeights the new vector of weights
         */
        void setWeights(vector<double> newWeights) {
            assert (newWeights.size() == moveList.size());
            total = 0;
            for (size_t i = 0; i < newWeights.size(); i++) {
                weights[i] = total;
                total += newWeights[i];
                
            }
        }

        /**
         * Set the weight of all moves to the new given weight
         * @param weight the new weight
         */
        void setAllWeights(double weight) {
            total = 0;
            for (size_t i = 0; i < weights.size(); i++) {
                weights[i] = total;
                total += weight;
                
            }
        }

        /**
         * Set the weight at a given index to a new weight
         * @param index the index of the move
         * @param newWeight the new weight
         */
        void setWeight(int index, double newWeight) {
            auto diff = newWeight - getWeight(index);
            total += diff;
            for (size_t i = index; i < weights.size() - 1; i++) {
                weights[i + 1] += diff;
            }
        }

        /**
         * set the weight of the last used move
         * @param newWeight the new weigth
         */
        void updateLastWeight(double newWeight) {
            if (lastIndex > -1 && lastIndex < ((int) weights.size())) setWeight(lastIndex, max(newWeight, 1.));
        }

        /**
         * get the weight of the move at the given index
         * @param index the index of the move in the move list
         * @return the weight of the move
         */
        double getWeight(int index) {
            if (index <= -1 || index >= ((int) weights.size())) return -1;
            if (index == (int) weights.size() - 1) 
                return total - weights[index];
            //else
            return weights[index + 1] - weights[index];
        }

        /**
         * get the weight of the last used move
         * @return the weight of the last used move
         */
        double getLastWeight() {
            return getWeight(lastIndex);
        }

        /**
         * @return a copy of the move list
         */
        virtual vector<Move> getMoveList() const {
            return moveList;
        }


        /**
         * select a move from the move list in accordance to their weights
         * @return the selected move
         */
        virtual Move* getMove() {
            double nextMove = ((double)rand() / RAND_MAX) * total;
            int index = binarySearch(weights, nextMove);
            return &moveList[index];
        }


};


#endif
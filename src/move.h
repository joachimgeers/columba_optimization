#ifndef MOVE_H
#define MOVE_H

#define NO_VALUE 0
#include <cstddef>
#include <map>
#include <iostream>
#include <algorithm>
#include "solution.h"
#include "sampling.h"
using namespace std;

enum StepSizeState {
    BETTER,
    EQUAL,
    WORSE
};

// an enum representing the possible moves
enum LocalMove {
    // NONE MOVE
    NONE,

    // WEIGHT MOVES
    RAISE_RANDOM, 
    LOWER_RANDOM, 
    EXCHANGE_RANDOM, 
    SWAP_RANDOM,
    RAISE_ALL, 
    LOWER_ALL,
    FLIP,
    COPY,
    RAISE_OUTER,
    LOWER_CENTER,
    UPDATE_WEIGHT_BY_SEEDS,
    NORMALIZE,
    MULTIPLY_ALL,

    // SEED MOVES
    MOVE_SEED_LEFT,
    MOVE_SEED_RIGHT,
    MOVE_ALL_SEEDS_LEFT,
    MOVE_ALL_SEEDS_RIGHT,
    
    MOVE_NEIGHBOURS_CLOSER,
    MOVE_NEIGHBOURS_FURTHER,

    MOVE_ALL_INWARD,
    MOVE_ALL_OUTWARD,

    RECENTER_SEED,
    INVERT_ALL_SEEDS,
};

class Move;

// a pointer to a move function
typedef Solution (Move::*MovePointer)(Solution solution) const;


// ============================================================================
// CLASS MOVE
// ============================================================================

class Move {
    private:
        // a pointer to the specified move
        MovePointer movePointer;
        double lower;
        double upper;
        StepSizeState formerState = BETTER;
        // the stepSize of the move (if applicable)
        double stepSize;
        // the enum corresponding with the move
        LocalMove id;
        // which weight seed to be updated (if applicable)
        int position;
        bool validStepSize;
        // a mapping of the possible moves to their functions
        map<LocalMove, MovePointer> moveMap {
            // NONE MOVE
            {NONE, &Move::none},
            // WEIGHT MOVES
            {RAISE_RANDOM, &Move::raiseRandom}, {LOWER_RANDOM, &Move::lowerRandom}, {SWAP_RANDOM, &Move::swapRandom},
            {RAISE_ALL, &Move::raiseAll}, {LOWER_ALL, &Move::lowerAll}, 
            {EXCHANGE_RANDOM, &Move::exchangeRandom}, {FLIP, &Move::flip}, {COPY, &Move::copyRandom}, 
            {RAISE_OUTER, &Move::raiseOuter}, {LOWER_CENTER, &Move::lowerCenter}, {UPDATE_WEIGHT_BY_SEEDS, &Move::updateWeightsBySeeds},
            {NORMALIZE, &Move::normalize}, {MULTIPLY_ALL, &Move::multiplyAll},
            // SEED MOVES
            {MOVE_SEED_LEFT, &Move::moveSeedLeft}, {MOVE_SEED_RIGHT, &Move::moveSeedRight}, 
            {MOVE_ALL_SEEDS_LEFT, &Move::moveAllSeedsLeft}, {MOVE_ALL_SEEDS_RIGHT, &Move::moveAllSeedsRight},
            {MOVE_NEIGHBOURS_CLOSER, &Move::moveNeighboursCloser}, {MOVE_NEIGHBOURS_FURTHER, &Move::moveNeighboursFurther},
            {MOVE_ALL_INWARD, &Move::moveAllInward}, {MOVE_ALL_OUTWARD, &Move::moveAllOutward},
            {RECENTER_SEED, &Move::recenterSeed}, {INVERT_ALL_SEEDS, &Move::invertSeeds}
        };

        // translation of moves to strings for printing purposes
        map<LocalMove, string> localMoveToStr {
            // NONE MOVE
            {NONE, "None"},
            // WEIGHT MOVES
            {RAISE_RANDOM, "RAISE WEIGHT"}, {LOWER_RANDOM, "LOWER WEIGHT"},  {SWAP_RANDOM, "SWAP WEIGHTS"},
            {RAISE_ALL, "RAISE ALL"}, {LOWER_ALL, "LOWER ALL"}, 
            {EXCHANGE_RANDOM, "EXCHANGE WEIGHTS"}, {FLIP, "FLIP"}, {COPY, "COPY"}, 
            {RAISE_OUTER, "RAISE OUTER"}, {LOWER_CENTER, "LOWER CENTER"}, {UPDATE_WEIGHT_BY_SEEDS, "UPDATE WEIGHTS BY SEEDS"},
            {NORMALIZE, "NORMALIZE"}, {MULTIPLY_ALL, "MULTIPLY ALL"},
            // SEED MOVES
            {MOVE_SEED_LEFT, "MOVE SEED LEFT"}, {MOVE_SEED_RIGHT, "MOVE SEED RIGHT"}, 
            {MOVE_ALL_SEEDS_LEFT, "MOVE ALL SEEDS LEFT"}, {MOVE_ALL_SEEDS_RIGHT, "MOVE ALL SEEDS RIGHT"},
            {MOVE_NEIGHBOURS_CLOSER, "MOVE NEIGHBOURS CLOSER"}, {MOVE_NEIGHBOURS_FURTHER, "MOVE NEIGHBOURS FURTHER"},
            {MOVE_ALL_INWARD, "MOVE ALL INWARD"}, {MOVE_ALL_OUTWARD, "MOVE ALL OUTWARD"},
            {RECENTER_SEED, "RE-CENTER SEED"}, {INVERT_ALL_SEEDS, "INVERT SEEDS"}
            };


        // WEIGHT MOVES


        /**
         * Increase/decrease a random weight in the given [solution] with a specified [amount]
         * @param solution The solution to be updated
         * @param amount the amount with which to update the weight
         **/
        Solution updateRandom(Solution solution, int amount) const;

        /**
         * Increase/decrease all weights in the given [solution] with a specified [amount]
         * @param solution The solution to be updated
         * @param amount the amount with which to update the weights
         **/
        Solution updateAll(Solution solution, int amount) const;
        

        /**
         * Increase/decrease a random weight in the given [solution] by the current weight * [factor]
         * @param solution The solution to be updated
         * @param factor the factor used to determine the amount used to update the weight (current weight * factor)
         **/
        Solution updateRandomRelative(Solution solution, double factor) const;
        
        /**
         * Increase/decrease all weights in the given [solution] by the current weight * [factor]
         * @param solution The solution to be updated
         * @param factor the factor used to determine the amount used to update the weights (current weight * factor)
         **/
        Solution updateAllRelative(Solution solution, double factor) const;


        //
        // the raise/decrease functions make use of either updaterandom or update random relative
        // but the amount/factor are inferred from the stepsize of the move
        //

        /**
         * Increase a random weight in the solution with the current stepsize
         * @param solution the solution to be updated
         **/
        Solution raiseRandom(Solution solution) const;

        /**
         * Increase all weights in the solution with the current stepsize
         * @param solution the solution to be updated
         **/
        Solution raiseAll(Solution solution) const;

        /**
         * Decrease a random weight in the solution with the current stepsize
         * @param solution the solution to be updated
         **/
        Solution lowerRandom(Solution solution) const;

        /**
         * Decrease all weights in the solution with the current stepsize
         * @param solution the solution to be updated
         **/
        Solution lowerAll(Solution solution) const;

        /**
         * Decrease a random weight, increase another (e.g. weights x, y -> x - stepSize, y + stepSize)
         * @param solution the solution to be updated
         **/
        Solution exchangeRandom(Solution solution) const;

        /**
         * swap two random weight values with each other
         * @param solution the solution to be updated
         **/
        Solution swapRandom(Solution solution) const;

        /**
         * flip the order of the weights (e.g. (a, b, ... y, z) -> (z, y, ..., b, a))
         * @param solution the solution to be updated
         **/
        Solution flip(Solution solution) const;

        /**
         * copy a random weight to another random weight (e.g. (a, b) -> (b, b))
         * @param solution the solution to be updated
         **/
        Solution copyRandom(Solution solution) const;

        /**
         * lower the values of all partitions except for the ones on the sides
         * @param solution the solution to be updated
         **/
        Solution lowerCenter(Solution solution) const;

        /**
         * Increase the weights of the outer partitions with the current stepsize
         * @param solution the solution to be updated
         **/
        Solution raiseOuter(Solution solution) const;


        /**
         * Increase/decrease the weights in accordance to their part sizes
         * @param solution the solution to be updated
         **/
        Solution updateWeightsBySeeds(Solution solution) const;

        /**
         * Divide all weights by the smallest weight
         * @param solution the solution to be updated
         **/
        Solution normalize(Solution solution) const;

        /**
         * Multiply all weights with the stepsize
         * @param solution the solution to be updated
         **/
        Solution multiplyAll(Solution solution) const;

        // SEED MOVES


        /**
         * Update a seed at a specified index with a specified amount
         * @param solution the solution to be updated
         **/
        Solution moveSeed(Solution solution, size_t index, double amount) const;

        /**
         * Decrease the value of a seed (move left)
         * @param solution the solution to be updated
         **/
        Solution moveSeedLeft(Solution solution) const;


        /**
         * Increase the value of a seed (move right)
         * @param solution the solution to be updated
         **/
        Solution moveSeedRight(Solution solution) const;

        /**
         * Decrease the value of all seeds (move left)
         * @param solution the solution to be updated
         **/
        Solution moveAllSeedsLeft(Solution solution) const;

        /**
         * Increase the value of all seeds (move right)
         * @param solution the solution to be updated
         **/
        Solution moveAllSeedsRight(Solution solution) const;

        /**
         * Bring two neighbouring seeds closer to each other (e.g. (a, b) -> (a + stepsize, b - stepsize))
         * @param solution the solution to be updated
         **/
        Solution moveNeighboursCloser(Solution solution) const;

        /**
         * Bring two neighbours further apart (e.g. (a, b) -> (a - stepsize, b + stepsize))
         * @param solution the solution to be updated
         **/
        Solution moveNeighboursFurther(Solution solution) const;

        /**
         * Move all seeds closer to the central seed (increase left of center, decrease right of center)
         * @param solution the solution to be updated
         **/
        Solution moveAllInward(Solution solution) const;


        /**
         * Move all seeds closer to the sides (decrease left of center, increase right of center)
         * @param solution the solution to be updated
         **/
        Solution moveAllOutward(Solution solution) const;

        /**
         * move a seed to a percentile (specified by stepSize) of its two neighbours
         * @param solution the solution to be updated
         **/
        Solution recenterSeed(Solution solution) const;

        /**
         * invert the proportions of the seeds (a, b -> 1-b, 1-a)
         * @param solution the solution to be updated
         **/
        Solution invertSeeds(Solution solution) const;

        // NONE MOVE

        /**
         * Do nothing
         * @param solution the solution to not be updated
         **/
        Solution none(Solution solution) const;
    
    public:
        // ----------------------------------------------------------------------------
        // CONSTRUCTOR
        // ----------------------------------------------------------------------------
        /**
         * @param move the enum representation of the move, the id
         * @param stepSize the step size for the move
         **/
        Move(enum LocalMove move, double stepSize, int position = -1, bool validStepSize=true);

        double getStepSize() const{
            return stepSize;
        }

        bool hasStepSize() const {
            return validStepSize;
        }

        void updateStepSize(enum StepSizeState state) {
            auto tmp = stepSize;
            if (upper < stepSize) upper = stepSize;

            if (lower > stepSize) lower = stepSize;
            switch(state) {
                case WORSE:
                    if (formerState == WORSE || upper <= lower) {
                        stepSize /= 1.8;
                    } else {
                        stepSize = (upper + lower) * 0.25;
                    }
                    upper = tmp;
                    break;
                case EQUAL:
                    if (formerState == EQUAL || upper <= lower) {
                        stepSize *= 1.8;
                    } else {
                        stepSize = (upper + lower) * 0.75;
                    }
                    lower = tmp;
                    break;
                case BETTER:
                    stepSize *= 1.01;
                    break;
                default:
                    break;
            }
            if (upper < lower) {
                tmp = lower;
                lower = upper * 0.99;
                upper = tmp * 1.01;
            }
            formerState = state;
        }
        // ----------------------------------------------------------------------------
        // MOVE FUNCTIONS
        // ----------------------------------------------------------------------------
        
        /**
         * apply the move on a solution (the move function corresponding to the ID will be used here)
         * @param solution the solution to be updated
         **/
        Solution apply(Solution solution) const;

        /**
         * convert a move to a string representation
         **/
        string getId(bool showStepSize=false) {
            std::ostringstream stringStream;
            string pos = position > -1 ? " " + to_string(position) : "";
            string step = showStepSize ? ", " + to_string(stepSize) : "";
            stringStream << "[" << localMoveToStr[this->id] << pos << step << "]";
            return stringStream.str();
        }
};

#endif
#ifndef ACCEPTOR
#define ACCEPTOR
#include "move.h"


// ACCEPTOR CLASS
//
// An acceptor represents a local search strategy
// this acceptor serves as an interface

class Acceptor {

    protected:
        Solution currentSolution;
        size_t currentCost;

    public:
        virtual ~Acceptor(){};
        Acceptor(Solution currentSolution) : currentSolution(currentSolution) {
            currentCost = currentSolution.getCost();
        }

        /**
         * Check if a current solution is accepted by the acceptor; 
         * if it is accepted the current solution will be overwritten by the new solution
         * @param solution the new solution to be compared
         * @return whether or not the solution was accepted
         **/
        virtual bool accept(Solution solution) = 0;
        
        /**
         * Reheats the search strategy; serves as a soft restart by increasing the current search bound
         **/
        virtual void reheat() {};

        /**
         * Spend an idle iteration; this makes it easier to work with outside criteria that proc an idling
         **/
        virtual void idle() {};

        /**
         * Checks whether an acceptor has met its criterium for reheating
         * @returns whether or not the acceptor should reheat
         **/
        virtual bool isReheat() const { return false; };

        /**
         * Return the current solution
         **/
        Solution getCurrentSolution() {
            return currentSolution;
        }

        /**
         * Return the cost of the current solution
         **/
        size_t getCurrentCost() {
            return currentCost;
        }
};

// SIMPLE DESCENT ACCEPTOR CLASS
//
// This class implements the simple descent/gradient descent local search heuristic
// basically just always accept when a solution is better
class SimpleDescentAcceptor: public Acceptor {
    public:
        /**
         * Check if a current solution is accepted by the acceptor; 
         * if it is accepted the current solution will be overwritten by the new solution
         * @param solution the new solution to be compared
         * @return whether or not the solution was accepted
         * The simple acceptor always accepts if the new solution is an improvement
         **/
        bool accept(Solution solution) {
            double cost = solution.getCost();
            if (cost < currentCost) {
                currentCost = cost;
                return true;
            }
            
            return false;
        };
};

class StepCountingHillClimbingAcceptor : public Acceptor {
    private:
        // the threshold
        double threshold;
        // the number of accepts before we update the threshold
        const int acceptsUntilUpdate;
        // the current number of accepts with our current threshold
        int accepts;
        // the maximum number of iterations we can idle before reheating
        const int maxIdleIterations;
        // the number of iterations we have idled inbetween accepts
        int idleIterations;
        // determines the strength of the reheat
        const double reheatFactor;

    public:

        StepCountingHillClimbingAcceptor(Solution currentSolution, int acceptsUntilUpdate, int maxIdleIterations, double reheatFactor) 
            : Acceptor(currentSolution), acceptsUntilUpdate(acceptsUntilUpdate), maxIdleIterations(maxIdleIterations), reheatFactor(reheatFactor) {
                threshold = currentCost;
                accepts = 0;
                idleIterations = 0;
        }

        /**
         * reheat when we overexceed the maximum idle iterations
         * @returns whether or not the reheat criterium is met
         **/
        bool isReheat() const {
            return idleIterations > maxIdleIterations;
        }


        /**
         * upon reheat, increase the threshold with the reheat factor (also resets idle iterations and accepts)
         **/
        void reheat() {
            idleIterations = 0;
            accepts = 0;
            threshold *= reheatFactor;
            //cout << "NEW THRESHOLD AFTER REHEATING: " << threshold << endl;
        }

        /**
         * Idling increases the number of idle iterations
         **/
        void idle() {
            idleIterations++;
        }


        /**
         * Check if a current solution is accepted by the acceptor; 
         * if it is accepted the current solution will be overwritten by the new solution
         * @param solution the new solution to be compared
         * @return whether or not the solution was accepted
         * The step counting hill climbing acceptor checks the current cost against a set threshold
         * after a set amount of accepts it lowers the threshold
         * if the acceptor does not accept, it idles
         **/
        bool accept(Solution solution) {
            size_t nextCost = solution.getCost();
            // we accept if our next cost is better than the threshold
            /*cout << "weights: ";
            for (auto weight: solution.getWeights()) {
                cout << weight << "; ";
            }
            cout << "seeds: ";
            for (auto seed: solution.getSeeds()) {
                cout << seed << "; ";
            }
            cout << "cost: " << nextCost << endl;*/
            if (threshold > nextCost) {
                accepts++;
                idleIterations = 0;
                currentSolution = solution;
                currentCost = nextCost;
                
                // update the threshold to the current cost when we have enough accepts
                // reset the accept count
                if (accepts > acceptsUntilUpdate) {
                    accepts = 0;
                    //cout << "New Threshold" << endl;
                    threshold = nextCost;
                }
                return true;
            }
            // else
            idle();
            return false;
        }
};


#endif
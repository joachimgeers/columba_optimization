#ifndef TRAININGSEARCHSTRATEGY_H
#define TRAININGSEARCHSTRATEGY_H
#include "searchstrategy.h"

using namespace std;

class TrainingSearchStrategy: public CustomSearchStrategy {
    private:
        vector<int> trainingWeights;
        vector<double> trainingSeeds;
        vector<double> trainingBegins;

        const std::vector<double> getSeedingPositions(const int& numParts, const int& maxScore) const override {
            assert(supportsMaxScore[maxScore - 1]);
            return trainingSeeds;
        };

        const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
            assert(supportsMaxScore[maxScore - 1]);
            return trainingWeights;
        };

        const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
            assert(supportsMaxScore[maxScore - 1]);
            for (auto t: trainingBegins) {
                cout << t << "; ";
            }
            cout << endl;
            return trainingBegins;
        };



        virtual const std::vector<double> getUniformSeeds(const int& maxScore) const {
            int numParts = this->calculateNumParts(maxScore);
            double u = 1.0 / (numParts - 1);
            std::vector<double> s;
            for (int i = 1; i < numParts - 1; i++) {
                s.push_back(i * u);
            }
            return s;
        }

        virtual const std::vector<double> getUniformBegins(const int& maxScore) const {
            int numParts = this->calculateNumParts(maxScore);
            double u = 1.0 / (numParts);
            std::vector<double> s;
            for (int i = 1; i < numParts; i++) {
                s.push_back(i * u);
            }
            return s;
        }

        virtual const std::vector<int> getUniformWeights(const int& maxScore) const {
            int numParts = this->calculateNumParts(maxScore);
            std::vector<int> w;
            for (int i = 0; i < numParts; i++) {
                w.push_back(1);
            }
            return w;
        }
    
    public:
        TrainingSearchStrategy(FMIndex& index, std::string pathToFolder,
                            PartitionStrategy p = DYNAMIC,
                            DistanceMetric metric = EDITOPTIMIZED,
                            bool verbose = false)
            : CustomSearchStrategy(index, pathToFolder, p, metric) {


            };

        void initializeUniformParameters(const int& maxScore) {
            this->trainingSeeds = getUniformSeeds(maxScore);
            this->trainingWeights = getUniformWeights(maxScore);
            this->trainingBegins = getUniformBegins(maxScore);
        };

        void setTrainingSeeds(vector<double> seeds) {
            this->trainingSeeds = seeds;
        };

        void setTrainingWeights(vector<int> weights) {
            this->trainingWeights = weights;
        };

        void setTrainingBegins(vector<double> begins) {
            this->trainingBegins = begins;
        };


        vector<double> getTrainingSeeds() {
            return this->trainingSeeds;
        };

        vector<int> getTrainingWeights() {
            return this->trainingWeights;
        };

        vector<double> getTrainingBegins() {
            return this->trainingBegins;
        }

        length_t getSumExactMatchRanges(const std::string& pattern, length_t maxED) const {
            Counters counters;
            vector<Substring> parts;
            int numParts = calculateNumParts(maxED);
            vector<SARangePair> exactMatchRanges(numParts);
            partition(pattern, parts, numParts, maxED, exactMatchRanges, counters);
            length_t sum = 0;
            for (auto range: exactMatchRanges) {
                sum += range.width();
        }

    return sum;
}


};

#endif
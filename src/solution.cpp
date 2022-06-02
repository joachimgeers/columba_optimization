#include "move.h"
#include "searchstrategy.h"
#include "solution.h"
#include <chrono>
#include <math.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <set>
#include <string.h>
#include "inputProcessing.h"

using namespace std;

/*size_t Solution::getCost() const{
    this->strategy->setTrainingSeeds(this->seeds);
    this->strategy->setTrainingWeights(this->weights);

    size_t targetChunk = chunkSize;
    FastQReader myReader(readsFile, "");
    myReader.startReaderThread(targetChunk, targetChunk * threads);

    OutputWriter myWriter("testReads.sam");
    myWriter.startTraining(threads);

    // start worker threads
    int ed = ED;
    vector<thread> worker(threads);
    for (size_t i = 0; i < worker.size(); i++) {
        worker[i] = thread(InputProcessing::threadEntry, std::ref(myReader), std::ref(myWriter),
                           strategy, ref(ed));
    }

    // wait for worker threads to finish
    for_each(worker.begin(), worker.end(), mem_fn(&thread::join));

    // wait for reader and writer to finish
    myReader.joinReaderThread();
    myWriter.joinWriterThread();
    return outputNodes;
}*/

size_t Solution::getCost() const{
    Counters counters;
    this->strategy->setTrainingSeeds(this->seeds);
    this->strategy->setTrainingBegins(this->seeds);
    this->strategy->setTrainingWeights(this->weights);

    //size_t totalNodes = 0;

    // just run the approximate pattern matching on our reads
    // the total nodes is our cost here
    // we could also do this with timings but that'd be more costly to measure
    for (unsigned int i = 0; i < reads.size(); i += 2) {
        const auto& p = reads[i];

        auto originalPos = p.first;
        string read = p.second;
        string revCompl = reads[i + 1].second;

        auto matches = strategy->matchApprox(read, ED, counters);
        //totalNodes += mapper->getNodes();

        // do the same for the reverse complement
        auto matchesRevCompl = strategy->matchApprox(revCompl, ED, counters);
        //totalNodes += mapper->getNodes();
    }

    return counters.nodeCounter;

}

TimeCost Solution::getTimeCost(int iterations) const {
    this->strategy->setTrainingSeeds(this->seeds);
    this->strategy->setTrainingWeights(this->weights);

    size_t cost = 0;
    double averageDuration = 0;
    vector<double> durations;
    double diff;
    for (int i = 0; i < iterations; i++) {
        auto startTime = chrono::high_resolution_clock::now();
        cost += this->getCost();
        auto stopTime = chrono::high_resolution_clock::now();
        diff = chrono::duration_cast<chrono::microseconds>(stopTime - startTime).count();
        durations.emplace_back(diff);
        averageDuration += diff;
    }

    averageDuration /= iterations;

    double var = 0;
    for (auto d : durations) {
        diff = d - averageDuration;
        var += diff * diff;
    }

    var /= (iterations - 1);

    TimeCost timeCost;
    timeCost.averageTime = averageDuration;
    timeCost.minimumTime = *min_element(durations.begin(), durations.end());
    timeCost.sd = sqrt(var);
    timeCost.cost = cost/iterations;

    return timeCost;
}


#include <iostream>
#include <iomanip>
#include <sstream>
#include <bits/stdc++.h>
#include "acceptor.h"
#include "searchstrategy.h"
#include "fmindex.h"
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include "inputProcessing.h"
#include "moveGenerator.h"


using namespace std;
typedef pair<string, string> readPair;

enum StepSizeState intToState(int i) {
    if (i < 0) return BETTER;
    if (i > 0) return WORSE;
    return EQUAL;
}

string getTimeString(int timeInSeconds) {
    int hours, minutes, seconds;
    seconds = timeInSeconds;
    
    hours = timeInSeconds / 3600;
    seconds -= hours * 3600;

    minutes = seconds/60;
    seconds -= minutes * 60;
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2) << seconds;
    return ss.str();
}

bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

void writeLine(ostream& outFile, Solution trainingSolution, Solution testSolution, string timeStamp) {

    auto resultTraining = trainingSolution.getCost();
    auto resultTest = testSolution.getCost();

    outFile << timeStamp << "; " //<< resultTraining.minimumTime << ", "
                                 //<< resultTraining.averageTime  << ", "
                                 //<< resultTraining.sd << ", "
                                 << resultTraining << "; "
                                 //<< resultTest.minimumTime << ", "
                                 //<< resultTest.averageTime << ", "
                                 //<< resultTest.sd << ", "
                                 << resultTest - resultTraining << "; ";

    for (auto w: trainingSolution.getWeights()) {
        outFile << w << " ";
    }
    outFile << "; ";

    for (auto s: trainingSolution.getSeeds()) {
        outFile << s << " ";
    }

    outFile << endl;
}


void showTrainingUsage() {
    cout << "this is supposed to print usage information" << endl;
}

int parseArguments(int argc, char* argv[], string& baseFile, string& readsFile,
                   string& searchscheme, PartitionStrategy& pStrat, int& ed, length_t& saSF,
                   DistanceMetric& metric, string& pathToCustomScheme, int& duration, 
                   int& sampleSize, int& threshold, string& outFileName, bool& intermediateWrites, int& seed, int& threads, 
                   int& maxWeight, int& reheatIterations, double& reheatFactor, int& acceptsUntilUpdate) {
        // process optional arguments
    int requiredArguments = 2;
    InputProcessing::parseBaseArguments(argc, argv, 
                baseFile, readsFile, 
                searchscheme, pStrat, ed, saSF, metric, pathToCustomScheme, 
                vector<string>({"--threshold", "--reheat-iterations", "--reheat-factor", "--accepts", "--max-weight", "-t", "--threads", "-d", "--duration", "-z", "--sample-size", "--out-file", "--intermediate-output", "--seed"}), showTrainingUsage);
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];
        if (arg == "-d" || arg == "--duration") {
            duration = stoi(argv[++i]);
        } else if (arg == "-z" || arg == "--sample-size") {
            sampleSize = stoi(argv[++i]);
        } else if (arg == "--out-file") {
            outFileName = argv[++i];
        } else if (arg == "--intermediate-output") {
            intermediateWrites = to_bool(argv[++i]);
        } else if (arg == "--seed") {
            seed = stoi(argv[++i]);
        } else if (arg == "-t" || arg == "--threads") {
            threads = stoi(argv[++i]);
        } else if (arg == "--max-weight") {
            maxWeight = stoi(argv[++i]);
        } else if (arg == "--reheat-factor") {
            reheatFactor = stod(argv[++i]);
        } else if (arg == "--reheat-iterations") {
            reheatIterations = stoi(argv[++i]);
        } else if (arg == "--accepts") {
            acceptsUntilUpdate = stoi(argv[++i]);
        } else if (arg == "--threshold") {
            threshold = stoi(argv[++i]);
        }
    }

    return requiredArguments;
}

void writeSample(string outFileName, vector<readPair> sample) {
    ofstream outFile(outFileName);
    if (outFile.is_open()) {
        for (auto s: sample) {
            std::string sepLine(s.second.size(), '#');
            outFile << "@" << s.first << endl;
            outFile << s.second << endl;
            outFile << "+" << endl;
            outFile << sepLine << endl;
        }
        outFile.close();
    } else {
        cout << "Couldn't open file " + outFileName + " for writing.";
        exit(1337);
    }
}

double getMinimumSeedDistance(size_t minimumPatternSize, int wordSize) {
    return (1. * wordSize)/minimumPatternSize;
}

size_t getShortestPatternSize(vector<readPair> reads) {
    size_t minPattern = reads[0].second.length();
    
    for (auto read: reads) {
        if (minPattern > read.second.length()) minPattern = read.second.length();
    }

    return minPattern;
}

vector<readPair> getReadSample(vector<readPair> reads, length_t threshold, int seed, length_t sampleSize, TrainingSearchStrategy* strategy, length_t maxED) {
    vector<readPair> readSample;
    vector<int> indices;
    
    for (length_t i = 0; i < reads.size(); i+=2) {
        const auto& p = reads[i];

        auto originalPos = p.first;
        string read = p.second;
        string revCompl = reads[i + 1].second;
        auto semr = max(strategy->getSumExactMatchRanges(read, maxED), strategy->getSumExactMatchRanges(revCompl, maxED));
        if (semr > threshold) {
            indices.push_back(i/2);
        }
    }


    //vector<int> indices(reads.size()/2); // generate n/2 indices (we consider forward/backward to be one entry)
    //iota(indices.begin(), indices.end(), 0); // initialize the indices

    shuffle(indices.begin(), indices.end(), std::default_random_engine(seed)); // shuffle indices

    // we only sample over half the sampleSize, because each entry is forward + backward read
    for (length_t sampleRound = 0; sampleRound < sampleSize/2; sampleRound++) {
        if (sampleRound < indices.size()){
            auto index = indices[sampleRound] * 2;
            readSample.push_back(reads[index]);
            readSample.push_back(reads[index+1]);
        }
    }

    return readSample;
}

void removeFile(string file) {
    if(remove(file.c_str()) != 0)
        perror("Could not clear training sample file");
    else
        puts("Training sample file successfully deleted");
}

void printResults(Solution training, size_t trainingCost, Solution test, size_t testCost, string name) {
    cout << name << " weights: ";

    for (int i: training.getWeights()) {
        cout << i << " ";
    }

    cout << " | " << name << " seeds: ";
    for (double i: training.getSeeds()) {
        cout << i << " ";
    }

    cout << " | cost on training sample: " << trainingCost;

    cout << " | cost on test sample: ";
    cout << testCost - trainingCost << " | total cost: " << testCost << endl;
}

void printResultLine(Solution training, int currentCost, int bestCost, int seconds) {
    cout << "\33[2K\r" << 
            "[" <<  getTimeString(seconds) << "] " << 
            "WEIGHTS: { ";
    
    for (auto weight: training.getWeights()) cout << weight << " ";
    
    cout << "} | SEEDS: { ";
    
    for (auto seed: training.getSeeds()) cout << seed << " ";
    
    cout << "}" << 
            " | CURRENT COST: " << currentCost << 
            " | BEST COST: " << bestCost;
    
    cout.flush();
}

void writeIntermediateResults(ofstream& outFile, Solution training, Solution test, int seconds) {
    auto time = getTimeString(seconds);
    cout << "\33[2K\r" << 
            "\u001b[31m" << 
            "[" <<  time << "] " << 
            "RUNTIME PAUSED: WRITING INTERMEDIATE RESULT TO FILE (THIS WILL TAKE A WHILE)" <<
            "\u001b[0m" <<
            "\r";
    cout.flush();
    
    if (outFile.is_open()) {
        writeLine(outFile, training, test, time);
    }
}

bool pairCompare(pair<string, double>& a, pair<string, double>& b){
    return a.second < b.second;
}

double minMap(map<string, double>& improvementMap, map<string, int>& moveMap) {
    // Declare vector of pairs
    vector<pair<string, double> > resultPairs;
  
    // Copy key-value pair from Map
    // to vector of pairs
    for (auto& it : improvementMap) {
        it.second /= moveMap[it.first];
        resultPairs.push_back(it);
    }
  
    // Sort using comparator function
    sort(resultPairs.begin(), resultPairs.end(), pairCompare);
    
    return resultPairs[0].second;
}

void printMap(map<string, double>& improvementMap, map<string, int>& moveMap){
  
    // Declare vector of pairs
    vector<pair<string, double> > resultPairs;
  
    // Copy key-value pair from Map
    // to vector of pairs
    for (auto& it : improvementMap) {
        resultPairs.push_back(it);
    }
  
    // Sort using comparator function
    sort(resultPairs.begin(), resultPairs.end(), pairCompare);
  
    // Print the sorted value
    for (auto& it : resultPairs) {
  
        cout << it.first << ' '
             << it.second << ' ' << moveMap[it.first] << endl;
    }
}


int main(int argc, char* argv[]) {

    //
    // PROCESS INPUT
    //
    InputProcessing::printWelcome("Weight finder for Columba");

    PartitionStrategy pStrat = STATIC;
    DistanceMetric metric = EDITOPTIMIZED;

    vector<readPair> reads;

    string outFileName = "optimization_out.txt", baseFile, readsFile, customFile = "",
           trainingSampleFile = "tmpSample.fastq";

    string searchscheme = "kuch1";

    length_t saSF;

    int ed, 
        duration = 60, // in seconds
        sampleSize = 2000, 
        seed = 1337, // seed to subdue sampling randomness
        threads = 1, maxWeight = 1000, reheatIterations=5000, acceptsUntilUpdate=2, threshold=0, chunkSize = 10000;
    
    double reheatFactor = 1.005;
    bool intermediateWrites = false; // whether or not we should write out results on the test set

    parseArguments(argc, argv, 
        baseFile, readsFile, 
        searchscheme, pStrat, ed, saSF, metric, customFile, 
        duration, sampleSize, threshold, 
        outFileName, intermediateWrites, 
        seed, threads, maxWeight, 
        reheatIterations, reheatFactor, acceptsUntilUpdate);

    reads = InputProcessing::readReads(readsFile);

    cout << "Start creation of BWT approximate matcher" << endl;
    FMIndex bwt = FMIndex(baseFile, saSF);
    TrainingSearchStrategy* strategy = InputProcessing::createTrainingStrategy(bwt, pStrat, metric, customFile);
    strategy->initializeUniformParameters(ed);

    //
    // LOCAL SEARCH OPTIMIZATION
    //

    // size of shortest pattern in reads
    size_t minPattern = getShortestPatternSize(reads);

    // Get move generator
    vector<double> seedStepSizes = {0.2};
    vector<double> weightStepSizes = {1};

    double maxMoveWeight = 1000.;
    double minWeight = 1.;
    vector <double> startSeeds;
    double minSeedDist;
    // initialize weights and seeds

    if (pStrat == STATIC) {
        minSeedDist = 1./minPattern;
        startSeeds = strategy->getTrainingBegins();
    } else {
        minSeedDist = getMinimumSeedDistance(minPattern, bwt.getWordSize()); // minimum required distance between seeds
        startSeeds = strategy->getTrainingSeeds();
    }
    
    
    vector<int> startWeights = strategy->getTrainingWeights();

    WeightedMoveGenerator *gen = new WeightedMoveGenerator(weightStepSizes, seedStepSizes, startSeeds.size(), startWeights.size(), maxMoveWeight);

    //MoveGenerator *gen = new UniformSeedMoveGenerator(seedStepSizes, startSeeds.size());
    cout << gen->getMoveList().size() << endl;
    // get training sample
    auto trainingSample = getReadSample(reads, threshold, seed, sampleSize, strategy, ed);
    writeSample(trainingSampleFile, trainingSample);

    // initialize acceptor
    
    Solution solution = Solution(startWeights, startSeeds, trainingSampleFile, &bwt, strategy, ed, minSeedDist, chunkSize, threads, maxWeight, trainingSample);

    Acceptor *acceptor = new StepCountingHillClimbingAcceptor(solution, acceptsUntilUpdate, reheatIterations, reheatFactor);
    std::unordered_set <Solution, SolutionHash> tabuList;
    
    cout << "\u001b[31m Running initial solution on both samples...\u001b[0m\r"; cout.flush();
    // run initial solution on training sample
    Solution bestSolution = acceptor->getCurrentSolution();
    size_t bestCost = acceptor->getCurrentCost();

    // run initial solution on test sample
    Solution testSolution = Solution(bestSolution.getWeights(), bestSolution.getSeeds(), readsFile, &bwt, strategy, ed, minSeedDist, chunkSize, threads, maxWeight, reads);
    
    printResults(bestSolution, bestCost, testSolution, testSolution.getCost(), "initial");
    
    auto startTime = chrono::high_resolution_clock::now();
    auto relativeTime = startTime-startTime;
    int64_t seconds = 0, previousPrintTime = 0, newTime = 0;
    int f1 = 5, fn = 5, tmp; // faux fibonacci sequence to determine when to print out values
    
    ofstream outFile (outFileName);

    if (intermediateWrites) {
        writeIntermediateResults(outFile, bestSolution, testSolution, 0);
    }
    map<string, double> improvementsPerMove;
    map<string, int> movesPerMove;
    //printMap(improvementsPerMove, movesPerMove);
    int moves = 0;
    do {
        startTime = chrono::high_resolution_clock::now();
        Move* move = gen->getMove();
        Solution next = move->apply(solution);

        // check for reheat criterium
        // clear the tabu list, and reheat the acceptor
        // could technically also phase shift but idk how complicated a mh is necessary
        //cout << endl << move.getId() << " " << improvementsPerMove[move.getId()] << endl;

        if (tabuList.find(next) == tabuList.end()) {
            // side effects update the current solution
            auto lastCost = acceptor->getCurrentCost();
            if (acceptor->accept(next)) {
                if (lastCost > acceptor->getCurrentCost()) {
                    gen->updateLastWeight(gen->getLastWeight() + 1);
                } else {
                    gen->updateLastWeight(gen->getLastWeight() - 1);
                }
                // keep track of best solution so far
                if (acceptor->getCurrentCost() < bestCost) {
                    bestSolution = next;
                    bestCost = acceptor->getCurrentCost();
                }
                
                // store current solution so we can apply move
                solution = acceptor->getCurrentSolution();
                move->updateStepSize(intToState(acceptor->getCurrentCost() - lastCost));
            } else {
                move->updateStepSize(WORSE);
                // idle
                gen->updateLastWeight(gen->getLastWeight() - 2);
            }
            improvementsPerMove[move->getId()] = improvementsPerMove[move->getId()] * .99 + (double) acceptor->getCurrentCost() - lastCost;
            movesPerMove[move->getId()] += 1;
            moves++;
        
        } else {
            acceptor->idle();
            if (move->hasStepSize()) {
                improvementsPerMove[move->getId()] = improvementsPerMove[move->getId()] * .98;
            } else {
                improvementsPerMove[move->getId()] = improvementsPerMove[move->getId()] * .98 - 10;
            }
            movesPerMove[move->getId()] += 1;
            moves++;
            move->updateStepSize(EQUAL);
            gen->updateLastWeight(gen->getLastWeight() - 1);
        }

        //cout << endl << move.getId() << " (" << gen->getLastWeight() << ") ";

        if (acceptor->isReheat()) {
            //gen->setAllWeights(maxMoveWeight);
            improvementsPerMove.clear();
            movesPerMove.clear();
            tabuList.clear();
            cout << endl << "\u001b[31m REHEATING \u001b[0m" << endl;
            acceptor->reheat();
        } else tabuList.insert(next);

        //cout << move->getId(true) << endl;
        
        // update time
        relativeTime += chrono::high_resolution_clock::now() - startTime; // current time relative to our starting time
        newTime = chrono::duration_cast<chrono::seconds>(relativeTime).count();

        // check how much time has passed
        if (newTime - previousPrintTime >= 1) {
            previousPrintTime = chrono::duration_cast<chrono::seconds>(relativeTime).count();
            // update our output line every second
            seconds = chrono::duration_cast<chrono::seconds>(relativeTime).count();
            printResultLine(acceptor->getCurrentSolution(), acceptor->getCurrentCost(), bestCost, seconds);
        }

        if (moves >= 100) {
            moves = 0;
            //previousUpdateTime = chrono::duration_cast<chrono::seconds>(relativeTime).count();
            map<string, double> tmp;
            auto minElement = minMap(improvementsPerMove, movesPerMove);
            auto s = 0;
            for (const auto &it : improvementsPerMove ) {
                tmp[it.first] = (it.second - minElement)/movesPerMove[it.first] + 1.;
                s += tmp[it.first];
            }

            // scale to maxweight
            for (const auto &it : tmp ) {
                tmp[it.first] *= (maxMoveWeight - minWeight);
                tmp[it.first] /= s;
                tmp[it.first] += minWeight;
            }

            //cout << endl;
            //printMap(tmp, movesPerMove);

            gen->setWeights(tmp);
        }

        if (newTime >= fn && intermediateWrites) {
            // write to output file at time FN (at most every 2 minutes)
            auto testSolution = Solution(bestSolution.getWeights(), bestSolution.getSeeds(), readsFile, &bwt, strategy, ed, minSeedDist, chunkSize, threads, maxWeight, reads);
            writeIntermediateResults(outFile, bestSolution, testSolution, newTime);
            // compute next time we need to write out
            tmp = f1;
            f1 = fn;
            fn += min(tmp, 120);
        }
    
    } while (newTime < duration);


    //printMap(improvementsPerMove, movesPerMove);

    cout << endl;
    cout << "\u001b[31m Running new solution on test sample...\u001b[0m \r";
    
    testSolution = Solution(bestSolution.getWeights(), bestSolution.getSeeds(), readsFile, &bwt, strategy, ed, minSeedDist, chunkSize, threads, maxWeight, reads);

    printResults(bestSolution, bestCost, testSolution, testSolution.getCost(), "Best found");
    if (intermediateWrites) writeIntermediateResults(outFile, bestSolution, testSolution, duration);

    delete strategy;
    delete gen;
    delete acceptor;
    //removeFile(trainingSampleFile);

    cout << "Bye...\n";
}

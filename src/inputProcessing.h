#ifndef INPUT_PROCESSING_H
#define INPUT_PROCESSING_H

#include "searchstrategy.h"
#include "trainingSearchstrategy.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <string.h>
#include <vector>

using namespace std;


/**
 * Class of mostly static functions to be reused by different main files
 */

class InputProcessing {
    private:
        static string getFileExt(const string& s) {

            size_t i = s.rfind('.', s.length());
            if (i != string::npos) {
                return (s.substr(i + 1, s.length() - i));
            }

            return ("");
        }

        static vector<pair<string, string>> getReads(const string& file) {
            vector<pair<string, string>> reads;
            reads.reserve(200000);

            const auto& extension = getFileExt(file);

            bool fasta =
                (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
            bool fastq = (extension == "fq") || (extension == "fastq");

            ifstream ifile(file.c_str());
            if (!ifile) {
                throw runtime_error("Cannot open file " + file);
            }
            if (!fasta && !fastq) {
                // this is a csv file
                if (extension != "csv") {
                    throw runtime_error("extension " + extension +
                                        " is not a valid extension for the readsfile");
                }
                string line;
                // get the first line we do not need this
                getline(ifile, line);

                while (getline(ifile, line)) {
                    istringstream iss{line};

                    vector<string> tokens;
                    string token;

                    while (getline(iss, token, ',')) {
                        tokens.push_back(token);
                    }
                    string position = tokens[1];
                    string read =
                        tokens[2]; // ED + 2 column contains a read with ED compared
                                // to the read at position position with length
                    string p = position;
                    reads.push_back(make_pair(p, read));
                }
            } else if (fasta) {
                // fasta file
                string read = "";
                string p = "";
                string line;
                while (getline(ifile, line)) {
                    if (!line.empty() && line[0] == '>') {

                        if (!read.empty()) {

                            reads.push_back(make_pair(p, read));
                            reads.push_back(
                                make_pair(p, Nucleotide::getRevCompl(read)));
                            read.clear();
                        }

                        p = (line.substr(1));

                    } else {
                        read += line;
                    }
                }
                if (!read.empty()) {

                    reads.push_back(make_pair(p, read));
                    reads.push_back(make_pair(p, Nucleotide::getRevCompl(read)));
                    read.clear();
                }
            } else {
                // fastQ
                string read = "";
                string id = "";
                string line;
                bool readLine = false;
                while (getline(ifile, line)) {
                    if (!line.empty() && line[0] == '@') {
                        if (!read.empty()) {

                            reads.push_back(make_pair(id, read));
                            reads.push_back(
                                make_pair(id, Nucleotide::getRevCompl(read)));
                            read.clear();
                        }
                        id = (line.substr(1));
                        readLine = true;
                    } else if (readLine) {
                        read = line;
                        readLine = false;
                    }
                }
                if (!read.empty()) {

                    reads.push_back(make_pair(id, read));
                    reads.push_back(make_pair(id, Nucleotide::getRevCompl(read)));
                    read.clear();
                }
            }

            return reads;
        }
    public:

        static TrainingSearchStrategy* createTrainingStrategy(FMIndex& fmindex,
                               PartitionStrategy& pStrat,
                               DistanceMetric& metric,
                               const string& pathToCustom) {
            
            return new TrainingSearchStrategy(fmindex, pathToCustom, pStrat, metric);
        };

        static SearchStrategy* createStrategy(const string& searchscheme, FMIndex& fmindex,
                               PartitionStrategy& pStrat,
                               DistanceMetric& metric,
                               const string& pathToCustom) {
            SearchStrategy* strategy;
            if (searchscheme == "kuch1") {
                strategy = new KucherovKplus1(fmindex, pStrat, metric);
            } else if (searchscheme == "kuch2") {
                strategy = new KucherovKplus2(fmindex, pStrat, metric);
            } else if (searchscheme == "kianfar") {
                strategy = new OptimalKianfar(fmindex, pStrat, metric);
            } else if (searchscheme == "manbest") {
                strategy = new ManBestStrategy(fmindex, pStrat, metric);
            } else if (searchscheme == "01*0") {
                strategy = new O1StarSearchStrategy(fmindex, pStrat, metric);
            } else if (searchscheme == "pigeon") {
                strategy = new PigeonHoleSearchStrategy(fmindex, pStrat, metric);
            } else if (searchscheme == "naive") {
                strategy = new NaiveBackTrackingStrategy(fmindex, pStrat, metric);
            } else if (searchscheme == "custom") {
                strategy =
                    new TrainingSearchStrategy(fmindex, pathToCustom, pStrat, metric);
            } else {
                // should not get here
                throw runtime_error(searchscheme +
                                    " is not on option as search scheme");
            }
            return strategy;
        };

        /*static void threadEntry(FastQReader& myReader, OutputWriter& myWriter,
                    SearchStrategy* strategy, int& maxED) {
            // local storage of reads
            vector<FastQRecord> input;
            OutputChunk output;

            size_t chunkID;
            while (myReader.getNextChunk(input, chunkID)) {

                processChunk(input, output, strategy, maxED);

                // mark output as ready
                output.set(chunkID);
                // add output to writer
                myWriter.commitChunk(output);
                output.clear();
            }

            myWriter.sendTermination(chunkID);
        };
        
        static void processChunk(vector<FastQRecord>& input, OutputChunk& output,
                    SearchStrategy* strategy, const int& maxED) {
            output.clear();
            for (const auto& record : input) {
                std::string id = "" + record.getSeqID();

                size_t nodes = 0;
                size_t matrixElements = 0;
                size_t totalOccNonRedundant = 0;

                // match along the forward strand
                vector<TextOccurrence> matches =
                    strategy->matchApprox(record.getRead(), maxED);
                nodes += strategy->getNodes();
                matrixElements += strategy->getMatrixElements();
                totalOccNonRedundant += strategy->getTotal();

                // match along the backward strand
                vector<TextOccurrence> matchesRevCompl = strategy->matchApprox(
                    Nucleotide::getRevCompl(record.getRead()), maxED);
                nodes += strategy->getNodes();
                matrixElements += strategy->getMatrixElements();
                totalOccNonRedundant += strategy->getTotal();

                OutputRecord out =
                    OutputRecord(matches, matchesRevCompl, nodes, matrixElements,
                                totalOccNonRedundant, id);

                output.addRecord(out);
            }
        };*/
        
        /**
         * Parse the base arguments
         * @param argc number of arguments
         * @param argv array of arguments
         * @param basefile file name of the base file (to be filled in by function)
         * @param readsFile file name of the reads file (to be filled in by function)
         * @param searchscheme name of the search scheme (to be filled in by function)
         * @param pStrat name of the used partitioning strategy (to be filled in by function)
         * @param ed the edit/hamming distance (to be filled in by function)
         * @param saSF sparseness factor of suffix array (to be filled in by function)
         * @param metric the used metric (to be filled in by function)
         * @param pathToCustomScheme file name of custom search scheme (if any) (to be filled in by function)
         * @param exclusions a vector of parameters to be ignored during parsing of arguments
         * @param showUsage display a message of how to use the command
         */
        static int parseBaseArguments(int argc, char* argv[], string& baseFile, string& readsFile,
                   string& searchscheme, PartitionStrategy& pStrat, int& ed, length_t& saSF,
                   DistanceMetric& metric, string& pathToCustomScheme, vector<string> exclusions, void (*showUsage)()) {
            
            vector<string> schemes = {"kuch1", "kuch2", "kianfar", "manbest", "pigeon", "01*0",  "naive",   "custom"};
            // process required arguments
            int requiredArguments = 2; // baseFile of files and file containing reads
            if (argc < requiredArguments) {
                cerr << "Insufficient number of arguments" << endl;
                showUsage();
                return EXIT_FAILURE;
            }
            if (argc == 2 && strcmp("help", argv[1]) == 0) {
                showUsage();
                return EXIT_SUCCESS;
            }

            baseFile = argv[argc - 2];
            readsFile = argv[argc - 1];

            // process optional arguments
            string saSparse = "1";
            string maxED = "0";
            searchscheme = "kuch1";

            pStrat = DYNAMIC;
            metric = EDITOPTIMIZED;

            map<string, enum PartitionStrategy> options = {
                {"uniform", UNIFORM},
                {"dynamic",DYNAMIC},
                {"static", STATIC},
            };

            for (int i = 1; i < argc - requiredArguments; i++) {
                const string& arg = argv[i];

                if (arg == "-p" || arg == "--partitioning") {
                    if (i + 1 < argc) {
                        string s = argv[++i];
                        if (options.find(s) != options.end()) {
                            pStrat = options[s];
                        } else {
                            std::stringstream ss;
                            for(auto it = options.begin(); it != options.end(); ++it) {
                                ss << it->first << ", ";
                            }
                            
                            throw runtime_error(
                                s + " is not a partitioning option\nOptions are: "
                                + ss.str());
                        }

                    } else {
                        throw runtime_error(arg + " takes 1 argument as input");
                    }
                } else if (arg == "-s" || arg == "--sa-sparseness") {
                    if (i + 1 < argc) {
                        saSparse = argv[++i];

                    } else {
                        throw runtime_error(arg + " takes 1 argument as input");
                    }
                } else if (arg == "-e" || arg == "--max-ed") {
                    if (i + 1 < argc) {
                        maxED = argv[++i];
                    } else {
                        throw runtime_error(arg + " takes 1 argument as input");
                    }
                } else if (arg == "-ss" || arg == "--search-scheme") {
                    if (i + 1 < argc) {
                        searchscheme = argv[++i];
                        if (find(schemes.begin(), schemes.end(), searchscheme) ==
                            schemes.end()) {
                            throw runtime_error(searchscheme +
                                                " is not an option as search scheme");
                        }

                        if (searchscheme == "custom") {
                            if (i + 1 < argc) {
                                pathToCustomScheme = argv[++i];
                            } else {
                                throw runtime_error(
                                    "custom search scheme takes a folder as argument");
                            }
                        }

                    } else {
                        throw runtime_error(arg + " takes 1 argument as input");
                    }

                } else if (arg == "-m" || arg == "-metric") {
                    if (i + 1 < argc) {
                        string s = argv[++i];
                        if (s == "editopt") {
                            metric = EDITOPTIMIZED;
                        } else if (s == "editnaive") {
                            metric = EDITNAIVE;
                        } else if (s == "hamming") {
                            metric = HAMMING;
                        } else {
                            throw runtime_error(s +
                                                " is not a metric option\nOptions are: "
                                                "editopt, editnaive, hamming");
                        }

                    } else {
                        throw runtime_error(arg + " takes 1 argument as input");
                    }
                }

                else {
                    // skip arguments unique to implementation
                    bool unknownArgument = true;
                    for (auto exclusion: exclusions) {
                        if (arg == exclusion) {
                            unknownArgument = false;
                            i++;
                            break;
                        }
                    }
                    if (unknownArgument) {
                        cerr << "Unknown argument: " << arg << " is not an option" << endl;
                        return false;
                    }
                }
            }

            ed = stoi(maxED);
            if (ed < 0 || ed > 4) {
                cerr << ed << " is not allowed as maxED should be in [0, 4]" << endl;

                return EXIT_FAILURE;
            }
            saSF = stoi(saSparse);
            if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
                cerr << saSF
                    << " is not allowed as sparse factor, should be in 2^[0, 8]"
                    << endl;
            }

            if (ed != 4 && searchscheme == "manbest") {
                throw runtime_error("manbest only supports 4 allowed errors");
            }

            return 2; // parsing succesfull
        };

        static void setParallel(bool p) {
            extern bool allowWrites;
            allowWrites = !p;
        }

        static void printWelcome(string programName) {
            cout << "Welcome to " << programName << "!\n";
        }

        static vector<pair<string, string>> readReads(string readsFile) {
            cout << "Reading in reads from " << readsFile << endl;

            try {
                return getReads(readsFile);
            } catch (const exception& e) {
                string er = e.what();
                er += " Did you provide a valid reads file?";
                throw runtime_error(er);
            }
        }
};

#endif
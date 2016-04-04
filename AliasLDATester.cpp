/*
 * GPUSPDPTester.cpp
 *
 *  Created on: 03/12/2012
 *      Author: polymorpher
 */
#include"ILDA.hpp"
#include"LDA.hpp"
//#include"SparseLDA_Space.hpp"
#include"SparseLDA_Speed.hpp"
#include"util.hpp"
#include<iostream>
#include<cstdio>
#include<json/json.h>
#include"FileReader.hpp"
#include"Logger.hpp"
#include"Dictionary.hpp"
#include"DataGroup.hpp"
#include<boost/regex.hpp>
#include<fstream>
#include<unordered_map>
#include<sys/time.h>
#include"PerplexityAnalyzer.hpp"
#include<map>
#include<glog/logging.h>
#include<boost/format.hpp>
#include<boost/algorithm/string.hpp>
//#include"AliasTable.hpp"
//#include "AliasLDA_Speed.hpp"
#include "AliasLDA_UltraSpeed.hpp"
//#include "AliasLDA_Space.hpp"
//#include"WordPDPLDA.hpp"
//#include"AliasWordPDPLDA.hpp"
//#include"GeneralizedStirlingNumber.hpp"
//#include"StructTester.hpp"
//#include"HDPLDA.hpp"
//#include"AliasHDPLDA.hpp"
namespace AliasLDATester {
bool __compare_reverseOrder(double a, double b) {
	return a > b;
}
class Tester {
public:
	int numIteration;
	AliasLDA::Logger logger;
	int V;
	int numTopicShown;
	int numWordsShown;
	Vector2DInt docs;
	Vector2DInt testdocs;
	DataGroup dataGroup;
	Dictionary dic;
	AliasLDA::ILDA *lda;
	std::unordered_map<int, int> internalDicMapping;
	std::unordered_map<int, int> rInternalDicMapping;
	std::string impl;
	std::string testImpl;
	double totalCPUElapsed; //Microseconds;
	Tester(std::string implementation) :Tester(implementation,"")
			{

	}
	Tester(std::string implementation, std::string testImplementation) :
			numTopicShown(0), numWordsShown(0), V(0), totalCPUElapsed(0), numIteration(
					0), logger("tester.log", true, true), impl(implementation),testImpl(testImplementation){
		lda=createInstance(impl);
	}
	AliasLDA::ILDA* createInstance(std::string implementation){
		if (implementation == "LDA") {
			return new AliasLDA::LDA();
//		} else if (implementation == "AliasLDA_Speed") {
//			return new AliasLDA::AliasLDA_Speed();
		} else if (implementation == "AliasLDA_UltraSpeed") {
			return new AliasLDA::AliasLDA_UltraSpeed();
//		} else if (implementation == "AliasLDA_Space") {
//			return new AliasLDA::AliasLDA_Space();
		} else if (implementation == "SparseLDA_Speed") {
			return new AliasLDA::SparseLDA_Speed();
//		} else if (implementation == "SparseLDA_Space") {
//			return new AliasLDA::SparseLDA_Space();
//		} else if (implementation == "WordPDPLDA") {
//			return new AliasLDA::WordPDPLDA();
//		} else if (implementation == "AliasWordPDPLDA") {
//			return new AliasLDA::AliasWordPDPLDA();
//		}else if (implementation == "StructTester") {
//			return new AliasLDA::StructTester();
//		}else if (implementation == "HDPLDA") {
//			return new AliasLDA::HDPLDA();
//		}else if (implementation == "AliasHDPLDA") {
//			return new AliasLDA::AliasHDPLDA();
		}
		return NULL;
	}
	bool parseTestFile() {
		const std::string &fileName = dataGroup.testFile;
		LOG(INFO) << (std::string("Parsing test file: ") + fileName);
		int linenum = 0;
		std::fstream file(fileName, std::fstream::in);
		if (file.good()) {
			std::string line = "0";
			while ((!line.empty()) && (!file.eof())) {
				std::getline(file, line);
				linenum++;
				std::vector<std::string> words;
				boost::regex re("\\s+");
				boost::sregex_token_iterator it(line.begin(), line.end(), re,
						-1);
				boost::sregex_token_iterator end;
				while (it != end)
					words.push_back(*it++);

				std::vector<int> thisdoc;
				for (size_t i = 0; i < words.size(); i++) {
					size_t id;
					__convert(words[i], id);
					id--;
					if (id < 0 || id >= dic.size()) {
						LOG(INFO) << format(
								"Test file contains invalid entry %1% at line %2%.") %
								words[i].c_str() % linenum;
						file.close();
						return false;
					}
					thisdoc.push_back(id);
				}
				if (!thisdoc.empty())
					testdocs.push_back(std::vector<int>(thisdoc));
			}
			file.close();
		} else {
			LOG(INFO) << "Cannot open test file: " << fileName.c_str();
			return false;
		}
		LOG(INFO) << ("Test file parsed successfully");
		return true;
	}

	bool parseTrainingFile() {
		const std::string &fileName = dataGroup.trainingFile;
		LOG(INFO) << std::string("Parsing training file: ") << fileName;
		int linenum = 0;
		std::fstream file(fileName, std::fstream::in);
		if (file.good()) {
			std::string line = "0";
			while ((!line.empty()) && (!file.eof())) {
				std::getline(file, line);
				linenum++;
				std::vector<std::string> words;
				//if(linenum>100)break;
//					std::regex r("\\s+");
//					std::regex_token_iterator<std::string::iterator> it (line.begin(),line.end(),r,-1);
//					std::regex_token_iterator<std::string::iterator> end;
				boost::regex re("\\s+");
				boost::sregex_token_iterator it(line.begin(), line.end(), re,
						-1);
				boost::sregex_token_iterator end;
				while (it != end)
					words.push_back(*it++);

				std::vector<int> thisdoc;
				for (size_t i = 0; i < words.size(); i++) {
					size_t id;
					__convert(words[i], id);
					id--;
					if (id < 0 || id >= dic.size()) {
						LOG(INFO) << format(
								"Training file contains invalid entry %1% at line %2%.") %
								words[i].c_str() % linenum;
						file.close();
						return false;
					}
					thisdoc.push_back(id);
				}
				if (!thisdoc.empty())
					docs.push_back(std::vector<int>(thisdoc));
			}
			file.close();
		} else {
			LOG(INFO) << "Cannot open trainning file: " << fileName.c_str();
			return false;
		}
		LOG(INFO) << "Training file parsed successfully";
		return true;
	}

	bool parseBagTrainingFile() {
		const std::string &fileName = dataGroup.trainingFile;
		LOG(INFO) << std::string("Parsing bag training file: ") + fileName;
		int linenum = 0;
		std::fstream file(fileName, std::fstream::in);
		if (file.good()) {
			std::string line = "0";
			std::getline(file, line);
			std::getline(file, line);
			std::getline(file, line);
			linenum += 3;
			std::vector<int> words;
			int curdoc = 0;
			while ((!line.empty()) && (!file.eof())) {
				line = "";
				std::getline(file, line);
				if (line.empty())break;
				linenum++;
				std::stringstream ss;
				ss << line;
				int docid = -1, tokenid = -1, numtoken = -1;
				ss >> docid >> tokenid >> numtoken;
				ss.clear();
				if (docid == -1 || tokenid == -1 || numtoken == -1) {
					LOG(INFO) << format(
							"Training file contains invalid entry %1% at line %2%.") %
							line.c_str() % linenum;
					file.close();
					return false;
				}
				if (!dic.check(tokenid)) {
					LOG(INFO) << format(
							"Training file contains invalid token %1% at line %2%.") %
							line.c_str() % linenum;
					file.close();
					return false;
				}
				if (docid != curdoc) {
					if (!words.empty()) {
						docs.push_back(std::vector<int>(words));
						words.clear();
					}
					curdoc = docid;
				}
				for (int i = 0; i < numtoken; i++) {
					words.push_back(tokenid);
				}
			}
			file.close();
		} else {
			LOG(INFO) << "Cannot open training file: " << fileName.c_str();
			return false;
		}
		LOG(INFO) << "Training file parsed successfully";
		return true;
	}

	bool parseBagTestFile() {
		const std::string &fileName = dataGroup.testFile;
		LOG(INFO) << std::string("Parsing bag test file: ") << fileName;
		int linenum = 0;
		std::fstream file(fileName, std::fstream::in);
		if (file.good()) {
			std::string line = "0";
			std::getline(file, line);
			std::getline(file, line);
			std::getline(file, line);
			linenum += 3;
			std::vector<int> words;
			int curdoc = 0;
			while ((!line.empty()) && (!file.eof())) {
				line = "";
				std::getline(file, line);
				if (line.empty())
					break;
				linenum++;
				std::stringstream ss;
				ss << line;
				int docid = -1, tokenid = -1, numtoken = -1;
				ss >> docid >> tokenid >> numtoken;
				ss.clear();
				if (docid == -1 || tokenid == -1 || numtoken == -1) {
					LOG(INFO) << format(
							"Test file contains invalid entry %1% at line %2%.") %
							line.c_str() % linenum;
					file.close();
					return false;
				}
				if (!dic.check(tokenid)) {
					LOG(INFO) << format(
							"Test file contains invalid token %1% at line %2%.") %
							line.c_str() % linenum;
					file.close();
					return false;
				}
				if (docid != curdoc) {
					if (!words.empty()) {
						testdocs.push_back(std::vector<int>(words));
						words.clear();
					}
					curdoc = docid;
				}
				for (int i = 0; i < numtoken; i++) {
					words.push_back(tokenid);
				}
			}
			file.close();
		} else {
			LOG(INFO) << "Cannot open test file: " << fileName.c_str();
			return false;
		}
		LOG(INFO) << "Test file parsed successfully";
		return true;
	}

	bool parseDataConfig(const char *path) {
		const char *jsonStr = AliasLDA::FileReader::readFile(path);
		if (jsonStr == NULL) {
			return false;
		}
		LOG(INFO) << jsonStr;
		Json::Value root;
		Json::Reader reader;
		bool parsingSuccessful = reader.parse(jsonStr, root);
		if (!parsingSuccessful) {
			LOG(INFO) << "Data configuration file syntax error.";
			return false;
		}
		int i = 0;
		LOG(INFO) << "Configuration read.";
		std::vector<std::string> m = root.getMemberNames();
		dataGroup.name = root[m[i]]["name"].asString();
		dataGroup.vocabularyFile = root[m[i]]["vocabularyFile"].asString();
		dataGroup.trainingFile = root[m[i]]["trainingFile"].asString();
		dataGroup.testFile = root[m[i]]["testFile"].asString();
		dataGroup.format = root[m[i]]["format"].asString();
		LOG(INFO) << "Group " << i << " settings read.";
		if (dataGroup.vocabularyFile.empty() || dataGroup.trainingFile.empty()
				|| dataGroup.testFile.empty()) {
			LOG(INFO) << "Some fields in data configuration file are invalid.";
			return false;
		}
		if (dataGroup.format == "bag") {
			if (!dic.parseBagVocabularyFile(dataGroup.vocabularyFile)) {
				LOG(INFO) << "Failed to parse bag vocabulary file";
				return false;
			}
			if (!parseBagTrainingFile()) {
				LOG(INFO) << "Failed to parse bag training file";
				return false;
			}
			if (!parseBagTestFile()) {
				LOG(INFO) << "Failed to parse bag testing file";
				return false;
			}
		} else {
			if (!dic.parseVocabularyFile(dataGroup.vocabularyFile)) {
				LOG(INFO) << "Failed to parse vocabulary file";
				return false;
			}
			if (!parseTrainingFile()) {
				LOG(INFO) << "Failed to parse training file";
				return false;
			}
			if (!parseTestFile()) {
				LOG(INFO) << "Failed to parse testing file";
				return false;
			}
		}
		LOG(INFO) << format("Group %1% vocabulary size %2%.") % i % dic.size();

		std::unordered_map<int, bool> actualVocab;
		for (int d = 0; d < docs.size(); d++) {
			VecInt &doc = docs[d];
			for (int l = 0; l < doc.size(); l++) {
				actualVocab[doc[l]] = true;
			}
		}
		for (int d = 0; d < testdocs.size(); d++) {
			VecInt &doc = testdocs[d];
			for (int l = 0; l < doc.size(); l++) {
				actualVocab[doc[l]] = true;
			}
		}
		dic.compress(actualVocab);
		LOG(INFO) << "Effective vocabulary size " << dic.size();
		V = (int) dic.size();
		int internalCount = 0;
		for (auto it = dic.dictionary.begin(); it != dic.dictionary.end();
			 ++it) {
			internalDicMapping[internalCount] = it->second;
			rInternalDicMapping[it->second] = internalCount;
			internalCount++;
		}
		for (int d = 0; d < docs.size(); d++) {
			VecInt &doc = docs[d];
			for (int l = 0; l < doc.size(); l++) {
				doc[l] = rInternalDicMapping[doc[l]];
			}
		}
		for (int d = 0; d < testdocs.size(); d++) {
			VecInt &doc = testdocs[d];
			for (int l = 0; l < doc.size(); l++) {
				doc[l] = rInternalDicMapping[doc[l]];
			}
		}
		return true;
	}

	void computePerplexity() {
		LOG(INFO) << "Computing test perplexity...";
		AliasLDA::ILDA *testlda;
		if (testImpl == "") {
			testlda = createInstance(impl);
			testlda->copyState(lda);
		} else {
			testlda = createInstance(testImpl);
			testlda->copyState((void *) lda, impl);
		}

		testlda->setDocuments(testdocs);
		testlda->initialiseAssignmentOnly();
		LOG(INFO) << "Initialised test instance. Start sampling...";
		for (int i = 0; i < 4; i++) {
			testlda->gibbsStep();
			LOG(INFO) << "Sampling iteration " << i << " complete";
		}
		testlda->computeTheta();
		const Vector2DDouble &theta = testlda->getTheta();
		lda->computePhi();
		const Vector2DDouble &phi = lda->getPhi();

		PerplexityAnalyzer *pa = new PerplexityAnalyzer(testdocs, phi, theta,
				testlda->getNumTopics(), V);
		double perp=pa->getPerplexity();

		delete pa;
		printLDAStatus(testlda);
		delete testlda;
		lda->clearPhi();
		LOG(INFO) << "Test Perplexity=" << perp;
		LOG(INFO) << "Test instance deleted.";

	}

	void printLDAStatus(AliasLDA::ILDA *_lda) {
		_lda->computeThetaSum();
		_lda->computePhi();
		_lda->computeTheta();

		const VecDouble &thetasum = _lda->getThetasum();
		const Vector2DDouble &phi = _lda->getPhi();
		const Vector2DDouble &theta = _lda->getTheta();
		numTopicShown = 10; //lda->getNumTopics();
		numWordsShown = 20;

		bool (*rFunc)(double, double) = __compare_reverseOrder;
		std::multimap<double, int, bool (*)(double, double)> topTopics(rFunc);
		//for(int i=0;i<I;i++){topTopics[i]=new std::multimap<float,int,bool(*)(float,float)>(rFunc);}
		std::vector<std::multimap<double, int, bool (*)(double, double)> *> topWords;
		//(I,std::vector(numTopicKeep,std::multimap<float,int,bool(*)(float,float)>(rFunc)));
		VecInt topicLookup;
		topicLookup.resize(numTopicShown);
		std::unordered_map<int, int> topicReverseLookup;
		topWords.resize(numTopicShown);
		for (size_t k = 0; k < numTopicShown; k++) {
			topWords[k] = new std::multimap<double, int,
					bool (*)(double, double)>(rFunc);
		}
		for (int k = 0; k < _lda->getNumTopics(); k++) {
			if (topTopics.size() < numTopicShown) {
				topTopics.insert(std::pair<double, int>(thetasum[k], (int) k));
			} else {
				double key = (--topTopics.end())->first;
				if (key < thetasum[k]) {
					topTopics.insert(std::pair<double, int>(thetasum[k], k));
					topTopics.erase(--(topTopics.end()));
				}
			}
		}
		int ind = 0;
		for (std::multimap<double, int, bool (*)(double, double)>::iterator it =
				topTopics.begin(); it != topTopics.end(); it++) {
			int k = it->second;
			LOG(INFO) << format("Topic %1%, Rank:%2% Probability:%3%") % k % (ind + 1) % it->first;
			topicLookup[ind] = k;
			topicReverseLookup[k] = ind;
			for (int w = 0; w < V; w++) {
				double phival = phi[w][k];
				if (topWords[ind]->size() < numWordsShown) {
					topWords[ind]->insert(std::pair<double, int>(phival, w));
				} else {
					double key = (--topWords[ind]->end())->first;
					if (key <= phival) {
						topWords[ind]->insert(
								std::pair<double, int>(phival, w));
						topWords[ind]->erase(--(topWords[ind]->end()));
					}
				}
			}
			std::stringstream ss;
			for (std::multimap<double, int, bool (*)(double, double)>::iterator it =
					topWords[ind]->begin(); it != topWords[ind]->end(); ++it) {
				ss << format("[%1%:%2%,%3%] ") % it->second %
						dic.rlookup(internalDicMapping[it->second]).c_str() %
						it->first;
			}
			LOG(INFO) << ss.str();
			LOG(INFO) << "--------------------------";
			ind++;
		}

		for (size_t k = 0; k < numTopicShown; k++) {
			delete topWords[k];
		}
//		if(_lda==lda){
//			printf("\nComputing Training-Perplexity (For convergence analysis)...\n");
//			PerplexityAnalyzer* pa = new PerplexityAnalyzer(docs, phi, theta,
//					lda->getNumTopics(), V);
//			printf("Training-Perplexity=%f\n", pa->getPerplexity());
//			delete pa;
//		}
		_lda->clearPhi();


	}

	void CPUStep() {
		struct timeval tv;
		struct timeval start_tv;
		double elapsed = 0.0;
		LOG(INFO) << format(
				"Gibbs CPU Sampling started, current total CPU time= %1% s, average=%2% s, numiteration=%3%)") %
				totalCPUElapsed % (totalCPUElapsed / (double) numIteration) % numIteration;
		gettimeofday(&start_tv, NULL);
		lda->gibbsStep();
		gettimeofday(&tv, NULL);
		numIteration++;
		elapsed = (tv.tv_sec - start_tv.tv_sec)
				+ (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
		totalCPUElapsed += elapsed;
		LOG(INFO) << format(
				"Gibbs CPU Sampling step completed, this iteration = %1% s, current total CPU time=%2% s, average=%3% s, numIteration=%4%") %
				elapsed % totalCPUElapsed %
				(totalCPUElapsed / (double) numIteration) % numIteration;
	}
};
}

int main(int argc, char *argv[]) {
	FLAGS_logtostderr=1;
	google::InitGoogleLogging(argv[0]);
	dup2(1, 2);
	LOG(INFO) << "instance started.";
	if (argc < 7) {
		LOG(INFO) << "Incorrect number of arguments";
		return 0;
	}
	std::string dataPath = std::string(argv[1]);
	std::string method = std::string(argv[2]);
	boost::algorithm::to_lower(method);
	double alpha = convertString<double>(std::string(argv[3]));
	double beta = convertString<double>(std::string(argv[4]));
	int numTopics = convertString<int>(std::string(argv[5]));
	int numIterations = convertString<int>(std::string(argv[6]));
//	numIterations = std::min(100, numIterations);
	AliasLDATester::Tester *tester;
	if (method == "lda") {
		tester = new AliasLDATester::Tester("LDA", "LDA");
	} else if (method == "sparselda") {
		tester = new AliasLDATester::Tester("SparseLDA_Speed", "SparseLDA_Speed");
	} else if (method == "aliaslda") {
		tester = new AliasLDATester::Tester("AliasLDA_UltraSpeed", "AliasLDA_UltraSpeed");
	} else {
		LOG(INFO) << "Unknown implementation "<<method;
		return 0;
	}
	LOG(INFO) << "initialised.";
	tester->parseDataConfig(dataPath.c_str());
	LOG(INFO) << "LDA starting...\n";
	tester->lda->setDocuments(tester->docs);
	tester->lda->setAlpha(alpha);
	tester->lda->setBeta(beta);
	tester->lda->setNumTopics(numTopics);
	tester->lda->setSizeVocabulary(tester->V);
	tester->lda->initialise();
	LOG(INFO) << "Starting instance. numDocs=" << tester->lda->getDocuments().size() << " K=" << tester->lda->getNumTopics() << " V=" << tester->V;
	LOG(INFO) << "Test docset numDocs=" << tester->testdocs.size();
//	for (int i = 0; i < 20; ++i) {
//		printf("%s ",tester.dic.rlookup(tester.internalDicMapping[tester.lda->getDocuments()[1][i]]).c_str());
//	}
//	return 0;
//	tester.unitTest_SparseLDA_oneStepSortInc();
//	tester.unitTest_SparseLDA_oneStepSortDec();
//	return 0;
	LOG(INFO) << "";
	for (int i = 0; i < numIterations; ++i) {
//		LOG(INFO)<<"iter="<<i;
		tester->CPUStep();
		if (i % 5 == 0) {
			tester->printLDAStatus(tester->lda);
			tester->computePerplexity();
		}
		//sleep(5);
	}
	exit(0);
	delete tester->lda;
	LOG(INFO) << ("instance deleted!");
	return 0;
}


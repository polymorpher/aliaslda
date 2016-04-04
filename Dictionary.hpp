/*
 * PDPDictionary.hpp
 *
 *  Created on: 03/12/2012
 *      Author: polymorpher
 */
//retval:
//0 success
//1 syntax error in vocabulary file
#ifndef PDPDICTIONARY_HPP_
#define PDPDICTIONARY_HPP_
#include<unordered_map>
#include<string>
#include<vector>
#include<fstream>
#include<regex>
#include<sstream>
#include<boost/regex.hpp>
template <class T, class R> void __convert(T &source, R& ret){
	std::stringstream s;
	s<<source;
	s>>ret;
}
namespace AliasLDATester{
	class Dictionary{
	public:
		std::unordered_map<std::string,int> dictionary;
		std::unordered_map<int,std::string> reverseDictionary;
		int index;
		Dictionary():index(0){
		}
		int addWord(std::string word){
			if(dictionary.count(word)==0){
				reverseDictionary[index]=word;
				dictionary[word]=index;
				index++;
				return index-1;
			}else{
				return dictionary[word];
			}
		}
		inline bool check(std::string word){
			return dictionary.find(word)!=dictionary.end();
		}
//		inline bool check(int id){
//			return reverseDictionary.find(id)!=reverseDictionary.end();
//		}
		inline bool check(int id){
			return id>0&&id<=dictionary.size();
		}
		inline int lookup(std::string word){
			return dictionary[word];
		}
		inline std::string& rlookup(int id){
			return reverseDictionary[id];
		}
		inline size_t size(){
			return dictionary.size();
		}
		bool parseVocabularyFile(std::string file){
			//logger.append("Parsing vocabulary file for group %d.\n",groupID);
			int maxID=0;
			std::fstream f(file,std::fstream::in);
			std::string line="0";
			while(!(f.eof()||line.empty())){
				std::getline(f,line);
				boost::regex r("(\\d+)(\\s+)(.+)");
				boost::cmatch m;
				if(!boost::regex_match(line.c_str(),m,r)){
					continue;
				}
				std::string word;
				int id;
				__convert(m[3],word);
				__convert(m[1],id);
				dictionary[word]=id-1;
				reverseDictionary[id-1]=word;
			}
			f.close();
			return true;
		}
		bool parseBagVocabularyFile(std::string file){
			//logger.append("Parsing vocabulary file for group %d.\n",groupID);
			int counter=1;
			std::fstream f(file,std::fstream::in);
			std::string line="0";
			while(!(f.eof()||line.empty())){
				std::getline(f,line);
				dictionary[line]=counter;
				reverseDictionary[counter]=line;
				counter++;
			}
			f.close();
			return true;
		}
		void compress(std::unordered_map<int,bool>& actualVocab){
			std::vector<int> removeListInt;
			std::vector<std::string> removeListString;
			for(auto it=dictionary.begin();it!=dictionary.end();++it){
				if(actualVocab.count(it->second)<=0){
					removeListInt.push_back(it->second);
					removeListString.push_back(it->first);
				}
			}
			printf("removing %d items\n",removeListInt.size());
			for(int i=0;i<removeListInt.size();i++){
				reverseDictionary.erase(removeListInt[i]);
			}
			for(int i=0;i<removeListString.size();i++){
				dictionary.erase(removeListString[i]);
			}
		}
	};
}


#endif /* PDPDICTIONARY_HPP_ */

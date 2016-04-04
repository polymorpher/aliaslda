/*
 * DataGroup.hpp
 *
 *  Created on: 03/12/2012
 *      Author: polymorpher
 */
#pragma once
#include<string>
namespace AliasLDATester {
class DataGroup {
public:
	std::string name;
	std::string vocabularyFile;
	std::string trainingFile;
	std::string testFile;
	std::string format;
};
}

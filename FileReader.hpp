/*
 * FileReader.hpp
 *
 *  Created on: 03/12/2012
 *      Author: polymorpher
 */

#pragma once
#include<cstdio>
#include<cstdlib>
namespace AliasLDA{
	class FileReader{
	public:
		static char* readFile(const char *filename);
	};
}


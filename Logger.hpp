/*
 * Logger.hpp
 *
 *  Created on: 03/12/2012
 *      Author: polymorpher
 */

#pragma once
#include<fstream>
#include<cstdio>
#include<iostream>
#include<string>
#include<cstdarg>
namespace AliasLDA{
	class Logger {
		std::fstream f;
	public:
		bool fileLogEnabled;
		bool printEnabled;
		void append(std::string s, ...){
			char outs[1024];
		    va_list argptr;
		    va_start(argptr, s);
		    vsnprintf(outs,1023,s.c_str(), argptr);
		    va_end(argptr);
		    if(printEnabled)printf("%s",outs);
			if(fileLogEnabled){f<<outs;f.flush();}
		}
		Logger(std::string fileName,bool fileLog, bool printLog):f(fileName),
				fileLogEnabled(fileLog),printEnabled(printLog){

		}
		Logger(){
			Logger("out.log",true,true);
		}
	};
}


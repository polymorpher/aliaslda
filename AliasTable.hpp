/*
 * AliasTable.hpp
 *
 *  Created on: 05/01/2014
 *      Author: polymorpher
 */

#pragma once
#include"util.hpp"
#include<climits>
namespace AliasLDA{
class AliasTable{
public:
	struct LabeledElement{
		int label;
		double prob;
		LabeledElement(int _label,double _prob):label(_label),prob(_prob){}
	};
	struct AliasElement{
		int label;
		int alias;
		double prob;
		AliasElement(int _label,int _alias, double _prob):label(_label),alias(_alias),prob(_prob){}
	};
protected:
	VecDouble proportions;
	std::vector<AliasElement> table;
	VecInt labels;	//main labels
	bool useLabels;
	std::mt19937_64 rgen;
	std::uniform_real_distribution<double> u01;
	int tableSize;
	int maxTableSize;
	MapIntInt labelPosLookup;
	bool locked;
public:
	const inline VecDouble& getProportions() const {
		return proportions;
	}
	AliasTable(size_t n):useLabels(false),u01(0,1),maxTableSize(n),locked(false){
		init();
	}
	AliasTable(size_t n, VecInt& _labels):labels(_labels),useLabels(true),u01(0,1),maxTableSize(n),locked(false){
		assert(_labels.size()==n);
		init();
		for(int i=0;i<labels.size();i++){
			labelPosLookup[labels[i]]=i;
		}
	}
	void init(){
		AliasElement ae(-1,-1,0);
		table.clear();
		table.resize(maxTableSize, ae);
		proportions.resize(maxTableSize);
		rgen.seed(time(NULL));
		labelPosLookup.set_empty_key(INT_MAX);
	}
	inline double rand01(){
		return u01(rgen);
	}
	inline int randInt(int limit){
		return (int)(rand01()*limit);
	}
	inline int getLabel(int pos){
		return useLabels?labels[pos]:pos;
	}
	void build(VecDouble& prop){
		int n=prop.size();
		tableSize=n;
		double sum=0;
		for(int i=0;i<prop.size();i++){
			sum+=prop[i];
		}
		std::vector<LabeledElement> small,large;
		proportions = prop;
		small.reserve(n);
		large.reserve(n);
		double mfactor=n/sum;
		for(int i=0;i<prop.size();i++){
			double val=prop[i]*mfactor;
			if(val>=1){
				large.push_back(LabeledElement(i,val));
			}else{
				small.push_back(LabeledElement(i,val));
			}
		}
		int largeCursor=0;
		int smallCursor=0;
		while(smallCursor<small.size() && largeCursor<large.size()){
			LabeledElement& sle=small[smallCursor];
			LabeledElement& lle=large[largeCursor];
			AliasElement& ae=table[sle.label];
			ae.label=getLabel(sle.label);
			ae.alias=getLabel(lle.label);
			ae.prob=sle.prob;
			lle.prob=lle.prob+sle.prob-1;
			if(lle.prob<1){
				small.push_back(lle);
				smallCursor++;
				largeCursor++;
			}else{
				smallCursor++;
			}
		}
		while(largeCursor<large.size()){
			LabeledElement& lle=large[largeCursor];
			AliasElement& ae=table[lle.label];
			ae.label=getLabel(lle.label);
			ae.alias=-1;
			ae.prob=1;
			largeCursor++;
		}
		while(smallCursor<small.size()){
//			printf("Warning: numerical instability. %d/%d\n", smallCursor,small.size());
			LabeledElement& sle=small[smallCursor];
			AliasElement& ae=table[sle.label];
			ae.label=getLabel(sle.label);
			ae.alias=-1;
			ae.prob=1;
			smallCursor++;
		}
	}
	int sample(){
		int roll=randInt(tableSize);
		AliasElement& ae=table[roll];
		double roll01=rand01();
		if(roll01<ae.prob){
			return ae.label;
		}else{
			return ae.alias;
		}
	}
	inline void lock(){
		while(locked){}
		locked=true;
	}
	inline bool try_lock(){
		if(locked){
			return false;
		}else{
			lock();
			return true;
		}
	}
	inline void unlock(){
		locked=false;
	}
	inline void force_lock(){
		locked=true;
	}
	inline bool isLocked(){
		return locked;
	}
	~AliasTable(){
		while(locked){printf("warning: delaying destruction of alias table\n");}
	}

};
typedef std::unique_ptr<AliasTable> ATableUPtr;
}

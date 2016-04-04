/*
 * ILDA.hpp
 *
 *  Created on: 24/12/2013
 *      Author: polymorpher
 */

#pragma once
#include"util.hpp"
namespace AliasLDA {

class ILDA{
public:
	virtual ~ILDA(){}
	//optional
	virtual void setGamma(double gamma){}
	virtual double getGamma(){return 0;}
	virtual void setDiscount(double a){}
	virtual double getDiscount(){return 0;}
	virtual void setConcentration(double b){}
	virtual double getConcentration(){return 0;}
	virtual void setRootConcentration(double b){}
	virtual double getRootConcentration(){return 0;}
	virtual void setBeta(double beta){}
	virtual double getBeta() const{return 0;}
	virtual void copyState(void* _lda, std::string impl){}
	//required
	virtual void setAlpha(double alpha){}
	virtual void setAlphaVec(const VecDouble& alphaVec){}
	virtual const VecDouble& getAlphaVec() const {return VecDouble();}
	virtual void initialise()=0;
	virtual void computeTheta()=0;
	virtual void computeThetaSum()=0;
	virtual void computePhi()=0;
	virtual void clearPhi()=0;
	virtual void clearTheta()=0;
	virtual void gibbsStep()=0;
	virtual const Vector2DDouble& getPhi() const=0;
	virtual const Vector2DDouble& getTheta() const=0;
	virtual const VecDouble& getThetasum() const=0;
	virtual void setDocuments(Vector2DInt& docs)=0;
	virtual const Vector2DInt& getDocuments() const=0;
	virtual void setNumTopics(int K)=0;
	virtual int getNumTopics() const=0;
	virtual void setSizeVocabulary(int V)=0;
	virtual void copyState(ILDA* _lda) =0;

	virtual void initialiseAssignmentOnly() =0;
};
} /* namespace AliasLDA */




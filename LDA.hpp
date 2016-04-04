/*
 * LDA.h
 *
 *  Created on: 17/12/2013
 *      Author: polymorpher
 */

#pragma once

#include"ILDA.hpp"
#include"util.hpp"
#include"SparseLDA_Speed.hpp"
//#include"AliasLDA_UltraSpeed.hpp"
namespace AliasLDA {

class LDA: public ILDA {
	VecInt nk;
	Vector2DInt docs;
	Vector2DInt z;	//Topic assignment
	//VecInt nkw;//Awesome SparseLDA special data structure
	Vector2DInt nwk;
	Vector2DInt ndk;
	//Uniform Prior
	double beta;
	double alpha;
	//Non-uniform Prior
	//VecDouble betaVec;
	VecDouble alphaVec;
	//Common between unifrom and non-unifrom prior
	double betasum;
	double alphasum;
	int K;	//numTopics
	int V;	//sizeVocabulary
	std::mt19937_64 rgen;
	std::uniform_real_distribution<double> u01;

	//results
//	Vector2DDouble phi;
//	Vector2DDouble theta;
	Vector2DDouble phi;
	Vector2DDouble theta;
	VecDouble thetasum;

	//temp buffer
	VecDouble topicProbabilities;

public:
	LDA():K(0),V(0),beta(0.1)/*,alpha(0.1)*/,betasum(0),alphasum(0),u01(0,1){
		rgen.seed(time(NULL));
	}

	inline double rand01(){
		return u01(rgen);
	}
	inline int randInt(int limit){
		return (int)(rand01()*limit);
	}

	inline void inc_ndk(int d,int k){
		ndk[d][k]++;
		//nd[d]++;
	}
	inline void dec_ndk(int d,int k){
		ndk[d][k]--;
		//nd[d]--;
	}
	inline void inc_nwk(int w,int k){
		nwk[w][k]++;
		nk[k]++;
	}
	inline void dec_nwk(int w,int k){
		nwk[w][k]--;
		nk[k]--;
	}
	void initialise(){
		initialiseBeta();
		initialiseAlpha();
		vecResize2D(nwk,V,K);
		initialiseAssignmentOnly();
	}
	virtual void initialiseAssignmentOnly(){
		vecResize2D(ndk,docs.size(),K);
		thetasum.resize(K);
		initialiseTopicsAndCounts();
		topicProbabilities.resize(K);
	}
	void initialiseTopicsAndCounts(){

		z.resize(docs.size());
		nk.resize(K);
		for(size_t d=0;d<docs.size();d++){
			VecInt& doc=docs[d];
			z[d].resize(doc.size());
			VecInt& zd=z[d];
			for(size_t l=0;l<doc.size();l++){
				int w=doc[l];
				int k=randInt(K);
				inc_ndk(d,k);
				inc_nwk(w,k);
				zd[l]=k;
			}
		}
	}
	void initialiseBeta(){
		betasum=beta*V;
	}
	void initialiseAlpha(){
		alphasum=0;
		if(alpha!=0){
			alphaVec.resize(K);
			for(size_t i=0;i<alphaVec.size();i++){
				alphaVec[i]=alpha;
			}
		}
		for(auto it=alphaVec.begin();it!=alphaVec.end();it++){
			alphasum+=*it;
		}
	}
	virtual void sampleWord(int d,int l){
		int w=docs[d][l];
		int k=z[d][l];
		dec_ndk(d,k);
		dec_nwk(w,k);
		double sump=0;
		for(int i=0;i<K;i++){
			double prob=(ndk[d][i]+alpha)*(beta+nwk[w][i])/(betasum+nk[i]);
			topicProbabilities[i]=prob;
			//printf("k=%d prob=%lf\n",i,topicProbabilities[i]);
			sump+=prob;
		}
		double r=rand01()*sump;
		int newTopic=-1;
		for(int i=0;i<K;i++){
			double curprob=topicProbabilities[i];
			if(curprob>=r){
				newTopic=i;
				break;
			}else{
				r-=curprob;
			}
		}
		if(abs(r)<=0.0001 && newTopic==-1){
			newTopic=K-1;
		}
		assert(newTopic!=-1);
		inc_ndk(d,newTopic);
		inc_nwk(w,newTopic);
		z[d][l]=newTopic;
	}
	virtual void gibbsStep(){
		for(size_t d=0;d<docs.size();d++){
			VecInt& doc=docs[d];
			for(size_t l=0;l<doc.size();l++){
				sampleWord(d,l);
			}
		}
	}
	virtual void computePhi(){
		phi.resize(V);
		for(int w=0;w<V;w++){
			phi[w].resize(K);
			for(int k=0;k<K;k++){
				phi[w][k]=(nwk[w][k]+beta)/(nk[k]+betasum);
			}
		}
	}
	virtual void computeTheta(){
		theta.resize(docs.size());
		for(size_t d=0;d<docs.size();d++){
			theta[d].resize(K);
			for(int k=0;k<K;k++){
				theta[d][k]=(ndk[d][k]+alpha)/(docs[d].size()+alphasum);
			}
		}
	}
	virtual void computeThetaSum(){
		for(int k=0;k<K;k++){
			thetasum[k]=0;
		}
		for(size_t d=0;d<docs.size();d++){
			for(int k=0;k<K;k++){
				double thetaval=(ndk[d][k]+alpha)/(docs[d].size()+alphasum);
				thetasum[k]+=thetaval/docs.size();
			}
		}
	}
	//getters and setters
	virtual void setDocuments(Vector2DInt& docs){
		this->docs.resize(docs.size());
		for(size_t i=0;i<docs.size();i++){
			auto& doc=this->docs[i];
			auto& thisdoc=docs[i];
			doc.resize(docs[i].size());
			for(size_t j=0;j<thisdoc.size();j++){
				doc[j]=thisdoc[j];
				//printf("%d ",thisdoc[j]);
			}
		}
		//this->docs=docs;
	}
//	double getAlpha() const {
//		return alpha;
//	}
//
//	void setAlpha(double alpha) {
//		this->alpha = alpha;
//	}

	virtual const VecDouble& getAlphaVec() const {
		return alphaVec;
	}

	virtual void setAlphaVec(const VecDouble& alphaVec) {
		this->alphaVec = alphaVec;
	}
	virtual void setAlpha(double alpha){
		this->alpha = alpha;
	}
	virtual double getBeta() const {
		return beta;
	}

	virtual void setBeta(double beta) {
		this->beta = beta;
	}
	virtual void setNumTopics(int K){
		this->K=K;
	}
	virtual void setSizeVocabulary(int V){
		this->V=V;
	}
//	virtual const Vector2DDouble& getPhi() const{
//		return phi;
//	}
//	virtual const Vector2DDouble& getTheta() const{
//		return theta;
//	}
	virtual const Vector2DDouble& getPhi() const{
		return phi;
	}
	virtual const Vector2DDouble& getTheta() const{
		return theta;
	}
	virtual const VecDouble& getThetasum() const{
		return thetasum;
	}
	virtual int getNumTopics() const{
		return K;
	}
	virtual const Vector2DInt& getDocuments() const{
		return docs;
	}

//	const VecDouble& getBetaVec() const {
//		return betaVec;
//	}
//
//	void setBetaVec(const VecDouble& betaVec) {
//		this->betaVec = betaVec;
//	}

	virtual void clearPhi(){
		phi.clear();
	}
	virtual void clearTheta(){
		theta.clear();
	}
	virtual void copyState(ILDA* _lda){
		LDA* lda=static_cast<LDA*> (_lda);
		K=lda->K;
		V=lda->V;
		nk=lda->nk;
		nwk=lda->nwk;
		beta=lda->beta;
		alpha=lda->alpha;
		alphaVec=lda->alphaVec;
		betasum=lda->betasum;
		alphasum=lda->alphasum;
	}
	virtual void copyState(void* _lda, std::string impl){
		if(impl=="SparseLDA_Speed"){
			SparseLDA_Speed* lda=static_cast<SparseLDA_Speed*> (_lda);
			K=lda->K;
			V=lda->V;
			nk=lda->nk;
			nwk.resize(V);
			for(int w=0;w<V;w++){
				nwk[w].resize(K);
				for(int k=0;k<K;k++){
					nwk[w][k]=lda->getnwk(w,k);
				}
			}
			beta=lda->beta;
			alpha=lda->alpha;
			alphaVec=lda->alphaVec;
			betasum=lda->betasum;
			alphasum=lda->alphasum;
//		}else if(impl=="AliasLDA_UltraSpeed"){
//			AliasLDA_UltraSpeed* lda=static_cast<AliasLDA_UltraSpeed*> (_lda);
//			K=lda->K;
//			V=lda->V;
//			nk=lda->nk;
//			nwk=lda->nwk;
//			beta=lda->beta;
//			alpha=lda->alpha;
//			alphaVec=lda->alphaVec;
//			betasum=lda->betasum;
//			alphasum=lda->alphasum;
		}else if(impl=="LDA"){
			LDA* lda=static_cast<LDA*> (_lda);
			K=lda->K;
			V=lda->V;
			nk=lda->nk;
			nwk=lda->nwk;
			beta=lda->beta;
			alpha=lda->alpha;
			alphaVec=lda->alphaVec;
			betasum=lda->betasum;
			alphasum=lda->alphasum;
		}
	}
};

} /* namespace AliasLDA */

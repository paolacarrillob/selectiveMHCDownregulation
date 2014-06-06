/*
 * Genes.h
 *
 *  Created on: Apr 20, 2011
 *      Author: paola
 */

#ifndef GENES_H_
#define GENES_H_

#include "Bitstring.h"
#include "math.h"
#include <fstream>
#include <sstream>

#include <list>
#include <map>
#include <utility>
#include <bitset>

class Gene {
public:
	Gene();
	virtual ~Gene(){};
	bool operator ==(Gene& rhs);
	void SetGeneFunctionality(bool functionality);
	void SetGeneExpression(bool expression);
	void SetGeneSpecificity(); //works DO I NEED IT?
	void SetGeneSpecificity(int s){specificity = s;};//DO I NEED IT?
	
	//void SetPseudogene(bool expression);
	//void SetPseudogene(int M_id);
	void SetPseudogene(int M_id_A, int M_id_B);
	
	int GetGeneMidA(){return mID_A;};
	int GetGeneMidB(){return mID_B;};
	
	int GetGeneSpecificity(){return specificity;};
	bool IsFunctional(){return functional;};
	bool IsExpressed(){return isExpressed;};
	bool IsPseudogene(){return pseudoGene;};
	Gene& Copy(Gene& rhsGene);//works;
	void SetGeneID(int _ID){geneID = _ID;};
	int GetGeneID(){return geneID;};

	void Mutate(double mutationRate); //works IS OLD
	void PointMutation();//works
	void PrintBits();//works
	void MutateSpecificity();//works
	bool IsGeneUnique(vector<Gene>& genePool, int counter);//works
	int BindMolecule(const Gene& anotherMolecule); //works NEW
	int CountHammingDistance(Gene& anotherMolecule);

	void SaveBackupGenes(fstream& backupFile);//works
	void SaveGenes (fstream& outfile);
	void PrintParameters();
	string RestoreGenes(stringstream& sgline);
	//pair <int,int> genePhenotype;

protected:
	bool functional;
	bool isExpressed;
	int geneID;
	int specificity;
	bool pseudoGene;
	//unsigned long int mID;
	unsigned long int mID_A;
	unsigned long int mID_B;


};


class MHCGene: public Gene
{
public:
	MHCGene();
	virtual ~ MHCGene(){};
};

class KIRGene: public Gene
{
public:
	KIRGene(){};
	KIRGene(int L);
	virtual ~ KIRGene(){};
	enum state{inhibitory, activating};

	bool operator ==(KIRGene& rhs);
	KIRGene& Copy(KIRGene& rhsGene);//works;
	void SaveBackupGenes(fstream& backupFile); //works
	void SaveGenes (fstream& outfile);//works
	string RestoreGenes(stringstream& sgline); //works
	int GetGeneType(){return geneType;}; //works
	void SetGeneType(int type);
	void MutateReceptorType();//works
	bool IsInhibitory();
	bool IsActivating();
	void PrintGenes();


protected:
	state geneType;
};


class GenePool
{
public:
	GenePool(){};
	virtual ~GenePool(){};
	float GetAlleleFreq(int alleleIndex); //works
	bool GeneAlreadyInPool(int geneID); //works
	unsigned int GetPoolSize(){return genes.size();}
	int RandomlyPickGene(bool distribution); //works
	vector<int>& GetGenes(){return genes;};
	void WriteOutGenes(const string& fileName);//works
	int ComparePools(GenePool& otherPool);

protected:
	vector<int> genes;//a vector of integers
	fstream genePoolFile;
};

class MHCGenePool: public GenePool
{
public:
	MHCGenePool(){};
	void FillMHCGenePool(int size);//works
	void FillMHCGenePoolWithSimilarMHCs(int size);
	virtual ~MHCGenePool(){};
protected:
	int poolSize;
};

class KIRGenePool: public GenePool
{
public:
	KIRGenePool(GenePool& mhcPool, bool dist, int specificity);
	virtual ~KIRGenePool(){};

};


class Map : public GenePool
{
public:
	Map(){};
	virtual ~Map(){};
	void FillMap(GenePool& mhcpoolA, GenePool& mhcpoolB, KIRGene& gene); //works now with more pools!
	bool IsGeneInMap(KIRGene& gene); //works
	unsigned int GetMapSize(){return mapGenes.size();} //works

	string RestoreMap(stringstream& smline); //works
	void SaveBackupMap(fstream& backupFile); //works
	multimap< pair< int,int>,  pair<int, pair<int,int> > >& GetMap() {return mapGenes;}; //works also with more pools
	
//protected:
	multimap <pair< int,int>, pair <int, pair<int,int> > > mapGenes; //the map has the design: (geneID,gene_type), (L, mID)
};


#endif /* GENES_H_ */

/*
 30.05. expanded with mID_A, mID_B, changed the maps accordignly. 
 *
 */

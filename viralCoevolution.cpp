//============================================================================
// Name        : CoevolutionKIR.cpp
// Author      : P.Carrillo-Bustamante
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string.h>
#include "World.h"
#include "Exceptions.h"
using namespace std;

int main(int argc, char*argv[])
{
	//testing the new Virus members
	/*Virus dummyVirus;
	Virus nastyVirus;
	nastyVirus.SetViralParameters(0.0, 0.0,0.0,0.0, 2, false, false, 0);
	cout << "dummyVirus"<<endl;
	dummyVirus.PrintParametersVirus();
	cout << endl<< "nastyVirus"<<endl;
	nastyVirus.PrintParametersVirus();
	dummyVirus.Copy(nastyVirus);
	cout<<endl<<"dummyVirus again:"<< endl;
	dummyVirus.PrintParametersVirus();
	bool bla_both = dummyVirus.IsDownregulatingMHC_Both();
	bool bla_b = dummyVirus.IsDownregulatingMHC_B();
	bool bla_a = dummyVirus.IsDownregulatingMHC_A();
	cout <<bla_both<<"|"<<bla_a<<"|"<<bla_b<<endl;
	exit(-1);
	/*
	//testing the assembly of the MHC genes
	World myWorld;
	myWorld.LoadParameterFile("Parameters.data");
	myWorld.Initialize();	
	
	myWorld.IntroduceVirus(dummyVirus,0.0, 2);
	
	list<Infection>::iterator inf;
	for(inf = myWorld.hosts.at(0).infections.begin(); inf!=myWorld.hosts.at(0).infections.end(); inf++) //check for every infection within one host, whether it's time to be cleared
	{
		myWorld.hosts.at(0).ClearInfection(0.0,*inf);
		cout << "host # 0"<<endl;
	}
	for(inf = myWorld.hosts.at(1).infections.begin(); inf!=myWorld.hosts.at(1).infections.end(); inf++) //check for every infection within one host, whether it's time to be cleared
	{
		myWorld.hosts.at(1).ClearInfection(0.0,*inf);
		cout << "host # 1"<<endl;
	}
	exit(-1); */
	
	//*/
	
	//testing the distribution of Lc,NKR<->MHC
	/*vector<int> Lcs;
	 
	 for(int i = 0; i<10000; i++)
	 {
	 MHCGenePool mhc;
	 mhc.FillMHCGenePool(30);
	 KIRGene dummy(2);
	 //mhc.FillMHCGenePoolWithSimilarMHCs(30);
	 
	 for(unsigned int j = 0; j <mhc.GetPoolSize(); j++)
	 {
	 Gene dummy_mhc;
	 dummy_mhc.SetGeneID(mhc.GetGenes().at(j));
	 int dummy_Lc = dummy.BindMolecule(dummy_mhc);
	 Lcs.push_back(dummy_Lc);
	 }
	 }
	 int size = 30;
	 for(int i = 0; i<=16; i++)
	 {
	 double mycount = count(Lcs.begin(), Lcs.end(), i);
	 mycount = mycount/10000;
	 double percentage_mhc = mycount/size;
	 
	 cout << i << " " << mycount <<" "<<percentage_mhc<<endl;
	 }
	 exit(-1);
	 //*/
	
	//testing the distribution of Lc, NKR<-> MHC with fixed L
	/*
	 for(int i = 1; i<=16; i++)
	 {
	 MHCGenePool mhc;
	 mhc.FillMHCGenePool(1000);
	 int size = mhc.GetPoolSize();
	 int counter = 0;
	 for(unsigned int j = 0; j <size; j++)
	 {
	 Gene dummy_mhc;
	 dummy_mhc.SetGeneID(mhc.GetGenes().at(j));
	 for(int k = 0; k<1000; k++)
	 {
	 KIRGene dummy(i);
	 int dummy_Lc = dummy.BindMolecule(dummy_mhc);
	 if(dummy_Lc >= dummy.GetGeneSpecificity())
	 {
	 counter ++;
	 }
	 }
	 }
	 double fraction = counter /(mhc.GetPoolSize());
	 double double_fraction = 30*fraction/1000;
	 
	 cout << i << " "<<fraction <<" " <<double_fraction<<endl;
	 }
	 
	 exit(-1);*/
	
	/*
	 //testing the overlap between two MHC pools
	 vector<int> overlap;
	 for (int i = 0; i<1000; i++)
	 {
	 
	 MHCGenePool poolA;
	 //poolA.FillMHCGenePool(20);
	 poolA.FillMHCGenePoolWithSimilarMHCs(30);
	 MHCGenePool poolB;
	 //poolB.FillMHCGenePool(20);
	 poolB.FillMHCGenePoolWithSimilarMHCs(30);
	 
	 int a = poolA.ComparePools(poolB);
	 overlap.push_back(a);
	 //cout <<"the overlap between poolA and pool B is: "<<a<<endl;
	 }
	 for(int i = 0; i<=30; i++)
	 {
	 double mycount = count(overlap.begin(), overlap.end(), i);
	 mycount = mycount/1000;
	 //double percentage_mhc = mycount/size;
	 
	 cout << i << " " << mycount <<endl;
	 }
	 
	 exit(-1);*/
	
	/*Map myMap;
	 KIRGene bla(1);
	 bla.SetGeneID(4758);
	 myMap.FillMap(poolA, poolB, bla);
	 bla.PrintGenes();
	 exit(-1);*/
	
	cerr << "Testing if stderr is redirected to stdoutput" << endl;
	if (argc<2)
	{
		cerr << "Usage: "<<argv[0] << " <Parameter file> <-b Loading Backup file> \n CoevolutionKIR simulates the evolution of the complex KIR system. For correct usage you have to indicate a parameter file and a backup file name. If you are not loading any backup, just give any random name"<< endl;
		exit(-1);
	}
	// delete the existing host_file
	char buffer[512];
	//Ouss: Use boost library here for boost:filesystem to delete files... etc
	cout <<"deleting old files ..."<<endl;
	//sprintf(buffer, "rm *.txt");
	sprintf(buffer, "rm *.log");
	system(buffer);

	string parameterFile(argv[1]);
	string backupFile;
	bool loadingBackup=false;
	for(int i=2; i< argc;i++)
	{
		if(strcmp(argv[i],"-b")==0 && argc > i)
		{
			loadingBackup = true;
			backupFile.assign(argv[++i]);
		}
	}

	World theWorld;
	try
	{
		theWorld.LoadParameterFile(parameterFile);
	}
	catch (OussException& e ) {
		cout << e.GetErrorData() << endl;
		exit(-1);
	}
	catch (...)
	{
		cout << "unknown exception thrown" <<endl;
		exit(-1);
	}
	
	//theWorld.LoadParameterFile(parameterFile);
	cout <<"welcome, hello" << endl;
	if(!loadingBackup)
	{
		// initialize host population
		theWorld.Initialize();
	}
	else
	{
		cout << "\n Loading data from backup file: "  << backupFile << endl;
		theWorld.LoadBackupFile(backupFile);
	}

	theWorld.Simulate();
	cout << "bye bye \n";

	return 0;
}


/*
 *24.10.2012 fixed bug in LoadBackupFile! -> fixed it again (07.10.2013)... learned about ss.clear() and getline :-)
 *
 */

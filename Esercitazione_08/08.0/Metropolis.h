#ifndef _Metropolis_h
#define _Metropolis_h

#include "random.h"
#include "funcvirtual.h"
#include <vector>

using namespace std;

class Metropolis{

private:

  Random _rnd;
  Function *_DistributionToBeSampled;    
  int _Dimensions;       
  double _StepWidth;   
  int _NumberOfSteps;    
  vector <vector <double> > _Positions;     
  vector <double> _NextPositions;            
  int _StepsDone;        
  int _StepAccepted;                               
  int _StepRefused;

  //Istogramma
  double _Extremes; //Estremi per il calcolo istogramma
	int _NumberOfBins; 
	vector<vector<double>> _Psi; //Serve per l'istogramma delle psi
	int _NumberOfIteractions; //Numero di istogrammi su cui mediare      

public:

  Metropolis(Function *DistributionToBeSampled, vector <double> InitialPositions, double StepWidth, int NumberOfSteps, int Dimensions, int MCSteps, int NumberOfBins);
  ~Metropolis();

  void Restart();
  void Equilibrate(int NumberOfEquilibrationSteps, const char * trial);
  void RunUniform(void);
  void RunGauss(void);
  void NextUniform(void);    //Distribuzione trial == uniforme
  void NextGauss(void);      //Distribuzione trial == gaussiana
  bool Acceptance(void);     //Probalità di accettazione
  void UniformStep(void);
  void GaussStep(void);
  void AcceptanceRate(void);
  void Print(const char*);
  double Integrate(Function *hamiltonian);
  void Histogram(int index); //Aggiornamento dell'istogramma della funzione d'onda di trial
	void AnalysisHistogram(int NumberOfBlocks, const char* OutputFile);//Data blocking sull' istogramma di psi
  
  vector<double> Rays();
  vector<vector <double>> *GetSampleWafefunction(void);
};

#endif

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

public:

  Metropolis(Function *DistributionToBeSampled, vector <double> InitialPositions, double StepWidth, int NumberOfSteps, int Dimensions);
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
  vector<double> Rays();

};

#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "DataBlocking.h"
#include "funcvirtual.h"
#include "Metropolis.h"

using namespace std;

int main (int argc, char *argv[]){

  int dimension = 1;
  double metropolisStep = 2.6;
  int numberOfSteps = 1000;
  int numberOfIntegrals = 10000;
  int numberOfEquilibratioSteps = 1000;
  int numberOfBlocks = 100;       //Inutile per questo punto
  int numberOfBins = 200;         //Inutile per questo punto

  WaveFunction* PsiTrial = new WaveFunction();
  Function* H = new Hamiltonian(PsiTrial);

	vector<double> initialPositions(dimension,0);
	initialPositions.at(0) = 0;

  Metropolis walker(PsiTrial, initialPositions, metropolisStep, numberOfSteps, dimension, numberOfIntegrals, numberOfBins);
  vector<double> resultsIntegration(numberOfIntegrals,0);
  ofstream OptimizationParameters;
  OptimizationParameters.open("../Results/08.0/optimization.parameter.prova.out");

  walker.Equilibrate(numberOfEquilibratioSteps, "uniform");
  for(double mi=0.795; mi<0.81; mi+=0.001){
    cout << "Mi: " <<mi << endl;
    PsiTrial->SetMi(mi);
    for(double sigma=0.61; sigma<0.62; sigma+=0.001){
      cout << "Sigma: " << sigma << endl;
      PsiTrial->SetSigma(sigma);
      for(int i=0; i<numberOfIntegrals; i++){
        walker.RunUniform();
        resultsIntegration.at(i) = walker.Integrate(H);
        walker.Restart();
      }
      OptimizationParameters << mi << "  " << sigma << "  " << Mean(resultsIntegration) << endl;
    }
  }
  OptimizationParameters.close();

  delete PsiTrial;
  delete H;
  return 0;
}

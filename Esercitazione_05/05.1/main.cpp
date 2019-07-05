#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"
#include "DataBlocking.h"
#include "funcvirtual.h"
#include "Metropolis.h"

using namespace std;

int main (int argc, char *argv[]){

  //DATI GENERALI:
  int dimensions = 3;                       //Dimensioni del problema in questione: sono punti nello spazio x,y,z
  int numberOfSamples = 1000000;
  int numberOfEquilibrationSteps = 1000;
  int numberOfBlocks = 100;

  cout << endl << "==========================================================" << endl;
  cout << "========================ORBITALE 1S=======================" << endl;
  cout << "==========================================================" << endl;
  cout << "Distribuzione di trial uniforme..." << endl;
  //Passo del metropolis in unità del raggio di Bohr
  double step_1s_uniform = 1.2;

  //Orbitale da campionare
  Function *orbitale_1s = new Uno_s(); 

  //Posizioni di partenza del Metropolis
	vector<double> initialPositions(dimensions,0);

	initialPositions.at(0) = 1.;      
	initialPositions.at(1) = 1.;
	initialPositions.at(2) = 1.;


  Metropolis walker_1s_uniform(orbitale_1s, initialPositions, step_1s_uniform, numberOfSamples, dimensions);

  walker_1s_uniform.Equilibrate(numberOfEquilibrationSteps, "uniform");
	walker_1s_uniform.RunUniform();
  walker_1s_uniform.AcceptanceRate();
	walker_1s_uniform.Print("../Results/position_1s_unif.out");

	vector<double> rays_1s_uniform = walker_1s_uniform.Rays();
  Analysis(numberOfBlocks, rays_1s_uniform, "../Results/rays_1s_unif.out");

  cout << endl << "----------------------------------------" << endl;
  cout << "Distribuzione di trial gaussiana..." << endl;

  //Passo del metropolis in unità del raggio di Bohr
  double step_1s_gauss = 0.75;

  Metropolis walker_1s_gauss(orbitale_1s, initialPositions, step_1s_gauss, numberOfSamples, dimensions);
  
  walker_1s_gauss.Equilibrate(numberOfEquilibrationSteps, "gauss");
	walker_1s_gauss.RunGauss();
  walker_1s_gauss.AcceptanceRate();
	walker_1s_gauss.Print("../Results/position_1s_gauss.out");

	vector<double> rays_1s_gauss = walker_1s_gauss.Rays();
  Analysis(numberOfBlocks, rays_1s_gauss, "../Results/rays_1s_gauss.out");

  cout << endl << "----------------------------------------" << endl;

  cout << "==========================================================" << endl;
  cout << "========================ORBITALE 2P=======================" << endl;
  cout << "==========================================================" << endl;
  cout << "Distribuzione di trial uniforme..." << endl;

  //Passo del metropolis in unità del raggio di Bohr
  double step_2p_uniform = 2.95;

  //Orbitale da campionare
  Function *orbitale_2p = new Due_p();
  //Posizioni di partenza del metropolis

  initialPositions.at(0) = 3.;
  initialPositions.at(1) = 3.;
  initialPositions.at(2) = 3.;

  Metropolis walker_2p_uniform(orbitale_2p, initialPositions, step_2p_uniform, numberOfSamples, dimensions);

	walker_2p_uniform.Equilibrate(numberOfEquilibrationSteps, "uniform");
	walker_2p_uniform.RunUniform();
  walker_2p_uniform.AcceptanceRate();
	walker_2p_uniform.Print("../Results/position_2p_unif.out");

	vector<double> rays_2p_uniform = walker_2p_uniform.Rays();
  Analysis(numberOfBlocks, rays_2p_uniform, "../Results/rays_2p_unif.out");

  cout << endl << "----------------------------------------" << endl;
  cout << "Distribuzione di trial gaussiana..." << endl;
  //Passo del metropolis in unità del raggio di Bohr
  double step_2p_gauss = 1.85;

  Metropolis walker_2p_gauss(orbitale_2p, initialPositions, step_2p_gauss, numberOfSamples, dimensions);
  
  walker_2p_gauss.Equilibrate(numberOfEquilibrationSteps, "gauss");
	walker_2p_gauss.RunGauss();
  walker_2p_gauss.AcceptanceRate();
	walker_2p_gauss.Print("../Results/position_2p_gauss.out");

	vector<double> rays_2p_gauss = walker_2p_gauss.Rays();
  Analysis(numberOfBlocks, rays_2p_gauss, "../Results/rays_2p_gauss.out");

  cout << endl << "----------------------------------------" << endl;

  delete orbitale_1s;
  delete orbitale_2p;

  return 0;

}

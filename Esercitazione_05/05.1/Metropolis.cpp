#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"
#include "funcvirtual.h"
#include "DataBlocking.h"
#include "Metropolis.h"

using namespace std;

//Costruttore
Metropolis::Metropolis(Function *DistributionToBeSampled, vector <double> InitialPositions, double StepWidth, int NumberOfSteps, int Dimensions): _NextPositions(Dimensions){

  SetRandomGenerator(&_rnd);

  _DistributionToBeSampled = DistributionToBeSampled;
  _StepWidth = StepWidth;
  _NumberOfSteps = NumberOfSteps;
  _Dimensions = Dimensions;
  _StepsDone = 0;
  _StepAccepted = 0;
  _StepRefused = 0;
  _Positions.push_back(InitialPositions);

}

//Distruttore
Metropolis::~Metropolis(){
  _rnd.SaveSeed();
}

//Metodo che riinizializza tutto a zero per riiniziare il metropolis dall'ultima posizione in memoria
void Metropolis::Restart(){

  _StepsDone = 0;
  _StepAccepted = 0;
  _StepRefused = 0;

  _Positions.erase(_Positions.begin(), _Positions.end() -1)  ;

}

//Metodo per equilibrare il sistema
void Metropolis::Equilibrate(int NumberOfEquilibrationSteps, const char * trial){

	int temp = _NumberOfSteps;  
	_NumberOfSteps = NumberOfEquilibrationSteps;

	//A seconda di quale distribuzione di trial sto usando utilizzo Run diversi
	if(strcmp(trial,"uniform") == 0){    
		RunUniform();  
	}
	else if(strcmp(trial,"gauss") == 0){
		RunGauss();   
	}

	Restart();
	_NumberOfSteps = temp;

}

//Possibile passo uniforme
void Metropolis::NextUniform(void){         

  for(int i=0; i<_Dimensions; i++){
    _NextPositions.at(i) = _rnd.Rannyu(-_StepWidth, _StepWidth) + _Positions[_StepsDone][i];
  }

}

//Possibile passo gaussiano
void Metropolis::NextGauss(void){          

  for(int i=0; i<_Dimensions; i++){
    _NextPositions.at(i) = _rnd.Gauss(0, _StepWidth) + _Positions[_StepsDone][i];     
  }

}

//Probabilità di accettazione
bool Metropolis::Acceptance(){

	//Distribuzione di probabilità valutata nel punto in cui mi trovo
	double actualProbability = _DistributionToBeSampled-> Eval(_Positions[_StepsDone]);
	//Distribuzione di probabilità valutata nel possibile nuovo punto
	double newProbability = _DistributionToBeSampled-> Eval(_NextPositions);
	//Rapporto tra le distribuizioni di probabilità
	double probabilityRatio = newProbability/actualProbability;
	//Probabilità di accettazione
	double acceptance = min(1., probabilityRatio);

	if(acceptance >= 1){
		return true;
	}

	else{
		double a = _rnd.Rannyu();
		if( a <= acceptance){
			return true;
		}
		return false;
	}

}

void Metropolis::UniformStep(){

	NextUniform();    
	bool acceptance = Acceptance();

	if(acceptance == true){

    _Positions.push_back(_NextPositions);
    _StepAccepted++;

	}
	else{
		vector<double> temp(_Dimensions);	//Vettore di appoggio
    	_StepRefused++;
		for(int i=0; i<_Dimensions; i++){
			temp.at(i) = _Positions[_StepsDone][i];
		}
		_Positions.push_back(temp);

	}

	//Incremento del numero totale di steps fatti
	_StepsDone ++;  

}

void Metropolis::GaussStep(){

	NextGauss();   
	bool acceptance = Acceptance();

	if(acceptance == true){

    _Positions.push_back(_NextPositions);
    _StepAccepted++;

	}
	else{
		vector<double> appoggio(_Dimensions);
    	_StepRefused++;
		for(int i=0; i<_Dimensions; i++){
			appoggio.at(i) = _Positions[_StepsDone][i];
		}
		_Positions.push_back(appoggio);

	}

	//Incremento del numero totale di steps fatti
	_StepsDone ++;

}

//UnifromStep ripetuto NumberOfSteps volte
void Metropolis::RunUniform(){

	for(int i=0; i < _NumberOfSteps; i++){
		UniformStep();
	}

}

//GaussStep ripetuto _StepAccepted volte
void Metropolis::RunGauss(){

	for(int i=0; i < _NumberOfSteps; i++){
		GaussStep();
	}

}

//Metodo per controllare la regola empirica del 50%
void Metropolis::AcceptanceRate(){
	cout << endl << "Accepted steps = " << _StepAccepted << endl;    
  	cout << endl << "Refused steps = " << _StepRefused << endl;
	cout << endl << "ACCEPTANCE RATE: " << ((double)_StepAccepted/(double)_NumberOfSteps) * 100. << "%" << endl;
}

//Metodo per stampare su file i punti campionati
void Metropolis::Print(const char* filename){

	ofstream Results;
	Results.open(filename);

	for(int j = 0; j < _NumberOfSteps+1; j++){
		for(int i=0; i <_Dimensions; i++){
			Results <<  _Positions[j][i] << " " ;
		}
		Results << endl;
	}

  Results.close();
}

//Metodo legato all'esercizio 05.1
//Metodo per ottenere i raggi dei punti (x,y,z) campionati
vector<double> Metropolis::Rays(){

	vector<double> rays(_NumberOfSteps+1,0);

	for(int j = 0; j < _NumberOfSteps+1; j++){
		for(int i=0; i < _Dimensions; i++){
			rays.at(j) += pow( _Positions[j][i], 2);
		}
	rays.at(j) = sqrt(rays.at(j));
	}

return rays;
}

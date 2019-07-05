#include "Integrate.h"
#include "random.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

//Costruttore
Integrate::Integrate(double LowerBound, double UpperBound, Function * function){

	_integrand = function;
	//ordino gli estremi in ordine crescente e se non dovessero già essere ordinati ne tengo conto
	_LowerBound = min(LowerBound,UpperBound);
	_UpperBound = max(LowerBound,UpperBound);
	if ( _LowerBound>= _UpperBound){
		_Sign=-1;
	}
	else _Sign=1;

	Start();

}

//Distruttore
Integrate::~Integrate() {
	_rnd.SaveSeed();
}

void Integrate::Start(){

	// preparo il generatore di numeri casuali
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");


	if (Primes.is_open()){
	   Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
	   while ( !input.eof() ){
		input >> property;
		if( property == "RANDOMSEED" ){
	  		 input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	  	 _rnd.SetRandom(seed,p1,p2);
	 	}
	     }
	    input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

}

//Calcolo dell'integrale montecarlo usando numeri distribuiti uniformemente tra [_LowerBound,_UpperBound)
double Integrate::Media(int NumberOfIntervals){

	_PartialSum = 0;
	double mean = 0;
	double interval = _UpperBound - _LowerBound;

	for (int i=1; i <= NumberOfIntervals; i++){
		double x = _rnd.Rannyu(_LowerBound,_UpperBound);
		double f = _integrand -> Eval(x);
		_PartialSum += f;
	}

	mean = (_PartialSum/double(NumberOfIntervals));
	_Integral = _Sign * mean * interval;

	return _Integral;

};

double Integrate::MediaRetta(int NumberOfIntervals){

	_PartialSum = 0;
	double mean = 0;
	double interval = _UpperBound - _LowerBound;;

	for (int i=1; i <= NumberOfIntervals; i++){
		double x = _rnd.Retta();
		double f = _integrand -> Eval(x);
		_PartialSum += f;
	}

	mean = (_PartialSum/double(NumberOfIntervals));
	_Integral = _Sign * mean * interval;

	return _Integral;

};

//restituisce il valore dell'integrale
double Integrate::GetResult() const{
	return _Integral;
};

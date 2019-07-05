#include "Option.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"
#include <string>

using namespace std;

//Costruttore
Option::Option(double InitialPrice, double StrikePrice, double RiskFreeRate, double Volatility){


  Start();

   _InitialPrice = InitialPrice;
   _StrikePrice = StrikePrice;
   _RiskFreeRate = RiskFreeRate;
   _Volatility = Volatility;
   _SpotPrice = InitialPrice;

}

//Distruttore
Option::~Option(){
  _rnd.SaveSeed();
}

void Option::Start(){

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

//Metodo che reimposta lo SpotPrice al valore iniziale
void Option::Restart(){
  _SpotPrice = _InitialPrice;
}

//Metodo per campionare valori da un moto browniano geometrico
void Option::GeometricBrownianMotion(double start, double end){
  _SpotPrice =  _SpotPrice*exp((_RiskFreeRate - 0.5*pow(_Volatility,2))*(end - start) + _Volatility*_rnd.Gauss(0,1)*sqrt(end - start));

}

//Pricing Put option con la possibilità di discretizzare il processo
double Option::PutPrice(double start, double end, int steps){
  double increase = (end - start)/steps;

  for(int i = 0; i< steps; i++) GeometricBrownianMotion(start + i*increase, start + (i+1)*increase);

  double actualization = exp(-_RiskFreeRate*(end-start));
  double profit = max(0., _StrikePrice - _SpotPrice);
  double PutPrice = actualization * profit;

  Restart();
  return PutPrice;

}

//Pricing Call option con la possibilità di discretizzare il processo
double Option::CallPrice(double start, double end, int steps){

  double increase = (end - start)/steps;
  for(int i = 0; i< steps; i++) GeometricBrownianMotion(start + i*increase, start + (i+1)*increase);

  double actualization = exp(-_RiskFreeRate*(end-start));
  double profit = max(0., _SpotPrice - _StrikePrice);
  double CallPrice = actualization * profit;
  
  Restart();
  return CallPrice;
}

#include <cmath>
#include "funcvirtual.h"

using namespace std;

//Implementazione delle funzioni da integrare

Function::~Function(){}

//Primo punto
double Coseno::Eval(double x) const {
	return (M_PI/2)*cos(x*(M_PI/2));
}
Coseno::~Coseno(){}

//Secondo punto
double gfunc::Eval(double x) const {
	return ((M_PI/2)*cos(x*(M_PI/2)))/(2*(1-x));
}

gfunc::~gfunc(){}
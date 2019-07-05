#ifndef _funcvirtual_h
#define _funcvirtual_h
#include <cmath>

// Classe contentente tutte le funzioni che dovrò intergrare

class Function {
	public:
		virtual double Eval(double x) const=0;
		virtual ~Function() = 0;
};


// Funzione da integrare nel primo punto dell' esercizio 02.1
class Coseno: public Function {

	public:
		virtual double Eval(double x) const;
		~Coseno();
};

//Funzione da integrare nel secondo punto dell'esercizio 02.1
class gfunc: public Function {				

	public:
		virtual double Eval(double x) const;
		~gfunc();
};

#endif

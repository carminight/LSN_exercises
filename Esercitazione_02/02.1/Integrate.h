#ifndef _Integrate_h
#define _Integrate_h
#include "random.h"
#include "funcvirtual.h"

class Integrate{

  private:

    Random _rnd; 
    //Utile per tenere conto di eventuali scambi negli estremi di integrazione
    int _Sign; 
    //Estremi di integrazione
    double _LowerBound, _UpperBound; 
    double _PartialSum; 
    double _Integral; 
    Function * _integrand; 

  public:

    //Costruttore
    Integrate(double LowerBound, double UpperBound, Function * function); 
    //Distruttore
  	~Integrate(); 
    //Calcola l'integrale montecarlo generando numeri uniformenete distribuiti tra _LowerBound e _UpperBound suddiviso in NumberOfIntervals intervalli con in metodo della media
  	double Media(int NumberOfIntervals);
    // Calcola l'integrale con numeri uniformemente distribuiti tra 0 e 1 con p(x) = 2(1-x)
    double MediaRetta(int NumberOfIntervals); 
    // Restituisce il risultato salvato in precedenza
  	double GetResult() const; 
    void Start();

};

#endif

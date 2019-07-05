#ifndef _Option_h
#define _Option_h
#include "random.h"

class Option{

private:

  Random _rnd;
  double _InitialPrice;
  double _StrikePrice;
  double _RiskFreeRate;
  double _Volatility;
  double _SpotPrice;


public:

  Option(double InitialPrice, double StrikePrice, double RiskFreeRate, double Volatility);
  ~Option();

  void Start();
  void Restart();
  void GeometricBrownianMotion( double start, double end);    //Moto browniano geometrico applicato allo SpotPrice
  double PutPrice(double start, double end, int steps);
  double CallPrice(double start, double end, int steps);

};

#endif

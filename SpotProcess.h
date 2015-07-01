#ifndef __SpotProcess_h__
#define __SpotProcess_h__

#include <vector>
#include <math.h>
using namespace std;


class SpotProcess{
public:
	//this class encapsulates the pay off of a european call option on spot
	//also allows retrieving the advection, difusion and source/sink terms
	//of the black scholes PDE of Spot process, with short term risk free rate 
	//and a cost of carry rate
	//the PDE is expressed in units of natural logarithm of Spot, since this way
	//the discretization error is lower than when using the PDE in Spot units
	SpotProcess()
	{
		Spot_m = 100; Strike_m = 100; shortTermRiskFreeRate_m = 1E-2; CostOfCarry_m = 1E-2; Sigma_m = 0.50; Expiry_m = 1.0;
	}
	
	SpotProcess(const double& Spot, const double& Strike, const double& ShortTermRiskFreeRate, const double& CostOfCarry, const double& Sigma, const double& Expiry)
	{
		Spot_m = Spot;
		Strike_m = Strike;
		shortTermRiskFreeRate_m = ShortTermRiskFreeRate;
		CostOfCarry_m = CostOfCarry;
		Sigma_m = Sigma;
		Expiry_m = Expiry;
	}
	
	inline double PdeAdvectionConvection(double S, double t) const 
	{
		return (shortTermRiskFreeRate_m - CostOfCarry_m - 0.5*Sigma_m*Sigma_m);
	}
	
	inline double PdeDifusion(double S, double t) const 
	{
		double val = 0.5*Sigma_m*Sigma_m;
		return val;
	}
	
	inline double PdeSourceSink(const double& S, const double& t) const 
	{
		return -shortTermRiskFreeRate_m;
	}
	
	virtual double Payoff(const double& x, const double& t) const 
	{
		return std::exp(x) - Strike_m > 0 ? std::exp(x) - Strike_m : 0;
	}
	
	double getSigma() const 
	{
		return Sigma_m;
	}
	
	double getShortTermRiskFreeRate() const 
	{
		return shortTermRiskFreeRate_m;
	}
	
	double getSpot() const 
	{
		return Spot_m;
	}

	void operator=(SpotProcess&) = delete;

	SpotProcess(const SpotProcess&) = delete;

private:
	double Spot_m;
	double Strike_m;
	double shortTermRiskFreeRate_m;
	double CostOfCarry_m;
	double Sigma_m;
	double Expiry_m;

};

#endif //__SpotProcess_h__




#ifndef MultiThreadPDE_Explicit_h
#define MultiThreadPDE_Explicit_h

#include <thread>
#include "ThreadPool.h"
#include "spline.h"

using namespace std;

/****************************************************************************
* Author: Horacio Aliaga	
* Copyright (C) 2015 Horacio Aliaga (horacio.aliaga at gmail.com)
* Description:  Library providing a multithreading Black Scholes PDE Solver
*               using the explicit method of order one.  
*				The algorithm in unconditionally stable provided the infinite norm
*				of the heat equation matrix is lower than 1
* Dependencies: this library is using ThreadPool, that can be found on:
*               https://github.com/progschj/ThreadPool
*
*				and is also using spline.h
* Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
*               
****************************************************************************/

template<class T>
class MultiThreadPDE_Explicit
{
public:

	MultiThreadPDE_Explicit(T &process) : process_m(process) {}

	MultiThreadPDE_Explicit(
		T &process,					     //process returs the diffusion, advection and sink/source terms of the parabolic equation
		int numberOfTimeNodes,           //number of partitions of the time domain
		double initialTime,              //lower boundary for time domain
		double finalTime,                //upper boundary for time domain
		int numberOfSpatialNodes,        //number of partitions of the spatial domain   
		double lowerBoundLogUnderlying,  //lower boundary for spatial domain      
		double upperBoundLogUnderlying   //upper boundary for time domain      
		) :
		process_m(process),
		totalTimePoints_m(numberOfTimeNodes),
		initialTime_m(initialTime),
		finalTime_m(finalTime),
		totalSpatialPoints_m(numberOfSpatialNodes),
		spatialLowerBound_m(lowerBoundLogUnderlying),
		spatialUpperBound_m(upperBoundLogUnderlying),
		deltaX_m((upperBoundLogUnderlying - lowerBoundLogUnderlying) / (numberOfSpatialNodes + 1.0)),
		deltaT_m((finalTime - initialTime) / (numberOfTimeNodes - 1.0)),
		spatialGrid_m(numberOfSpatialNodes + 2),
		currentSolution_m(numberOfSpatialNodes + 2),
		previousSolution_m(numberOfSpatialNodes + 2)
	{
		//initializes lower and upper boundaries of spatialGrid_m and currentSolution_m
		spatialGrid_m.front() = spatialLowerBound_m;
		currentSolution_m.front() = process_m.Payoff(spatialGrid_m.front(), initialTime_m);
		spatialGrid_m.back() = spatialUpperBound_m;
		currentSolution_m.back() = process_m.Payoff(spatialGrid_m.back(), initialTime_m);

		//initializes non-boundaries values of spatialGrid_m and currentSolution_m
		for (int r = 1; r <= totalSpatialPoints_m; r++)
		{
			spatialGrid_m[r] = spatialLowerBound_m + r*deltaX_m;
			currentSolution_m[r] = process_m.Payoff(spatialGrid_m[r], initialTime_m);
		}
	}

	virtual ~MultiThreadPDE_Explicit() {}

	//propagates solution from time i+1 to time i
	//using the kolmogorov backward equation
	//the propagation is using order one on explicit scheme
	//which subdivides the spatial nodes in different groups (depending on hardware concurrency
	//the tasks are allocated on working threads by ThreadPool
	//confirmation is a dummy variable used foor synchronization of each iteration
	//multithreading modifies the currentSolution_m[1,totalSpatialPoints_m] elements
	//the update of 0 and totalSpatialPoints_m of currenSolution_m (boundaries) is performed
	//separately by the main thread
	void propagate()
	{
		int num_threads = std::thread::hardware_concurrency();
		int blockSize = int(totalSpatialPoints_m / num_threads);
		ThreadPool pool(num_threads);
		for (int n = 0; n<(totalTimePoints_m - 1); n++)
		{
			std::vector<std::future<bool> >confirmation;
			double currentTime = initialTime_m + n*deltaT_m;
			previousSolution_m.swap(currentSolution_m);
			int block_start = 1;
			for (int i = 0; i < int(num_threads-1); ++i)
			{
				int block_end = block_start;
				block_end += blockSize;
				block_end = std::min(block_end, totalSpatialPoints_m + 1);
				confirmation.push_back(pool.enqueue(std::bind(&MultiThreadPDE_Explicit<T>::propagateOneTimeStep, this, block_start, block_end, currentTime)));
				block_start = block_end + 1;
			}
			int block_end = block_start;
			block_end += blockSize;
			block_end = std::min(block_end, totalSpatialPoints_m);
			confirmation.push_back(pool.enqueue(std::bind(&MultiThreadPDE_Explicit<T>::propagateOneTimeStep, this, block_start, block_end, currentTime)));
			
			//wait for all threads to end iteration
			for (int k = 0; k < num_threads; ++k) confirmation[k].get();
			
			//updates lower and upper boundary of currentSolution_m
			currentSolution_m[0] = 0.0;
			currentSolution_m[totalSpatialPoints_m + 1] = process_m.Payoff(spatialGrid_m.back(), currentTime);
		}
	}

	//gets the Present Value of the payoff specified
	const double getPV()
	{
		tk::spline spl;
		spl.set_points(spatialGrid_m, currentSolution_m);
		double interp = spl.operator()(std::log(process_m.getSpot()));
		return interp;
	}
	//Stability is based on calculation of infinite Norm of transition matrix
	//for the particular case of Black Scholes with constant parameters,
	//this matrix is independent of Underlying and time
	//the function below is true only for this special case
	// x = 0.0 (underlying at spot, and time t=1.0  are completely arbitrary parameters
	//that allows quick determination of infinite norm
	bool isStable()
	{
		double x = 0;
		double t = 1.0;
		double advection = process_m.PdeAdvectionConvection(x, t) * 0.5 / deltaX_m;
		double diffusion = process_m.PdeDifusion(x, t) / deltaX_m / deltaX_m;
		double shortRate = process_m.PdeSourceSink(x, t); 
		double L = deltaT_m*(diffusion - advection);
		double D = (1.0 - deltaT_m*(2.0*diffusion - shortRate));
		double U = deltaT_m*(diffusion + advection);
		double infiniteNorm = std::abs(L) + std::abs(D) + std::abs(U);
		return (infiniteNorm < 1.0 ? true : false);
	}

	MultiThreadPDE_Explicit(const MultiThreadPDE_Explicit&) = delete;
	void operator=(const MultiThreadPDE_Explicit&) = delete;

private:

	//mission is to update currentSolution_m
	//minimum start is ONE
	//maximum end is totalSpatialPoints_m
	bool propagateOneTimeStep(const int& begin, const int& end, const double& t)
	{
		std::thread::id id = std::this_thread::get_id();
		// minimum begin needs to be 1
		if (begin < 1)
			throw std::invalid_argument("propagateOneTimeStep:: received index array less than 1");
		if (end > totalSpatialPoints_m)
			throw std::invalid_argument("propagateOneTimeStep:: received index array greater than than totalSpatialPoints_m");
		// maximum end needs to be totalSpatialPoints_m
		double tempInfiniteNorm = -1E+12;
		for (int r = begin; r <= end; r++)
		{
			double x = spatialGrid_m[r];
			//advection, diffusion and shortRate are constant
			//on Black Scholes model with constant parameters
			//the unoptimized generalized calculation below
			//is left for a future extension to non constant parameters
			double advection = process_m.PdeAdvectionConvection(x, t) * 0.5/ deltaX_m;		
			double diffusion = process_m.PdeDifusion(x, t)/ deltaX_m / deltaX_m;					
			double shortRate = process_m.PdeSourceSink(x, t);	

			double L = deltaT_m*(diffusion - advection);
			double D = (1.0 - deltaT_m*(2.0*diffusion - shortRate));
			double U = deltaT_m*(diffusion + advection);

			currentSolution_m[r] = previousSolution_m[r - 1] * L;
			currentSolution_m[r] +=    previousSolution_m[r] * D;
			currentSolution_m[r] +=previousSolution_m[r + 1] * U;
		}
		return true;
	}

	double initialTime_m;
	double finalTime_m;
	int totalTimePoints_m;
	double spatialLowerBound_m;
	double spatialUpperBound_m;
	int totalSpatialPoints_m;
	double deltaX_m, deltaT_m;

	T& process_m;
	std::vector<double> spatialGrid_m;              
	std::vector<double> currentSolution_m;         
	std::vector<double> previousSolution_m;         

};


#endif

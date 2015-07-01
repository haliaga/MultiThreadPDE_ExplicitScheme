# MultiThreadPDE_ExplicitScheme
Multithreaded 1D Black Scholes PDE Solver C++11
* Author: Horacio Aliaga	
* Copyright (C) 2015 Horacio Aliaga (horacio.aliaga at gmail.com)
* Description:  Library providing a multithreading Black Scholes PDE Solver
*               using the explicit method of order one.
*               The model has constant parameters. Multithreading is achieved 
*               by means of C++11 functionality
*				        The algorithm in unconditionally stable provided the infinite norm
*				        of the heat equation matrix is lower than 1
* 

Dependencies: this library is using:
*               ThreadPool, that can be found on:
*               https://github.com/progschj/ThreadPool
*
*				        spline.h
* Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)

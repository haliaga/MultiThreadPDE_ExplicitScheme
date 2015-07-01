
#include "stdafx.h"
#include "SpotProcess.h"
#include "MultiThreadPDE_Explicit.h"

using namespace std;
using namespace tk;

int main()
{
	SpotProcess sp;
	double maxS = std::log(400);
	double minS = std::log(0.1);
	MultiThreadPDE_Explicit<SpotProcess> mtpde(sp, 100, 0.0, 1.0, 50, minS, maxS);
	double pv = 0.0;
	if (mtpde.isStable())
	{
		mtpde.propagate();
		pv = mtpde.getPV();
	}
	return 0;
}

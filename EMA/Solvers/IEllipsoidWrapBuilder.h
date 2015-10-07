#pragma once
#include <vector>
#include "../Drawing.h"
#include "IEllipsoidWrapper.h"

class IEllipsoidWrapBuilder
{
public:
	virtual IEllipsoidWrapper* GetEllipseWrap(Ellipsoid* el) = 0;
	virtual IEllipsoidWrapper* GetEllipseWrap(const vector<Point> & points, const vector<int> & selectedIndexes) = 0;
	virtual double CalcFitness(const vector<Point>& points, IEllipsoidWrapper* ellipsoid) = 0;
};
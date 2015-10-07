#pragma once
#include "IEllipsoidBuilder.h"
#include "IEllipsoidWrapBuilder.h"
#include "Support/Petunin.h"

class PetuninBuilder : public IEllipsoidBuilder, public IEllipsoidWrapBuilder
{
public:
	Ellipsoid* Exec(const vector<Point>& points, Window* window);
	Ellipsoid2D* Exec(const vector<Point2D>& points, Window* window);
	IEllipsoidWrapper* GetEllipseWrap(Ellipsoid* el);
	IEllipsoidWrapper* GetEllipseWrap(const vector<Point>& points, const vector<int>& selectedIndexes);
	double CalcFitness(const vector<Point>& points, IEllipsoidWrapper* ellipsoid);
};
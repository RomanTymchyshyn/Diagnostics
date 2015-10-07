#pragma once
#include <vector>
#include "../GeomFigures.h"

class IEllipsoidWrapper: Cloneable<IEllipsoidWrapper>
{
public:
	bool CheckPercents;
	double PercentPointsIn(const vector<Point>& points)
	{
		int count = 0;
		int size = points.size();
		for (int i = 0; i < size; ++i)
			count += _ellipsoid->Inside(points[i]);
		return count * 1.0 / size;
	}
	virtual double Volume() { return M_PI * _ellipsoid->Axes()[0] * _ellipsoid->Axes()[1]; };
	virtual double CalcFitness(const vector<Point> & points) = 0;
	IEllipsoidWrapper* Clone() = 0;
	double fitness;
	double probability;
	Ellipsoid* _ellipsoid;
};
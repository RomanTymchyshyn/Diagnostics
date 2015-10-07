#pragma once
#include "ConvHull.h"
#include "../IEllipsoidWrapper.h"

//added as temporary solution
///TODO: refactor this (do better architecture, such that all algorithms working with same data structures)
class PetuninEllipsoid :public IEllipsoidWrapper
{
private:
	Ellipsoid2D* ellipsoid;
public:
	PetuninEllipsoid()
	{
		this->_ellipsoid = NULL;
		this->ellipsoid = NULL;
		this->CheckPercents = true;
	}
	PetuninEllipsoid(Ellipsoid* ellipsoid)
	{
		this->_ellipsoid = ellipsoid;
		this->ellipsoid = new Ellipsoid2D();
		*this->ellipsoid = *ellipsoid;
		this->CheckPercents = true;
	}
	PetuninEllipsoid(Ellipsoid2D* ellipsoid)
	{
		this->_ellipsoid = new Ellipsoid(2);
		this->ellipsoid = ellipsoid;
		*this->_ellipsoid = *ellipsoid;
	}
	double CalcFitness(const vector<Point> & points);
	IEllipsoidWrapper* Clone();
	~PetuninEllipsoid();
};
Ellipsoid2D* PetuninAlgo(const vector<Point2D> & points, Window* window = NULL);
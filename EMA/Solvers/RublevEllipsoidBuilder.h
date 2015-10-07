#pragma once
#include "Support/RublevSup.h"
#include "IEllipsoidBuilder.h"

class RublevEllipsoidBuilder :public IEllipsoidBuilder
{
public:
	Ellipsoid* Exec(const vector<Point>& points, Window* window);
	Ellipsoid2D* Exec(const vector<Point2D> & points, Window* window = NULL);
};
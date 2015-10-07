#pragma once
#include "Support/KhachiyanSup.h"
#include "IEllipsoidBuilder.h"

class KhachiyanEllipsoidBuilder : public IEllipsoidBuilder
{
public:
	Ellipsoid* Exec(const vector<Point>& points, Window* window);
	Ellipsoid2D* Exec(const vector<Point2D> & points, Window* window = NULL);
};
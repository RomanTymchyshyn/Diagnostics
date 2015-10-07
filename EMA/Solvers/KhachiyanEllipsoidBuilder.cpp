#include "KhachiyanEllipsoidBuilder.h"

Ellipsoid* KhachiyanEllipsoidBuilder::Exec(const vector<Point>& points, Window* window)
{
	return KhachiyanAlgo(points, 0.01);
}

Ellipsoid2D* KhachiyanEllipsoidBuilder::Exec(const vector<Point2D> & points, Window* window)
{
	vector<Point> vertexes(points.size(), Point());
	for (int i = 0; i < points.size(); ++i)
		vertexes[i] = points.at(i);
	Ellipsoid* el = KhachiyanAlgo(vertexes, 0.01);
	Ellipsoid2D* res = new Ellipsoid2D();
	(*res) = (*el);
	delete el;
	el = NULL;
	return res;
}
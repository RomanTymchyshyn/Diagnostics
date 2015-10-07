#include "RublevEllipsoidBuilder.h"

Ellipsoid* RublevEllipsoidBuilder::Exec(const vector<Point>& points, Window* window)
{
	vector<Point2D> vertexes(points.size(), Point2D());
	for (int i = 0; i < points.size(); ++i)
		vertexes[i] = Point2D(points.at(i).Coord(0), points.at(i).Coord(1));
	Ellipsoid2D* el = RublevAlg(vertexes, window);
	vertexes.clear();
	Point centre;
	centre = el->Centre();
	Ellipsoid* res = new Ellipsoid(centre, el->Eigenvectors(), el->Axes());

	delete el; el = NULL;
	return res;
}

Ellipsoid2D* RublevEllipsoidBuilder::Exec(const vector<Point2D> & points, Window* window)
{
	return RublevAlg(points, window);
}
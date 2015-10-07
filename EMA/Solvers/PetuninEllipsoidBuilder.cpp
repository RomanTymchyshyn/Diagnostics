#include "PetuninEllipsoidBuilder.h"

IEllipsoidWrapper* PetuninBuilder::GetEllipseWrap(Ellipsoid* el)
{
	return new PetuninEllipsoid(el);
}

IEllipsoidWrapper* PetuninBuilder::GetEllipseWrap(const vector<Point>& points, const vector<int>& selectedIndexes)
{
	int pointsSize = points.size();
	int selectedSize = selectedIndexes.size();
	int dim = 2;

	vector<Point2D> vertexes;
	for (int i = 0; i < selectedSize; ++i)
	{
		Point2D p(points[selectedIndexes[i]].Coord(0), points[selectedIndexes[i]].Coord(1));
		vertexes.push_back(p);
	}

	PetuninEllipsoid* ellipsoid = new PetuninEllipsoid(PetuninAlgo(vertexes));

	ellipsoid->fitness = ellipsoid->CalcFitness(points);
	return ellipsoid;
}

double PetuninBuilder::CalcFitness(const vector<Point>& points, IEllipsoidWrapper* ellipsoid)
{
	PetuninEllipsoid* el = new PetuninEllipsoid(ellipsoid->_ellipsoid);
	double fitness = el->CalcFitness(points);
	delete el; el = NULL;
	return fitness;
}

Ellipsoid* PetuninBuilder::Exec(const vector<Point>& points, Window* window)
{
	vector<Point2D> vertexes(points.size(), Point2D());
	for (int i = 0; i < points.size(); ++i)
		vertexes[i] = Point2D(points.at(i).Coord(0), points.at(i).Coord(1));
	Ellipsoid2D* el = PetuninAlgo(vertexes, window);
	vertexes.clear();
	Point centre;
	centre = el->Centre();
	Ellipsoid* res = new Ellipsoid(centre, el->Eigenvectors(), el->Axes());

	delete el; el = NULL;
	return res;
}

Ellipsoid2D* PetuninBuilder::Exec(const vector<Point2D>& points, Window* window)
{
	return PetuninAlgo(points, window);
}
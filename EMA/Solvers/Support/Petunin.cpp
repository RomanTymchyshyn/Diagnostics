#include "Petunin.h"

vector<Point2D> ConvertAnyTo2D(const vector<Point> & points)
{
	int n = points.size();
	vector<Point2D> res;
	for (int i = 0; i < n; ++i)
	{
		Point2D p = new Point2D(points[i].Coord(0), points[i].Coord(1));
		res.push_back(p);
	}
	return res;
}

void FindFarthestPair(const vector<Point2D> & points, Point2D & diag1, Point2D & diag2, double & distance)
{
	int n = points.size();
	distance = -1;
	for (int i = 0; i < n; ++i)
	{
		Point2D p1 = points[i];
		for (int j = i + 1; j < n; ++j)
		{
			Point2D p2 = points[j];
			double dist = distancePToP(p1, p2);
			if (dist > distance)
			{
				diag1 = p1;
				diag2 = p2;
				distance = dist;
			}
		}
	}
	return;
}

void FindFarthestFromLine(Line2D l, const vector<Point2D> & points, Point2D & top, Point2D & bottom, double & maxTop, double & maxBot)
{
	int n = points.size();
	maxTop = -1, maxBot = -1;
	for (int i = 0; i < n; ++i)
	{
		Point2D p = points[i];
		double dist = l.Distance(p);
		int deviation = l.Deviation(p);
		if (deviation == -1)
		{
			if (dist > maxTop)
			{
				top = p;
				maxTop = dist;
			}
		}
		else if (deviation == 1)
		{
			if (dist > maxBot)
			{
				bottom = p;
				maxBot = dist;
			}
		}
	}
	return;
}

vector<Point2D> RotatePoints(const vector<Point2D> & points, const double & angle)
{
	vector<Point2D> res;
	for (Point2D p : points)
		res.push_back(rotate(p, angle));
	return res;
}

vector<Point2D> TranslatePoints(const vector<Point2D> & points, const double & xOffset, const double & yOffset)
{
	vector<Point2D> res;
	for (Point2D p : points)
		res.push_back(p.Translation(xOffset, yOffset));
	return res;
}

vector<Point2D> ScalePoints(const vector<Point2D> & points, const double & byX, const double & byY)
{
	vector<Point2D> res;
	for (Point2D p : points)
		res.push_back(p.Scale(byX, byY));
	return res;
}

real_2d_array GetEigenVectors()
{
	real_2d_array eVectors;
	eVectors.setlength(2, 2);
	eVectors[0][0] = 1;
	eVectors[0][1] = 0;
	eVectors[1][0] = 0;
	eVectors[1][1] = 1;
	return eVectors;
}

double GetRadius(const vector<Point2D> & points, const Point2D & centre)
{
	double max = -1;
	for (Point2D p : points)
	{
		double dist = distancePToP(p, centre);
		if (dist > max) max = dist;
	}
	return max;
}

Rectangle2D GetRectangle(const vector<Point2D> & points, Window* window = NULL)
{
	Point2D left, right;
	double distance = 0.0;
	FindFarthestPair(points, left, right, distance);

	DrawEvent* draw = NULL;
	if (window != NULL)
	{
		draw = new DrawEvent();
		left.SetColor(Color(1.0, 0.0, 0.0));
		left.SetSize(8.0);
		right.SetColor(Color(1.0, 0.0, 0.0));
		right.SetSize(8.0);
		draw->add(2, &left, &right);
		window->AddToLastScene(draw);
	}

	Line2D L(left, right);
	Point2D top(left), bottom(right);
	double maxTopDist, maxBotDist;
	FindFarthestFromLine(L, points, top, bottom, maxTopDist, maxBotDist);
	if (top == left)
		top = left.Y() > right.Y() ? left : right;
	if (bottom == right)
		bottom = left.Y() < right.Y() ? left : right;

	if (window != NULL)
	{
		draw->clear();
		top.SetColor(Color(1.0, 0.0, 0.0));
		top.SetSize(8.0);
		bottom.SetColor(Color(1.0, 0.0, 0.0));
		bottom.SetSize(8.0);
		draw->add(3, &top, &bottom, &L);
		window->AddToLastScene(draw);
	}


 	Line2D L1 = L.GetParallel(top);
	L1.SetColor(Color(1.0, 0.0, 0.0));
	Line2D L2 = L.GetParallel(bottom);
	L2.SetColor(Color(0.0, 1.0, 0.0));
	Line2D L3 = L.GetPerpendicular(left);
	L3.SetColor(Color(0.0, 0.0, 1.0));
	Line2D L4 = L.GetPerpendicular(right);
	L4.SetColor(Color(1.0, 1.0, 1.0));

	if (window != NULL)
	{
		draw->clear();
		draw->add(4, &L1, &L2, &L3, &L4);
		window->AddToLastScene(draw);
	}

	Point2D p1 = L1.Intersection(L3);
	Point2D p2 = L1.Intersection(L4);
	Point2D p3 = L2.Intersection(L4);
	Point2D p4 = L2.Intersection(L3);

	if (window != NULL)
	{
		delete draw;
		draw = NULL;
	}

	return Rectangle2D(p1, p2, p3, p4);
}

void GetParameters(Rectangle2D & rect, double* offset, double* rectParams, double & angle)
{
	int leftBottomIndex = rect.LeftBottomIndex();
	offset[0] = -rect.Vertex(leftBottomIndex).X();
	offset[1] = -rect.Vertex(leftBottomIndex).Y();

	int neighbour1 = leftBottomIndex == 0 ? 3 : leftBottomIndex - 1;
	int neighbour2 = leftBottomIndex == 3 ? 0 : leftBottomIndex + 1;
	int secondPoint = rect.Vertex(neighbour1).X() > rect.Vertex(neighbour2).X() ? neighbour1 : neighbour2;
	vector<double> v1 = { rect.Vertex(secondPoint).X() - rect.Vertex(leftBottomIndex).X(), 
		rect.Vertex(secondPoint).Y() - rect.Vertex(leftBottomIndex).Y() };
	angle = acos(v1[0] / sqrt(v1[0] * v1[0] + v1[1] * v1[1]));
	if (v1[1] < 0) angle *= -1.0;

	rectParams[0] = distancePToP(rect.Vertex(secondPoint), rect.Vertex(leftBottomIndex));
	rectParams[1] = distancePToP(rect.Vertex(secondPoint == neighbour1 ? neighbour2 : neighbour1), 
		rect.Vertex(leftBottomIndex));
}

Ellipsoid2D* PetuninAlgo(const vector<Point2D> & points, Window* window)
{
	////Drawing
	if (window != NULL)
		window->clear();
	////Drawing

	//Convex hull
	Point2D left, right, up, bot;
	vector<Point2D>* conv = ConvHull(points, left, right, up, bot, NULL);
	
	//Getting Rectangle
	Rectangle2D rect = GetRectangle(*conv, window);

	DrawEvent* draw = NULL;
	////Drawing
	if (window != NULL)
	{
		draw = new DrawEvent();
		draw->add(1, &rect);
		window->AddNewDrawEvent(draw);
	}
	////Drawing

	double offset[2];
	double rectParams[2];
	double angle;
	GetParameters(rect, offset, rectParams, angle);

	//Translating to origin
	
	Rectangle2D modified = rect.Translation(offset[0], offset[1]);
	vector<Point2D> modifiedConv = TranslatePoints(*conv, offset[0], offset[1]);

	////Drawing
	Color green;
	if (window != NULL)
	{
		green = Color(0.0, 1.0, 0.0);
		draw->clear();
		draw->add(&modifiedConv, green);
		draw->add(1, &modified);
		window->AddNewDrawEvent(draw);
	}
	////Drawing

	modified = modified.Rotate(angle);
	modifiedConv = RotatePoints(modifiedConv, angle);

	////Drawing
	if (window != NULL)
	{
		draw->clear();
		draw->add(&modifiedConv, green);
		draw->add(1, &modified);
		window->AddNewDrawEvent(draw);
	}
	////Drawing
	
	double alpha = rectParams[1] / rectParams[0];

	modified = modified.Scale(alpha, 1.0);
	modifiedConv = ScalePoints(modifiedConv, alpha, 1.0);

	////Drawing
	if (window != NULL)
	{
		draw->clear();
		draw->add(&modifiedConv, green);
		draw->add(1, &modified);
		window->AddNewDrawEvent(draw);
	}
	////Drawing

	Point2D centre(rectParams[1] / 2.0, rectParams[1] / 2.0);
	double R = GetRadius(modifiedConv, centre);

	real_2d_array eVectors = GetEigenVectors();
	real_1d_array axes;
	axes.setlength(2);
	axes[0] = axes[1] = R;
	Ellipsoid2D* petEl = new Ellipsoid2D(centre, eVectors, axes);

	////Drawing
	if (window != NULL)
	{
		draw->clear();
		draw->add(1, petEl);
		window->AddToLastScene(draw);
	}
	////Drawing

	(*petEl) = petEl->Scale(1.0 / alpha, 1.0);
	modified = modified.Scale(1.0 / alpha, 1.0);
	modifiedConv = ScalePoints(modifiedConv, 1.0 / alpha, 1.0);

	////Drawing
	if (window != NULL)
	{
		draw->clear();
		draw->add(2, petEl, &modified);
		draw->add(&modifiedConv, green);
		window->AddNewDrawEvent(draw);
	}
	////Drawing

	(*petEl) = petEl->Rotate(-angle);
	modified = modified.Rotate(-angle);
	modifiedConv = RotatePoints(modifiedConv, -angle);

	////Drawing
	if (window != NULL)
	{
		draw->clear();
		draw->add(2, petEl, &modified);
		draw->add(&modifiedConv, green);
		window->AddNewDrawEvent(draw);
	}
	////Drawing

	(*petEl) = petEl->Translation(-offset[0], -offset[1]);
	modified = modified.Translation(-offset[0], -offset[1]);
	modifiedConv = TranslatePoints(modifiedConv, -offset[0], -offset[1]);

	////Drawing
	if (window != NULL)
	{
		draw->clear();
		draw->add(2, petEl, &modified);
		draw->add(&modifiedConv, green);
		window->AddNewDrawEvent(draw);
		delete draw;
		draw = NULL;
	}
	////Drawing

	delete conv;
	conv = NULL;
	return petEl;
}

double PetuninEllipsoid::CalcFitness(const vector<Point>& points)
{
	double fitness = 0.0;
	double volume = ellipsoid->Area();
	double percentOfPoints = PercentPointsIn(points);
	if (percentOfPoints >= 0.95)
		fitness = (1 / (volume + 1))*percentOfPoints;
	else fitness = (1 / (volume + 1))*(percentOfPoints > 0.5 ? (1.0 / 1000) : 0);
	return fitness;
}

IEllipsoidWrapper* PetuninEllipsoid::Clone()
{
	IEllipsoidWrapper* copy = new PetuninEllipsoid(_ellipsoid);
	copy->fitness = fitness;
	copy->probability = probability;
	return copy;
}

PetuninEllipsoid::~PetuninEllipsoid()
{
	delete ellipsoid;
	ellipsoid = NULL;
	delete _ellipsoid;
	_ellipsoid = NULL;
}
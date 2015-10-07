#pragma once
#include "IEllipsoidBuilder.h"
#include "IEllipsoidWrapBuilder.h"
#include "Support/MatrixAlgo.h"

class MatrixBuilder : public IEllipsoidBuilder, public IEllipsoidWrapBuilder
{
	alglib::real_2d_array GetPointSet(const vector<Point> & points, const vector<int> & selectedIndexes);
	alglib::real_1d_array Mean(const alglib::real_2d_array & data, const int & num_row, const int & num_col);
	Ellipsoid* GetEllipse(const alglib::real_1d_array & mean, const alglib::real_2d_array & covmatrix, const int & dim);
public:
	IEllipsoidWrapper* GetEllipseWrap(Ellipsoid* el);
	IEllipsoidWrapper* GetEllipseWrap(const vector<Point>& points, const vector<int>& selectedIndexes);
	double CalcFitness(const vector<Point>& points, IEllipsoidWrapper* ellipsoid);
	Ellipsoid2D* Exec(const vector<Point2D>& points, Window* window = NULL);
	Ellipsoid* Exec(const vector<Point>& points, Window* window);
};
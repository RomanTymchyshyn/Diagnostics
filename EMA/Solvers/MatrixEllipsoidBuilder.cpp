#include "MatrixEllipsoidBuilder.h"

IEllipsoidWrapper* MatrixBuilder::GetEllipseWrap(Ellipsoid* el)
{
	MatrixEllipsoid* ellipsoid = new MatrixEllipsoid(el);
	return ellipsoid;
}

IEllipsoidWrapper* MatrixBuilder::GetEllipseWrap(const vector<Point>& points, const vector<int>& selectedIndexes)
{
	int pointsSize = points.size();
	int selectedSize = selectedIndexes.size();
	int dim = points[0].Dim();

	//get set of selected points
	alglib::real_2d_array pointSet = GetPointSet(points, selectedIndexes);

	//calculation of mean-vector
	alglib::real_1d_array mean = Mean(pointSet, selectedSize, dim);

	//calculating of covariance matrix
	alglib::real_2d_array covmatrix;
	alglib::covm(pointSet, covmatrix);

	MatrixEllipsoid* ellipsoid = new MatrixEllipsoid(GetEllipse(mean, covmatrix, dim));

	ellipsoid->fitness = ellipsoid->CalcFitness(points);
	return ellipsoid;
}

alglib::real_2d_array MatrixBuilder::GetPointSet(const vector<Point> & points, const vector<int> & selectedIndexes)
{
	int dim = points[0].Dim();
	int N = points.size();
	int k = selectedIndexes.size();
	alglib::real_2d_array set;
	set.setlength(k, dim);
	for (int j = 0; j < k; ++j)
	{
		for (int l = 0; l < dim; ++l)
		{
			set(j, l) = points[selectedIndexes[j]].Coord(l);
		}
	}
	return set;
}

alglib::real_1d_array MatrixBuilder::Mean(const alglib::real_2d_array & data, const int & num_row, const int & num_col)
{
	alglib::real_1d_array mean;
	mean.setlength(num_col);

	for (int i = 0; i < num_col; ++i)
	{
		mean[i] = 0.0;
		for (int j = 0; j < num_row; ++j)
			mean[i] += data[j][i];
		mean[i] /= num_row;
	}

	return mean;
}

Ellipsoid* MatrixBuilder::GetEllipse(const alglib::real_1d_array & mean, const alglib::real_2d_array & covmatrix, const int & dim)
{
	//getting shape matrix of the ellipse
	alglib::real_2d_array shapeMatrix(covmatrix);
	alglib::ae_int_t info;
	alglib::matinvreport rep;

	alglib::rmatrixinverse(shapeMatrix, info, rep);

	alglib::real_2d_array eVectors;
	eVectors.setlength(dim, dim);
	alglib::real_1d_array eValues;
	eValues.setlength(dim);

	alglib::smatrixevd(shapeMatrix, dim, 1, false, eValues, eVectors);

	for (int i = 0; i < dim; ++i)
		eValues[i] = 1 / sqrt(eValues[i]);
	Ellipsoid* el = new Ellipsoid(Point(mean), eVectors, eValues);
	return el;
}

double MatrixBuilder::CalcFitness(const vector<Point>& points, IEllipsoidWrapper* ellipsoid)
{
	alglib::real_2d_array CovMatrix;
	alglib::real_1d_array mean;
	MatrixEllipsoid* el = new MatrixEllipsoid(ellipsoid->_ellipsoid);
	double fitness = el->CalcFitness(points);
	delete el;
	el = NULL;
	return fitness;
}

Ellipsoid2D* MatrixBuilder::Exec(const vector<Point2D>& points, Window* window)
{
	vector<Point> vertexes(points.size(), Point());
	for (int i = 0; i < points.size(); ++i)
		vertexes[i] = points.at(i);
	Ellipsoid* el = MatrixAlgo(vertexes);
	Ellipsoid2D* res = new Ellipsoid2D();
	(*res) = (*el);
	delete el;
	el = NULL;
	return res;
}

Ellipsoid* MatrixBuilder::Exec(const vector<Point>& points, Window* window)
{
	return MatrixAlgo(points);
}
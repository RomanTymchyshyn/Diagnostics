#include "MatrixAlgo.h"

double MatrixEllipsoid::Volume(const alglib::real_2d_array & covMatrix, const alglib::real_1d_array & mean, const vector<Point> & points)
{
	int n = points.size();
	int dim = points[0].Dim();
	int h = (n + dim + 1.0) / 2.0;
	double* MahDist = new double[n];
	alglib::real_2d_array vector;
	vector.setlength(1, dim);

	alglib::real_2d_array shapeMatrix(covMatrix);
	alglib::ae_int_t info;
	alglib::matinvreport rep;

	double det = alglib::rmatrixdet(covMatrix);
	det = sqrt(det);
	alglib::rmatrixinverse(shapeMatrix, info, rep);

	for (int i = 0; i < n; ++i)
	{
		alglib::real_2d_array res1;
		alglib::real_2d_array res2;
		res1.setlength(1, dim);
		res2.setlength(1, 1);
		for (int j = 0; j < dim; ++j)
			vector(0, j) = points[i].Coord(j) - mean(j);
		alglib::rmatrixgemm(1, dim, dim, 1, vector, 0, 0, 0, shapeMatrix, 0, 0, 0, 0, res1, 0, 0);
		alglib::rmatrixgemm(1, 1, dim, 1, res1, 0, 0, 0, vector, 0, 0, 1, 0, res2, 0, 0);
		MahDist[i] = res2(0, 0);
	}
	sort(MahDist, MahDist + n);
	double volume = det * pow(MahDist[h], dim);
	volume *= (volume < 0 ? -1 : 1);
	return volume;
}

double MatrixEllipsoid::CalcFitness(const vector<Point> & points)
{
	alglib::real_2d_array CovMatrix;
	alglib::real_1d_array mean;
	GetCovAndMean(CovMatrix, mean);
	double fitness = 0.0;
	double volume = Volume(CovMatrix, mean, points);
	double percentOfPoints = PercentPointsIn(points);
	if (percentOfPoints >= 0.95)
		fitness = (1 / (volume + 1))*percentOfPoints;
	else fitness = (1 / (volume + 1))*(1.0 / 10);
	return fitness;
}

void MatrixEllipsoid::GetCovAndMean(alglib::real_2d_array & covMatrix, alglib::real_1d_array & mean)
{
	covMatrix.setlength(_ellipsoid->Dim(), _ellipsoid->Dim());
	mean = _ellipsoid->Centre().Coords();

	alglib::real_2d_array eVecs = _ellipsoid->Eigenvectors();
	alglib::real_2d_array eVals;
	eVals.setlength(_ellipsoid->Dim(), _ellipsoid->Dim());

	for (int i = 0; i < _ellipsoid->Dim(); ++i)
	{
		for (int j = 0; j < _ellipsoid->Dim(); ++j)
		{
			if (i == j)
			{
				eVals(i, j) = 1 / (_ellipsoid->Axes()[i] * _ellipsoid->Axes()[i]);
			}
			else eVals(i, j) = 0.0;
		}
	}

	alglib::real_2d_array shape;
	shape.setlength(_ellipsoid->Dim(), _ellipsoid->Dim());
	alglib::real_2d_array res;
	res.setlength(_ellipsoid->Dim(), _ellipsoid->Dim());
	rmatrixgemm(_ellipsoid->Dim(), _ellipsoid->Dim(), _ellipsoid->Dim(), 1, eVecs, 0, 0, 0, eVals, 0, 0, 0, 0, res, 0, 0);
	rmatrixgemm(_ellipsoid->Dim(), _ellipsoid->Dim(), _ellipsoid->Dim(), 1, res, 0, 0, 0, eVecs, 0, 0, 1, 0, covMatrix, 0, 0);

	alglib::ae_int_t info;
	alglib::matinvreport rep;

	alglib::rmatrixinverse(covMatrix, info, rep);
}

IEllipsoidWrapper* MatrixEllipsoid::Clone()
{
	IEllipsoidWrapper* copy = new MatrixEllipsoid(_ellipsoid);
	copy->fitness = fitness;
	copy->probability = probability;
	return copy;
}

MatrixEllipsoid::~MatrixEllipsoid()
{
	delete _ellipsoid;
	_ellipsoid = NULL;
}

alglib::real_2d_array GetPointSet(const vector<Point> & points, const vector<int> & selectedIndexes)
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

alglib::real_1d_array Mean(const alglib::real_2d_array & data, const int & num_row, const int & num_col)
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

Ellipsoid* GetEllipse(const alglib::real_1d_array & mean, const alglib::real_2d_array & covmatrix, const int & dim)
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

MatrixEllipsoid* GetEllipseWrap(const vector<Point> & points, const vector<int> & selectedIndexes)
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

Ellipsoid* MatrixAlgo(const vector<Point>& points)
{
	int n = points.size();
	vector<int> selectedIndexes;
	for (int i = 0; i < n; ++i)
		selectedIndexes.push_back(i);
	MatrixEllipsoid* elWrap = GetEllipseWrap(points, selectedIndexes);
	Ellipsoid* el = new Ellipsoid(*elWrap->_ellipsoid);
	delete elWrap;
	return el;
}
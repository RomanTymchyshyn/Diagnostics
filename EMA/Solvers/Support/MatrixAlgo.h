#pragma once
#include "../../GeomFigures.h"
#include "../../Vendor/cpp/src/linalg.h"
#include "../../Vendor/cpp/src/statistics.h"
#include <time.h>
#include "../../Drawing.h"
#include "../IEllipsoidWrapper.h"

//added as temporary solution
///TODO: refactor this (do better architecture, such that all algorithms working with same data structures)
class MatrixEllipsoid:public IEllipsoidWrapper
{
public:
		MatrixEllipsoid()
		{
			this->_ellipsoid = NULL;
			this->CheckPercents = false;
		}
		MatrixEllipsoid(Ellipsoid* ellipsoid)
		{
			this->_ellipsoid = ellipsoid;
			this->CheckPercents = false;
		}
		double CalcFitness(const vector<Point> & points);
		double Volume(const alglib::real_2d_array & CovMatrix, const alglib::real_1d_array & mean, const vector<Point> & points);
		void GetCovAndMean(alglib::real_2d_array & CovMatrix, alglib::real_1d_array & mean);
		IEllipsoidWrapper* Clone();
		~MatrixEllipsoid();
};

Ellipsoid* MatrixAlgo(const vector<Point>& points);
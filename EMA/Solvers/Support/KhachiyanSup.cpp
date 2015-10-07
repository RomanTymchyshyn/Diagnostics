#include "KhachiyanSup.h"


alglib::real_2d_array shape(const alglib::real_1d_array & p, const alglib::real_2d_array & data)
{
	alglib::real_2d_array res;
	int n = data.rows();
	res.setlength(n,n);
	alglib::rmatrixgemm(n,n,1,p[0],data,0,0,0,data,0,0,1,0,res,0,0);
	for(int i = 1; i < data.cols(); ++i)
		alglib::rmatrixgemm(n,n,1,p[i],data,0,i,0,data,0,i,1,1,res,0,0);
	return res;
}

alglib::real_1d_array findCentre(const alglib::real_2d_array& data,const alglib::real_1d_array& P)
{
	alglib::real_1d_array res;
	res.setlength(data.rows());
	for (int i = 0; i < data.rows(); i++) {
		res[i] = P[0]*data(i,0);
	}
	for (int j = 1; j < data.cols(); j++) {
		for (int i = 0; i < data.rows(); i++) {
			res[i] += P[j]*data(i,j);
		}
	}
	return res;
}

alglib::real_1d_array computeDecentered(const alglib::real_2d_array& data, const double& eps)
{
	alglib::real_1d_array P;
	int n = data.rows(), m = data.cols();
	P.setlength(m);
	for(int i = 0; i < m ; ++i)
		P[i] = 1./m;
	alglib::real_2d_array inv;
	alglib::ae_int_t info;
	alglib::matinvreport rep;

	double K = n+1;
	while(K > (1+eps)*n || K < (1-eps)*n)
	{
		K = 0;
		alglib::real_2d_array tempK;
		alglib::real_2d_array temp;
		tempK.setlength(1,1);
		temp.setlength(1,n);
		int j = 0;
		inv.setlength(n,n);
		for (int i = 0; i < m; i++) {
			inv = shape(P,data);
			alglib::rmatrixinverse(inv,info,rep);
			alglib::rmatrixgemm(1,n,n,1,data,0,i,1,inv,0,0,0,0,temp,0,0);
			alglib::rmatrixgemm(1,1,n,1,temp,0,0,0,data,0,i,0,0,tempK,0,0);
			if(tempK(0,0) > K)
			{
				K = tempK(0,0);
				j = i;
			}
		}
		double Beta = (K - n)/(n*(K-1));
		for (int i = 0; i < m; i++) {
			P[i] = (1 - Beta)*P[i];
			if(i==j)
				P[i] += Beta;
		}
	}
	return P;
}

alglib::real_2d_array computeMatrixForm(const alglib::real_2d_array& data, const double& eps, alglib::real_1d_array & centre)
{
	alglib::real_2d_array dataNew;
	dataNew.setlength(data.rows() + 1, data.cols());
	for (int j = 0; j < data.cols(); j++) {
		for (int i = 0; i < data.rows(); i++) {
			dataNew(i,j) = data(i,j);
		}
		dataNew(data.rows(),j) = 1;
	}
	alglib::real_1d_array Ptemp = computeDecentered(dataNew,eps);
	alglib::real_1d_array P;
	P.setlength(data.cols());
	for (int i = 0; i < P.length(); i++) {
		P[i] = Ptemp[i];// + Ptemp[i + data.cols()];
	}
	centre = findCentre(data,P);
	int n = data.rows(), m = data.cols();
	alglib::real_2d_array c;
	c.setlength(n,1);
	for (int i = 0; i < n; i++) {
		c(i,0) = centre[i];
	}
	alglib::real_2d_array inv = shape(P,data);
	alglib::real_2d_array temp;
	temp.setlength(n,n);
	alglib::rmatrixgemm(n,n,1,1,c,0,0,0,c,0,0,1,0,temp,0,0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			inv(i,j) -= temp(i,j);
		}
	}
	alglib::ae_int_t info;
	alglib::matinvreport rep;
	alglib::rmatrixinverse(inv,info,rep);
	for (int i = 0; i < inv.rows(); i++) {
		for (int j = 0; j < inv.cols(); j++) {
			inv(i,j) /= (1+eps)*n;
		}
	}
	return inv;
}

Ellipsoid* KhachiyanAlgo(const vector<Point> & points, const double & eps)
{
	if (points.size() < 3) return new Ellipsoid();
	int dim = points[0].Dim();
	int n = points.size();
	real_2d_array data;
	data.setlength(dim, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < dim; ++j)
			data(j,i) = points[i].Coord(j);
	real_1d_array centre;
	centre.setlength(dim);
	real_2d_array form = computeMatrixForm(data, eps, centre);

	bool isUpper = false;
	int ZNeeded = 1;
	real_1d_array eValues;
	eValues.setlength(dim);
	real_2d_array eVectors;
	eVectors.setlength(dim, dim);
	smatrixevd(form, dim, ZNeeded, isUpper, eValues, eVectors);

	Point cent(centre);

	for (int i = 0; i < eValues.length(); ++i)
		eValues[i] = 1.0/sqrt(eValues[i]);

	Ellipsoid* min_el = new Ellipsoid(cent, eVectors, eValues);

	return min_el;
}

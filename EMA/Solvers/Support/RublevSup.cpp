#include "RublevSup.h"
#include "ConvHull.h"
#define DOUBLE_MIN numeric_limits<double>::min()
#ifndef M_PI
#define M_PI (3.141592653589793238462643383279502884)
#endif
#define M_2PI (2.*M_PI)

//---------------------------------------------------------------------------------------------------------
// The ellipse is (ellipsePoint.x/axes[0])^2 + (ellipsePoint.y/axes[1])^2 = 1 with axes[0] >= axes[1](!!!).
// The query point is 'point' with point.x >= 0 and point.y >= 0. The function returns the distance from
// the query point to the ellipse. It also computes the ellipse point 'ellipsePoint'
// in the first quadrant that is closest to 'point'.
//---------------------------------------------------------------------------------------------------------
double DistancePointEllipseSpecial(const double axes[2], const Point2D & point, Point2D & ellipsePoint)
{
	double distance = 0.0;
	double x[2] = { ellipsePoint.X(), ellipsePoint.Y() };
	if (point.Y() > 0.0)
	{
		if (point.X() > 0.0)
		{
			// Bisect to compute the root of F(t) for t >= -e1*e1.
			double esqr[2] = { axes[0]*axes[0], axes[1]*axes[1] };
			double ey[2] = { axes[0]*point.X(), axes[1]*point.Y() };
			double t0 = -esqr[1] + ey[1];
			double t1 = -esqr[1] + sqrt(ey[0]*ey[0] + ey[1]*ey[1]);
			double t = t0;

			const int imax = 2*numeric_limits<double>::max_exponent;

			for (int i = 0; i < imax; ++i)
			{
				t = 0.5*(t0 + t1);
				if (t == t0 || t == t1)
				{
					break;
				}
				double r[2] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]) };
				double f = r[0]*r[0] + r[1]*r[1] - 1.0;
				if (f > 0.0)
				{
					t0 = t;
				}
				else if (f < 0.0)
				{
					t1 = t;
				}
				else
				{
					break;
				}
			}

			x[0] = esqr[0]*point.X()/(t + esqr[0]);
			x[1] = esqr[1]*point.Y()/(t + esqr[1]);
			double d[2] = { x[0] - point.X(), x[1] - point.Y() };
			distance = sqrt(d[0]*d[0] + d[1]*d[1]);
		}
		else //point.X() = 0.0
		{
			x[0] = 0.0;
			x[1] = axes[1];
			distance = fabs(point.Y() - axes[1]);
		}
	}
	else // point.Y() = 0.0
	{
		double denom0 = axes[0]*axes[0] - axes[1]*axes[1];
		double e0y0 = axes[0]*point.X();
		if (e0y0 < denom0)
		{
			// y0 is inside the subinterval.
			double x0de0 = e0y0 / denom0;
			double x0de0sqr = x0de0 * x0de0;
			x[0] = axes[0]*x0de0;
			x[1] = axes[1]*sqrt(fabs(1.0 - x0de0sqr));
			double d0 = x[0] - point.X();
			distance = sqrt(d0*d0 + x[1]*x[1]);
		}
		else
		{
			// y0 is outside the subinterval. The closest ellipse point has
			// x1 == 0 and is on the domain-boundary interval (x0/e0)^2 = 1.
			x[0] = axes[0];
			x[1] = 0.0;
			distance = fabs(point.X() - axes[0]);
		}
	}
	ellipsePoint = Point2D(x[0], x[1]);
	return distance;
};

/* The ellipse is (ellipsePoint.x/axes[0])^2 + (ellipsePoint.y/axes[1])^2 = 1.
   The query point is 'point'. The function returns the distance from
   the query point to the ellipse. It also computes the ellipse point 'ellipsePoint'
   in the first quadrant that is closest to 'point'.
*/
double DistancePointEllipseCA (const double axes[2], const Point2D & point, Point2D & ellipsePoint)
{
	// Determine reflections for y to the first quadrant.
	bool reflect[2] = { false, false };
	int i = 0, j = 0;
	double y[2] = { point.X(), point.Y() };
	double x[2] = { ellipsePoint.X(), ellipsePoint.Y() };
	for (i = 0; i < 2; ++i)
	{
		reflect[i] = (y[i] < 0.0);
	}
	// Determine the axis order for decreasing extents.
	int permute[2] = { 0, 0 };
	if (axes[0] < axes[1])
	{
		permute[0] = 1; 
		permute[1] = 0;
	}
	else
	{
		permute[0] = 0; 
		permute[1] = 1;
	}
	int invpermute[2] = { 0, 0 };
	for (i = 0; i < 2; ++i)
	{
		invpermute[permute[i]] = i;
	}
	double locE[2] = { 0.0, 0.0 };
	double locY[2] = { 0.0, 0.0 };
	for (i = 0; i < 2; ++i)
	{
		j = permute[i];
		locE[i] = axes[j];
		locY[i] = y[j];
		if (reflect[j])
		{
			locY[i] = -locY[i];
		}
	}
	double locX[2] = { 0.0, 0.0 };
	Point2D translatedPoint(locY[0], locY[1]);
	double distance = DistancePointEllipseSpecial(locE, translatedPoint, ellipsePoint);
	locX[0] = ellipsePoint.X();
	locX[1] = ellipsePoint.Y();
	// Restore the axis order and reflections.
	for (i = 0; i < 2; ++i)
	{
		j = invpermute[i];
		if (reflect[i])
		{
			locX[j] = -locX[j];
		}
		x[i] = locX[j];
	}
	ellipsePoint = Point2D(x[0], x[1]);
	return distance;
}

/* The ellipse is non-centered and non aligned.
   The query point is 'point'. The function returns the distance from
   the query point to the ellipse. It also computes the ellipse point 'ellipsePoint'
   in the first quadrant that is closest to 'point'.
*/
double DistancePointEllipse(const Ellipsoid2D & e, const Point2D & p, Point2D & ellipsePoint)
{
	double angle = acos(e.Eigenvectors(0, 0) / sqrt(e.Eigenvectors(0,0) * e.Eigenvectors(0,0)
												  + e.Eigenvectors(1,0) * e.Eigenvectors(1,0)));//angle between (1,0)
																								//and eigenvector which 
																								//corresponds axe a
	if (e.Eigenvectors(1,0) < 0) angle*=(-1.0);
	Point2D centre(0.0, 0.0);
	real_2d_array eigenvectors(e.Eigenvectors());
	real_1d_array axes(e.Axes());
	Ellipsoid2D el(centre, eigenvectors, axes);
	el = rotate(el, angle);
	Point2D p1(p.X() - e.Centre().X(), p.Y() - e.Centre().Y());
	p1 = rotate(p1, angle);
	double ax[2] = { el.Axes(0), el.Axes(1) };
	double distance = DistancePointEllipseCA(ax, p1, ellipsePoint);
	Point2D temp(ellipsePoint);
	ellipsePoint = rotate(temp, (-1.0) * angle);
	temp = ellipsePoint;
	ellipsePoint = Point2D(temp.X() + e.Centre().X(), temp.Y() + e.Centre().Y());
	return distance;
}

//finding the farthest to ellipse point
int Farthest(const Ellipsoid2D & e, vector<Point2D> points)
{
	//at first we must to check which points lies outside the ellipse and for them find the distance
	double max_dist = 0.0;
	double angle = acos(e.Eigenvectors(0, 0) / sqrt(e.Eigenvectors(0,0) * e.Eigenvectors(0,0)
												  + e.Eigenvectors(1,0) * e.Eigenvectors(1,0)));//angle between (1,0)
																								//and eigenvector which 
																								//corresponds axe a
	if (e.Eigenvectors(1,0) < 0) angle*=(-1.0);

	int max_point_ind = -1;
	for (int i = 0; i < points.size(); ++i)
	{
		Point2D current(points[i]);
		Point2D temp(current);
		current = Point2D(temp.X() - e.Centre().X(), temp.Y() - e.Centre().Y());
		temp = current;
		current = rotate(temp, angle);
		double dist = 0.0;
		Point2D elp;
		if (((current.X()*current.X())/(e.Axes(0)*e.Axes(0)) + (current.Y()*current.Y())/(e.Axes(1)*e.Axes(1))) > 1.0)
		{
			dist = DistancePointEllipse(e, points[i], elp);
			if (dist > max_dist)
			{
				max_dist = dist;
				max_point_ind = i;
			}
		}
	}
	return max_point_ind;
}

template<class T>
void Swap(T & a, T & b)
{
	T temp = a;
	a = b;
	b = temp;
	return;
}

//finding determinant
double Det(double A11, double A12, double A21, double A22)
{
	return (A11*A22 - A21*A12);
}

#define EPS 0.000000001
inline bool betw (double l, double r, double x)
{
	return min(l,r) <= x + EPS && x <= max(l,r) + EPS;
}
 
inline bool intersect_1d (double a, double b, double c, double d)
{
	if (a > b)  Swap (a, b);
	if (c > d)  Swap (c, d);
	return max (a, c) <= min (b, d) + EPS;
}

bool operator<(const Point2D & p1, const Point2D & p2)
{
	return p1.X() < p2.X() - EPS || fabs(p1.X()-p2.X()) < EPS && p1.Y() < p2.Y() - EPS;
}

inline Point2D max(const Point2D & p1, const Point2D & p2)
{
	return (p1 < p2) ? p2 : p1;
}

//AB & CD
//return left and right ends of the segment which is common for two segments
//when this is the only point that is returned left and right will be the same
bool intersect (Point2D A, Point2D B, Point2D C, Point2D D, Point2D & left, Point2D & right)
{
	if (! intersect_1d (A.X(), B.X(), C.X(), D.X()) || ! intersect_1d (A.Y(), B.Y(), C.Y(), D.Y()))
		return false;
	double A1 = B.Y() - A.Y();
	double B1 = A.X() - B.X();
	double C1 = B.X()*A.Y() - A.X()*B.Y();
	double z = sqrt (A1*A1 + B1*B1);
	if (fabs(z) > EPS) A1 /= z,  B1 /= z,  C1 /= z;
	double A2 = D.Y() - C.Y();
	double B2 = C.X() - D.X();
	double C2 = D.X()*C.Y() - C.X()*D.Y();
	z = sqrt (A2*A2 + B2*B2);
	if (fabs(z) > EPS) A2 /= z,  B2 /= z,  C2 /= z;
	double zn = Det (A1, B1, A2, B2);
	if (fabs (zn) < EPS)
	{
		if (fabs (A1*C.X() + B1*C.Y() + C1) > EPS || fabs (A2*A.X() + B2*A.Y() + C2) > EPS)
			return false;
		if (B < A) Swap (A, B);
		if (D < C) Swap(C, D);
		left = max (A, C);
		right = min (B, D);
		return true;
	}
	else
	{
		double leftX, rightX, leftY, rightY;
		leftX = rightX = - Det (C1, B1, C2, B2) / zn;
		leftY = rightY = - Det (A1, C1, A2, C2) / zn;
		left = Point2D(leftX, leftY);
		right = Point2D(rightX, rightY);
		return betw (A.X(), B.X(), left.X())
			&& betw (A.Y(), B.Y(), left.Y())
			&& betw (C.X(), D.X(), left.X())
			&& betw (C.Y(), D.Y(), left.Y());
	}
}

bool Contains(const Ellipsoid2D & e, const vector<Point2D> & points)
{
	real_2d_array form;
	form.setlength(2,2);
	real_2d_array V;
	V.setlength(2,2);
	V(0,0) = e.Eigenvectors(0,0);
	V(1,0) = e.Eigenvectors(1,0);
	V(0,1) = e.Eigenvectors(0,1);
	V(1,1) = e.Eigenvectors(1,1);
	real_2d_array prom_res;
	prom_res.setlength(2,2);
	real_2d_array lambda = "[[0.0,0.0],[0.0,0.0]]";
	lambda(0,0) = 1.0/(e.Axes(0)*e.Axes(0));
	lambda(1,1) = 1.0/(e.Axes(1)*e.Axes(1));

	rmatrixgemm(2,2,2,1,V,0,0,0,lambda,0,0,0,0,prom_res,0,0);
	rmatrixgemm(2,2,2,1,prom_res,0,0,0,V,0,0,1,0,form,0,0);

	for (int i = 0; i < points.size(); ++i)
	{
		real_2d_array x;
		x.setlength(2,1);
		x(0,0) = points[i].X() - e.Centre().X();
		x(1,0) = points[i].Y() - e.Centre().Y();
		real_2d_array prom_res;
		prom_res.setlength(1, 2);
		real_2d_array res;
		res.setlength(1,1);

		rmatrixgemm(1,2,2,1,x,0,0,1,form,0,0,0,0,prom_res,0,0);
		rmatrixgemm(1,1,2,1,prom_res,0,0,0,x,0,0,0,0,res,0,0);
		double result = res(0,0);
		if (result > 1.0 + 0.01)
		{
			return false;
		}
	}
	return true;
}

//finding afin changing f=Mx on 3 points
vector<vector<double>> afin_change(const Point2D & p1, const Point2D & p2, const Point2D & p3, const Point2D & res1, const Point2D & res2, const Point2D & res3)
{
	double syst[3][4] = { { p1.X(), p1.Y(), 1, res1.X() }, { p2.X(), p2.Y(), 1, res2.X() }, { p3.X(), p3.Y(), 1, res3.X() } };
	
	double det = syst[0][0]*syst[1][1]*syst[2][2] + syst[0][1]*syst[1][2]*syst[2][0] + syst[0][2]*syst[1][0]*syst[2][1]
				- syst[0][2]*syst[1][1]*syst[2][0] - syst[0][0]*syst[1][2]*syst[2][1] - syst[0][1]*syst[1][0]*syst[2][2];
	double delta1 = syst[0][3]*syst[1][1]*syst[2][2] + syst[0][1]*syst[1][2]*syst[2][3] + syst[0][2]*syst[1][3]*syst[2][1]
				- syst[0][2]*syst[1][1]*syst[2][3] - syst[0][3]*syst[1][2]*syst[2][1] - syst[0][1]*syst[1][3]*syst[2][2];
	double delta2 = syst[0][0]*syst[1][3]*syst[2][2] + syst[0][3]*syst[1][2]*syst[2][0] + syst[0][2]*syst[1][0]*syst[2][3]
				- syst[0][2]*syst[1][3]*syst[2][0] - syst[0][0]*syst[1][2]*syst[2][3] - syst[0][3]*syst[1][0]*syst[2][2];
	double delta3 = syst[0][0]*syst[1][1]*syst[2][3] + syst[0][1]*syst[1][3]*syst[2][0] + syst[0][3]*syst[1][0]*syst[2][1]
				- syst[0][3]*syst[1][1]*syst[2][0] - syst[0][0]*syst[1][3]*syst[2][1] - syst[0][1]*syst[1][0]*syst[2][3];
	
	double m11 = delta1 / det;
	double m12 = delta2 / det;
	double v1 = delta3 / det;
	syst[0][3] = res1.Y();
	syst[1][3] = res2.Y();
	syst[2][3] = res3.Y();
	delta1 = syst[0][3]*syst[1][1]*syst[2][2] + syst[0][1]*syst[1][2]*syst[2][3] + syst[0][2]*syst[1][3]*syst[2][1]
			- syst[0][2]*syst[1][1]*syst[2][3] - syst[0][3]*syst[1][2]*syst[2][1] - syst[0][1]*syst[1][3]*syst[2][2];
	delta2 = syst[0][0]*syst[1][3]*syst[2][2] + syst[0][3]*syst[1][2]*syst[2][0] + syst[0][2]*syst[1][0]*syst[2][3]
			- syst[0][2]*syst[1][3]*syst[2][0] - syst[0][0]*syst[1][2]*syst[2][3] - syst[0][3]*syst[1][0]*syst[2][2];
	delta3 = syst[0][0]*syst[1][1]*syst[2][3] + syst[0][1]*syst[1][3]*syst[2][0] + syst[0][3]*syst[1][0]*syst[2][1]
			- syst[0][3]*syst[1][1]*syst[2][0] - syst[0][0]*syst[1][3]*syst[2][1] - syst[0][1]*syst[1][0]*syst[2][3];
	double m21 = delta1 / det;
	double m22 = delta2 / det;
	double v2 = delta3 / det;
	vector<vector<double>> res;
	vector<double> temp;
	temp.push_back(m11);
	temp.push_back(m12);
	temp.push_back(v1);
	res.push_back(temp);
	temp.clear();
	temp.push_back(m21);
	temp.push_back(m22);
	temp.push_back(v2);
	res.push_back(temp);
	return res;
}

/* Cubic equation solution. Real coefficients case.

   int Cubic(double *x,double a,double b,double c);
   Parameters:
   x - solution array (size 3). On output:
       3 real roots -> then x is filled with them;
       1 real + 2 complex -> x[0] is real, x[1] is real part of 
                             complex roots, x[2] - non-negative 
                             imaginary part.
   a, b, c - coefficients, as described 
   Returns: 3 - 3 real roots;
            1 - 1 real root + 2 complex;
            2 - 1 real root + complex roots imaginary part is zero 
                (i.e. 2 real roots). 
*/
int Cubic (double *x, double a, double b, double c) 
{
	double q, r, r2, q3;
	q = (a*a - 3.0*b) / 9.0;
	r = ( a * (2.0*a*a - 9.0*b) + 27.0*c) / 54.0;
	r2 = r*r;
	q3 = q*q*q;
	if(r2 < q3) 
	{
		double t = acos(r / sqrt(q3));
		a /= 3.0; q = -2.0*sqrt(q);
		x[0] = q*cos( t/3.0) - a;
		if (fabs(x[0]) < 1e-10) x[0] = 0.0;
		x[1] = q*cos((t + M_2PI) / 3.0) - a;
		if (fabs(x[1]) < 1e-10) x[1] = 0.0;
		x[2] = q*cos((t - M_2PI) / 3.0) - a;
		if (fabs(x[2]) < 1e-10) x[2] = 0.0;
		return(3);
	}
	else 
	{
		double aa,bb;
		if(r <= 0.0) r = -r;
		aa = -pow(r + sqrt(r2 - q3), 1.0/3.0); 
		if(aa != 0.0) bb = q / aa;
		else bb = 0.;
		a /= 3.0;
		q = aa + bb; 
		r = aa - bb; 
		x[0] = q - a;
		x[1] = (-0.5)*q - a;
		x[2] = (sqrt(3.0)*0.5)*fabs(r);
		if(!(fabs(x[2]) > DOUBLE_MIN)) return(2);
		return(1);
	}
}

//solve the system of linear algebraic equations
vector<double> SLAEsolver(const vector<vector<double>> & matrix)
{
	vector<vector<double>> matr(matrix);
	int n = matr.size();
	if (matr[0].size() > n + 1) return vector<double>();
	for (int col = 0; col < n; ++col)
	{
		int row=0;
		for (row = col; row < n && matr[row][col] == 0.0; ++row);
		if (row == n) return vector<double>();
		swap(matr[row],matr[col]);

		double x = matr[col][col];
		for (int j = 0; j < n + 1; ++j)
			matr[col][j] = matr[col][j]/x;

		vector<double> temp_x;
		for (int i = col; i < n; ++i)
		{
			if (i != col)
			{
				temp_x = matr[col];
				x = matr[i][col];
				for (int k = 0; k < n+1; ++k)
				{
					temp_x[k] = temp_x[k]*x;
					matr[i][k] = matr[i][k] - temp_x[k];
				}
			}
		}
	}

	//Зворотній хід методу Гауса
	for (int i = n - 1; i > 0; --i)
	{
		for (int k = i - 1; k >= 0; --k)
		{
			double x = matr[k][i];
			vector<double> temp_x = matr[i];
			for (int j = 0; j <= n; ++j)
			{
				temp_x[j]=temp_x[j]*x;
				matr[k][j] = matr[k][j] - temp_x[j];
			}
		}
	}

	vector<double> ans(n,0);
	double temp = 0.0;
	ans[n-1] = matr[n-1][n];
	for (int i = n-2; i >= 0; --i)
		ans[i] = matr[i][n];
	return ans;
}

//for that function the first parametr must be the convex polygon, the second parametr - fixed vertex
vector<int>* BaseTri(const vector<Point2D> & conv_hull, const int & num_of_vertex, Window* window)
{
	int size = conv_hull.size();
	if (size < 3) return new vector<int>();
	if (num_of_vertex > size) return new vector<int>();
	if (size == 3)
	{
		int temp[3] = {0, 1, 2};
		return new vector<int>(temp, temp + 3);
	}

	int A = num_of_vertex, B = 0, C = 0;
	int num3 = 0;
	if (num_of_vertex >= size - 2)
	{
		if (num_of_vertex >= size - 1)
		{
			B = 0;
			num3 = B + 1;
		}
		else
		{
			B = num_of_vertex + 1;
			num3 = 0;
		}
	}
	else
	{
		B = num_of_vertex + 1;
		num3 = B + 1;
	}
	C = num3;

	Triangle current_tri(conv_hull[A], conv_hull[B], conv_hull[C]);
	double S0 = current_tri.Area();

	//
	DrawEvent* draw = new DrawEvent();
	DrawEvent* context = NULL;
	if (window != NULL)
	{
		context = window->LastDrawEvent();
		draw->add(1, &current_tri);
		window->AddToLastScene(draw);
	}
	//
	
	for (int i = ((num3 == size - 1) ? 0 : num3 + 1) ; i != num_of_vertex; ++i)
	{
		Triangle temp_tri(conv_hull[A], conv_hull[B], conv_hull[i]);
		//
			if (window != NULL)
			{
				draw->clear();
				draw->add(1, &temp_tri);
				draw->add(context);
				window->AddNewDrawEvent(draw);
			}
			//
		if (temp_tri.Area() > S0)
		{
			num3 = i;
			S0 = temp_tri.Area();
		}
		if (i == size - 1) i = -1;
	}

	C = num3;

	int k = B + 1;

	while(k != num_of_vertex)
	{
		double current_area = 0.0;
		int temp_num3 = num3;
		for (int i = num3; i != num_of_vertex; ++i)
		{
			Triangle temp_tri(conv_hull[A], conv_hull[k], conv_hull[i]);
			//
			if (window != NULL)
			{
				draw->clear();
				draw->add(1, &temp_tri);
				draw->add(context);
				window->AddNewDrawEvent(draw);
			}
			//
			if (temp_tri.Area() > current_area)
			{
				temp_num3 = i;
				current_area = temp_tri.Area();
			}
			if (i == size - 1) i = -1;
		}

		if (current_area > S0)
		{
			S0 = current_area;
			num3 = temp_num3;
			B = k;
			C = num3;
			//
			if (window != NULL)
			{
				draw->clear();
				Triangle* triangle = new Triangle(conv_hull[A], conv_hull[B], conv_hull[C]);
				draw->add(1, triangle);
				delete triangle;
				draw->add(context);
				window->AddNewDrawEvent(draw);
			}
			//

		}

		if (k == size - 1) k = 0;
		else ++k;

		if (k == num3)
			if (num3 == size - 1) num3 = 0;
			else ++num3;
	}

	Triangle max_tri(conv_hull[A], conv_hull[B], conv_hull[C]);
	vector<int>* max_tri_index = new vector<int>();
	max_tri_index->push_back(A);
	max_tri_index->push_back(B);
	max_tri_index->push_back(C);

	//
	if (window != NULL)
	{
		draw->clear();
		draw->add(context);
		draw->add(1, &max_tri);
		window->AddNewDrawEvent(draw);
		draw->del(draw->NumberOfElements() - 1);
		window->AddNewDrawEvent(draw);
		draw->clear();
		window->AddNewDrawEvent(context);
		delete context;
	}
	delete draw;
	//

	return max_tri_index;
}

//for that function the first parametr must be the convex polygon
Triangle MaxTri(const vector<Point2D> & conv_hull, Window* window)
{
	int size = conv_hull.size();
	if (size < 3) return Triangle();
	if (size == 3)
		return Triangle(conv_hull[0], conv_hull[1], conv_hull[2]);
	
	vector<int>* base_tri_index = BaseTri(conv_hull, 0, window);//indexes of vertexes in convex hull of base triangle of vertex 0

	int A1 = base_tri_index->at(0);
	int A2 = base_tri_index->at(1);
	int A3 = base_tri_index->at(2);

	delete base_tri_index;
	base_tri_index = NULL;

	Triangle base_tri(conv_hull[A1], 
					  conv_hull[A2], 
					  conv_hull[A3]);
	double S0 = base_tri.Area();

	//
	DrawEvent* draw = new DrawEvent();
	DrawEvent* context = NULL;
	if (window != NULL)
	{
		context = window->LastDrawEvent();
		Point2D* p1 = new Point2D(base_tri.P1());
		Point2D* p2 = new Point2D(base_tri.P2());
		Point2D* p3 = new Point2D(base_tri.P3());
		p1->SetColor(Color(1.0f, 0.0f, 0.0f));
		p2->SetColor(Color(1.0f, 0.0f, 0.0f));
		p3->SetColor(Color(1.0f, 0.0f, 0.0f));
		context->add(3, p1, p2, p3);
		draw->add(1, &base_tri);
		draw->add(context);
		window->AddNewDrawEvent(draw);
		delete p1;
		delete p2;
		delete p3;
	}
	//
	
	int Ak = A1 + 1;

	int An2 = A2;
	int An3 = A3;

	int Ak2 = A2;
	int Ak3 = A3;

	int max_tri_indexes[3] = {A1, A2, A3};

	double current_area = 0.0;
	for ( ; Ak != ((A2 == size - 1) ? 0 : A2 + 1); ++Ak)
	{
		for (Ak2 = ((An2 == Ak) ? ((An2 < size-1)?An2+1:0) : An2) ; Ak2 != ((A3 == size - 1) ? 0 : A3 + 1); ++Ak2)
		{
			for (Ak3 = ((An3 == Ak2) ? ((An3 < size-1)?An3+1:0) : An3) ; Ak3 != ((A1 == size - 1) ? 0 : A1 + 1); ++Ak3)
			{
				Triangle temp_tri(conv_hull[Ak], conv_hull[Ak2], conv_hull[Ak3]);

				//
				if( window != NULL)
				{
					draw->clear();
					draw->add(1, &temp_tri);
					draw->add(context);
					window->AddNewDrawEvent(draw);
				}
				//
				double temp_area = 0.0;
				temp_area = temp_tri.Area();
				if (temp_area > current_area)
				{
					An2 = Ak2;
					An3 = Ak3;
					current_area = temp_area;
				}
				if (Ak3 == size - 1) Ak3 = -1;
			}

			if (Ak2 == size - 1) Ak2 = -1;
		}
		if (current_area > S0)
		{
			S0 = current_area;
			max_tri_indexes[0] = Ak;
			max_tri_indexes[1] = An2;
			max_tri_indexes[2] = An3;

			//
			if (window != NULL)
			{
				Triangle* tri = new Triangle(conv_hull[Ak], conv_hull[An2], conv_hull[An3]);
				draw->clear();
				draw->add(1, tri);
				draw->add(context);
				window->AddNewDrawEvent(draw);
				delete tri;
			}
			//
		}
	}
	
	Triangle max_tri(conv_hull[max_tri_indexes[0]], 
					 conv_hull[max_tri_indexes[1]], 
					 conv_hull[max_tri_indexes[2]]);

	//
	if (window != NULL)
	{
		draw->clear();
		Point2D * p1 = new Point2D(max_tri.P1());
		Point2D * p2 = new Point2D(max_tri.P2());
		Point2D * p3 = new Point2D(max_tri.P3());
		p1->SetColor(Color(0.0f, 0.0f, 1.0f));
		p2->SetColor(Color(0.0f, 0.0f, 1.0f));
		p3->SetColor(Color(0.0f, 0.0f, 1.0f));
		draw->add(4, p1, p2, p3, &max_tri);
		draw->add(context);
		window->AddNewDrawEvent(draw);
		delete p1;
		delete p2;
		delete p3;
		window->AddNewDrawEvent(context);
		delete context;
		context = NULL;
	}
	delete draw;
	draw = NULL;
	//
	
	return max_tri;
}

///TODO: rewrite next functions in better manner

Ellipsoid2D* ellipse3(const Triangle & tri,Window* window)//ellipse circumscribing the triangle
{
	double phi = acos((tri.P2().X() - tri.P1().X()) / sqrt((tri.P2().X() - tri.P1().X()) * (tri.P2().X() - tri.P1().X())
															+ (tri.P2().Y() - tri.P1().Y())*(tri.P2().Y() - tri.P1().Y())));
	if (tri.P2().Y() < tri.P1().Y()) phi*=-1.0;
	
	//At first we must find afin triangle for our triangle which is placed:
	//AB is placed on OX, and the middle of AB is placed at origin,
	//and the vertex C is in upper half-plane

	//rotate triangle by phi
	Point2D A = rotate(tri.P1(), phi);
	Point2D B = rotate(tri.P2(), phi);
	Point2D C = rotate(tri.P3(), phi);

	//
	DrawEvent* draw = new DrawEvent();
	DrawEvent* context = NULL;
	if (window!=NULL)
	{
		context = window->LastDrawEvent();
		Triangle* rot = new Triangle(A, B, C);
		rot->SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->clear();
		draw->add(context);
		draw->add(1, rot);
		window->AddNewDrawEvent(draw);
		delete rot;
	}
	//

	//shift triangle into the origin
	double kx = (-1) * (B.X() + A.X()) / 2.0;
	double ky = (-1) * (B.Y() + A.Y()) / 2.0;
	
	double x = A.X() + kx;
	double y = A.Y() + ky;
	A = Point2D(x, y);

	x = B.X() + kx;
	y = B.Y() + ky;
	B = Point2D(x, y);

	x = C.X() + kx;
	y = C.Y() + ky;
	C = Point2D(x, y);

	//
	if (window!=NULL)
	{
		Triangle* shifted = new Triangle(A, B, C);
		shifted->SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->clear();
		draw->add(context);
		draw->add(1, shifted);
		window->AddNewDrawEvent(draw);
		delete shifted;
	}
	//

	bool mirror = false;
	//C must be in upper half-plane
	if (C.Y() < 0)
	{
		mirror = true;
		double x = A.X();
		double y = A.Y() * (-1.0);
		A = Point2D(x, y);

		x = B.X();
		y = B.Y() * (-1.0);
		B = Point2D(x, y);

		x = C.X();
		y = C.Y() * (-1.0);
		C = Point2D(x, y);
	}

	//
	if (window!=NULL)
	{
		Triangle* mirror = new Triangle(A, B, C);
		mirror->SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->clear();
		draw->add(context);
		draw->add(1, mirror);
		window->AddToLastScene(draw);
		delete mirror;
	}
	//

	Triangle afin_tri(A, B, C);

	//
	if (window!=NULL)
	{
		draw->clear();
		draw->add(context);
		afin_tri.SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->add(1, &afin_tri);
		window->AddNewDrawEvent(draw);
	}
	//

	double a = afin_tri.P1().X();
	double b = afin_tri.P3().X();
	double c = afin_tri.P3().Y();
	double a11 = c * c;
	double a22 = b * b + 3 * a * a;
	double a12 = (-1) * b * c;
	double a13 = 0;
	double a23 = (-1) * a * a * c;
	double a33 = a23 * c;
	Ellipsoid2D* min_el = new Ellipsoid2D(a11, a22, a12, a13, a23, a33);//цей еліпс потрібно ще повернути, в обереній послідовності до трикутника
	
	//
	if (window!=NULL)
	{
		draw->clear();
		draw->add(context);
		min_el->SetColor(Color(1.0f, 0.0f, 0.0f));
		draw->add(2, &afin_tri, min_el);
		window->AddNewDrawEvent(draw);
	}
	//

	//And now we must put ellipse into his place

	//At first mirror reflection
	if (mirror)
	{
		Point2D centre(min_el->Centre().X(), (-1.0)*min_el->Centre().Y());
		real_2d_array VR(min_el->Eigenvectors());
		VR(1,0)*=(-1.0);
		VR(1,1)*=(-1.0);
		real_1d_array axes(min_el->Axes());
		delete min_el;
		min_el = new Ellipsoid2D(centre, VR, axes);
	}

	//
	if (window!=NULL)
	{
		draw->clear();
		draw->add(context);
		min_el->SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->add(1, min_el);
		window->AddNewDrawEvent(draw);
	}
	//

	//shifting Ellipse
	x = min_el->Centre().X() - kx;
	y = min_el->Centre().Y() - ky;
	Point2D centre(x, y);
	real_1d_array* axes = new real_1d_array(min_el->Axes());
	real_2d_array* eigenvectors = new real_2d_array(min_el->Eigenvectors());
	delete min_el;
	min_el = new Ellipsoid2D(centre, *eigenvectors, *axes);
	delete eigenvectors;
	delete axes;

	//
	if (window!=NULL)
	{
		draw->clear();
		draw->add(context);
		min_el->SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->add(1, min_el);
		window->AddNewDrawEvent(draw);
	}
	//

	//rotating Ellipse
	double a1 = cos((-1)*phi);
	double a2 = sin((-1)*phi);
	real_2d_array rotate_matrix = "[ [0, 0], [0, 0] ]";
	rotate_matrix(0,0) = a1;
	rotate_matrix(0,1) = a2;
	rotate_matrix(1,0) = (-1)*a2;
	rotate_matrix(1,1) = a1;

	//rotating centre of ellipse
	real_2d_array point = "[ [0], [0] ]";
	point(0,0) = min_el->Centre().X();
	point(1,0) = min_el->Centre().Y();
	real_2d_array centre_res = "[ [0], [0] ]";
	int M = 2, N = 1, K = 2;
	rmatrixgemm(M, N, K, 1, rotate_matrix, 0, 0, 0, point, 0, 0, 0, 0, centre_res, 0, 0);
	Point2D cent = Point2D(centre_res(0,0), centre_res(1,0));
	
	//rotating eigenvectors
	real_2d_array eigenv_res = "[ [0, 0], [0, 0] ]";
	axes = new real_1d_array(min_el->Axes());
	eigenvectors = new real_2d_array(min_el->Eigenvectors());
	M=2, N = 2, K =2;
	rmatrixgemm(M, N, K, 1, rotate_matrix, 0, 0, 0, *eigenvectors, 0, 0, 0, 0, eigenv_res, 0, 0);
	delete min_el;

	min_el = new Ellipsoid2D(cent, eigenv_res, *axes);

	//
	if (window!=NULL)
	{
		draw->clear();
		draw->add(context);
		min_el->SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->add(1, min_el);
		window->AddNewDrawEvent(draw);
		window->AddNewDrawEvent(context);
		delete context;
		context = NULL;
	}
	delete draw;
	draw = NULL;
	//
	
	delete eigenvectors;
	eigenvectors = NULL;
	delete axes;
	axes = NULL;

	return min_el;
}

//for that function quadrangle should not have self-intersections
Ellipsoid2D* ellipse4(const Point2D & p1, const Point2D & p2, const Point2D & p3, const Point2D & p4)
{
	//
	Point2D A(p1);
	Point2D B(p2);
	Point2D C(p3);
	Point2D D(p4);
	//

	Point2D left;
	Point2D right;
	intersect(A, C, B, D,left,right);
	Point2D O(left);

	double AO = distancePToP(A, O);
	double OC = distancePToP(O, C);
	double BO = distancePToP(B, O);
	double OD = distancePToP(O, D);
	double k = sqrt(OC/AO);
	double n = sqrt(BO/OD);

	if (k > 1)
	{
		Swap(A,C);
		AO = distancePToP(A, O);
		OC = distancePToP(O, C);
		k = sqrt(OC/AO);
	}
	if (n > 1)
	{
		Swap(B,D);
		BO = distancePToP(B, O);
		OD = distancePToP(O, D);
		n = sqrt(BO/OD);
	}
	if (k > n)
	{
		Swap(A,B);
		Swap(C,D);
		AO = distancePToP(A, O);
		OC = distancePToP(O, C);
		BO = distancePToP(B, O);
		OD = distancePToP(O, D);
		k = sqrt(OC/AO);
		n = sqrt(BO/OD);
		if (k > 1)
		{
			Swap(A,C);
			AO = distancePToP(A, O);
			OC = distancePToP(O, C);
			k = sqrt(OC/AO);
		}
		if (n > 1)
		{
			Swap(B,D);
			BO = distancePToP(B, O);
			OD = distancePToP(O, D);
			n = sqrt(BO/OD);
		}
	}
				
	double solution[3] = { };
	double t = 0.0;
	double I = (n*n*k*k + 1)*(n*n + k*k);
	double J = 2*n*k*(1 - n*n)*(1 - k*k);
	double L = 4*n*n*k*k;

	Ellipsoid2D* el = NULL;
	if (n >= k && k*k >= (n*n)/(1 + n*n + n*n*n*n))
	{
		int cubic = Cubic(solution, 2*J/L, (2*L - 3*I)/L, J/L);
		if (cubic == 1) t = solution[0];
		else if (cubic == 2)
			{
				if (solution[0] >= 0.0 && solution[0] <= 1.0) t = solution[0];
				if (solution[1] >= 0.0 && solution[1] <= 1.0) t = solution[1];
			} else
				{
					if (solution[0] >= 0.0 && solution[0] <= 1.0) t = solution[0];
					if (solution[1] >= 0.0 && solution[1] <= 1.0) t = solution[1];
					if (solution[2] >= 0.0 && solution[2] <= 1.0) t = solution[2];
				}
		double phi = acos(t);
		double tsin = sin(phi);

		Point2D F((1.0 - k*k)/2.0, (k*(1.0 - n*n) - (1.0 - k*k)*n*t)/(2.0*n*tsin));
		double R = (I - J*t - L*t*t)/(4*n*n*tsin*tsin);

		Point2D afA(1.0,0.0);
		Point2D afB(-k*n*t, -k*n*tsin);

		Point2D afC(-k*k, 0.0);

		vector<vector<double>> af = afin_change(A, B, C, afA, afB, afC);

		double AA = af[0][0]*af[0][0] + af[1][0]*af[1][0];
		double BB = af[0][0]*af[0][1] + af[1][0]*af[1][1];
		double CC = af[0][1]*af[0][1] + af[1][1]*af[1][1];
		double EE = af[0][1]*(af[0][2] - F.X()) + af[1][1]*(af[1][2] - F.Y());
		double DD = af[0][0]*(af[0][2] - F.X()) + af[1][0]*(af[1][2] - F.Y());
		double FF = (af[0][2] - F.X())*(af[0][2] - F.X()) + (af[1][2] - F.Y())*(af[1][2] - F.Y()) - R;

		el = new Ellipsoid2D(AA, CC, BB, DD, EE, FF);
	}
	return el;
}

//for that function the first parametr must be the convex polygon
Ellipsoid2D* points4(const vector<Point2D> & conv, const int & fixed, Window* window)
{
	const int size = conv.size();
	Ellipsoid2D* el = NULL;
	if (size < 4) return el;
	if (fixed > size - 1) return el;
	int A, B, C, D = fixed;
	if (fixed < size - 3) { A = fixed + 3, B = fixed + 2, C = fixed + 1; }
	else if (fixed == size - 3) { A = 0, B = fixed + 2, C = fixed + 1; }
	else if (fixed == size - 2) { A = 1, B = 0, C = fixed + 1; }
	else if (fixed == size - 1) { A = 2, B = 1, C = 0; }

	bool finded = false;

	//
	DrawEvent* context = NULL;
	if (window != NULL)
	{
		context = window->LastDrawEvent();
	}
	//
	
	double min_area = numeric_limits<double>::max();
	Ellipsoid2D *min_el = NULL;

	int stop = 0;
	if (D > 2) stop = D - 3;
	if (D == 0) stop = size - 3;
	if (D == 1) stop = size - 2;
	if (D == 2) stop = size - 1;

	for (int d1 = D; d1 != stop; ++d1)
	{
		int d = d1;
		if (d < size - 3) { A = d + 3, B = d + 2, C = d + 1; }
		else if (d == size - 3) { A = 0, B = d + 2, C = d + 1; }
		else if (d == size - 2) { A = 1, B = 0, C = d + 1; }
		else if (d == size - 1) { A = 2, B = 1, C = 0; }
		for (int a1 = A; a1 != D; ++a1)
		{
			int a = a1;
			for (int b1 = B; b1 != a1; ++b1)
			{
				int b = b1;
				for (int c1 = C; c1 != b1; ++c1)
				{
					int c = c1;
					Point2D left;
					Point2D right;
					if (intersect(conv[a], conv[b], conv[c], conv[d],left,right)) Swap(a,c);
					else if (intersect(conv[a],conv[d],conv[b],conv[c],left,right)) Swap(a,b);
					el = ellipse4(conv[a], conv[b], conv[c], conv[d]);

					//
					if (window != NULL)
					{
						DrawEvent* draw = new DrawEvent();
						Point2D* p1 = new Point2D(conv[a]);
						p1->SetColor(Color(1.0f,0.0f,1.0f));
						p1->SetSize(8.0f);
						Point2D* p2 = new Point2D(conv[b]);
						p2->SetColor(Color(1.0f,0.0f,1.0f));
						p2->SetSize(8.0f);
						Point2D* p3 = new Point2D(conv[c]);
						p3->SetColor(Color(1.0f,0.0f,1.0f));
						p3->SetSize(8.0f);
						Point2D* p4 = new Point2D(conv[d]);
						p4->SetColor(Color(1.0f,0.0f,1.0f));
						p4->SetSize(8.0f);
						draw->add(context);
						draw->add(5, p1, p2, p3, p4, el);
						window->AddNewDrawEvent(draw);
						delete p1;
						delete p2;
						delete p3;
						delete p4;
						delete draw;
					}
					//

					if(el != NULL && Contains(*el, conv) && el->Area() < min_area)
					{
						min_area = el->Area();
						delete min_el;
						min_el = new Ellipsoid2D(*el);
					}

					delete el;
					el = NULL;

					if (c1 == size - 1) c1 = -1;
				}
				if (b1 == size - 1) b1 = -1;
			}
			if (a1 == size - 1) a1 = -1;
		}
		if (d1 == size - 1) d1 = -1;
	}

	//
	if (window != NULL)
	{
		window->AddNewDrawEvent(context);
		delete context;
		context = NULL;
	}
	//

	return min_el;
}

//for that function pentagon should not have self-intersections
Ellipsoid2D* ellipse5(const Point2D & p1, const Point2D & p2, const Point2D & p3, const Point2D & p4, const Point2D & p5)
{
	vector<vector<double>> SLAE;
	vector<double>* temp = new vector<double>();
	temp->push_back(2*p1.X()*p1.Y());
	temp->push_back(p1.Y()*p1.Y());
	temp->push_back(2*p1.X());
	temp->push_back(2*p1.Y());
	temp->push_back(1);
	temp->push_back((-1)*p1.X()*p1.X());
	SLAE.push_back(*temp);
	temp->clear();

	temp->push_back(2*p2.X()*p2.Y());
	temp->push_back(p2.Y()*p2.Y());
	temp->push_back(2*p2.X());
	temp->push_back(2*p2.Y());
	temp->push_back(1);
	temp->push_back((-1)*p2.X()*p2.X());
	SLAE.push_back(*temp);
	temp->clear();

	temp->push_back(2*p3.X()*p3.Y());
	temp->push_back(p3.Y()*p3.Y());
	temp->push_back(2*p3.X());
	temp->push_back(2*p3.Y());
	temp->push_back(1);
	temp->push_back((-1)*p3.X()*p3.X());
	SLAE.push_back(*temp);
	temp->clear();

	temp->push_back(2*p4.X()*p4.Y());
	temp->push_back(p4.Y()*p4.Y());
	temp->push_back(2*p4.X());
	temp->push_back(2*p4.Y());
	temp->push_back(1);
	temp->push_back((-1)*p4.X()*p4.X());
	SLAE.push_back(*temp);
	temp->clear();

	temp->push_back(2*p5.X()*p5.Y());
	temp->push_back(p5.Y()*p5.Y());
	temp->push_back(2*p5.X());
	temp->push_back(2*p5.Y());
	temp->push_back(1);
	temp->push_back((-1)*p5.X()*p5.X());
	SLAE.push_back(*temp);
	temp->clear();

	delete temp;
	temp = NULL;

	vector<double> sol = SLAEsolver(SLAE);
	if (sol.size() < 5) return new Ellipsoid2D();
	double A = 1;
	double B = sol[0];
	double C = sol[1];
	double D = sol[2];
	double E = sol[3];
	double F = sol[4];
	return new Ellipsoid2D(A, C, B, D, E, F);
}

//for that function the first parametr must be the convex polygon
Ellipsoid2D* points5(const vector<Point2D> & conv, const int & fixed, Window* window)
{
	const int size = conv.size();
	Ellipsoid2D* el = NULL;
	if (size < 5) return el;
	if (fixed > size - 1) return el;
	int A, B, C, D, E = fixed;
	if (fixed < size - 4) { A = fixed + 4, B = fixed + 3, C = fixed + 2; D = fixed + 1; }
	else if (fixed == size - 4) { A = 0, B = fixed + 3, C = fixed + 2, D = fixed + 1; }
	else if (fixed == size - 3) { A = 1, B = 0, C = fixed + 2, D = fixed + 1; }
	else if (fixed == size - 2) { A = 2, B = 1, C = 0, D = fixed + 1; }
	else if (fixed == size - 1) { A = 3, B = 2, C = 1, D = 0; }

	//
	DrawEvent* context = NULL;
	if (window != NULL)
	{
		context = window->LastDrawEvent();
	}
	//

	bool finded = false;
	double min_area = numeric_limits<double>::max();
	Ellipsoid2D *min_el = NULL;

	int stop = 0;
	if (E > 3) stop = E - 4;
	if (E == 0) stop = size - 4;
	if (E == 1) stop = size - 3;
	if (E == 2) stop = size - 2;
	if (E == 3) stop = size - 1;
	
	for (int e1 = E; e1 != stop; ++e1)
	{
		int e = e1;
		if (e < size - 4) { A = e + 4, B = e + 3, C = e + 2; D = e + 1; }
		else if (e == size - 4) { A = 0, B = e + 3, C = e + 2, D = e + 1; }
		else if (e == size - 3) { A = 1, B = 0, C = e + 2, D = e + 1; }
		else if (e == size - 2) { A = 2, B = 1, C = 0, D = e + 1; }
		else if (e == size - 1) { A = 3, B = 2, C = 1, D = 0; }
		for (int a1 = A; a1 != E; ++a1)
		{
			int a = a1;
			for (int b1 = B; b1 != a1; ++b1)
			{
				int b = b1;
				for (int c1 = C; c1 != b1; ++c1)
				{
					int c = c1;
					for (int d1 = D; d1 != c1; ++d1)
					{
						int d = d1;
					
						el = ellipse5(conv[a], conv[b], conv[c], conv[d], conv[e]);

						if (window != NULL)
						{
							DrawEvent* draw = new DrawEvent();
							Point2D* p1 = new Point2D(conv[a]);
							p1->SetColor(Color(1.0f,0.0f,1.0f));
							p1->SetSize(8.0f);
							Point2D* p2 = new Point2D(conv[b]);
							p2->SetColor(Color(1.0f,0.0f,1.0f));
							p2->SetSize(8.0f);
							Point2D* p3 = new Point2D(conv[c]);
							p3->SetColor(Color(1.0f,0.0f,1.0f));
							p3->SetSize(8.0f);
							Point2D* p4 = new Point2D(conv[d]);
							p4->SetColor(Color(1.0f,0.0f,1.0f));
							p4->SetSize(8.0f);
							Point2D* p5 = new Point2D(conv[e]);
							p5->SetColor(Color(1.0f,0.0f,1.0f));
							p5->SetSize(8.0f);
							draw->add(context);
							draw->add(6, p1, p2, p3, p4, p5, el);
							window->AddNewDrawEvent(draw);
							delete p1;
							delete p2;
							delete p3;
							delete p4;
							delete p5;
							delete draw;
						}

						if(el != NULL && Contains(*el, conv) && el->Area() < min_area)
						{
							min_area = el->Area();
							delete min_el;
							min_el = new Ellipsoid2D(*el);
						}

						delete el;
						el = NULL;

						if (d1 == size - 1) d1 = -1;
					}
					if (c1 == size - 1) c1 = -1;
				}
				if (b1 == size - 1) b1 = -1;
			}
			if (a1 == size - 1) a1 = -1;
		}
		if (e1 == size - 1) e1 = -1;
	}
	//
	if (window != NULL)
	{
		window->AddNewDrawEvent(context);
		delete context;
		context = NULL;
	}
	//

	return min_el;
}

//for this function points must be connvex hull
/* This function returns ellipse for elementary polygon
 */
Ellipsoid2D* elementaryPolygon(const vector<Point2D> & points, Window* window = NULL)
{
	int n = points.size();
	if (n > 6 || n < 3) return NULL;
	
	Triangle max_tri(MaxTri(points, window));

	//
	DrawEvent* context = NULL;
	if (window != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		context = window->LastDrawEvent();
		draw->add(context);
		draw->add(1, & max_tri);
		window->AddNewDrawEvent(draw);
		delete draw;
		draw = NULL;
	}
	//
	
	Ellipsoid2D* el = ellipse3(max_tri, window);

	//
	if (window != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		draw->add(1, el);
		window->AddToLastScene(draw);
		delete draw; 
		draw = NULL;
	}
	//

	if (el != NULL && Contains(*el, points)) return el;
	int fixed = Farthest(*el, points);
	delete el;

	int fixedA, fixedB, fixedC, fixedD;
	for(int i = 0; i < points.size(); ++i)
	{
		if (points[i] == max_tri.P1()) fixedA = i;
		if (points[i] == max_tri.P2()) fixedB = i;
		if (points[i] == max_tri.P3()) fixedC = i;
	}
	int arr[5] = {fixedA, fixedB, fixedC, fixed};
	sort(arr,arr+4);
	fixedA = arr[0],fixedB = arr[1], fixedC = arr[2], fixedD = arr[3];

	el = ellipse4(points[fixedA], points[fixedB], points[fixedC], points[fixedD]);
	if (el != NULL && Contains(*el, points)) return el;
	else
	{
		delete el;
		el = points4(points, fixed, window);
	}

	//
	if (window != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		Point2D* p = new Point2D(points[fixed]);
		p->SetColor(Color(1.0f,0.0f,0.0f));
		p->SetSize(8.0f);
		draw->add(1, p);
		window->AddToLastScene(draw);
		delete p;
		p = NULL;
		delete draw;
		draw = NULL;
	}
	//

	//
	if (window != NULL && el != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		draw->add(context);
		draw->add(1, el);
		window->AddNewDrawEvent(draw);
		delete draw;
		draw = NULL;
	}
	//

	if (el != NULL)
	{
		if (Contains(*el, points)) return el;
		else
		{
			fixed = Farthest(*el, points);
			arr[4] = fixed;
			sort(arr,arr+5);
			int fixedE;
			fixedA = arr[0], fixedB = arr[1], fixedC = arr[2], fixedD = arr[3], fixedE = arr[4];

			//
			if (window != NULL)
			{
				DrawEvent* draw = new DrawEvent();
				Point2D* p = new Point2D(points[fixed]);
				p->SetColor(Color(1.0f,0.0f,0.0f));
				p->SetSize(8.0f);
				draw->add(1, p);
				window->AddToLastScene(draw);
				delete p;
				p = NULL;
				delete draw;
				draw = NULL;
			}
			//

			delete el;
	
			el = ellipse5(points[fixedA], points[fixedB], points[fixedC], points[fixedD], points[fixedE]);
		}
	}
	else if (n == 5)
			{
				el = ellipse5(points[0], points[1], points[2], points[3], points[4]);
			}

	if (el != NULL && Contains(*el, points)) return el;
	else { if (el != NULL) fixed = Farthest(*el, points);
		   else fixed = 0;
		 }

	delete el;
	el = points5(points, fixed, window);

	return el;
}

Ellipsoid2D* RublevAlg(const vector<Point2D> & points, Window* window)
{
	if (points.size() < 3) return NULL;
	vector<Point2D>* v = ConvHull(points);
	vector<Point2D> conv(*v);
	delete v;
	for (int i = 0; i < conv.size() - 1; ++i)
	{
		int j = i + 1;
		while (j < conv.size() && conv[i] == conv[j])
			conv.erase(conv.begin() + j - 1);
	}
	
	//
	if (window != NULL)
	{
		DrawEvent* draw_event = new DrawEvent();
		for (int i = 0; i < (int)conv.size(); ++i)
		{
			conv.at(i).SetColor(Color(0.0f,1.0f,0.0f));
			conv.at(i).SetSize(8.0f);
			draw_event->add(1, &conv.at(i));
		}
		window->AddToLastScene(draw_event);
		draw_event->clear();
		for (int i = 1; i < (int)conv.size(); ++i)
		{
			Segment* s = new Segment(conv.at(i-1),conv.at(i));
			draw_event->add(1, s);
			delete s;
		}
		Segment* s = new Segment(conv.at(conv.size() - 1), conv.at(0));
		draw_event->add(1, s);
		delete s;
		window->AddToLastScene(draw_event);
		delete draw_event;
		draw_event = NULL;
	}
	//

	Ellipsoid2D* el = new Ellipsoid2D();
	if (conv.size() <= 6)
	{
		delete el;
		el = elementaryPolygon(conv, window);
		return el;
	}

	Triangle max_tri(MaxTri(conv, window));

	//
	DrawEvent* context = NULL;
	if (window != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		context = window->LastDrawEvent();
		draw->add(context);
		draw->add(1, & max_tri);
		window->AddNewDrawEvent(draw);
		delete draw;
		draw = NULL;
	}
	//
	
	delete el;
	el = ellipse3(max_tri, window);

	//
	if (window != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		el->SetColor(Color(0.0f,1.0f,0.0f));
		draw->add(context);
		draw->add(1, el);
		window->AddNewDrawEvent(draw);
		window->AddNewDrawEvent(context);
		delete draw; 
		draw = NULL;
	}
	//

	if (Contains(*el, conv)) return el;
	int fixed = Farthest(*el, conv);
	delete el;

	el = points4(conv, fixed, window);

	//
	if (window != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		Point2D* p = new Point2D(conv[fixed]);
		p->SetColor(Color(1.0f,0.0f,0.0f));
		p->SetSize(8.0f);
		draw->add(context);
		draw->add(1, p);
		window->AddNewDrawEvent(draw);
		delete p;
		p = NULL;
		delete draw;
		draw = NULL;
	}
	//

	//
	if (window != NULL && el != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		draw->add(context);
		el->SetColor(Color(0.0f,1.0f,0.0f));
		draw->add(1, el);
		window->AddNewDrawEvent(draw);
		delete draw;
		draw = NULL;
	}
	//

	if (el != NULL)
	{
		if (Contains(*el, conv)) return el;
		else fixed = Farthest(*el, conv);
	}

	//
	if (window != NULL)
	{
		DrawEvent* draw = new DrawEvent();
		Point2D* p = new Point2D(conv[fixed]);
		p->SetColor(Color(1.0f,0.0f,0.0f));
		p->SetSize(8.0f);
		draw->add(context);
		draw->add(1, p);
		window->AddNewDrawEvent(draw);
		delete p;
		p = NULL;
		delete draw;
		draw = NULL;
	}
	//

	delete el;
	el = points5(conv, fixed, window);

	return el;
}
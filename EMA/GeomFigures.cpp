#include "GeomFigures.h"
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif
#define DOUBLE_MIN numeric_limits<double>::min()
#include <algorithm>

Color Color::operator=(const Color & color)
{
	this->RED = color.RED;
	this->GREEN = color.GREEN;
	this->BLUE = color.BLUE;
	return (*this);
}

Point2D::Point2D(const double & x, const double & y):
	Drawable()
{
	_x = x;
	_y = y;
}

Point2D::Point2D(const Point2D &p):
	Drawable(p._color, p._size)
{
	_x = p._x;
	_y = p._y;
};

Point2D::Point2D(Point2D * p):
	Drawable(p->_color, p->_size)
{
	if (p != NULL)
	{
		_x = p->_x;
		_y = p->_y;
	}
	else
	{
		_x = 0.0;
		_y = 0.0;
	}
};

Point2D Point2D::Scale(const double& xCoef, const double& yCoef)
{
	return Point2D(this->_x * xCoef, this->_y * yCoef);
}

Point2D Point2D::Translation(const double& byX, const double& byY)
{
	return Point2D(this->X() + byX, this->Y() + byY);
}

void Point2D::Draw()
{
	glColor3f(_color.RED, _color.GREEN, _color.BLUE);
	glPointSize(_size);
	glBegin(GL_POINTS);
	glVertex2d(_x, _y);
	glEnd();
};

Point2D* Point2D::Clone()
{
	Point2D* point = new Point2D(*this);
	return point;
}

Point2D::~Point2D()
{
	_x = 0;
	_y = 0;
};

Point2D& Point2D::operator=(const Point2D & p)
{
	_color = p._color;
	_size = p._size;
	_x = p._x;
	_y = p._y;
	return (*this);
};

bool operator==(const Point2D & p1, const Point2D & p2)
{
	return (p1.X() == p2.X() && p1.Y() == p2.Y());
}

bool operator!=(const Point2D & p1, const Point2D & p2)
{
	return (p1.X() != p2.X() || p1.Y() != p2.Y());
}

Point::Point(const int & dim)
{
	_dim = dim;
	_coords.setlength(dim);
	for (int i = 0; i < dim; ++i)
		_coords[i] = 0.0;
}

Point::Point(const vector<double> & coords)
{
	_dim = coords.size();
	_coords.setlength(_dim);
	for (int i = 0; i < _dim; ++i)
		_coords[i] = coords[i];
}

Point::Point(const real_1d_array & coords)
{
	_dim = coords.length();
	_coords.setlength(_dim);
	for (int i = 0; i < _dim; ++i)
		_coords[i] = coords[i];
}

Point::Point(const int & dim, const double & coord ...)
{
	_dim = dim;
	_coords.setlength(dim);
	va_list v1;
	va_start(v1, coord);
	double current = coord;
	for (int i = 0; i < _dim; ++i)
	{
		_coords[i] = current;
		current = va_arg(v1, double);
	}
	va_end(v1);
	return;
}

Point::Point(const Point & point)
{
	_dim = point.Dim();
	_coords = point._coords;
}

Point& Point::operator=(const Point & p)
{
	_dim = p.Dim();
	_coords = p._coords;
	return (*this);
}

Point& Point::operator=(const Point2D & p2d)
{
	_dim = 2;
	_coords.setlength(_dim);
	_coords[0] = p2d.X();
	_coords[1] = p2d.Y();
	return (*this);
}

bool Compare_Points_by_X(const Point2D & p1, const Point2D & p2)
{
	return p1.X() < p2.X();
};

bool Compare_Points_by_Y(const Point2D & p1, const Point2D & p2)
{
	return p1.Y() < p2.Y();
};

//*//linear and afin changes (two-dimensional case)
Point2D rotate(const Point2D & p, const double & angle)
{
	double a1 = cos(angle);
	double a2 = sin(angle);
	real_2d_array matrix = "[ [0, 0], [0, 0] ]";
	matrix(0,0) = a1;
	matrix(0,1) = a2;
	matrix(1,0) = (-1)*a2;
	matrix(1,1) = a1;
	real_2d_array point = "[ [0], [0] ]";
	point(0,0) = p.X();
	point(1,0) = p.Y();
	real_2d_array result = "[ [0], [0] ]";
	int M = 2, N = 1, K = 2;
	rmatrixgemm(M, N, K, 1, matrix, 0, 0, 0, point, 0, 0, 0, 0, result, 0, 0);
	Point2D res(result(0,0), result(1,0));
	return res;
}

vector<double> rotate(vector<double> a, const double & phi)
{
	double a11 = cos(phi);
	double a12 = sin(phi);
	vector< vector<double> > matr;
	vector< double > result;
	result.push_back(a11);
	result.push_back(a12);
	matr.push_back(result);
	result.clear();
	result.push_back((-1.0) * a12);
	result.push_back(a11);
	matr.push_back(result);
	result.clear();
	result = matr * a;
	return result;
}

//rotation of ellipse relative to the center of ellipse
Ellipsoid2D rotate(const Ellipsoid2D & e, const double & phi)
{
	vector<double> e1;
	e1.push_back(e.Eigenvectors(0, 0));
	e1.push_back(e.Eigenvectors(1, 0));
	vector<double> e2;
	e2.push_back(e.Eigenvectors(0, 1));
	e2.push_back(e.Eigenvectors(1, 1));
	Point2D centre(e.Centre());
	vector<double> e1_new = rotate(e1, phi);
	vector<double> e2_new = rotate(e2, phi);
	real_2d_array eigenvectors = " [ [0, 0], [0, 0] ]";
	eigenvectors(0, 0) = e1_new[0];
	eigenvectors(1, 0) = e1_new[1];
	eigenvectors(0, 1) = e2_new[0];
	eigenvectors(1, 1) = e2_new[1];
	real_1d_array axes(e.Axes());
	Ellipsoid2D el(centre, eigenvectors, axes);
	return el;
}
//*//

Segment::Segment():
	Drawable()
{
	_p1 = new Point2D(0.0, 0.0);
	_p2 = new Point2D(0.0, 0.0);
};

Segment::Segment(const Point2D & p1, const Point2D & p2):
	Drawable()
{
	_p1 = new Point2D(p1);
	_p2 = new Point2D(p2);
};

Segment::Segment(const Segment &s):
	Drawable(s._color, s._size)
{
	if (s._p1 != NULL) _p1 = new Point2D(*s._p1);
	else _p1 = new Point2D();
	if (s._p2 != NULL) _p2 = new Point2D(*s._p2);
	else _p2 = new Point2D();
};

Segment::Segment(Segment * s):
	Drawable(s->_color, s->_size)
{
	if (s == NULL || s->_p1 == NULL || s->_p2 == NULL)
	{
		_p1 = new Point2D();
		_p2 = new Point2D();
	}
	else
	{
		_p1 = new Point2D(*s->_p1);
		_p2 = new Point2D(*s->_p2);
	}
};

void Segment::Draw()
{
	if (_p1 == NULL || _p2 == NULL) return;
	glColor3f(_color.RED, _color.GREEN, _color.BLUE);
	glLineWidth(_size);
	glBegin(GL_LINES);
	glVertex2d(_p1->X(), _p1->Y());
	glVertex2d(_p2->X(), _p2->Y());
	glEnd();
};

Segment* Segment::Clone()
{
	Segment* segment = new Segment(*this);
	return segment;
}

Segment::~Segment()
{
	if (_p1 != NULL)
	{
		delete _p1;
		_p1 = NULL;
	}
	if (_p2 != NULL)
	{
		delete _p2;
		_p2 = NULL;
	}
};

Segment& Segment::operator=(const Segment & s)
{
	_color = s._color;
	_size = s._size;
	_p1 = new Point2D(*s._p1);
	_p1 = new Point2D(*s._p1);
	return (*this);
};

int Deviation(const Point2D & p, const Segment & edge)
{
	double C = edge.P1().X() * edge.P2().Y() - edge.P2().X() * edge.P1().Y();
	double dev = (p.Y() * (edge.P2().X() - edge.P1().X()) + 
			  p.X() * (edge.P1().Y() - edge.P2().Y()) + C) * (C > 0 ? -1 : 1);
	if (dev > 0.0) return 1; 
	else
		if (dev < 0) return -1;
		else return 0;
}

double distancePToEdge(const Point2D & p, const Segment & edge)
{
	double A, B, C;
	A = edge.P2().Y() - edge.P1().Y();
	B = edge.P1().X() - edge.P2().X();
	C = edge.P1().Y() * edge.P2().X() - edge.P2().Y() * edge.P1().X();
	double distance = fabs(p.X() * A + p.Y() * B + C)/sqrt(A*A + B*B);
	return distance;
}

double distancePToP(const Point2D & p1, const Point2D & p2)
{
	return sqrt((p2.X() - p1.X()) * (p2.X() - p1.X()) + (p2.Y() - p1.Y()) * (p2.Y() - p1.Y()));
}

vector<Point2D>* Sort_by_X(const vector<Point2D> & points)
{
	vector<Point2D>* sorted_points = new vector<Point2D>(points);
	sort(sorted_points->begin(), sorted_points->end(), Compare_Points_by_X);
	return sorted_points;
}

vector<Point2D>* Sort_by_Y(const vector<Point2D> & points)
{
	vector<Point2D>* sorted_points = new vector<Point2D>(points);
	sort(sorted_points->begin(), sorted_points->end(), Compare_Points_by_Y);
	return sorted_points;
}

Triangle::Triangle():
	Drawable()
{
	_p1 = new Point2D(0.0, 0.0);
	_p2 = new Point2D(0.0, 0.0);
	_p3 = new Point2D(0.0, 0.0);
};

Triangle::Triangle(const Point2D &p1, const Point2D & p2, const Point2D & p3):
	Drawable()
{
	_p1 = new Point2D(p1);
	_p2 = new Point2D(p2);
	_p3 = new Point2D(p3);
};

Triangle::Triangle(const Triangle & tr):
	Drawable(tr._color, tr._size)
{
	_p1 = new Point2D(tr._p1);
	_p2 = new Point2D(tr._p2);
	_p3 = new Point2D(tr._p3);
};

double Triangle::Area()
{
	Segment* s = new Segment(_p1, _p2);
	double height = distancePToEdge(_p3, *s);
	return 0.5 * s->length() * height;
}

void Triangle::Draw()
{
	if (_p1 == NULL || _p2 == NULL || _p3 == NULL) return;
	glColor3f(_color.RED, _color.GREEN, _color.BLUE);
	glLineWidth(_size);
	glBegin(GL_LINE_LOOP);
	glVertex2d(_p1->X(), _p1->Y());
	glVertex2d(_p2->X(), _p2->Y());
	glVertex2d(_p3->X(), _p3->Y());
	glEnd();
};

void Triangle::Fill()
{
	if (_p1 == NULL || _p2 == NULL || _p3 == NULL) return;
	glColor3f(_color.RED, _color.GREEN, _color.BLUE);
	glLineWidth(_size);
	glBegin(GL_TRIANGLES);
	glVertex2d(_p1->X(), _p1->Y());
	glVertex2d(_p2->X(), _p2->Y());
	glVertex2d(_p3->X(), _p3->Y());
	glEnd();
};

Triangle* Triangle::Clone()
{
	Triangle* triangle = new Triangle(*this);
	return triangle;
}

Triangle::~Triangle()
{
	if (_p1 != NULL)
	{
		delete _p1;
		_p1 = NULL;
	}
	if (_p2 != NULL)
	{
		delete _p2;
		_p2 = NULL;
	}
	if (_p3 != NULL)
	{
		delete _p3;
		_p3 = NULL;
	}
}

Triangle& Triangle::operator=(const Triangle & rhs)
{
	_color = rhs._color;
	_size = rhs._size;
	_p1 = new Point2D(rhs.P1());
	_p2 = new Point2D(rhs.P2());
	_p3 = new Point2D(rhs.P3());
	return (*this);
}

bool Ellipsoid2D::CheckVectorsCollinearity(const real_1d_array& v1, const real_1d_array& v2)
{
	return fabs(v1[0] * v2[1] - v1[1] * v2[0]) < 0.0001;
}

Ellipsoid2D::Ellipsoid2D()
{
	_eigenvectors.setlength(2, 2);
	_axes.setlength(2);
	_centre = new Point2D(675, 350);
}

Ellipsoid2D::Ellipsoid2D(const Ellipsoid2D & e)
{
	_eigenvectors = e._eigenvectors;
	_axes = e._axes;
	_color = e._color;
	_size = e._size;
	if (e._centre != NULL) _centre = new Point2D(*e._centre);
	else _centre = new Point2D();
}

Ellipsoid2D::Ellipsoid2D(const Point2D & centre, const real_2d_array & eigenvectors, const real_1d_array & axes)
{
	_centre = new Point2D(centre);
	_eigenvectors = eigenvectors;
	_axes = axes;
}

Ellipsoid2D::Ellipsoid2D(double a11, double a22, double a12, double a13, double a23, double a33)//setting an ellipse using the general equation of the curve of the second order
{
	_axes.setlength(2);
	_eigenvectors.setlength(2, 2);
	double D = a11 * a22 - a12 * a12;
	double delta = a11 * a22 * a33 + 2 * a12 * a23 * a13 - a13 * a22 * a13 - a12 * a12 * a33 - a23 * a23 * a11;
	double I = a11 + a22;
	if (D > 0 && (delta * I) < 0)
	{
		double phi = 0.0;
		if (fabs(a12) > DOUBLE_MIN) phi = atan(2.0 * a12 / (a11 - a22 + 0.000000001)) / 2.0;
		if (phi < 0) phi+=M_PI/2.0;//phi must be >0 and <M_PI/2

		double alfa1 = a11 * cos(phi) * cos(phi) + a12 * sin(2 * phi) + a22 * sin(phi) * sin(phi);
		double alfa2 = a11 * sin(phi) * sin(phi) - a12 * sin(2 * phi) + a22 * cos(phi) * cos(phi);
		double z1 = a13 * cos(phi) + a23 * sin(phi);
		double z2 = a23 * cos(phi) - a13 * sin(phi);
		double right_part = z1 * z1 / alfa1 + z2 * z2 / alfa2 - a33;

		double xc = (a12 * a23 - a13 * a22) / D; //centre of ellipsoid
		double yc = (a13 * a12 - a11 * a23) / D; //centre of ellipsoid
		_centre = new Point2D(xc, yc);

		double lambda1 = right_part / alfa1;//eigenvalue 1
		double lambda2 = right_part / alfa2;//eigenvalue 2
		
		_axes[0] = sqrt(lambda1);
		_axes[1] = sqrt(lambda2);

		//first eigenvector is the vector (1, 0) rotated phi, and the second - the first rotated M_PI / 2
		//e1 <-> _axes[0], e2 <-> _axes[1]
		double e1x = cos(phi);
		double e1y = sin(phi);
		double e2x = (-1) * e1y;
		double e2y = e1x;

		double e_matrix[4]; //eigenvectors
		e_matrix[0] = e1x;
		e_matrix[1] = e2x;
		e_matrix[2] = e1y;
		e_matrix[3] = e2y;

		_eigenvectors.setcontent(2, 2, e_matrix);
	}
	else
	{
		_centre = new Point2D(0, 0);
	}
}

double Ellipsoid2D::Eigenvectors(const int & i, const int & j) const
{
	if (i < 0 || i >= _eigenvectors.rows() || j < 0 || j >= _eigenvectors.cols()) return 0.0; 
	return _eigenvectors(i,j);
}

double Ellipsoid2D::Axes(const int & i) const
{
	if (i < 0 || i >= _axes.length()) return -1.0;
	return _axes[i];
}

Ellipsoid2D* Ellipsoid2D::Clone()
{
	Ellipsoid2D* e = new Ellipsoid2D(*this);
	return e;
}

Ellipsoid2D Ellipsoid2D::Rotate(const double& angle)
{
	Ellipsoid2D rotated = rotate(*this, angle);
	Point2D rotatedCentre = rotate(this->_centre, angle);
	delete rotated._centre;
	rotated._centre = NULL;
	rotated._centre = new Point2D(rotatedCentre);
	return rotated;
}

Ellipsoid2D Ellipsoid2D::Scale(const double& xCoef, const double& yCoef)
{
	Point2D centre = this->Centre().Scale(xCoef, yCoef);
	real_1d_array axes;
	axes.setlength(2);
	real_1d_array e1;
	e1.setlength(2);
	real_1d_array xAxe;
	xAxe.setlength(2);
	xAxe[0] = 1;
	xAxe[1] = 0;
	e1[0] = this->Eigenvectors(0, 0);
	e1[1] = this->Eigenvectors(1, 0);
	double x = 0.0, y = 0.0;
	x = this->Axes()[0];
	y = this->Axes()[1];
	if (CheckVectorsCollinearity(e1, xAxe))
	{
		axes[0] = x * xCoef;
		axes[1] = y * yCoef;
	}
	else
	{
		axes[0] = x * yCoef;
		axes[1] = y * xCoef;
	}
	return Ellipsoid2D(centre, this->Eigenvectors(), axes);
}

Ellipsoid2D Ellipsoid2D::Translation(const double& byX, const double& byY)
{
	Point2D centre = this->Centre().Translation(byX, byY);
	return Ellipsoid2D(centre, this->Eigenvectors(), this->Axes());
}

void Ellipsoid2D::Draw()
{
	double angle = acos(_eigenvectors(0, 0) / sqrt(_eigenvectors(0,0) * _eigenvectors(0,0)
												+ _eigenvectors(1,0) * _eigenvectors(1,0)));//angle between (1,0)
																							//and eigenvector which 
																							//corresponds axe a
	if (_eigenvectors(1,0) < 0) angle*=(-1.0);

	glColor3f(1.0f, 0.0f, 0.0f);
	glPointSize(8.0);
	glBegin(GL_POINTS);
	glVertex2d(_centre->X(), _centre->Y());
	glEnd();
	glColor3f(_color.RED, _color.GREEN, _color.BLUE);
	glLineWidth(_size);
	glBegin(GL_LINE_STRIP);
	for (double t = 0.0; t <= 2 * M_PI + 0.001; t+=0.001)
		glVertex2d(_axes[0] * cos(t) * cos(angle) - _axes[1] * sin(t) * sin(angle) + _centre->X(), _axes[0] * cos(t) * sin(angle) + _axes[1] * sin(t) * cos(angle) + _centre->Y());
	glEnd();
}

Ellipsoid2D::~Ellipsoid2D()
{
	if (_centre != NULL) 
	{
		delete _centre;
		_centre = NULL;
	}
}

Ellipsoid2D& Ellipsoid2D::operator=(const Ellipsoid2D & rhs)
{
	_color = rhs._color;
	_size = rhs._size;
	_centre = new Point2D(rhs.Centre());
	_eigenvectors = rhs.Eigenvectors();
	_axes = rhs.Axes();
	return (*this);
}

Ellipsoid2D& Ellipsoid2D::operator=(const Ellipsoid & rhs)
{
	if (rhs.Dim() != 2) return (*this);
	_color = Color();
	_size = 3.0;
	_centre = new Point2D(rhs.Centre().Coord(0), rhs.Centre().Coord(1));
	_eigenvectors = rhs.Eigenvectors();
	_axes = rhs.Axes();
	return (*this);
}

Ellipsoid::Ellipsoid(const Point & point, const real_2d_array & eigenvectors, const real_1d_array & axes)
{
	_dimension = point.Dim();
	_centre = Point(point);
	_eigenvectors = eigenvectors;
	_axes = axes;
}

Ellipsoid::Ellipsoid(const Ellipsoid & e)
{
	_dimension = e.Dim();
	_centre = Point(e.Centre());
	_eigenvectors = e.Eigenvectors();
	_axes = e.Axes();
}

Ellipsoid& Ellipsoid::operator=(const Ellipsoid & e)
{
	_dimension = e.Dim();
	_centre = Point(e.Centre());
	_eigenvectors = e.Eigenvectors();
	_axes = e.Axes();
	return (*this);
}

Ellipsoid& Ellipsoid::operator=(const Ellipsoid2D & e)
{
	_dimension = 2;
	_centre = e.Centre();
	_eigenvectors = e.Eigenvectors();
	_axes = e.Axes();
	return (*this);
}

bool Ellipsoid2D::Inside(const Point2D & p)
{
	double angle = acos(this->Eigenvectors(0, 0) / sqrt(this->Eigenvectors(0,0) * this->Eigenvectors(0,0)
												  + this->Eigenvectors(1,0) * this->Eigenvectors(1,0)));//angle between (1,0)
																										//and eigenvector which 
																										//corresponds axe a
	if (this->Eigenvectors(1,0) < 0) angle*=(-1.0);

	Point2D temp = Point2D(p.X() - this->Centre().X(), p.Y() - this->Centre().Y());
	Point2D current = rotate(temp, angle);

	if (((current.X()*current.X())/(this->Axes(0)*this->Axes(0)) + (current.Y()*current.Y())/(this->Axes(1)*this->Axes(1))) > 1.0 + 0.001) return false;
	return true;
}

bool Ellipsoid::Inside(const Point & p)
{
	int dim = p.Dim();
	alglib::real_2d_array eVecs(_eigenvectors);
	alglib::real_2d_array eVals;
	eVals.setlength(dim, dim);
	alglib::real_2d_array vec;
	vec.setlength(1, dim);
	for (int i = 0; i < dim; ++i)
	{
		vec(0, i) = p.Coord(i) - _centre.Coord(i);
		for (int j = 0; j < dim; ++j)
		{
			if (i == j)
			{
				double temp = _axes[i] * _axes[i];
				temp = 1 / temp;
				eVals(i, j) = temp;
			}
			else eVals(i, j) = 0.0;
		}
	}
	alglib::real_2d_array shape;
	shape.setlength(dim, dim);
	alglib::real_2d_array res;
	res.setlength(dim, dim);
	rmatrixgemm(dim, dim, dim, 1, eVecs, 0, 0, 0, eVals, 0, 0, 0, 0, res, 0, 0);
	rmatrixgemm(dim, dim, dim, 1, res, 0, 0, 0, eVecs, 0, 0, 1, 0, shape, 0, 0);

	res.setlength(1, dim);
	rmatrixgemm(1, dim, dim, 1, vec, 0, 0, 0, shape, 0, 0, 0, 0, res, 0, 0);
	alglib::real_2d_array temporary;
	temporary.setlength(1, 1);
	rmatrixgemm(1, 1, dim, 1, res, 0, 0, 0, vec, 0, 0, 1, 0, temporary, 0, 0);
	double value = temporary(0, 0);

	if (value <= 1 + 0.0001) return true;
	return false;
}

void Line2D::BuildSegment()
{
	double x1, y1, x2, y2;
	if (B != 0)
	{
		x1 = -10, x2 = 1400;
		y1 = (-A * x1 - C) / B;
		y2 = (-A * x2 - C) / B;
	}
	else
	{
		y1 = -10, y2 = 800;
		x1 = (-B * y1 - C) / A;
		x2 = (-B * y2 - C) / A;
	}
	Point2D p1(x1, y1);
	Point2D p2(x2, y2);
	_s = new Segment(p1, p2);
}

Line2D::Line2D()
{
	A = 0, B = 1, C = 0;
	BuildSegment();
}

Line2D::Line2D(const double& a, const double& b, const double& c)
{
	if (A == 0 && B == 0) throw;
	A = a, B = b, C = c;
	BuildSegment();
}

Line2D::Line2D(Point2D& p1, Point2D& p2)
{
	A = p2.Y() - p1.Y();
	B = p1.X() - p2.X();
	C = p1.Y() * p2.X() - p2.Y() * p1.X();
	BuildSegment();
}

Line2D::Line2D(const Line2D& l)
{
	this->A = l.A;
	this->B = l.B;
	this->C = l.C;
	this->_s = new Segment(l._s);
}

void Line2D::Draw()
{
	_s->Draw();
}

Line2D* Line2D::Clone()
{
	Line2D* copy = new Line2D(*this);
	return copy;
}

Line2D& Line2D::operator=(const Line2D& l)
{
	this->A = l.A;
	this->B = l.B;
	this->C = l.C;
	this->_s = new Segment(l._s);
	return (*this);
}

int Line2D::Deviation(const Point2D& p)
{
	double dev = (A * p.X() + B * p.Y() + C) * (C > 0 ? -1 : 1);
	if (dev > 0.0) return 1;
	if (dev < 0) return -1;
	return 0;
}

double Line2D::Distance(const Point2D& p)
{
	double distance = fabs(p.X() * A + p.Y() * B + C) / sqrt(A*A + B*B);
	return distance;
}

Line2D Line2D::GetParallel(const Point2D& p)
{
	double c = -A * p.X() - B * p.Y();
	Line2D parallel(A, B, c);
	return parallel;
}

Line2D Line2D::GetPerpendicular(const Point2D& p)
{
	double c = -B * p.X() + A * p.Y();
	Line2D perpendicular(B, -A, c);
	return perpendicular;
}

Point2D Line2D::Intersection(const Line2D& l)
{
	double delta = A * l.B - l.A * B;
	double Px = (l.C * B - C * l.B) / delta;
	double Py = (l.A * C - A * l.C) / delta;
	return Point2D(Px, Py);
}

Line2D::~Line2D()
{
	delete _s;
	_s = NULL;
}

Rectangle2D::Rectangle2D()
{
	_vertexes[0] = Point2D();
	_vertexes[1] = Point2D();
	_vertexes[2] = Point2D();
	_vertexes[3] = Point2D();
}

Rectangle2D::Rectangle2D(const Point2D& p1, const Point2D& p2, const Point2D& p3, const Point2D& p4)
{
	_vertexes[0] = p1;
	_vertexes[1] = p2;
	_vertexes[2] = p3;
	_vertexes[3] = p4;
}

Rectangle2D::Rectangle2D(const Rectangle2D& r)
{
	_vertexes[0] = r._vertexes[0];
	_vertexes[1] = r._vertexes[1];
	_vertexes[2] = r._vertexes[2];
	_vertexes[3] = r._vertexes[3];
}

Rectangle2D* Rectangle2D::Clone()
{
	Rectangle2D* copy = new Rectangle2D(*this);
	return copy;
}

void Rectangle2D::Draw()
{
	Segment s1 = S1();
	Segment s2 = S2();
	Segment s3 = S3();
	Segment s4 = S4();
	s1.SetColor(_color);
	s2.SetColor(_color);
	s3.SetColor(_color);
	s4.SetColor(_color);
	s1.Draw();
	s2.Draw();
	s3.Draw();
	s4.Draw();
	_vertexes[0].Draw();
	_vertexes[1].Draw();
	_vertexes[2].Draw();
	_vertexes[3].Draw();
}

int Rectangle2D::LeftBottomIndex()
{
	vector<Point2D>* byY = Sort_by_Y({ _vertexes[0], _vertexes[1], _vertexes[2], _vertexes[3] });
	Point2D& leftBottom = byY->at(0);
	if (byY->at(0).Y() == byY->at(1).Y())
		leftBottom = byY->at(0).X() < byY->at(1).X() ? byY->at(0) : byY->at(1);
	int leftBottomIndex = 0;
	for (int i = 0; i < 4; ++i)
		if (leftBottom == _vertexes[i]) leftBottomIndex = i;
	return leftBottomIndex;
}

Rectangle2D Rectangle2D::Scale(const double& xCoef, const double& yCoef)
{
	return Rectangle2D(_vertexes[0].Scale(xCoef, yCoef), _vertexes[1].Scale(xCoef, yCoef),
		_vertexes[2].Scale(xCoef, yCoef), _vertexes[3].Scale(xCoef, yCoef));
}

Rectangle2D Rectangle2D::Translation(const double& byX, const double& byY)
{
	return Rectangle2D(_vertexes[0].Translation(byX, byY), _vertexes[1].Translation(byX, byY),
		_vertexes[2].Translation(byX, byY), _vertexes[3].Translation(byX, byY));
}

Rectangle2D Rectangle2D::Rotate(const double& angle)
{
	Point2D rotatedP1 = rotate(_vertexes[0], angle);
	Point2D rotatedP2 = rotate(_vertexes[1], angle);
	Point2D rotatedP3 = rotate(_vertexes[2], angle);
	Point2D rotatedP4 = rotate(_vertexes[3], angle);
	return Rectangle2D(rotatedP1, rotatedP2, rotatedP3, rotatedP4);
}
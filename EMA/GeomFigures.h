#ifndef __GEOM_FIGURES__
#define __GEOM_FIGURES__

#include <cstdlib>
#include <cmath>
#include <gl\glut.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <vector>
#include <stdarg.h>
#include "Vendor/cpp/src/linalg.h"
#include <fstream>
#define M_PI 3.141592653589793238462643383279502884

using namespace std;
using namespace alglib;

struct Color
{
	float RED;
	float GREEN;
	float BLUE;
	Color(float red = 0.0f, float green = 0.0f, float blue = 0.0f)
	{
		RED = red;
		GREEN = green;
		BLUE = blue;
	}
	Color(Color & color)
	{
		RED = color.RED;
		GREEN = color.GREEN;
		BLUE = color.BLUE;
	}
	Color operator=(Color const & color);
};

template<class T>
class Cloneable
{
	public:
	virtual T* Clone() = 0;
};

class Drawable: public Cloneable<Drawable>
{
	protected:
		Color _color;
		float _size;
	public:
		Drawable () { _color = Color(); _size = 3.0f; };
		Drawable (const Color & color, const float & size = 3.0f) { _color = color; _size = size; };
		virtual void Draw() = 0;
		virtual Drawable* Clone() = 0;
		void SetColor(const Color & color) { _color = color; };
		void SetSize(const float & size) { _size = size; };
};

class Point2D: public Drawable, public Cloneable<Point2D>
{
	private:
		double _x;
		double _y;
	public:
		Point2D(): Drawable(), _x(0.0), _y(0.0) { };
		Point2D(const double & x, const double & y);
		Point2D(const Point2D &p);
		Point2D(Point2D * p);
		double X() const {return _x;};
		double Y() const {return _y;};
		Point2D Scale(const double & xCoef, const double & yCoef);
		Point2D Translation(const double & byX, const double & byY);
		void Draw();
		Point2D* Clone();
		int Dim() const { return 2; }
		double Coord(const int & i) const { return i == 0 ? _x : _y; }
		real_1d_array Coords() const
		{
			real_1d_array coords;
			coords.setlength(2);
			coords[0] = _x;
			coords[1] = _y;
			return coords;
		}
		~Point2D();
		Point2D& operator=(const Point2D & p);
};

class Point
{
	private:
		real_1d_array _coords;
		int _dim;
	public:
		Point(const int & dim = 2);
		Point(const vector<double> & coords);
		Point(const real_1d_array & coords);
		Point(const int & dim, const double & coord ...);
		Point(const Point & point);
		int Dim() const {return _dim;}
		double Coord(const int & i) const {return _coords[i];}
		real_1d_array Coords() const {return _coords;}
		Point& operator=(const Point & p);
		Point& operator=(const Point2D & p2d);
};

class Segment;

bool operator==(const Point2D & p1, const Point2D & p2);

bool operator!=(const Point2D & p1, const Point2D & p2);

bool Compare_Points_by_X(const Point2D & p1, const Point2D & p2);

bool Compare_Points_by_Y(const Point2D & p1, const Point2D & p2);

int Deviation(const Point2D & p, const Segment & edge);

double distancePToEdge(const Point2D & p, const Segment & edge);

double distancePToP(const Point2D & p1, const Point2D & p2);

class Segment: public Drawable, public Cloneable<Segment>
{
	private:
		Point2D *_p1;
		Point2D *_p2;
	public:
		Segment();
		Segment(const Point2D & p1, const Point2D & p2);
		Segment(const Segment &s);
		Segment(Segment * s);
		double length() {return distancePToP(_p1, _p2);};
		Point2D P1() const { if (_p1 != NULL) return *_p1; else return Point2D();};
		Point2D P2() const { if (_p2 != NULL) return *_p2; else return Point2D();};
		void Draw();
		Segment* Clone();
		~Segment();
		Segment& operator=(const Segment & s);
};

class Triangle: public Drawable, public Cloneable<Triangle>
{
	private:
		Point2D *_p1;
		Point2D *_p2;
		Point2D *_p3;
	public:
		Triangle();
		Triangle(const Point2D & p1, const Point2D & p2, const Point2D & p3);
		Triangle(const Triangle & tr);
		Point2D P1() const {if (_p1 != NULL) return *_p1; else return Point2D();};
		Point2D P2() const {if (_p2 != NULL) return *_p2; else return Point2D();};
		Point2D P3() const {if (_p3 != NULL) return *_p3; else return Point2D();};
		double Area();
		void Draw();
		void Fill();
		Triangle* Clone();
		~Triangle();
		Triangle& operator=(const Triangle & rhs);
};

class Ellipsoid;

class Ellipsoid2D: public Drawable, public Cloneable<Ellipsoid2D>
{
	private:
		Point2D* _centre;
		real_2d_array _eigenvectors;//column i corresponds to the eigenvector number i
		real_1d_array _axes;//axe i corresponds to the eigenvector number i
		bool CheckVectorsCollinearity(const real_1d_array & v1, const real_1d_array & v2);
	public:
		Ellipsoid2D();
		Ellipsoid2D(double a11, double a22, double a12, double a13, double a23, double a33);
		Ellipsoid2D(const Ellipsoid2D & e);
		Ellipsoid2D(const Point2D & centre, const real_2d_array & eigenvectors, const real_1d_array & axes);
		Point2D Centre() const { if (_centre != NULL) return *_centre; else return Point2D(); };
		real_2d_array Eigenvectors() const { return _eigenvectors; };
		real_1d_array Axes() const { return _axes; };
		double Eigenvectors(const int & i, const int & j) const;
		double Axes(const int & i) const;
		double Area() const { return M_PI*_axes[0]*_axes[1]; }
		int Dim() const { return 2; }
		void Draw();
		Ellipsoid2D* Clone();
		Ellipsoid2D Rotate(const double & angle);
		Ellipsoid2D Scale(const double & xCoef, const double & yCoef);
		Ellipsoid2D Translation(const double & byX, const double & byY);
		~Ellipsoid2D();
		Ellipsoid2D& operator=(const Ellipsoid2D & e);
		Ellipsoid2D& operator=(const Ellipsoid & e);
		bool Inside(const Point2D & p);//returns true if the point is inside ellipse
};

template<class T>
vector<T> operator*(const vector<vector<T>> b, const vector<T> a)
{
	if (a.size() == 0) return vector <T>();
	if (b.size() == 0 || b.size() != int(a.size())) return vector<T>();
	vector<T> result(b.size(), (T)0);
	for (int i = 0; i < b.size(); ++i)
	{
		T temp = 0;
		for (int j = 0; j < b.size(); ++j)
		{
			temp = b[i][j] * a[j];
			result[i] = result[i] + temp;
			temp = (T)0;
		}
	}
	return result;
}

template<class T>
vector< vector<T> > operator*(const vector< vector<T> > b, const vector< vector<T> > a)
{
	int n = b.size();
	if (a.size() == 0) return Matrix<T>();
	if (n == 0 || n != int(a.size())) return Matrix<T>();
	vector<vector<T>> result(n);
	for (int i = 0; i < n; ++i)
	{
		T temp = 0;
		for (int j = 0, k = 0; k < n; ++j)
		{
			if (j == n){ j = 0; ++k; }
			if (k == n) break;
			temp = b[i][j] * a[j][k];
			result[i][k] = result[i][k] + temp;
			temp = (T)0;
		}
		temp = (T)0;
	}
	return result;
}

Point2D rotate(const Point2D & p, const double & angle);
vector<double> rotate(vector<double> a, const double & phi);
Ellipsoid2D rotate(const Ellipsoid2D & e, const double & phi);
vector<Point2D>* Sort_by_X(const vector<Point2D> & points);

vector<Point2D>* Sort_by_Y(const vector<Point2D> & points);

class Line2D : public Drawable, public Cloneable<Line2D>
{
	private:
		double A;
		double B;
		double C;
		Segment* _s;//for drawing
		void BuildSegment();
	public:
		Line2D();
		Line2D(const double & a, const double & b, const double & c);
		Line2D(Point2D & p1, Point2D & p2);
		Line2D(const Line2D & l);
		void Draw();
		Line2D* Clone();
		Line2D& operator=(const Line2D & e);
		int Deviation(const Point2D & p);
		double Distance(const Point2D & p);
		Line2D GetParallel(const Point2D & p);
		Line2D GetPerpendicular(const Point2D & p);
		Point2D Intersection(const Line2D & l);
		void SetColor(const Color & color) { _s->SetColor(color); }
		~Line2D();
};

class Rectangle2D : public Drawable, public Cloneable<Rectangle2D>
{
	private:
		Point2D _vertexes[4];
	public:
		Rectangle2D();
		Rectangle2D(const Point2D & p1, const Point2D & p2, const Point2D & p3, const Point2D & p4);
		Rectangle2D(const Rectangle2D & r);
		Point2D P1() { return _vertexes[0]; }
		Point2D P2() { return _vertexes[1]; }
		Point2D P3() { return _vertexes[2]; }
		Point2D P4() { return _vertexes[3]; }
		double Height() { return distancePToP(_vertexes[0], _vertexes[1]); }
		double Weight() { return distancePToP(_vertexes[1], _vertexes[2]); }
		Rectangle2D* Clone();
		void Draw();
		Segment S1() { return Segment(_vertexes[0], _vertexes[1]); }
		Segment S2() { return Segment(_vertexes[1], _vertexes[2]); }
		Segment S3() { return Segment(_vertexes[2], _vertexes[3]); }
		Segment S4() { return Segment(_vertexes[3], _vertexes[0]); }
		Line2D L1() { return Line2D(_vertexes[0], _vertexes[1]); }
		Line2D L2() { return Line2D(_vertexes[1], _vertexes[2]); }
		Line2D L3() { return Line2D(_vertexes[2], _vertexes[3]); }
		Line2D L4() { return Line2D(_vertexes[3], _vertexes[0]); }
		Point2D Vertex(const int & i) { return _vertexes[i]; }
		Point2D* Vertexes() { return _vertexes; }
		int LeftBottomIndex();
		Rectangle2D Scale(const double & xCoef, const double & yCoef);
		Rectangle2D Translation(const double & byX, const double & byY);
		Rectangle2D Rotate(const double & angle);
		~Rectangle2D(){}
};

class Ellipsoid
{
	private:
		int _dimension;
		Point _centre;
		real_2d_array _eigenvectors;//column i corresponds to the eigenvector number i
		real_1d_array _axes;//axe i corresponds to the eigenvector number i
	public:
		Ellipsoid(int dim = 2): _dimension(dim), _centre(dim) { _eigenvectors.setlength(dim, dim); _axes.setlength(dim); }
		Ellipsoid(const Point & point, const real_2d_array & eigenvectors, const real_1d_array & axes);
		Ellipsoid(const Ellipsoid & e);
		Point Centre() const { return _centre; }
		int Dim() const { return _dimension; }
		real_2d_array Eigenvectors() const { return _eigenvectors; }
		real_1d_array Axes() const { return _axes; }
		~Ellipsoid() { };
		Ellipsoid& operator=(const Ellipsoid & e);
		Ellipsoid& operator=(const Ellipsoid2D & e);
		bool Inside(const Point & p);//returns true if the point is inside ellipse
};



#endif
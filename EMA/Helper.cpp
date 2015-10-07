#include "Helper.h"

vector<Point2D>* Helper::ReadPoints(string NameOfFile)
{
	ifstream fi(NameOfFile);
	vector<Point2D>* points = new vector<Point2D>();
	if (fi.fail()) return NULL;
	double x, y;
	while (fi >> x)
	{
		if (fi >> y)
		{
			Point2D p(x, y);
			points->push_back(p);
		}
	}
	fi.close();
	return points;
};

vector<Point>* Helper::Read(string NameOfFile)
{
	ifstream fi(NameOfFile);
	int dim = 0;
	vector<Point>* points = new vector<Point>();
	if (fi.fail()) return NULL;
	fi >> dim;
	string str;
	fi >> str;
	while (!fi.eof())
	{
		double coord;
		vector<double> coords;
		for (int i = 0; i < dim; ++i)
			if (fi >> coord) coords.push_back(coord);
		if (coords.size() == dim)
		{
			Point p(coords);
			points->push_back(p);
		}
	}
	return points;
}

void Helper::WritePoints(const vector<Point2D> & points, string NameOfFile)
{
	ofstream fo(NameOfFile, ios::app);
	fo.precision(7);
	fo << "\n=================================================\nPoints:\n";
	for (int i = 0; i < points.size(); ++i)
	{
		fo.width(15);
		fo << points[i].X();
		fo << "\t";
		fo.width(15);
		fo << points[i].Y();
		fo << endl;
	}
	fo << "\n=================================================\n";
	fo.close();
	return;
}

void Helper::WritePoints(const vector<Point> & points, string NameOfFile)
{
	ofstream fo(NameOfFile, ios::app);
	if (!points.empty()) fo << "Dimension: " << points[0].Dim();
	fo << "\n=================================================\nPoints:\n";
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = 0; j < points[0].Dim(); ++j)
		{
			fo.width(15);
			fo << points[i].Coord(j) << "\t";
		}
		fo << endl;
	}
	fo << "\n=================================================\n";
	fo.close();
	return;
}

vector<Point2D>* Helper::GenRandVertexes(int number)
{
	Window* window = Window::Create();
	double x = 0.0, y = 0.0;
	vector<Point2D>* points = new vector<Point2D>();
	for (int i = 0; i < number; ++i)
	{
		x = (double)(rand() % (int)(window->Width - window->mainMenu->Width() - 100.0) + 50.0);
		y = (double)(rand() % (int)(window->Height - 100.0) + 50.0);
		Point2D p(x, y);
		points->push_back(p);
	}
	return points;
}

vector<Point2D>* Helper::GenRandVertexes()
{
	Window* window = Window::Create();
	int number_of_vertexes = 0;
	int r = rand() % 3;
	switch (r)
	{
	case 0:
		number_of_vertexes = rand() % 2000;
		break;
	case 1:
		number_of_vertexes = rand() % 1000;
		break;
	case 2:
		number_of_vertexes = rand() % 500;
		break;
	}
	double x = 0.0, y = 0.0;
	vector<Point2D>* points = new vector<Point2D>();
	for (int i = 0; i < number_of_vertexes; ++i)
	{
		x = (double)(rand() % (int)(window->Width - WIDTH - 100.0) + 50.0);
		y = (double)(rand() % (int)(window->Height - 100.0) + 50.0);
		Point2D p(x, y);
		points->push_back(p);
	}
	return points;
};
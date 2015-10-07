#include <fstream>
#include <vector>
#include "Drawing.h"

using namespace std;

class Helper
{
public:
	static vector<Point2D>* ReadPoints(string NameOfFile);//Rublev input
	static void WritePoints(const vector<Point2D> & points, string NameOfFile);
	static void WritePoints(const vector<Point> & points, string NameOfFile);
	static vector<Point2D>* GenRandVertexes(int number);
	static vector<Point2D>* GenRandVertexes();
	static vector<Point>* Read(string NameOfFile);//Multidimensional input
};
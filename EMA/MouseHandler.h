#pragma once
#include <vector>
#include <fstream>
#include "Drawing.h"
#include "Solvers/EllipsoidBuilders.h"
#include "Helper.h"

using namespace std;

class MouseHandler
{
	static IEllipsoidBuilder* GetCorrectBuilder(string nameOfAlgo);
	static void Build2DEllipse(string nameOfAlgo, IEllipsoidBuilder* builder, bool draw);
	static void BuildNDEllipsoid(string nameOfAlgo, IEllipsoidBuilder* builder);
	static void WriteEllipse(ofstream & fo, string nameOfAlgo, const alglib::real_2d_array & eVectors, const alglib::real_1d_array & eValues, const int & dimension, const alglib::real_1d_array & centre);
	static const string Rublev;
	static const string Khachiyan;
	static const string GeneticMatrix;
	static const string GeneticPetunin;
	static const string Petunin;
	static void AddVertex(Point2D & p);
	static void ReadVertexesFromFile();
	static void GenerateRandomVertexes();
	static void ProcessExit();
	static void ReadMultidimensionalVertexes();
	static void WriteResultsToFile(string nameOfFile, string nameOfAlgo, Ellipsoid2D* el, bool multiDim = false);
	static void WriteResultsToFile(string nameOfFile, string nameOfAlgo, Ellipsoid* el, bool multiDim = false);
	static void BuildEllipsoid(string nameOfAlgo, bool multiDim = false, bool draw = false);
public:
	static void HandleMouseClick(int button, int state, int x, int y);
};
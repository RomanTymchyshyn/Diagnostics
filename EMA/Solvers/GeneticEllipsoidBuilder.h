#pragma once
#include "../GeomFigures.h"
#include "../Drawing.h"
#include "IEllipsoidWrapBuilder.h"
#include "IEllipsoidBuilder.h"

struct MatrixEllipsoid;

class GeneticBuilder : public IEllipsoidBuilder
{
		double _eps;
		IEllipsoidWrapBuilder* _builder;

		Ellipsoid* GeneticAlgo(const vector<Point> & points, const double & eps);
		
		void CountProbability(vector<IEllipsoidWrapper*> & population);
		
		IEllipsoidWrapper* GetFittestWrap(const vector<IEllipsoidWrapper*> & population);
		
		void CalculateParams(int numberOfPoints, 
			int & numberOfWorkPoints, 
			int & populationSize, 
			int & selectionNumber,
			int & crossoverNumber, 
			int & mutationNumber, 
			int & transNumber);
		
		vector<IEllipsoidWrapper*> GenerateNewPopulation(const vector<Point> & points, const int & numOfWorkPoints, const int & populationSize);
		
		vector<int> SelectRandomPoints(const int & numOfAllPoints, const int & numOfSelectedPoints);
		
		vector<IEllipsoidWrapper*> NextGeneration(const vector<Point> & points, 
			const vector<IEllipsoidWrapper*> & currentGeneration,
			const int & numOfAllPoints, 
			const int & numOfWorkPoints, 
			const int & populationSize, 
			const int & selectionSize,
			const int & crossoverSize, 
			const int & mutationSize, 
			const int & transplantingSize);
		
		void CalculateFitnesses(vector<IEllipsoidWrapper*> & population, const vector<Point> & points);
		
		vector<IEllipsoidWrapper*> Select(const vector<IEllipsoidWrapper*> & population, const int & selectionSize);
		
		vector<IEllipsoidWrapper*> Crossover(const vector<IEllipsoidWrapper*> & population, const int & crossoverSize);
		
		vector<IEllipsoidWrapper*> Mutation(const vector<IEllipsoidWrapper*> & population, const int & mutationSize);
		
		vector<IEllipsoidWrapper*> Transplanting(const vector<Point> & points, const int & transplSize, const int & numOfWorkPoints);
	public:
		
		GeneticBuilder(IEllipsoidWrapBuilder* builder, const double & eps);
		
		Ellipsoid2D* Exec(const vector<Point2D> & points, Window* window = NULL);
		
		Ellipsoid* Exec(const vector<Point>& points, Window* window);
		
		~GeneticBuilder();
};
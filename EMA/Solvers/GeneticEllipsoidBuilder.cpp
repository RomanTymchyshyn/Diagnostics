#include "GeneticEllipsoidBuilder.h"
#include <iostream>

using namespace std;

Ellipsoid* Mutate(Ellipsoid* el)
{
	alglib::real_1d_array eVals(el->Axes());
	alglib::real_2d_array eVecs(el->Eigenvectors());
	alglib::real_1d_array means(el->Centre().Coords());
	double change = (rand() % 50) / 100.0;
	switch (rand() % 3)
	{
	case 0:
		change += 0.75;
		for (int i = 0; i < eVals.length(); i++) {
			eVals[i] *= change;
		}
		break;
	case 1:
		change += 0.75;
		for (int i = 0; i < means.length(); i++) {
			means[i] *= change;
		}
		break;
	case 2:
		change -= 0.25;
		//Givens rotation
		int k = 0, l = 0;
		while (l == 0)
		{
			l = rand() % eVecs.rows();
		}
		k = rand() % l;
		double s = sqrt(1 - change*change);
		double* giv = new double[eVecs.rows() * eVecs.cols()];
		for (int i = 0; i < eVecs.rows(); i++) {
			for (int j = 0; j < eVecs.cols(); j++) {
				if (i == j)
				{
					if (i == k || i == l)
						*(giv + i*eVecs.cols() + j) = change;
					else
						*(giv + i*eVecs.cols() + j) = 1;
				}
				else if (i == l && j == k)
				{
					*(giv + i*eVecs.cols() + j) = s;
				}
				else if (i == k && j == l)
				{
					*(giv + i*eVecs.cols() + j) = -s;
				}
				else
				{
					*(giv + i*eVecs.cols() + j) = 0;
				}
			}
		}
		alglib::real_2d_array G;
		G.setcontent(eVecs.rows(), eVecs.cols(), giv);
		alglib::rmatrixgemm(el->Dim(), el->Dim(), el->Dim(), 1, eVecs, 0, 0, 0, G, 0, 0, 0, 0, eVecs, 0, 0);
		delete[] giv;
		giv = NULL;
		break;
	}
	Ellipsoid* mutated = new Ellipsoid(Point(means), eVecs, eVals);
	return mutated;
}

bool Contains(const vector<int> & indexes, const int & index)
{
	int n = indexes.size();
	for (int i = 0; i < n; ++i)
	{
		if (indexes[i] == index) return true;
	}
	return false;
}

IEllipsoidWrapper* GetFittestWrap(const vector<IEllipsoidWrapper*> & population)
{
	int fittest = 0;
	int size = population.size();
	for (int i = 1; i < size; ++i)
		if (population[i]->fitness > population[fittest]->fitness)
			fittest = i;
	return population[fittest];
}

void CalculateImprovement(const vector<IEllipsoidWrapper*> & currentPopulation, 
	const vector<IEllipsoidWrapper*> & nextPopulation, 
	double & bestImprovement, 
	double & averageImprovement)
{
	double currentFitness = 0.0;
	double nextFitness = 0.0;
	int currentSize = currentPopulation.size();
	int nextSize = nextPopulation.size();
	for (int i = 0; i < nextSize; ++i)
		currentFitness += currentPopulation[i]->fitness;
	currentFitness /= currentSize;
	for (int i = 0; i < nextSize; ++i)
		nextFitness += nextPopulation[i]->fitness;
	nextFitness /= nextSize;
	averageImprovement = fabs(nextFitness - currentFitness);
	bestImprovement = fabs(GetFittestWrap(currentPopulation)->fitness - GetFittestWrap(nextPopulation)->fitness);
}

void ClearVector(vector<IEllipsoidWrapper*> & v)
{
	int n = v.size();
	for (int i = 0; i < n; ++i)
	{
		delete v[i];
		v[i] = NULL;
	}
	v.clear();
}

Ellipsoid* GeneticBuilder::GeneticAlgo(const vector<Point>& points, const double& eps)
{
	int numberOfPoints = points.size();//p

	int numberOfWorkPoints = 0,//k
		populationSize = 0,//N
		selectionNumber = 0,//s
		crossoverNumber = 0,//c
		mutationNumber = 0,//m
		transplantingNumber = 0;//N-s-c-m

	CalculateParams(numberOfPoints, numberOfWorkPoints, populationSize, selectionNumber, crossoverNumber, mutationNumber, transplantingNumber);

	int h = (points.size() + points[0].Dim() + 1.0) / 2.0;
	vector<IEllipsoidWrapper*> initPopulation = GenerateNewPopulation(points, numberOfWorkPoints, populationSize);
	vector<IEllipsoidWrapper*> nextGeneration = initPopulation;
	initPopulation.clear();

	double bestImprovement = 1;
	double averageImprovement = 1;
	IEllipsoidWrapper* resWrap = NULL;
	bool toDo;
	do
	{
		vector<IEllipsoidWrapper*> currentGeneration = nextGeneration;
		CountProbability(currentGeneration);
		nextGeneration = NextGeneration(points, currentGeneration, numberOfPoints, numberOfWorkPoints, populationSize,
			selectionNumber, crossoverNumber, mutationNumber, transplantingNumber);
		CalculateImprovement(currentGeneration, nextGeneration, bestImprovement, averageImprovement);
		resWrap = GetFittestWrap(nextGeneration);
		toDo = resWrap->CheckPercents ? bestImprovement > eps || resWrap->PercentPointsIn(points) < 0.95 : bestImprovement > eps;
		ClearVector(currentGeneration);
	} while (toDo);

	Ellipsoid* res = new Ellipsoid(*resWrap->_ellipsoid);
	ClearVector(nextGeneration);
	resWrap = NULL;

	return res;
}

void GeneticBuilder::CountProbability(vector<IEllipsoidWrapper*>& population)
{
	int size = population.size();
	double sumOfFitnesses = 0;
	for (int i = 0; i < size; ++i)
		sumOfFitnesses += population[i]->fitness;
	for (int i = 0; i < size; i++)
		population[i]->probability = population[i]->fitness / (sumOfFitnesses == 0 ? 1 : sumOfFitnesses);
}

IEllipsoidWrapper* GeneticBuilder::GetFittestWrap(const vector<IEllipsoidWrapper*>& population)
{
	int fittest = 0;
	int size = population.size();
	for (int i = 1; i < size; ++i)
		if (population[i]->fitness > population[fittest]->fitness)
			fittest = i;
	return population[fittest];
}

void GeneticBuilder::CalculateParams(int numberOfPoints, 
	int& numberOfWorkPoints, 
	int& populationSize, 
	int& selectionNumber, 
	int& crossoverNumber, 
	int& mutationNumber, 
	int& transNumber)
{
	int minNOfWorkPoints = 23;
	int minPopulationSize = 100;

	minNOfWorkPoints = numberOfPoints < minNOfWorkPoints ? numberOfPoints : minNOfWorkPoints;
	minPopulationSize = numberOfPoints < minPopulationSize ? numberOfPoints : minPopulationSize;

	numberOfWorkPoints = numberOfPoints * 0.039;//k
	numberOfWorkPoints = numberOfWorkPoints < minNOfWorkPoints ? minNOfWorkPoints : numberOfWorkPoints;

	populationSize = numberOfPoints * 0.1;//N
	populationSize = populationSize < minPopulationSize ? minPopulationSize : populationSize;

	selectionNumber = populationSize * 0.45;//s
	crossoverNumber = populationSize * 0.25;//c
	mutationNumber = populationSize * 0.15;//m
	transNumber = populationSize - selectionNumber - crossoverNumber - mutationNumber;
}

vector<IEllipsoidWrapper*> GeneticBuilder::GenerateNewPopulation(const vector<Point>& points, 
	const int& numOfWorkPoints, 
	const int& populationSize)
{
	int numOfAllPoints = points.size();
	vector<IEllipsoidWrapper*> result(populationSize);
	for (int i = 0; i < populationSize; ++i)
	{
		vector<int> selectedIndexes = SelectRandomPoints(numOfAllPoints, numOfWorkPoints);
		IEllipsoidWrapper* el = _builder->GetEllipseWrap(points, selectedIndexes);
		result[i] = el->Clone();
		delete el;
	}
	CountProbability(result);
	return result;
}

vector<int> GeneticBuilder::SelectRandomPoints(const int& numOfAllPoints, 
	const int& numOfSelectedPoints)
{
	vector<int> res(numOfSelectedPoints, -1);
	for (int i = 0; i < numOfSelectedPoints; ++i)
	{
		int r = rand() % numOfAllPoints;
		if (Contains(res, r)) --i;
		else res[i] = r;
	}
	return res;
}

vector<IEllipsoidWrapper*> GeneticBuilder::NextGeneration(const vector<Point>& points, 
	const vector<IEllipsoidWrapper*>& currentGeneration, 
	const int& numOfAllPoints, 
	const int& numOfWorkPoints, 
	const int& populationSize, 
	const int& selectionSize, 
	const int& crossoverSize, 
	const int& mutationSize, 
	const int& transplantingSize)
{
	vector<IEllipsoidWrapper*> result;

	vector<IEllipsoidWrapper*> selected = Select(currentGeneration, selectionSize);
	CalculateFitnesses(selected, points);
	vector<IEllipsoidWrapper*> crossed = Crossover(currentGeneration, crossoverSize);
	CalculateFitnesses(crossed, points);
	vector<IEllipsoidWrapper*> mutated = Mutation(currentGeneration, mutationSize);
	CalculateFitnesses(mutated, points);
	vector<IEllipsoidWrapper*> transplanted = Transplanting(points, transplantingSize, numOfWorkPoints);
	CalculateFitnesses(transplanted, points);

	result.insert(result.end(), selected.begin(), selected.end());
	result.insert(result.end(), crossed.begin(), crossed.end());
	result.insert(result.end(), mutated.begin(), mutated.end());
	result.insert(result.end(), transplanted.begin(), transplanted.end());
	selected.clear();
	crossed.clear();
	mutated.clear();
	transplanted.clear();

	return result;
}

void GeneticBuilder::CalculateFitnesses(vector<IEllipsoidWrapper*>& population, 
	const vector<Point>& points)
{
	int size = population.size();
	for (int i = 0; i < size; ++i)
		population[i]->fitness = population[i]->CalcFitness(points);
}

vector<IEllipsoidWrapper*> GeneticBuilder::Select(const vector<IEllipsoidWrapper*>& population, 
	const int& selectionSize)
{
	vector<IEllipsoidWrapper*> res(selectionSize);
	vector<int> selectedIndexes(selectionSize, -1);
	int populationSize = population.size();
	for (int i = 0; i < selectionSize; ++i)
	{
		double r = ((double)rand() / (RAND_MAX));
		double left = 0.0;
		double right = population[0]->probability;
		int j = 0;
		while (j < populationSize - 1 && !(r >= left && r < right))
		{
			left = right;
			right += population[j + 1]->probability;
			++j;
		}
		if (!Contains(selectedIndexes, j))
		{
			selectedIndexes[i] = j;
			res[i] = population[j]->Clone();
		}
		else --i;
	}
	return res;
}

vector<IEllipsoidWrapper*> GeneticBuilder::Crossover(const vector<IEllipsoidWrapper*>& population, 
	const int& crossoverSize)
{
	vector<IEllipsoidWrapper*> result;

	for (int i = 0; i < crossoverSize; ++i)
	{
		vector<IEllipsoidWrapper*> temporary = Select(population, 3);
		Ellipsoid* el = new Ellipsoid(temporary[0]->_ellipsoid->Centre(), temporary[1]->_ellipsoid->Eigenvectors(), temporary[2]->_ellipsoid->Axes());
		IEllipsoidWrapper* newEllipsoid = _builder->GetEllipseWrap(el);
		result.push_back(newEllipsoid);
		ClearVector(temporary);
	}
	return result;
}

vector<IEllipsoidWrapper*> GeneticBuilder::Mutation(const vector<IEllipsoidWrapper*>& population, 
	const int& mutationSize)
{
	vector<IEllipsoidWrapper*> result;

	for (int i = 0; i < mutationSize; ++i)
	{
		vector<IEllipsoidWrapper*> mutant = Select(population, 1);
		Ellipsoid* el = Mutate(mutant[0]->_ellipsoid);
		IEllipsoidWrapper* newEllipsoid = _builder->GetEllipseWrap(el);
		result.push_back(newEllipsoid);
		ClearVector(mutant);
	}
	return result;
}

vector<IEllipsoidWrapper*> GeneticBuilder::Transplanting(const vector<Point>& points, 
	const int& transplSize, 
	const int& numOfWorkPoints)
{
	return GenerateNewPopulation(points, numOfWorkPoints, transplSize);
}

GeneticBuilder::GeneticBuilder(IEllipsoidWrapBuilder* builder, const double& eps)
{
	_builder = builder;
	_eps = eps;
}

Ellipsoid2D* GeneticBuilder::Exec(const vector<Point2D>& points, Window* window)
{
	vector<Point> vertexes(points.size(), Point());
	for (int i = 0; i < points.size(); ++i)
		vertexes[i] = points.at(i);
	Ellipsoid* el = GeneticAlgo(vertexes, _eps);
	Ellipsoid2D* res = new Ellipsoid2D();
	(*res) = (*el);
	delete el; el = NULL;
	return res;
}

Ellipsoid* GeneticBuilder::Exec(const vector<Point>& points, Window* window)
{
	return GeneticAlgo(points, _eps);
}

GeneticBuilder::~GeneticBuilder()
{
	delete _builder;
	_builder = NULL;
}
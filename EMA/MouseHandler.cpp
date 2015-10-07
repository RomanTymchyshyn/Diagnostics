#include "MouseHandler.h"

const string MouseHandler::Rublev = "Rublev";
const string MouseHandler::Khachiyan = "Khachiyan";
const string MouseHandler::GeneticMatrix = "GeneticMatrix";
const string MouseHandler::GeneticPetunin = "GeneticPetunin";
const string MouseHandler::Petunin = "Petunin";

IEllipsoidBuilder* MouseHandler::GetCorrectBuilder(string nameOfAlgo)
{
	if (nameOfAlgo == Rublev)
		return new RublevEllipsoidBuilder();
	if (nameOfAlgo == Khachiyan)
		return new KhachiyanEllipsoidBuilder();
	if (nameOfAlgo == GeneticMatrix)
	{
		IEllipsoidWrapBuilder* matrixBuilder = new MatrixBuilder();
		return new GeneticBuilder(matrixBuilder, 0.0000001);
	}
	if (nameOfAlgo == GeneticPetunin)
	{
		IEllipsoidWrapBuilder* petuninBuilder = new PetuninBuilder();
		return new GeneticBuilder(petuninBuilder, 0.0001);
	}
	if (nameOfAlgo == Petunin)
		return new PetuninBuilder();
}

void MouseHandler::Build2DEllipse(string nameOfAlgo, IEllipsoidBuilder* builder, bool toDraw)
{
	Window* window = Window::Create();
	if (toDraw == true && (nameOfAlgo != Rublev && nameOfAlgo != Petunin)) toDraw = false;

	if (window->vertexes->size() < 3) return;
	Ellipsoid2D* el;

	vector<Point2D>* conv = NULL;
	if (nameOfAlgo == GeneticMatrix || nameOfAlgo == Khachiyan) conv = ConvHull(*window->vertexes);

	vector<Point2D> points;
	if (conv != NULL) points = (*conv);
	else points = *window->vertexes;

	if (toDraw) el = builder->Exec(points, window);
	else el = builder->Exec(points);

	DrawEvent* draw = new DrawEvent();
	draw->add(window->vertexes);
	if (el != NULL)
	{
		el->SetColor(Color(0.0f, 1.0f, 0.0f));
		draw->add(1, el);
	}
	window->AddNewDrawEvent(draw);
	delete draw; draw = NULL;

	WriteResultsToFile("Resources/Output Data.txt", nameOfAlgo, el);

	delete el; el = NULL;
}

void MouseHandler::BuildNDEllipsoid(string nameOfAlgo, IEllipsoidBuilder* builder)
{
	Window* window = Window::Create();
	auto size = window->MultiDimVert->size();
	if (size <= window->MultiDimVert->at(0).Dim()) return;
	Ellipsoid* el = builder->Exec(*window->MultiDimVert);
	MouseHandler::WriteResultsToFile("Resources/Output Data.txt", nameOfAlgo, el);

	delete el; el = NULL;
}

void MouseHandler::WriteEllipse(ofstream & fo, string nameOfAlgo, const alglib::real_2d_array& eVectors, const alglib::real_1d_array& eValues, const int& dimension, const alglib::real_1d_array& centre)
{
	fo << ">>>" << nameOfAlgo << "<<<<\n";
	fo << "Parameters of the ellipsoid:\n";
	fo << "Centre of the ellipsoid:\n";
	fo << "(";
	for (int i = 0; i < dimension; ++i)
	{
		fo << centre[i];
		if (i != dimension - 1) fo << ", ";
	}
	fo << ")\n";
	fo << "Eigenvectors and Semi-axes:\n";
	for (int i = 0; i < dimension; ++i)
	{
		fo << "Eigenvector#" << i + 1 << ": (";
		for (int j = 0; j < dimension; ++j)
		{
			fo << eVectors(i, j);
			if (j != dimension - 1) fo << ", ";
		}
		fo << ") And the semi-axe which corresponds to this eigenvector: " << eValues[i] << "\n";
	}
}

void MouseHandler::AddVertex(Point2D & p)
{
	Window* window = Window::Create();
	window->vertexes->push_back(p);
	DrawEvent* draw = new DrawEvent();
	draw->add(1, &p);
	window->AddNext(draw);
	delete draw;
	draw = NULL;
}

void MouseHandler::ReadVertexesFromFile()
{
	Window* window = Window::Create();
	if (window->vertexes != NULL)
	{
		if (!window->vertexes->empty()) window->vertexes->clear();
		delete window->vertexes;
		window->vertexes = NULL;
	}
	window->clear();
	window->vertexes = Helper::ReadPoints("Resources/Input Data.txt");
	DrawEvent* draw = new DrawEvent();
	draw->add(window->vertexes);
	window->AddNext(draw);
	delete draw;
	draw = NULL;
}

void MouseHandler::GenerateRandomVertexes()
{
	Window* window = Window::Create();
	if (window->vertexes != NULL)
	{
		if (!window->vertexes->empty()) window->vertexes->clear();
		delete window->vertexes;
		window->vertexes = NULL;
	}
	window->clear();
	window->vertexes = Helper::GenRandVertexes();
}

void MouseHandler::ProcessExit()
{
	Window* window = Window::Create();
	delete window;
	window = NULL;
	exit(0);
}

void MouseHandler::ReadMultidimensionalVertexes()
{
	Window *window = Window::Create();
	if (window->vertexes != NULL)
	{
		if (!window->vertexes->empty()) window->vertexes->clear();
		delete window->vertexes;
		window->vertexes = NULL;
	}
	window->clear();
	window->MultiDimVert = Helper::Read("Resources/Multidimensional Input.txt");
}

void MouseHandler::WriteResultsToFile(string nameOfFile, string nameOfAlgo, Ellipsoid2D* el, bool multiDim)
{
	Window* window = Window::Create();
	if (multiDim)
		Helper::WritePoints(*window->MultiDimVert, nameOfFile);
	else
		Helper::WritePoints(*window->vertexes, nameOfFile);

	ofstream fo(nameOfFile, ios::app);
	WriteEllipse(fo, nameOfAlgo, el->Eigenvectors(), el->Axes(), el->Dim(), el->Centre().Coords());
	if (!multiDim) fo << "Volume:\t" << el->Area() << endl;
	fo.close();
}

void MouseHandler::WriteResultsToFile(string nameOfFile, string nameOfAlgo, Ellipsoid* el, bool multiDim)
{
	Window* window = Window::Create();
	if (multiDim)
		Helper::WritePoints(*window->MultiDimVert, nameOfFile);
	else
		Helper::WritePoints(*window->vertexes, nameOfFile);

	ofstream fo(nameOfFile, ios::app);
	WriteEllipse(fo, nameOfAlgo, el->Eigenvectors(), el->Axes(), el->Dim(), el->Centre().Coords());
	fo.close();
}

void MouseHandler::BuildEllipsoid(string nameOfAlgo, bool multiDim, bool draw)
{
	window->clear();
	IEllipsoidBuilder* builder = MouseHandler::GetCorrectBuilder(nameOfAlgo);
	if (multiDim) BuildNDEllipsoid(nameOfAlgo, builder);
	else Build2DEllipse(nameOfAlgo, builder, draw);
	delete builder; builder = NULL;
}

void MouseHandler::HandleMouseClick(int button, int state, int x, int y)
{
	Window* window = Window::Create();
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && window != NULL)
	{
		if (window->vertexes == NULL) window->vertexes = new vector<Point2D>();
		Point2D p((double)x - window->StartX, window->Height - (double)y - window->StartY);
		if (!(p.X() >= XSTART || p.Y() > window->Height))
		{
			if (window->EnableMouseInput)
				MouseHandler::AddVertex(p);
		}
		else
		{
			if (window->mainMenu != NULL)
			{
				if (window->mainMenu->at(Window::MainMenu::PutByMouse)->In(x, window->Height - y))
					window->EnableMouseInput = !window->EnableMouseInput;

				if (window->mainMenu->at(Window::MainMenu::RenderingInc)->In(x, window->Height - y))
					if (window->rendering > 100)
						window->rendering -= 100;
					else window->rendering -= 20;

					if (window->mainMenu->at(Window::MainMenu::RenderingDec)->In(x, window->Height - y))
						window->rendering += 50;

					if (window->mainMenu->at(Window::MainMenu::Exit)->In(x, window->Height - y))
						MouseHandler::ProcessExit();

					if (window->mainMenu->at(Window::MainMenu::Clear)->In(x, window->Height - y))
					{
						window->clear();
						window->vertexes->clear();
					}

					if (window->mainMenu->at(Window::MainMenu::GenRand)->In(x, window->Height - y))
						MouseHandler::GenerateRandomVertexes();

					if (window->mainMenu->at(Window::MainMenu::ReadFromFile)->In(x, window->Height - y))
						MouseHandler::ReadVertexesFromFile();

					if (window->mainMenu->at(Window::MainMenu::ReadMultidimVertexes)->In(x, window->Height - y))
						MouseHandler::ReadMultidimensionalVertexes();

					if (window->mainMenu->at(Window::MainMenu::RublevResult)->In(x, window->Height - y))
						if (window->vertexes != NULL)
							MouseHandler::BuildEllipsoid(MouseHandler::Rublev);

					if (window->mainMenu->at(Window::MainMenu::RublevSteps)->In(x, window->Height - y))
						if (window->vertexes != NULL)
							MouseHandler::BuildEllipsoid(MouseHandler::Rublev, false, true);

					if (window->mainMenu->at(Window::MainMenu::Khachijan)->In(x, window->Height - y))
					{
						if (window->vertexes != NULL && !window->vertexes->empty())
							MouseHandler::BuildEllipsoid(MouseHandler::Khachiyan);

						if (window->MultiDimVert != NULL && !window->MultiDimVert->empty())
							MouseHandler::BuildEllipsoid(MouseHandler::Khachiyan, true, false);
					}
					if (window->mainMenu->at(Window::MainMenu::GeneticMatrix)->In(x, window->Height - y))
					{
						if (window->vertexes != NULL && !window->vertexes->empty())
							MouseHandler::BuildEllipsoid(MouseHandler::GeneticMatrix);

						if (window->MultiDimVert != NULL && !window->MultiDimVert->empty())
							MouseHandler::BuildEllipsoid(MouseHandler::GeneticMatrix, true, false);
					}
					if (window->mainMenu->at(Window::MainMenu::GeneticPetunin)->In(x, window->Height - y))
						if (window->vertexes != NULL && !window->vertexes->empty())
							MouseHandler::BuildEllipsoid(MouseHandler::GeneticPetunin);

					if (window->mainMenu->at(Window::MainMenu::PetuninResult)->In(x, window->Height - y))
						if (window->vertexes != NULL && !window->vertexes->empty())
							MouseHandler::BuildEllipsoid(MouseHandler::Petunin);

					if (window->mainMenu->at(Window::MainMenu::PetuninSteps)->In(x, window->Height - y))
						if (window->vertexes != NULL && !window->vertexes->empty())
							MouseHandler::BuildEllipsoid(MouseHandler::Petunin, false, true);
			}
		}
	}
	return;
}
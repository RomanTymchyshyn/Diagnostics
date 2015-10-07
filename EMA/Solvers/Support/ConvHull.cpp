#include "ConvHull.h"

//*//finding Convex Hull by fast recursive method

vector<Point2D>* Recursia(const vector<Point2D> & vertexes, const Point2D & left, const Point2D & right, Point2D & top, Window* window)
{
	vector<Point2D>* points = new vector<Point2D>(vertexes);
	int points_size = points->size();
	if (points_size < 2) return points;
	top = points->at(0);
	Segment* base_edge = new Segment(left, right);

	//
	DrawEvent* draw = new DrawEvent();
	if (window != NULL)
	{
		draw->add(1, base_edge);
		Point2D* p = new Point2D(left);
		draw->add(1, p);
		delete p;
		p = new Point2D(right);
		draw->add(1, p);
		delete p;
		for (int i = 0; i < (int)points->size(); ++i)
		{
			draw->add(1, &points->at(i));
		}
		window->AddNewDrawEvent(draw);
	}
	//

	double dist = distancePToEdge(top, *base_edge);
	for (int i = 1; i < points_size; ++i)
	{
		double temp_dist = distancePToEdge(points->at(i), *base_edge);
		if (temp_dist > dist)
		{
			top = points->at(i);
			dist = temp_dist;
		}
		else
			if (temp_dist == dist)
			{
				Segment * s1 = new Segment(left, top);
				if (Deviation(points->at(i), *s1) != Deviation(right, *s1))
				{
					top = points->at(i);
					dist = temp_dist;
				}
				delete s1;
			}
	}
	delete base_edge;

	vector<Point2D>* S1 = new vector<Point2D>();
	vector<Point2D>* S2 = new vector<Point2D>();
	Segment* edge1 = new Segment(left, top);
	Segment* edge2 = new Segment(right, top);

	//
	if (window != NULL)
	{
		draw->clear();
		Point2D* temp = new Point2D(top);
		temp->SetSize(8.0f);
		draw->add(1, temp);
		delete temp;
		draw->add(1, edge1);
		draw->add(1, edge2);
		window->AddToLastScene(draw);
	}
	//

	int dev = Deviation(right, *edge1);
	for (int i = 0; i < points_size; ++i)
		if (points->at(i) != top && Deviation(points->at(i), *edge1) != dev)
			S1->push_back(points->at(i));

	//
	if (window != NULL)
	{
		draw->clear();
		for (int i = 0; i < (int)S1->size(); ++i)
			draw->add(1, &S1->at(i));
		window->AddToLastScene(draw);
	}
	//

	dev = Deviation(left, *edge2);
	for (int i = 0; i < points_size; ++i)
		if (points->at(i) != top && Deviation(points->at(i), *edge2) != dev)
			S2->push_back(points->at(i));

	//
	if (window != NULL)
	{
		draw->clear();
		for (int i = 0; i < (int)S2->size(); ++i)
			draw->add(1, &S2->at(i));
		window->AddToLastScene(draw);
	}
	delete draw;
	draw = NULL;
	//

	vector<Point2D>* result = Recursia(*S1, left, top, Point2D(), window);
	result->push_back(top);
	vector<Point2D>* temp_res = Recursia(*S2, top, right, Point2D(), window);
	for (int i = 0; i < (int)temp_res->size(); ++i)
		result->push_back(temp_res->at(i));
	points->clear();
	delete points;
	S1->clear();
	delete S1;
	S2->clear();
	delete S2;
	delete edge1;
	delete edge2;
	temp_res->clear();
	delete temp_res;
	return result;
}

vector<Point2D>* ConvHull(const vector<Point2D> & points, Point2D & left, Point2D & right, Point2D & up, Point2D & bot, Window* window)
{
	vector<Point2D>* conv = new vector<Point2D>();
	vector<Point2D>* sorted = Sort_by_X(points);

	if (points.size() < 3) return conv;
	int size_of_sorted = sorted->size();
	left = sorted->at(0);
	right = sorted->at(size_of_sorted - 1);
	vector<Point2D>* top = new vector<Point2D>();
	vector<Point2D>* bottom = new vector<Point2D>();
	Segment* base_edge = new Segment(left, right);

	//
	DrawEvent* draw = new DrawEvent();
	if (window != NULL)
	{
		draw->add(1, base_edge);
		window->AddNewDrawEvent(draw);
	}
	//

	for (int i = 0; i < size_of_sorted; ++i)
	{
		if (Deviation(sorted->at(i), *base_edge) == -1) bottom->push_back(sorted->at(i));
		else
			if (Deviation(sorted->at(i), *base_edge) == 1) top->push_back(sorted->at(i));
	}

	//
	if (window != NULL)
	{
		draw->clear();
		for (int i = 0; i < (int)top->size(); ++i)
			draw->add(1, &top->at(i));
		draw->add(1, base_edge);
		window->AddNewDrawEvent(draw);
		draw->clear();
		for (int i = 0; i < (int)bottom->size(); ++i)
			draw->add(1, &bottom->at(i));
		draw->add(1, base_edge);
		window->AddNewDrawEvent(draw);
	}
	//
	delete base_edge;

	conv->push_back(left);

	vector<Point2D>* conv_part1 = new vector<Point2D>(*Recursia(*top, left, right, up, window));
	for (int i = 0; i < (int)conv_part1->size(); ++i)
	{
		conv->push_back(conv_part1->at(i));
	}
	conv->push_back(right);

	vector<Point2D>* conv_part2 = new vector<Point2D>(*Recursia(*bottom, left, right, bot, window));
	reverse(conv_part2->begin(), conv_part2->end());
	for (int i = 0; i < (int)conv_part2->size(); ++i)
	{
		conv->push_back(conv_part2->at(i));
	}

	//
	if (window != NULL)
	{
		draw->clear();
		for (int i = 0; i < (int)sorted->size(); ++i)
			draw->add(1, &sorted->at(i));
		for (int i = 0; i < (int)conv->size(); ++i)
		{
			conv->at(i).SetSize(8.0f);
			conv->at(i).SetColor(Color(0.0f, 1.0f, 0.0f));
			draw->add(1, &conv->at(i));
		}
		window->AddNewDrawEvent(draw);
	}
	delete draw;
	draw = NULL;
	//

	sorted->clear();
	delete sorted;
	top->clear();
	delete top;
	bottom->clear();
	delete bottom;
	conv_part1->clear();
	delete conv_part1;
	conv_part2->clear();
	delete conv_part2;

	return conv;
}
//*//
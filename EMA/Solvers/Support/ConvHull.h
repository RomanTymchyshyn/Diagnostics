#pragma once
#include "../../Drawing.h"

vector<Point2D>* ConvHull(const vector<Point2D> & conv, Point2D & left = Point2D(), Point2D & right = Point2D(), Point2D & up = Point2D(), Point2D & bot = Point2D(), Window* window = NULL);
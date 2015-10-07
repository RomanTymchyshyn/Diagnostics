#include "Drawing.h"
#define RENDERING 200

DrawEvent::DrawEvent(vector<Drawable*>* obj)
{
	_objects = new vector<Drawable*>();
	if (obj != NULL)
	{
		for (int i = 0; i < (int)obj->size(); ++i)
		{
			if (obj->at(i) != NULL)
			{
				Drawable* drawable_obj = obj->at(i)->Clone();
				_objects->push_back(drawable_obj);
			}
		}
	}
}

DrawEvent::DrawEvent(DrawEvent & draw)
{
	_objects = new vector<Drawable*>();
	if (draw._objects != NULL)
	{
		for (int i = 0; i < (int)draw._objects->size(); ++i)
		{
			if (draw._objects->at(i) != NULL) 
			{
				Drawable* drawable_obj = draw._objects->at(i)->Clone();
				_objects->push_back(drawable_obj);
			}
		}
	}
}

void DrawEvent::Draw()
{
	if (_objects == NULL) return;
	for (int i = 0; i < (int)_objects->size(); ++i)
		if (_objects->at(i) != NULL) _objects->at(i)->Draw();
	return;
};

void DrawEvent::add(int number, Drawable * drawable_obj ...)
{
	if (_objects == NULL)_objects = new vector<Drawable*>();
	if (drawable_obj == NULL){ _objects = new vector<Drawable*>(); return;}
	va_list v1;
	va_start(v1, drawable_obj);
	Drawable* current = NULL;
	if (drawable_obj != NULL) current = drawable_obj->Clone();
	for (int i = 0; i < number; ++i)
	{
		if (current != NULL) _objects->push_back(current->Clone());
		current = va_arg(v1,Drawable*);
	}
	va_end(v1);
	return;
}

void DrawEvent::add(vector<Drawable*>* objects)
{
	if (_objects == NULL) _objects = new vector<Drawable*>();
	if (objects == NULL) return;
	for (int i = 0; i < (int)objects->size(); ++i)
		if (objects->at(i) != NULL)
		{
			Drawable* drawable = objects->at(i)->Clone();
			_objects->push_back(drawable);
		}
	return;
}

void DrawEvent::add(vector<Point2D>* objects)
{
	if (_objects == NULL) _objects = new vector<Drawable*>();
	if (objects == NULL) return;
	for (int i = 0; i < (int)objects->size(); ++i)
	{
			Drawable* drawable = objects->at(i).Clone();
			_objects->push_back(drawable);
	}
	return;
}

void DrawEvent::add(vector<Point2D>* objects, Color & color)
{
	if (_objects == NULL) _objects = new vector<Drawable*>();
	if (objects == NULL) return;
	for (int i = 0; i < (int)objects->size(); ++i)
	{
		Drawable* clone = objects->at(i).Clone();
		clone->SetColor(color);
		Drawable* drawable = clone;
		_objects->push_back(drawable);
	}
	return;
}

void DrawEvent::add(DrawEvent* drawing)
{
	if (_objects == NULL) _objects = new vector<Drawable*>();
	if (drawing == NULL || drawing->_objects == NULL) return;
	for (int i = 0; i < drawing->NumberOfElements(); ++i)
		if (drawing->_objects->at(i) != NULL)
		{
			Drawable* drawable = drawing->_objects->at(i)->Clone();
			_objects->push_back(drawable);
		}
	return;
}

Drawable* DrawEvent::last()
{
	if (_objects != NULL) return _objects->back();
	else return NULL;
}

void DrawEvent::del(int index)
{
	if (index < 0 || index >= (int)_objects->size()) return;
	if (_objects->at(index) != NULL) delete _objects->at(index);
	_objects->erase(_objects->begin() + index);
}

void DrawEvent::clear()
{
	if (_objects == NULL) return;
	while(!_objects->empty())
	{
		if (_objects->at(_objects->size() - 1) != NULL) delete _objects->at(_objects->size() - 1);
		_objects->pop_back();
	}
	return;
}

DrawEvent::~DrawEvent()
{
	if (_objects != NULL)
	{
		while(!_objects->empty())
		{
			if (_objects->at(_objects->size() - 1) != NULL) delete _objects->at(_objects->size() - 1);
			_objects->pop_back();
		}
		delete _objects;
		_objects = NULL;
	}
}

Window* Window::_instance = NULL;

Window::Window()
{
	StartX = 0;
	StartY = 0;
	Width = 1350;
	Height = 700;
	DrawEvents = new queue<DrawEvent*>();
	CurrentWindowContent = new DrawEvent();
	BackColor = Color(0.0f, 153.0f/255.0f, 1.0f);
	pause = false;
	EnableMouseInput = false;
	rendering = RENDERING;
	mainMenu = NULL;
	vertexes = NULL;
	MultiDimVert = NULL;
}

Window::Window(int width, int height)
{
	StartX = 0;
	StartY = 0;
	Width = width;
	Height = height;
	DrawEvents = new queue<DrawEvent*>();
	CurrentWindowContent = new DrawEvent();
	BackColor = Color(0.0f, 153.0f/255.0f, 1.0f);
	pause = false;
	rendering = RENDERING;
	EnableMouseInput = false;
	mainMenu = NULL;
	vertexes = NULL;
	MultiDimVert = NULL;
}

Window::Window(int width, int height, Color back_color)
{
	StartX = 0;
	StartY = 0;
	Width = width;
	Height = height;
	DrawEvents = new queue<DrawEvent*>();
	CurrentWindowContent = new DrawEvent();
	BackColor = Color(back_color);
	pause = false;
	rendering = RENDERING;
	EnableMouseInput = false;
	mainMenu = NULL;
	vertexes = NULL;
	MultiDimVert = NULL;
}

Window* Window::Create()
{
	if (_instance == NULL)
	{
		_instance = new Window();
	}
	return _instance;
};

Window* Window::Create(int width, int height)
{
	if (_instance == NULL)
	{
		_instance = new Window(width, height);
	}
	return _instance;
}

Window* Window::Create(int width, int height, Color back_color)
{
	if (_instance == NULL)
	{
		_instance = new Window(width, height, back_color);
	}
	return _instance;
}

void Window::SetStandartMainMenu()
{
	Menu main_menu;
	main_menu.AddToStandartOrder("Put points by mouse");
	main_menu.AddToStandartOrder("Read vertexes from file");
	main_menu.AddToStandartOrder("Read (multidimensional)");
	main_menu.AddToStandartOrder("Generate random vertexes");
	main_menu.AddToStandartOrder("Rublev Algo (Step-by-Step)");
	main_menu.AddToStandartOrder("Rublev Algo (Only Result)");
	main_menu.AddToStandartOrder("Khachiyan Algo (Only Result)");
	main_menu.AddToStandartOrder("Genetic Algo Matrix Variant (Only Result)");
	main_menu.AddToStandartOrder("Genetic Algo Petunin Variant (Only Result)");
	main_menu.AddToStandartOrder("Petunin Algo (Only Result)");
	main_menu.AddToStandartOrder("Petunin Algo (Steps)");
	main_menu.AddToStandartOrder("Clear");
	main_menu.AddToStandartOrder("EXIT");
	mainMenu = new Menu(main_menu);
	MenuItem* mItem = new MenuItem(window->Width - mainMenu->Width() + 20.0,
								window->Height - 50.0, mainMenu->Width()/3.0 - 20.0, 40, "REND+");
	mainMenu->AddMenuItem(mItem, RenderingInc);
	delete mItem;
	mItem = NULL;
	mItem = new MenuItem(window->Width - mainMenu->Width() + (mainMenu->Width()/3.0) + 15.0,
								window->Height - 50.0, mainMenu->Width()/3.0 - 20.0, 40, "REND");
	mainMenu->AddMenuItem(mItem, Rendering);
	delete mItem;
	mItem = NULL;
	mItem = new MenuItem(window->Width - mainMenu->Width() + 2.0*(mainMenu->Width()/3.0) + 15.0,
								window->Height - 50.0, mainMenu->Width()/3.0 - 20.0, 40, "REND-");
	mainMenu->AddMenuItem(mItem, RenderingDec);
	delete mItem;
	mItem = NULL;
	return;
}

void Window::AddNewDrawEvent(DrawEvent* draw_event)
{
	if (draw_event != NULL)
	{
		DrawEvent* copy = new DrawEvent(*draw_event);
		DrawEvents->push(copy);
	}
	return;
};

void Window::AddToLastScene(DrawEvent* draw_event)
{
	if (draw_event == NULL) return;
	if (DrawEvents->size() == 0)
	{
		this->AddNewDrawEvent(draw_event);
		return;
	}
	DrawEvent* last_scene = new DrawEvent(*LastDrawEvent());
	last_scene->add(draw_event);
	DrawEvents->push(last_scene);
	return;
};

void Window::AddNext(DrawEvent * draw_event)
{
	if (draw_event == NULL) return;
	if (DrawEvents->size() == 0)
	{
		this->AddNewDrawEvent(draw_event);
		return;
	}
	if (CurrentWindowContent != NULL)
	{
		CurrentWindowContent->add(draw_event);
	}
	DrawEvents->front()->add(draw_event);
	return;
}

DrawEvent* Window::LastDrawEvent()
{
	if (DrawEvents != NULL && !DrawEvents->empty()) return new DrawEvent(*DrawEvents->back());
	else return new DrawEvent();
};

DrawEvent* Window::FirstDrawEvent()
{
	if (DrawEvents != NULL && !DrawEvents->empty()) return new DrawEvent(*DrawEvents->front());
	else return new DrawEvent();
};

void Window::PopDrawEvent()
{
	if (DrawEvents != NULL)
	{ 
		if (DrawEvents->front() != NULL)
		{
			DrawEvent* draw = DrawEvents->front();
			delete draw;
			draw = NULL;
		}
		DrawEvents->pop();
	}
}

void Window::Translate(const double & x, const double & y)
{
	glTranslated(x, y, 0.0);
	StartX = x;
	StartY = y;
}

void Window::clear()
{
	if (CurrentWindowContent != NULL)
	{
		CurrentWindowContent->clear();
	}
	if (DrawEvents != NULL)
	{
		while (!DrawEvents->empty())
		{
			if (DrawEvents->front() != NULL)
			{
				DrawEvent* draw = DrawEvents->front();
				delete draw;
				draw = NULL;
			}
			DrawEvents->pop();
		}
	}
}

Window::~Window()
{
	if (MultiDimVert != NULL)
	{
		MultiDimVert->clear();
		delete MultiDimVert;
		MultiDimVert = NULL;
	}
	if (vertexes != NULL)
	{
		vertexes->clear();
		delete vertexes;
		vertexes = NULL;
	}
	if (CurrentWindowContent != NULL)
	{
		delete CurrentWindowContent;
		CurrentWindowContent = NULL;
	}
	if (mainMenu != NULL)
	{
		delete mainMenu;
		mainMenu = NULL;
	}
	if (DrawEvents != NULL)
	{
		while (!DrawEvents->empty())
		{
			if (DrawEvents->front() != NULL)
			{
				DrawEvent* draw = DrawEvents->front();
				delete draw;
				draw = NULL;
			}
			DrawEvents->pop();
		}
		delete DrawEvents;
		DrawEvents = NULL;
	}
	_instance = NULL;	
}

MenuItem::MenuItem():
	Drawable(Color(0.0f,(float)(204.0/255.0),1.0f))
{
	_xStart = 0.0;
	_yStart = 0.0;
	_width = 0.0;
	_height = 0.0;
	_text = "";
}

MenuItem::MenuItem(const double & xCoord, const double & yCoord, const double & width, const double & height, Color color):
	Drawable(color)
{
	_xStart = xCoord;
	_yStart = yCoord;
	_width = width;
	_height = height;
	_text = "";
}

MenuItem::MenuItem(const double & xCoord, const double & yCoord, const double & width, const double & height, const string & text, Color color):
	Drawable(color)
{
	_xStart = xCoord;
	_yStart = yCoord;
	_width = width;
	_height = height;
	_text = text;
}

MenuItem::MenuItem(const MenuItem & mItem):
	Drawable(mItem._color, mItem._size)
{
	_xStart = mItem._xStart;
	_yStart = mItem._yStart;
	_width = mItem._width;
	_height = mItem._height;
	_text = mItem._text;
}

void MenuItem::Draw()
{
	glClearColor(0,0,0,0);
	glColor3f(_color.RED, _color.GREEN, _color.BLUE);
	glLineWidth(_size);
	glEnd();
	glPushMatrix();
	glLoadIdentity();
	glRectd(_xStart,_yStart, _xStart + _width, _yStart + _height);
	glColor3f(0.0f,0.0f,0.0f);
	glBegin(GL_LINE_LOOP);
	glVertex2d(_xStart, _yStart);
	glVertex2d(_xStart, _yStart + _height);
	glVertex2d(_xStart + _width, _yStart + _height);
	glVertex2d(_xStart + _width, _yStart);
	glVertex2d(_xStart, _yStart);
	glEnd();
	double length = glutStrokeLength(GLUT_STROKE_ROMAN, (unsigned char *)_text.c_str());
	double coefX = (_width - 40.0) / length;
	double coefY = (_height - 5.0)*1.0/152.38;
	glTranslatef(_xStart + (_width - length*coefX)/2.0, _yStart + (_height - 152.38*coefY) / 2.0 , 0);
	glScalef(coefX, coefY, 0.0f);
	glColor3f(1.0f,0.0f,0.0f);
	for (int i = 0; i < _text.length(); ++i)
		glutStrokeCharacter(GLUT_STROKE_ROMAN, _text[i]);
	glPopMatrix();
	glViewport(0, 0, window->Width, window->Height);
}

MenuItem* MenuItem::Clone()
{
	MenuItem* mItem = new MenuItem(*this);
	return mItem;
}

bool MenuItem::In(const double & xCoord, const double & yCoord)
{
	double x = xCoord;
	double y = yCoord;
	if (_xStart <= x && x <= (_xStart + _width) && y >= _yStart && y <= (_yStart + _height))
		return true;
	else return false;
}

void MenuItem::AddText(const string & text)
{
	_text = text;
}

MenuItem::~MenuItem()
{
	_xStart = 0.0;
	_yStart = 0.0;
	_width = 0.0;
	_height = 0.0;
	_text.clear();
}

Menu::Menu(const double & width, const double & height, const double & xStart, const double & yStart, Color color):
	Drawable(color)
{
	_width = width;
	_height = height;
	_xStart = xStart;
	_yStart = yStart;
}

Menu::Menu(vector<MenuItem*> mItems, const double & width, const double & height, const double & xStart, const double & yStart, Color color):
	Drawable(color)
{
	for (int i = 0; i < mItems.size(); ++i)
	{
		if (mItems[i] != NULL)
		{
			MenuItem* mItem = new MenuItem(*mItems[i]);
			_mItems.push_back(mItem);
		}
	}
	_width = width;
	_height = height;
	_xStart = xStart;
	_yStart = yStart;
}

Menu::Menu(const Menu & menu):
	Drawable(menu._color)
{
	for (int i = 0; i < menu._mItems.size(); ++i)
	{
		if (menu._mItems[i] != NULL)
		{
			MenuItem* mItem = new MenuItem(*(menu._mItems[i]));
			_mItems.push_back(mItem);
		}
	}
	_width = menu._width;
	_height = menu._height;
	_xStart = menu._xStart;
	_yStart = menu._yStart;
}

void Menu::AddMenuItem(const int & number, MenuItem* mItem ...)
{
	if (mItem == NULL) return;
	va_list v1;
	va_start(v1, mItem);
	MenuItem* current = NULL;
	if (mItem != NULL) current = mItem->Clone();
	for (int i = 0; i < number; ++i)
	{
		if (current != NULL) _mItems.push_back(current->Clone());
		current = va_arg(v1,MenuItem*);
	}
	va_end(v1);
	return;
}

void Menu::AddMenuItem(MenuItem* mItem, const int & ordinal_number)
{
	_mItems.insert(_mItems.begin() + ordinal_number,1,mItem->Clone());
	return;
}

void Menu::AddMenuItem(const int & ordinal_number, string text)
{
	MenuItem* mItem = new MenuItem();
	mItem->AddText(text);
	_mItems.insert(_mItems.begin() + ordinal_number,1,mItem->Clone());
	delete mItem;
	mItem = NULL;
	return;
}

void Menu::StandartOrder()
{
	int n = _mItems.size();
	double height = _height/(double)n - 25.0;
	if (height > 80.0) height = 80.0;
	double width = _width - 40.0;
	double xStart = _xStart + 20.0;
	double yStart = _yStart + 20.0;
	for(int i = 0; i < n; ++i)
	{
		_mItems[i]->SetWidth(width);
		_mItems[i]->SetHeight(height);
		_mItems[i]->SetCoords(xStart, yStart + (n-1-i)*height + (n-1-i)*20.0);
	}
	return;
}

void Menu::AddToStandartOrder(const string & text, const int & ordinal_number)
{
	if (ordinal_number < 0 && ordinal_number != -1) return;
	MenuItem* mItem = new MenuItem();
	mItem->AddText(text);
	if (ordinal_number == -1 || ordinal_number >= _mItems.size()) _mItems.push_back(mItem->Clone());
	else _mItems.insert(_mItems.begin() + ordinal_number, 1, mItem->Clone());
	delete mItem;
	this->StandartOrder();
}

void Menu::clear()
{
	if (_mItems.size() == 0) return;
	while(!_mItems.empty())
	{
		if (_mItems[_mItems.size() - 1] != NULL) delete _mItems[_mItems.size() - 1];
		_mItems.pop_back();
	}
	return;
}

void Menu::Draw()
{
	glClearColor(0,0,0,0);
	glColor3f(_color.RED, _color.GREEN, _color.BLUE);
	glLineWidth(_size);
	glPushMatrix();
	glLoadIdentity();
	glRectd(_xStart,_yStart, _xStart + _width, _yStart + _height);
	glPopMatrix();
	if (_mItems.size() == 0) return;
	for (int i = 0; i < (int)_mItems.size(); ++i)
		if (_mItems[i] != NULL) _mItems[i]->Draw();
	return;
}

Menu* Menu::Clone()
{
	Menu* menu = new Menu(*this);
	return menu;
}

Menu::~Menu()
{
	if (_mItems.size() != 0)
	{
		while(!_mItems.empty())
		{
			if (_mItems[_mItems.size() - 1] != NULL) delete _mItems[_mItems.size() - 1];
			_mItems.pop_back();
		}
	}
}

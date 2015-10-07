#ifndef __DRAWING__
#define __DRAWING__

#include <queue>
#include <string.h>
#include "GeomFigures.h"

using namespace std;

class Menu;

class DrawEvent
{
	private:
		vector<Drawable*>* _objects;
	public:
		DrawEvent() { _objects = NULL; };
		DrawEvent(vector<Drawable*>* objects);
		DrawEvent(DrawEvent & draw);
		int NumberOfElements() {if (_objects != NULL) return _objects->size(); else return 0;};
		void Draw();
		void add(int number, Drawable * drawable_obj ...);
		void add(vector<Drawable*>* objects);
		void add(vector<Point2D>* objects);
		void add(vector<Point2D>* objects, Color& color);
		void add(DrawEvent* drawing);
		Drawable* last();
		void del(int index);
		void clear();
		~DrawEvent();
};

class Window
{
	private:
		Window();
		Window(int width, int height);
		Window(int width, int height, Color back_color);
		static Window* _instance;
	public:
		int Width;
		int Height;
		int StartX;//origin
		int StartY;//origin
		bool DrawVertexes = true;
		bool pause;
		int rendering;
		bool EnableMouseInput;
		Color BackColor;

		enum MainMenu{PutByMouse = 0, ReadFromFile, ReadMultidimVertexes, GenRand, RublevSteps, RublevResult,
					  Khachijan, GeneticMatrix, GeneticPetunin, PetuninResult, PetuninSteps, Clear, Exit,
					  RenderingInc, Rendering, RenderingDec};
		Menu* mainMenu;

		queue<DrawEvent*>* DrawEvents;
		DrawEvent* CurrentWindowContent;
		vector<Point2D>* vertexes;
		vector<Point>* MultiDimVert;//n-dimensional points

		static Window* Create();
		static Window* Create(int width, int height);
		static Window* Create(int width, int height, Color back_color);

		void SetStandartMainMenu();
		void AddNewDrawEvent(DrawEvent* draw_event);
		void AddToLastScene(DrawEvent* draw_event);
		void AddNext(DrawEvent * draw_event);
		DrawEvent* LastDrawEvent();
		DrawEvent* FirstDrawEvent();
		int NumberOfDrawEvents() {if (DrawEvents != NULL) return DrawEvents->size(); else return 0; };
		void PopDrawEvent();
		void Translate(const double & x, const double & y);
		void clear();
		~Window();
};

extern Window* window;

#define WIDTH 400.0
#define HEIGHT (window->Height)
#define XSTART (-window->StartX + window->Width - WIDTH)
#define YSTART (window->StartY)
#define YELLOW (Color(1.0f,1.0f,(float)(204.0/255.0)))
#define SKYBLUE (Color(0.0f,(float)(204.0/255.0),1.0f))

class MenuItem: public Drawable
{
	private:
		double _xStart;//left-bottom
		double _yStart;//
		double _width;
		double _height;
		string _text;

	public:
		MenuItem();
		MenuItem(const double & xCoord, const double & yCoord, const double & width, const double & height, Color color = SKYBLUE);
		MenuItem(const double & xCoord, const double & yCoord, const double & width, const double & height, const string & text = "", Color color = SKYBLUE);
		MenuItem(const MenuItem & mItem);
		void Draw();
		MenuItem* Clone();
		bool In(const double & xCoord, const double & yCoord);
		void AddText(const string & text);
		void SetWidth(const double & width) { _width = width; }
		void SetHeight(const double & height) { _height = height; }
		void SetCoords(const double & x, const double & y) { _xStart = x, _yStart = y; }
		~MenuItem();
};

class Menu: public Drawable
{
	private:
		vector<MenuItem*> _mItems;
		double _width;
		double _height;
		double _xStart;//left-bottom
		double _yStart;//
	public:
		Menu(const double & width = WIDTH, const double & height = HEIGHT, const double & xStart = XSTART, const double & yStart = YSTART, Color color = YELLOW);
		Menu(vector<MenuItem*> mItems, const double & width = WIDTH, const double & height = HEIGHT, const double & xStart = XSTART, const double & yStart = YSTART, Color color = YELLOW);
		Menu(const Menu & menu);
		void AddMenuItem(const int & number, MenuItem* m_item ...);//add some amount of menu items
		void AddMenuItem(MenuItem* mItem, const int & ordinal_number);//add menu item in special ordinal number
		void AddMenuItem(const int & ordinal_number, string text);//add menu item with special ordinal number and text
		void AddToStandartOrder(const string & text, const int & ordinal_number = -1);
		MenuItem* operator[](const int & i) { return _mItems[i]; }
		MenuItem* at(const int & i) { return _mItems[i]; }
		void clear();
		void Draw();
		Menu* Clone();
		void StandartOrder();
		double Width() const { return _width; }
		double Height() const {return _height; }
		double XStart() const {return _xStart; }
		double yStart() const { return _yStart; }
		~Menu();
};

#endif
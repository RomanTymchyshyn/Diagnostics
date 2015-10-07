#include <fstream>
#include <time.h>
#include <Windows.h>
#include "Diagnostics.h"
#include <iostream>
#include "Solvers/EllipsoidBuilders.h"
#include "MouseHandler.h"

Window* window;

void glut_init(int argc, char *argv[], int window_width, int window_height, int window_pos_x, int window_pos_y);
void tmf(int value);
void Display();
void Reshape(int w, int h);
void Keyboard(unsigned char key, int x, int y);

void main(int argc, char *argv[])
{
	srand(time(NULL));
	Helper::Read("Resources/Multidimensional Input.txt");
	glut_init(argc, argv, 1350, 700, 0, 0);
	glutCreateWindow("Ellipse of minimal volume");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutTimerFunc(window->rendering, tmf, 0);
	glutKeyboardFunc(Keyboard);
	glutMouseFunc(MouseHandler::HandleMouseClick);
	glutMainLoop();
}

void glut_init(int argc, char *argv[], int window_width, int window_height, int window_pos_x, int window_pos_y)
{
	window = Window::Create(window_width, window_height);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(window_width, window_height);
	glutInitWindowPosition(window_pos_x, window_pos_y);
	ofstream fo("Resources/Output Data.txt",ios::trunc);
	fo<<"<<<<<<<<<Minimal volume ellipsoid>>>>>>>>>>>>>>\n";
	fo.close();
	window->SetStandartMainMenu();
	return;
}

void tmf(int value) // Timer function
{
	if (window == NULL) return;
	if (!window->pause)
	{	
		if (window->DrawEvents != NULL) 
		{
			if(!window->DrawEvents->empty())
			{
				delete window->CurrentWindowContent;
				window->CurrentWindowContent = window->FirstDrawEvent();
			}
			if (window->NumberOfDrawEvents() > 1)
				window->PopDrawEvent();
		}
	}
	glutPostRedisplay();// Redraw windows
	glutTimerFunc(window->rendering, tmf, 0);// Setup next timer
}

void Display()
{
	glClearColor(window->BackColor.RED, window->BackColor.GREEN, window->BackColor.BLUE, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	glLineWidth(1.0);
	if (window != NULL)
	{
		if (window->CurrentWindowContent != NULL) window->CurrentWindowContent->Draw();
		if (window->vertexes != NULL && window->DrawVertexes)
		{
			for (int i = 0; i < window->vertexes->size(); ++i)
				window->vertexes->at(i).Draw();
		}
		if (window->mainMenu != NULL) window->mainMenu->Draw();
	}
	glFlush();
	glutSwapBuffers();
};

void Reshape(int w, int h)
{
	if (window == NULL) window = Window::Create(w, h);
	window->Width = w;
	window->Height = h;

	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, -1.0, 1.0);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
};

void Keyboard(unsigned char key, int x, int y)
{
	#define ESCAPE 27
	#define SPACE 32
	#define ENTER 13
	
	switch (key)
	{
		case ESCAPE:
			if (window != NULL) 
			{
				delete window;
				window = NULL;
			}
			exit(0);
			break;
		case SPACE://pause
			if (window != NULL) window->pause = !window->pause;
			break;
		case ENTER://next move
			if (window->DrawEvents != NULL) 
			{
				if(!window->DrawEvents->empty())
				{
					delete window->CurrentWindowContent;
					window->CurrentWindowContent = window->FirstDrawEvent();
				}
				if (window->NumberOfDrawEvents() > 1)
					window->PopDrawEvent();
			}
			break;
	}
};
/*
    Copyright 2012 thedukeofzill

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _ENGINE
#define _ENGINE

#include <stdio.h>
#include <stdlib.h>
#include <curses.h>
#include <locale.h>
#include <unistd.h>
#include <math.h>

int row,col,r,c,ccount;
short** data;
double** zbuffer;
WINDOW * mainwin;
int ch;
int wireframe;

int stagesize; //how much we want visible in map
int dist;  //eye point
int zval;  //image plane
int maxDepth;

int camx;
int camy;
int camz;
int camrx;
int camry;
int camrz;

double framerate;

typedef struct {
    double x;
    double y;
} vertex;

vertex point(double xx,double yy);

typedef struct {
    double x;
    double y;
    double z;
} vertex3;

vertex3 point3(double x,double y,double z);

typedef struct {
    int size;
    short color;
    vertex* vertices;
} poly;

poly pol(vertex* v, short cr, int s);

typedef struct {
    vertex3* v3;
    vertex3 pcenter;
    vertex* v2;
    int length;
    int size;
    short color;
} shape;

typedef struct{
    vertex3* v3;
    short color;
} triangle;

triangle tri(vertex3* v3, short color);

shape makeshape(vertex3* v3, int len, int s, short color);

void swap(int* a,int* b);

void drawLine(double xx1,double yy1,double xx2,double yy2, short color);

void drawPoly(poly p);

short** make2dshortarray(short arraySizeX, short arraySizeY);

double** make2ddoublearray(int arraySizeX, int arraySizeY);

void resetColors();

void resetZBuffer();
 
void setcolor(int bg, int fg);

void drawBox(double x,double y,double l,double w,short color);

void drawSolidBox(double x,double y,double l,double w,short color);

void startengine();

void killengine();

void drawTriangle(triangle t);

void drawscreen();

vertex3* cubevertices(double x,double y, double z, int size);

shape cube(double x,double y, double z, int size, short color);

int visible(vertex3* array, int size);

void translatecube(shape* sh);

void drawCube(shape* sh,short wall);

void rotateaxis(vertex3* v, int size, double xr, double yr, double zr);

void rotateinplace(shape* sh, double xr, double yr, double zr);

int tritest(double x0,double y0,double x1,double y1,double x2,double y2,double px, double py);

int test2(double px,double py,double m,double b);

int test1(double px,double py,double m,double b,double lx,double ly);


#endif

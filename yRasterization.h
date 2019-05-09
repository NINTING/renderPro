#pragma once
#include<iostream>
#include<algorithm>
#include<fstream>
const int Width = 526, Height = 526;
//#define vertex2(x,y) vertex(x,y,0)

#define R 1
#define T 2
#define L 4
#define B 8
const int Lx = 0, Rx = Width, Ty = Height, By = 0;
typedef
struct Vector{
	float x, y,z;
	Vector() {};
	Vector(float x, float y,float rz=0):x(x), y(y) {
		z = rz;
	};

}vertex,Vector;

//typedef
//struct vector2 {
//	float x, y;
//	vector2() {};
//	vector2(float x, float y) :x(x), y(y) {};
//
//}vertex2,vector2;


Vector VecCross(const Vector& a, const Vector& b);

int getcode(const vertex &v);
void clip(int& code, vertex& v, const vertex &u);
void lineclip(vertex& a, vertex& b);

void setPixel(int x, int y,unsigned char  color = 0);
void fillTriangle1(vertex v0,vertex v1,vertex v2);

void fillTriangle2(const vertex& v0, const vertex& v1,const vertex& v2);
Vector barycentric(const vertex& v0, const vertex &v1, const vertex& v2, const vertex& v);

void bresenham(const vertex &begin, const vertex& end);
void bresenham(int x1, int y1, int x2, int y2);
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

typedef
struct Vector4 {
	float _x, _y, _z,_w;
	Vector4() {};
	Vector4(float x, float y, float z,float w) :_x(x), _y(y), _z(z),_w(w) {};



}vertex4, Vector4;



class Matrix {
public:
	Matrix(float x11,float x12,float x13,float x14,
			float x21,float x22 ,float x23,float x24,
			float x31,float x32,float x33,float x34,
			float x41,float x42,float x43,float x44);

	float& operator () (int i, int j) { return value[i][j]; }
	float operator () (int i, int j)const { return value[i][j]; }

	


	float value[4][4];

};
Matrix operator * (const Matrix &lrh, const Matrix &rsh);
Matrix& operator *= (Matrix &lrh, Matrix &rsh);

Vector4 operator * (const Vector4 &lrh, const Matrix &rsh);
Vector4& operator *= (vertex4 &lrh, Matrix &rsh);


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
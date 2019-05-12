#pragma once
#include<iostream>
#include<algorithm>
#include<fstream>
#include<cstring>
const int Width = 600, Height = 600;
//#define vertex2(x,y) Vector4(x,y,0)

#define R 1
#define T 2
#define L 4
#define B 8
const int Lx = 0, Rx = Width, Ty = Height, By = 0;

//---------------------------------
//       向量声明				  |
//---------------------------------

#define indeiceBuffer int

typedef
struct Vector4 {
	float _x, _y, _z,_w;
	Vector4() {
		_x = 0.0f, _y = 0.0f, _z = 0.0f,_w = 0.0f;
	};
	Vector4(float x, float y, float z,float w) :_x(x), _y(y), _z(z),_w(w) {};


}vertex4, Vector4vertexBuffer;


struct Color
{
	unsigned char _r, _g, _b;
	Color(unsigned char r, unsigned char g, unsigned char b):_r(r),_g(g),_b(b){}
};
const Color Black(0, 0, 0);
//---------------------------------
//       矩阵声明				  |
//---------------------------------

class Matrix {
public:
	Matrix(float x11,float x12,float x13,float x14,
			float x21,float x22 ,float x23,float x24,
			float x31,float x32,float x33,float x34,
			float x41,float x42,float x43,float x44);
	Matrix() {
		memset(value, 0.0f, sizeof(value));
	};

	float& operator () (int i, int j) { return value[i][j]; }
	float operator () (int i, int j)const { return value[i][j]; }
	
	float value[4][4];

};
Matrix getIdentity();

Matrix operator * (const Matrix &lrh, const Matrix &rsh);
Matrix& operator *= (Matrix &lrh, Matrix &rsh);

Vector4 operator * (const Vector4 &lrh, const Matrix &rsh);
Vector4& operator *= (vertex4 &lrh, Matrix &rsh);

void VecCross(Vector4 *out,const Vector4& a, const Vector4& b);

//typedef
//struct vector2 {
//	float x, y;
//	vector2() {};
//	vector2(float x, float y) :x(x), y(y) {};
//
//}vertex2,vector2;



//---------------------------------
//       	直线裁剪声明		  |
//---------------------------------


int getcode(const Vector4 &v);
void clip(int& code, Vector4& v, const Vector4 &u);
void lineclip(Vector4 a, Vector4 b);



void setPixel(int x, int y,Color color);

//---------------------------------
//       光栅化三角形			  |
//---------------------------------

void fillTriangle1(Vector4 v0,Vector4 v1,Vector4 v2);

//void fillTriangle2(const Vector4& v0, const Vector4& v1,const Vector4& v2);
//Vector4 barycentric(const Vector4& v0, const Vector4 &v1, const Vector4& v2, const Vector4& v);

//---------------------------------
//       bresenham画线			  |
//---------------------------------

void bresenham(const Vector4 &begin, const Vector4& end);
void bresenham(int x1, int y1, int x2, int y2);
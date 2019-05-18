#pragma once
#include<iostream>
#include<algorithm>
#include<fstream>
#include<cstring>
#include<vector>
const int Width =600, Height = 600;
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

struct Color
{
	unsigned char _r, _g, _b;
	Color() = default;
	Color(unsigned char r, unsigned char g, unsigned char b) :_r(r), _g(g), _b(b) {}
};
Color operator +(const Color& a, const Color&b);
Color& operator +=(Color& a, const Color&b);
Color operator *(const Color& a, const Color&b);
Color& operator *=(Color& a, const Color&b);
Color operator *(const Color& a, float t);
Color& operator *=(Color& a, float t);


const Color Black(0, 0, 0);
const Color White(255, 255, 255);
const Color Red(255, 0, 0);
const Color Yellow(0, 255, 0);
const Color Blue(0, 0, 255);


typedef
struct Vector4 {
	float _x, _y, _z,_w;
	Color _c;
	float rhw;
	Vector4() {
		_x = 0.0f, _y = 0.0f, _z = 0.0f,_w = 0.0f;
		_c = Black; rhw = _z;
	};
	
	Vector4(float x, float y, float z, float w) :_x(x), _y(y), _z(z), _w(w) { _c = Black; rhw = _z; };
	Vector4(float x, float y, float z, float w, Color c) :_x(x), _y(y), _z(z), _w(w), _c(c) { rhw = _z; };
	Vector4(float x, float y, float z, float w, Color c,float rhw) :_x(x), _y(y), _z(z), _w(w), _c(c),rhw(rhw) { };

}vertex4, Vector4;

Vector4 operator - (const Vector4& lsh, const Vector4& rsh);

Vector4& operator -= (Vector4& lsh, const Vector4& rsh);

Vector4 operator + (const Vector4& lsh, const Vector4& rsh);

Vector4& operator += (Vector4& lsh, const Vector4& rsh);

float operator * (const Vector4& lsh, const Vector4& rsh );

void VecNormalize(Vector4 & out,const Vector4& in);




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
	const float& operator () (int i, int j)const { return value[i][j]; }
	
	float value[4][4];

};
Matrix getIdentity();

Matrix operator * (const Matrix &lrh, const Matrix &rsh);
Matrix& operator *= (Matrix &lrh, const Matrix &rsh);

Vector4 operator * (const Vector4 &lrh, const Matrix &rsh);
Vector4& operator *= (vertex4 &lrh, Matrix &rsh);

void VecCross(Vector4 *out,const Vector4& a, const Vector4& b);

void MatrxIdentity(Matrix &out);

//旋转矩阵
void RotationAxis(Matrix&out, const Vector4& axis, float theta);
void RotationAxisX(Matrix& out, float theta);
void RotationAxisY(Matrix& out, float theta);
void RotationAxisZ(Matrix& out, float theta);

//投影矩阵
void PerspectiveMatrix(Matrix * Pm);		//不好使
void PerspectiveMatrix(Matrix * Pm, float fov, float Aspect, float near, float far);
//视图变换矩阵

void viewMatrix(Matrix&out,const Vector4& at,const Vector4& view,const Vector4& up);


//映射矩阵
void toScreen(Vector4 *p);

//矩阵应用
void MatrixApply(Vector4 *out, const Vector4 &v, const Matrix &m);


//背面消隐
bool backCull(const vertex4& a, const vertex4& b,const vertex4& c);

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
void lineclip(Vector4 &a, Vector4 &b);



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


//---------------------------------
//       线性插值				  |
//---------------------------------

float interp(float a, float b, float t);
void vertexInterp(vertex4 *out, const vertex4 &v0, const vertex4 &v1, float t);
//void Zslerp(const vertexArr *va1, const vertexArr *va2, vertexArr *out, float t);

//---------------------------------
//       流水线处理 			  |
//---------------------------------

void IndeicesProcessPipeline(std::vector<vertex4>* outlist, const indeiceBuffer* ib, const std::vector<vertex4>&vb);

void vertexProcessPipeline(std::vector<vertex4>* outlist, const std::vector<vertex4>& list);
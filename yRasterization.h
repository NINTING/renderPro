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
extern char img[(Width + 10)*(Height + 10) * 3];


//顶点格式
#define FVFcolor 1
#define FVFtexture 2
#define FVFwireframe 4


#define indeiceBuffer int



//---------------------------------
//      颜色声明				  |
//---------------------------------

struct Color
{
	float _r, _g, _b;
	Color() = default;
	Color(float r, float g, float b) :_r(r), _g(g), _b(b) {}
};
Color operator +(const Color& a, const Color&b);
Color& operator +=(Color& a, const Color&b);
Color operator /(const Color& a, const Color&b);
Color& operator /=(Color& a, const Color&b);
Color operator *(const Color& a, float t);
Color operator *(float t,const Color& a);
Color& operator *=(Color& a, float t);
Color operator *(const Color& a,const Color &b);
Color& operator *=(Color& a, const Color b );


const Color Black(0, 0, 0);
const Color White(1, 1, 1);
const Color Red(1, 0, 0);
const Color Green(0, 1, 0);
const Color Blue(0, 0, 1);
//---------------------------------
//       向量声明				  |
//---------------------------------

#define TS_WORLD 1
#define TS_VIEW 2
#define TS_PROJECTION 4

typedef
struct Vector4
{
	float _x, _y, _z,_w;
	Vector4() = default;
	Vector4(float x, float y, float z,float w) :_x(x), _y(y), _z(z),_w(w) {};

}Vertex4, Vector4;

typedef
struct VertexAtrr {
	Vector4 _v;
	//float _x, _y, _z,_w;

	//vertex attribute
	Color _c;
	float rhw;		// 1/z
	float _tu, _tv;
	float nx, ny, nz;
	VertexAtrr() = default;
	//VertexAtrr(): _v(Vector4(0, 0, 0, 0)){};
	
	VertexAtrr(float x, float y, float z, float w) :_v(Vector4(x,y,z,w)) { _c = White; };
	VertexAtrr(float x, float y, float z, float w, Color c) :_v(Vector4(x, y, z, w)), _c(c) { };
	VertexAtrr(float x, float y, float z, float w, Color c, float u, float v): _v(Vector4(x, y, z, w)), _c(c), _tu(u), _tv(v) { };
	VertexAtrr(float x, float y, float z, float w,float u,float v) :_v(Vector4(x, y, z, w)),_tu(u),_tv(v) { _c = White; };

};

Vector4 operator - (const Vector4& lsh, const Vector4& rsh);
Vector4 operator - (const Vector4& rsh);
Vector4& operator -= (Vector4& lsh, const Vector4& rsh);

Vector4 operator + (const Vector4& lsh, const Vector4& rsh);

Vector4& operator += (Vector4& lsh, const Vector4& rsh);

float operator * (const Vector4& lsh, const Vector4& rsh );

void VecNormalize(Vector4 & out,const Vector4& in);

void getNormal(Vector4 *out, const Vertex4& a, const Vertex4& b, const Vertex4& c);


//---------------------------------
//       材质声明				  |
//---------------------------------

struct Matreial {
	Color _Ambient;		//环境光
	Color _Specular;		//镜面光
	Color _Diffuse;		//漫反射
	float _power;		//控制镜面反射的光斑
	//Color emition
	Matreial() {};
	Matreial(Color Ambient,Color Specular,Color Diffuse,float power):
		_Ambient(Ambient),_Specular(Specular),_Diffuse(Diffuse),_power(power){}

};

void setMatreial(Matreial *Mtr);

//---------------------------------
//       三角类					  |
//---------------------------------



struct Triangle {
	VertexAtrr v0,v1,v2;
	Vector4 normal;
	Matreial _matreial;
	Triangle(VertexAtrr a, VertexAtrr b, VertexAtrr c) :v0(a), v1(b), v2(c) {
		getNormal(&normal, v0._v, v1._v, v2._v);
	}

};

void getNormal(Triangle * tri);

//---------------------------------
//       光照相关				  |
//---------------------------------


enum LightType
{
	LightSpot = 1,
	LightPoint = 2,
	LightDirectional = 3,
};

struct Light {
	Light() {};
	Vector4 position;		//point ,spot light
	Vector4 direction;
	Color _Ambient;		    //环境光
	Color _Specular;		//镜面光
	Color _Diffuse;		    //漫反射
	LightType type;
	float Attenuation0;
	float Attenuation1;
	float Attenuation2;
	float range;
};

Light* Lightarr();
void setLight(int num,Light*light);
Light* getLight(int index);
void releaseLight();

Light initDirectionalLight(const Vector4& direction,const Color &color);
Light initPointLight(Vector4 direction, Vector4);
Light initSpotLight(Vector4 direction, Vector4);


void getSpecuC(Color *out, const Matreial& M,const Light &Light,const Vector4& view,const Vector4 &normal);
void getAmbientC(Color *out, const Matreial& M, const Light &Light);
void getdiffuseC(Color *out, const Matreial& M, const Light &Light, const Vector4 &normal);

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
Vector4& operator *= (Vector4 &lrh, Matrix &rsh);

void VecCross(Vector4 *out,const Vector4& a, const Vector4& b);

void MatrixIdentity(Matrix &out);
void MatrixZero(Matrix &out);
//平移矩阵

void MatrixTranslation(Matrix& out,float x,float y,float z);


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
bool backCull(const Vertex4& a, const Vertex4& b,const Vertex4& c);
bool backCull(const Triangle& tri);
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


int getcode(const Vertex4 &v);
void clip(int& code, VertexAtrr& va, const VertexAtrr & ua);
void lineclip(VertexAtrr &a, VertexAtrr &b);



void setPixel(int x, int y,Color color);

//---------------------------------
//       光栅化三角形			  |
//---------------------------------

void fillTriangle1(VertexAtrr v0, VertexAtrr v1, VertexAtrr v2);
void fillTriangle1(Triangle& tri);
void setPixel(int x, int y, Color color);
void setPixel(const Vertex4& v);

//void fillTriangle2(const Vector4& v0, const Vector4& v1,const Vector4& v2);
//Vector4 barycentric(const Vector4& v0, const Vector4 &v1, const Vector4& v2, const Vector4& v);

//---------------------------------
//       bresenham画线			  |
//---------------------------------

void bresenham(const VertexAtrr &begin, const VertexAtrr& end);
void bresenham(int x1, int y1, int x2, int y2);


//---------------------------------
//       线性插值				  |
//---------------------------------

float interp(float a, float b, float t);
void vertexInterp(VertexAtrr *out, const VertexAtrr &v0, const VertexAtrr &v1, float t);
//void Zslerp(const vertexArr *va1, const vertexArr *va2, vertexArr *out, float t);

//---------------------------------
//       流水线处理 			  |
//---------------------------------

void IndeicesProcessPipeline(std::vector<Vertex4>* outlist, const indeiceBuffer* ib, const std::vector<Vertex4>&vb);

void IndeicesProcessPipeline(std::vector<Triangle>* outlist, const indeiceBuffer* ib, const std::vector<Vertex4>&vb);

void vertexProcessPipeline(std::vector<Vertex4>* outlist, const std::vector<Vertex4>& list);



//---------------------------------
//       文件处理	 			  |
//---------------------------------
void getTexture(const char* file);

void draw(const char *name);


//---------------------------------
//			纹理设置			   |
//---------------------------------

struct Texture {
	int width, height;
	float *tex;
	~Texture() {
		delete[]tex;
	}
};
void getTexPixel(const Texture *Texobj, Color* c, float x, float y);

void Texture_set(Texture *texobj, void *texImg, int w, int h, int pitch);

void RsetFVF(int FvF);
void RsetTex(Texture *tex);

void init_Texture(Texture* texobj, int w, int h);


//---------------------------------
//			Z-buffer			   |
//---------------------------------



void setZbuffer(int w,int h);
void releaseZbuffer();
void writeZbuffer(int x, int y,float z);
float getZbuffer(int x, int y);



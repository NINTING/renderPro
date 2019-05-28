#include"yRasterization.h"
#include<vector>
#include<iostream>
#include<cstdio>
#define PI 3.14159265359f


using namespace std;
void releaseall();

int LCount;

Matrix tranMatrix;
Matrix ProjectMatr, ViewMatr, WorldMatr;

void setMatrix(int flag, Matrix * m) {
	if (flag & TS_WORLD)
		WorldMatr = *m;
	if (flag & TS_VIEW)
		ViewMatr = *m;
	if (flag & TS_PROJECTION)
		ProjectMatr = *m;
}

void init() {
	MatrixIdentity(ProjectMatr);
	MatrixIdentity(ViewMatr);
	MatrixIdentity(WorldMatr);
}



void setup() {
	//��������
	init();
	memset(img,0, sizeof(img));
	Texture tex;
	init_Texture(&tex, 255, 255);
	//RsetFVF(FVFtexture);
	RsetFVF(FVFtexture);
	RsetTex(&tex);
	PresentTex();
	Matrix Vm, Pm;
	PerspectiveMatrix(&Pm, PI * 0.5, Width / Height, 1, 500);
	Vector4 at(-3.0f, 0.0f, 1.0f, 1), view(0.0f, 0.0f, 0.0f, 1), up(0, 1,0 , 0);
	viewMatrix(Vm, at, view, up);
	setMatrix(TS_VIEW, &Vm);
	setMatrix(TS_PROJECTION, &Pm);
	
	//Zbuffer��ʼ��
	setZbuffer(Width, Height);
	
	//����
	Vector4 direction(1, 1, 0, 0);
	Color Lc(1, 1, 1);
	Light direL = initDirectionalLight(direction, Lc);
	setLight(0, &direL);
	LCount++;


	//����
	Color whiteM(1, 1, 1);
	Matreial matr1(whiteM,whiteM,whiteM,1.0f);
	setMatreial(&matr1);

}

void fillMesh(std::vector<VertexAtrr>& list) {
	size_t size =list.size();
	for (int i = 0; i < size; i+=3) {
		fillTriangle1(list[i], list[i+1], list[i+2]);
	}
}

void fillMesh(std::vector<Triangle>& list) {
	size_t size = list.size();
	for (int i = 0; i < size; i ++) {
		fillTriangle1(list[i]);
	}
}

void wireFrame(std::vector<VertexAtrr> list) {
	for (int i = 0; i < list.size(); i +=3) {
		lineclip(list[i], list[i + 1]);
		bresenham(list[i], list[i + 1]);
		
		lineclip(list[i], list[i + 2]);
		bresenham(list[i], list[i + 2]);
		
		lineclip(list[i+1], list[i + 2]);
		bresenham(list[i+1], list[i + 2]);
	}
}



////Ĭ��ͶӰ�����ԭ��Ϊ1 ��d = 1��
//void PerspectiveMatrix(Matrix * Pm) {
//	float Right = Width*0.5f, Left = -Right;	
//	float Top = Height*0.5f, Bottom = -Top;
//	const float near = 1.0f, far = 500.0f;
//	float a = (float)far / (float)(far - near), b = -(float)(far*near) / (float)(far - near);
//	 (*Pm)(0, 0) = 2 * near / (Right - Left);
//	 (*Pm)(1, 1) = 2 * near / (Top - Bottom);
//	 (*Pm)(2, 0) = (float)(Left + Right) / (float)(Left - Right)*1.0f;
//	 (*Pm)(2, 1) = (float)(Bottom + Top) / (float)(Bottom - Top)*1.0f;
//	 (*Pm)(2, 2) = a;
//	 (*Pm)(3, 2) = b;
//	 (*Pm)(2, 3) = 1.0f;
//}


void PerspectiveMatrix(Matrix * Pm, float fov,float Aspect,float near,float far) {
	float cot = 1 / (float)tanf(fov*0.5f);

	(*Pm)(0, 0) = cot / Aspect;
	(*Pm)(1, 1) = cot;
	(*Pm)(2, 2) = far / (far - near);
	(*Pm)(3, 2) = -far*near / (far - near);
	(*Pm)(2, 3) = 1.0f;
}


void toScreen(Vector4 *p){
	
	int X = 0, Y = 0;
	float dw = 1.0f / p->_w;
	p->_x = (p->_x*dw + 1.0f)*Width *0.5f+X;
	p->_y = (-p->_y*dw + 1.0f)*Height*0.5f+Y;
	p->_z = p->_z*dw;
	p->_w = 1.0f;
}


void MatrixApply(Vector4 *out,const Vector4 &v,const Matrix &m) {
	*out = v * m;
}

void applyForIndeices(int** varr, size_t size) {
	(*varr) = new int[size];
}



template<class t>
void release(t * arr) {
	delete []arr;
	arr = 0;
}



//�������˳��--���������*�ӽǾ���*ͶӰ����



void getTranMatrix() {
	tranMatrix = WorldMatr * ViewMatr;
}

void IndeicesProcessPipeline(vector<Triangle>* outlist,const indeiceBuffer* ib, const vector<VertexAtrr>&vb,int TriangleNum) {

	getTranMatrix();
	//createTriMesh()
	size_t size = TriangleNum*3;
	
	VertexAtrr a, b, c;
	
	for (int i = 0; i < LCount; i++) {
		MatrixApply(&(getLight(i)->direction), getLight(i)->direction, ViewMatr);
	}
		
	for (int i = 0; i <size; i += 3) {
		
		MatrixApply(&a._v, vb[ib[i]]._v, tranMatrix);
		
		MatrixApply(&b._v, vb[ib[i + 1]]._v, tranMatrix);
		
		MatrixApply(&c._v, vb[ib[i + 2]]._v, tranMatrix);
	
		//��������
		Triangle tri(a, b, c);
		getNormal(&tri);
		setTriNormal(tri,tri.normal);
		
		float rhw;
		if (backCull(tri)) {
			VertexShader(tri.v0,LCount);
			VertexShader(tri.v1,LCount);
			VertexShader(tri.v2,LCount);

			MatrixApply(&tri.v0._v , tri.v0._v, ProjectMatr);
			tri.v0.rhw = 1.0f / tri.v0._v._w;
			AttrMulRhw(tri.v0, vb[ib[i]], tri.v0.rhw);
			
			MatrixApply(&tri.v1._v, tri.v1._v, ProjectMatr);
			tri.v1.rhw = 1.0f / tri.v1._v._w;
			AttrMulRhw(tri.v1, vb[ib[i+1]], tri.v1.rhw);
			
			MatrixApply(&tri.v2._v, tri.v2._v, ProjectMatr);
			tri.v2.rhw = 1.0f / tri.v1._v._w;
			AttrMulRhw(tri.v2, vb[ib[i+2]], tri.v2.rhw);
			
			outlist->push_back(tri);

		}
	}


	size = outlist->size();
	for (int i = 0; i < size; i++) {
		toScreen(&(*outlist)[i].v0._v);
		toScreen(&(*outlist)[i].v1._v);
		toScreen(&(*outlist)[i].v2._v);
	}
	
}

//void vertexProcessPipeline(vector<VertexAtrr>* outlist,const vector<VertexAtrr>& list) {
//	Matrix Pm,Vm , tranMatrix;
//	
//	PerspectiveMatrix(&Pm);
//	MatrxIdentity(Vm);
//	Vector4 at(0, 0,-150.0f, 1), view(0, 0, 1, 0), up(0, 1,0 , 0);
//
//	viewMatrix(Vm, at, view, up);
//	tranMatrix =  Vm*Pm;
//	(*outlist).resize(list.size());
//	for (int i = 0; i < list.size(); i += 3) {
//		MatrixApply(&(*outlist)[i], list[i], tranMatrix);
//		MatrixApply(&(*outlist)[i+1], list[i + 1], tranMatrix);
//		MatrixApply(&(*outlist)[i+2], list[i + 2], tranMatrix);
//	}
//	for (int i = 0; i < list.size(); i++) {
//		toScreen(&(*outlist)[i]);
//	}
//}

void triangleTest() {
	
	VertexAtrr a(5.0f, 0, 0.0f, 1,Blue,0,0), b(0.0f, 5.0f, 0.0f, 1, Red,0,1), c(0.0f, -5.0f, 0.0f, 1, Green,1,0);
	VertexAtrr d(-5.0f, 0, 0.0f, 1, Blue,1,1);
	indeiceBuffer* ib = 0;
	vector<VertexAtrr>list{ a,b,c ,d};
	applyForIndeices(&ib, 2*3);
	ib[0] = 0, ib[1] = 1, ib[2] = 2;
	ib[3] = 2, ib[4] = 1, ib[5] = 3;
	vector<Triangle> outlist;

	IndeicesProcessPipeline(&outlist, ib, list, 2);
	//wireFrame(outlist);
	fillMesh(outlist);
	release(ib);
	draw("color10.ppm");
}
	
void pyramid()
{
	VertexAtrr a(1.0f, 0, 0.0f, 1, Red), b(0.0f, 0.0f, 1.0f, 1, Green), c(-1.0f, 0.0f, 0.0f, 1, Blue),
		d(0.0f, 0.0f, -1.0f, 1, Red), e(0.0f, 1.0f, 0.0f, 1.0f, White);

	indeiceBuffer* ib = 0;
	vector<VertexAtrr>list{ a,b,c,d,e };
	applyForIndeices(&ib, 4 * 3);
	
	//front
	ib[0] = 0, ib[1] = 1, ib[2] = 4;

	ib[3] = 1, ib[4] = 2, ib[5] = 4;
	
	ib[6] = 2, ib[7] = 3, ib[8] = 4;
	
	ib[9] = 3, ib[10] = 0, ib[11] = 4;

	vector<Triangle> outlist;
	
	IndeicesProcessPipeline(&outlist, ib, list, 4);
	//wireFrame(outlist);
	fillMesh(outlist);
	release(ib);
	draw("color2.ppm");

}


void Cubetest() {
	VertexAtrr v1(-1.0f, 1.0, -1.0f, 1, Red), v2(1.0f, 1.0, -1.0f, 1, Red),
		v3(1.0f, -1.0, -1.0f, 1, Red), v4(-1.0f, -1.0f, -1.0f, 1, Red);
	
	VertexAtrr v5(-1.0f, 1.0, 1.0f, 1.0f, Blue), v6(1.0f, 1.0, 1.0f, 1, Blue),
		v7(1.0f, -1.0, 1.0f, 1, Blue), v8(-1.0f, -1.0, 1.0f, 1, Blue);

	indeiceBuffer* ib = 0;
	vector<VertexAtrr>list{ v1,v2,v3,v4,v5,v6,v7,v8 };
	applyForIndeices(&ib, 12 * 3);
	//front
	ib[0] = 0, ib[1] = 3, ib[2] = 1;
	ib[3] = 3, ib[4] = 2, ib[5] = 1;

	//left
	ib[6] = 0, ib[7] = 4, ib[8] = 3;
	ib[9] = 3, ib[10] = 4, ib[11] = 7;

	//right
	ib[12] = 2, ib[13] = 5, ib[14] = 1;
	ib[15] = 2, ib[16] = 6, ib[17] = 5;

	//back
	ib[18] = 4, ib[19] = 5, ib[20] = 7;
	ib[21] = 7, ib[22] = 5, ib[23] = 6;

	//top
	ib[24] = 0, ib[25] = 1, ib[26] = 4;
	ib[27] = 1, ib[28] = 5, ib[29] = 4;

	//bottom
	ib[30] = 7, ib[31] = 6, ib[32] = 3;
	ib[33] = 6, ib[34] = 2, ib[35] = 3;


	vector<Triangle> outlist;																													
	IndeicesProcessPipeline(&outlist, ib, list, 12);
	//wireFrame(outlist);
	fillMesh(outlist);
	release(ib);
		//fillTriangle2(a, b, c);
		//wareFrame(list);
	draw("color7.ppm");
}

void textureCube() {

	setup();
	VertexAtrr v0(-1.0f, 1.0, -1.0f, 1,0,0), v1(1.0f, 1.0, -1.0f, 1, 0,1), 
		v2(1.0f, -1.0, -1.0f, 1, 1,1), v3(-1.0f, -1.0f, -1.0f, 1, 1,0);

	VertexAtrr v4(-1.0f, 1.0, 1.0f, 1.0f, 0,1),v5(1.0f,1.0, 1.0f,1,0,0),
		v6(1.0f, -1.0, 1.0f,1.0f,1.0f,0.0f), v7(-1.0f, -1.0, 1.0f, 1, 1.0f,1.0f);

	indeiceBuffer* ib = 0;
	vector<VertexAtrr>list{ v0,v1,v2,v3,v4,v5,v6,v7 };
	
	applyForIndeices(&ib, 12 * 3);
	//front
	ib[0] = 0, ib[1] = 3, ib[2] = 1;
	ib[3] = 3, ib[4] = 2, ib[5] = 1;

	//left
	ib[6] = 0, ib[7] = 4, ib[8] = 3;
	ib[9] = 3, ib[10] = 4, ib[11] = 7;

	//right
	ib[12] = 2, ib[13] = 5, ib[14] = 1;
	ib[15] = 2, ib[16] = 6, ib[17] = 5;

	//back
	ib[18] = 4, ib[19] = 5, ib[20] = 7;
	ib[21] = 7, ib[22] = 5, ib[23] = 6;

	//top
	ib[24] = 0, ib[25] = 1, ib[26] = 4;
	ib[27] = 1, ib[28] = 5, ib[29] = 4;

	//bottom
	ib[30] = 7, ib[31] = 6, ib[32] = 3;
	ib[33] = 6, ib[34] = 2, ib[35] = 3;

	
	vector<Triangle> outlist;
	IndeicesProcessPipeline(&outlist, ib, list, 12);
	//wireFrame(outlist);
	fillMesh(outlist);
	release(ib);
	releaseZbuffer();
	releaseLight();
	//fillTriangle2(a, b, c);
	//wareFrame(list);
	
}

void LightCube() {

	setup();
	VertexAtrr v0(-1.0f, 1.0, -1.0f, 1, 0, 0), v1(1.0f, 1.0, -1.0f, 1, 0, 1),
		v2(1.0f, -1.0, -1.0f, 1, 1, 1), v3(-1.0f, -1.0f, -1.0f, 1, 1, 0);

	VertexAtrr v4(-1.0f, 1.0, 1.0f, 1.0f, 0, 1), v5(1.0f, 1.0, 1.0f, 1, 0, 0),
		v6(1.0f, -1.0, 1.0f, 1.0f, 1.0f, 0.0f), v7(-1.0f, -1.0, 1.0f, 1, 1.0f, 1.0f);

	indeiceBuffer* ib = 0;
	vector<VertexAtrr>list{ v0,v1,v2,v3,v4,v5,v6,v7 };
	
	applyForIndeices(&ib, 12 * 3);
	//front
	ib[0] = 0, ib[1] = 3, ib[2] = 1;
	ib[3] = 3, ib[4] = 2, ib[5] = 1;

	//left
	ib[6] = 0, ib[7] = 4, ib[8] = 3;
	ib[9] = 3, ib[10] = 4, ib[11] = 7;

	//right
	ib[12] = 2, ib[13] = 5, ib[14] = 1;
	ib[15] = 2, ib[16] = 6, ib[17] = 5;

	//back
	ib[18] = 4, ib[19] = 5, ib[20] = 7;
	ib[21] = 7, ib[22] = 5, ib[23] = 6;

	//top
	ib[24] = 0, ib[25] = 1, ib[26] = 4;
	ib[27] = 1, ib[28] = 5, ib[29] = 4;

	//bottom
	ib[30] = 7, ib[31] = 6, ib[32] = 3;
	ib[33] = 6, ib[34] = 2, ib[35] = 3;
	vector<Triangle> triArr;

	vector<Triangle> outlist;
	IndeicesProcessPipeline(&outlist, ib, list, 12);
	//wireFrame(outlist);
	fillMesh(outlist);
	release(ib);
	releaseall();
	//fillTriangle2(a, b, c);
	//wareFrame(list);

}


void releaseall() {
	
	releaseZbuffer();
	releaseLight();
	releaseCM();
	releasePreTex();
}



int main() {
	

	//pyramid();
	
	//Cubetest();
	
	LightCube();
	/*Matrix Tm;
	MatrixTranslation(Tm, 0.5f, 0.0f, 2.0f);
	setMatrix(TS_WORLD, &Tm);
	textureCube();*/
	
	draw("color5.ppm");
	//triangleTest();
	//getTexture("1.jpg");

}




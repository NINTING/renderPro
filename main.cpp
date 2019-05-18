#include"yRasterization.h"
#include<vector>
#define PI 3.14159265359f


using namespace std;


unsigned char img[(Width+10)*(Height+10)*3];


void draw(const char *name) {
	std::ofstream fout(name);
	fout << "P3\n" << Width << " " << Height << "\n255\n";

	for (int i = 0; i <=Height ; i++) {
		unsigned char *p = img + (i * Width) * 3;
		for (int j = 0; j < Width; j++,p+=3) {

			fout << (int)p[0] << " " << (int)p[1] << " " << (int)p[2] << "\n";
		}
	}

	/*unsigned char *p = img;
	for (int i = Height - 1; i >= 0; i--) {
		for (int j = 0; j < Width; j++) {
			p += 3;
			fout << (int)p[0] << " " << (int)p[1] << " " << (int)p[2] << "\n";
		}
	}*/
	fout.close();

}

void setPixel(int x, int y,Color color) {
	unsigned char *p = img + (y*Width + x) * 3;
	p[0] = color._r, p[1] = color._g, p[2] = color._b;
}
void setPixel(const vertex4& v) {
	int y = v._y,	x = v._x;
	unsigned char *p = img + (y*Width + x) * 3;
	p[0] = v._c._r, p[1] = v._c._g, p[2] = v._c._b;
}

void fillMesh(std::vector<Vector4>& list) {
	for (int i = 0; i < list.size(); i+=3) {
		fillTriangle1(list[i], list[i+1], list[i+2]);
	}
}

void wireFrame(std::vector<Vector4> list) {
	for (int i = 0; i < list.size(); i +=3) {
		lineclip(list[i], list[i + 1]);
		bresenham(list[i], list[i + 1]);
		
		lineclip(list[i], list[i + 2]);
		bresenham(list[i], list[i + 2]);
		
		lineclip(list[i+1], list[i + 2]);
		bresenham(list[i+1], list[i + 2]);
	}
}



////默认投影面距离原点为1 即d = 1；
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

Matrix tranMatrix;

void setMatrix(const Matrix& m) {
	tranMatrix *= m;
}

void IndeicesProcessPipeline(vector<vertex4>* outlist,const indeiceBuffer* ib, const vector<vertex4>&vb,int TriangleNum) {
	Matrix Pm, Vm;
	size_t size = TriangleNum * 3;

	PerspectiveMatrix(&Pm,PI*0.5,Width/Height,1,500);
	//(*outlist).resize(size);
	MatrxIdentity(Vm);
	
	//Vector4 at(0,0.0f,-5.0f, 1), view(0.0f, 0.0f,0.0f, 1), up(0,1, 0, 0);
	Vector4 at(3, 3.0f, -3.0f, 1), view(0.0f, 0.0f, 0.0f, 1), up(0, 1, 0, 0);

	
	viewMatrix(Vm, at, view, up);
	tranMatrix = Vm * Pm;
	Vector4 a, b, c;
	for (int i = 0; i <size; i += 3) {
		
		MatrixApply(&a, vb[ib[i]], tranMatrix);
		a._c = vb[ib[i]]._c; a.rhw = a._z;
		
		MatrixApply(&b, vb[ib[i + 1]], tranMatrix);
		b._c = vb[ib[i + 1]]._c; b.rhw = b._z;
		
		MatrixApply(&c, vb[ib[i + 2]],tranMatrix);
		c._c = vb[ib[i + 2]]._c; c.rhw =c._z;
	
		//背面消隐
		if (backCull(a, b, c)) {
			outlist->push_back(a);
			outlist->push_back(b);
			outlist->push_back(c);
		}
	}



	//for (int i = 0; i < size; i += 3) {
	//	MatrixApply(&(*outlist)[i], vb[ib[i]], Vm);
	//	MatrixApply(&(*outlist)[i], (*outlist)[i], Pm);
	//	MatrixApply(&(*outlist)[i + 1], vb[ib[i + 1]],Vm);
	//	MatrixApply(&(*outlist)[i + 1], (*outlist)[i + 1], Pm);
	//	MatrixApply(&(*outlist)[i + 2], vb[ib[i + 2]], Vm);
	//	MatrixApply(&(*outlist)[i + 2], (*outlist)[i + 2], Pm);
	//}
	size = outlist->size();
	for (int i = 0; i < size; i++) {
		toScreen(&(*outlist)[i]);
	}
	
}

//void vertexProcessPipeline(vector<vertex4>* outlist,const vector<vertex4>& list) {
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
	
	vertex4 a(1.0f, 0, 0.0f, 1, Red), b(0.0f, 0.50f, 0.0f, 1, Red), c(0.0f, -1.0f, 0.0f, 1, Red);
	indeiceBuffer* ib = 0;
	vector<Vector4>list{ a,b,c };
	applyForIndeices(&ib, 3);

	ib[0] = 0, ib[1] = 1, ib[2] = 2;

	vector<vertex4> outlist;
	Matrix Pm, Vm;
	MatrxIdentity(tranMatrix);
	PerspectiveMatrix(&Pm, PI*0.5, Width / Height, 1, 500);
	MatrxIdentity(Vm);
	Vector4 at(3, 3.0f, -3.0f, 1), view(0.0f, 0.0f, 0.0f, 1), up(0, 1, 0, 0);
	viewMatrix(Vm, at, view, up);

	setMatrix(Pm);
	IndeicesProcessPipeline(&outlist, ib, list, 1);
	//wireFrame(outlist);
	fillMesh(outlist);
	draw("color1.ppm");
}

void pyramid()
{
	vertex4 a(1.0f, 0, 0.0f, 1, Red), b(0.0f, 0.0f, 1.0f, 1, Yellow), c(-1.0f, 0.0f, 0.0f, 1, Blue), d(0.0f, 0.0f, -1.0f, 1, Red), e(0.0f, 1.0f, 0.0f, 1.0f, White);

	indeiceBuffer* ib = 0;
	vector<Vector4>list{ a,b,c,d,e };
	applyForIndeices(&ib, 4 * 3);

	ib[0] = 0, ib[1] = 1, ib[2] = 4;
	ib[3] = 1, ib[4] = 2, ib[5] = 4;
	ib[6] = 2, ib[7] = 3, ib[8] = 4;
	ib[9] = 3, ib[10] = 0, ib[11] = 4;

	vector<vertex4> outlist;
	IndeicesProcessPipeline(&outlist, ib, list, 4);
	//wireFrame(outlist);
	fillMesh(outlist);

	draw("color2.ppm");

}


void planetest() {
	vertex4 a(6.0f, 0, 0.0f, 1, Red), b(0.0f, 0.0f, 6.0f, 1, Red), c(-6.0f, 0.0f, 0.0f, 1, Red), d(0.0f, 0.0f, -6.0f, 1, Red), e(0.0f, 0.0f, 0.0f, 1.0f, Red);
	indeiceBuffer* ib = 0;
	vector<Vector4>list{ a,b,c,d,e };
	applyForIndeices(&ib, 4 * 3);
	ib[0] = 0, ib[1] = 1, ib[2] = 4;
	ib[3] = 1, ib[4] = 2, ib[5] = 4;
	ib[6] = 2, ib[7] = 3, ib[8] = 4;
	ib[9] = 3, ib[10] = 0, ib[11] = 4;

	vector<vertex4> outlist;
	IndeicesProcessPipeline(&outlist, ib, list, 4);
	//wireFrame(outlist);
	fillMesh(outlist);
	//	release(ib);
		//fillTriangle2(a, b, c);
		//wareFrame(list);
	draw("color3.ppm");
}

int main() {
	memset(img, 255, sizeof(img));
	//float cx = Width * 0.5f - 0.5f, cy = Height * 0.5f - 0.5f;
	vertex4 a(1.0f, 0, 0.0f, 1,Red), b(0.0f, 0.0f, 1.0f, 1,Yellow), c(-1.0f, 0.0f,0.0f, 1,Blue),d(0.0f,0.0f,-1.0f,1,Red	),e(0.0f,1.0f,0.0f,1.0f,White);
	//
	
	/*vertex4 a(1.0f, 0, 0.0f, 1, Red), b(0.0f, 0.50f, 0.0f, 1, Red), c(0.0f, -1.0f, 0.0f, 1, Red);
	vector<Vector4>list{ a,b,c };*/
	indeiceBuffer* ib = 0;
	vector<Vector4>list{a,b,c,d,e};
	applyForIndeices(&ib, 4*3);
	/*ib[0] = 0, ib[1] = 1, ib[2] = 2;*/
	

	ib[0] = 0, ib[1] = 1, ib[2] = 4;
	ib[3] = 1, ib[4] = 2, ib[5] = 4;
	ib[6] = 2, ib[7] = 3, ib[8] = 4;
	ib[9] = 3, ib[10] = 0, ib[11] = 4;

	vector<vertex4> outlist;	
	IndeicesProcessPipeline(&outlist, ib, list,4);
	//wireFrame(outlist);
	fillMesh(outlist);
//	release(ib);
	//fillTriangle2(a, b, c);
	//wareFrame(list);
	draw("color3.ppm");

	
}
















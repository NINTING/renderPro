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


void fillMesh(std::vector<Vector4>& list) {
	for (int i = 0; i < list.size(); i+=3) {
		fillTriangle1(list[i], list[i+1], list[i+2]);
	}
}

void wireFrame(std::vector<Vector4>& list) {
	for (int i = 0; i < list.size(); i +=3) {
		lineclip(list[i], list[i + 1]);
		lineclip(list[i], list[i + 2]);
		lineclip(list[i+1], list[i + 2]);
	}
}

//默认投影面距离原点为1 即d = 1；
void PerspectiveMatrix(Matrix * Pm) {
	float Right = Width*0.5f, Left = -Right;
	float Top = Height*0.5f, Bottom = -Top;
	const float near = 1.0f, far = 500.0f;
	float a = (float)far / (float)(far - near), b = (float)(far*near) / (float)(near - far);
	 (*Pm)(0, 0) = 2 * near / (Right - Left);
	 (*Pm)(1, 1) = 2 * near / (Top - Bottom);
	 (*Pm)(2, 2) = a;
	 (*Pm)(3, 2) = b;
	 (*Pm)(2, 3) = 1.0f;
}

void toScreen(Vector4 *p){
	float dw = 1.0f / p->_w;
	p->_x = (p->_x*dw + 1.0f)*Width *0.5f;
	p->_y = (-p->_y*dw + 1.0f)*Height*0.5f;
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


void IndeicesProcessPipeline(vector<vertex4>* outlist,const indeiceBuffer* ib, const vector<vertex4>&vb) {
	Matrix Pm;
	size_t size = vb.size() * 3;
	PerspectiveMatrix(&Pm);
	(*outlist).resize(size);
	for (int i = 0; i <size; i += 3) {
		MatrixApply(&(*outlist)[i], vb[ib[i]], Pm);
		MatrixApply(&(*outlist)[i + 1], vb[ib[i + 1]], Pm);
		MatrixApply(&(*outlist)[i + 2], vb[ib[i + 2]], Pm);
	}
	for (int i = 0; i < size; i++) {
		toScreen(&(*outlist)[i]);
	}
	
}

void vertexProcessPipeline(vector<vertex4>* outlist,const vector<vertex4>& list) {
	Matrix Pm;
	PerspectiveMatrix(&Pm);
	(*outlist).resize(list.size());
	for (int i = 0; i < list.size(); i += 3) {
		MatrixApply(&(*outlist)[i], list[i], Pm);
		MatrixApply(&(*outlist)[i+1], list[i + 1], Pm);
		MatrixApply(&(*outlist)[i+2], list[i + 2], Pm);
	}
	for (int i = 0; i < list.size(); i++) {
		toScreen(&(*outlist)[i]);
	}
}

int main() {
	memset(img, 255, sizeof(img));
	//float cx = Width * 0.5f - 0.5f, cy = Height * 0.5f - 0.5f;
	vertex4 a(-300, 0, 2, 1), b(300, 0, 2, 1), c(0, 300, 2, 1),d(0,0,300,1);
	
	indeiceBuffer* ib = 0;
	vector<Vector4>list{a,b,c,d};
	applyForIndeices(&ib, list.size()*3);
	
	ib[0] = 0, ib[1] = 1, ib[2] = 2;
	ib[3] = 1, ib[4] = 3, ib[5] = 2;
	ib[6] = 0, ib[7] = 3, ib[8] = 2;
	ib[9] = 0, ib[10] = 1, ib[11] = 3;
	vector<vertex4> outlist;	
	IndeicesProcessPipeline(&outlist, ib, list);
	wireFrame(outlist);

	release(ib);
	//fillTriangle2(a, b, c);
	//wareFrame(list);
	draw("color1.ppm");

	
}
















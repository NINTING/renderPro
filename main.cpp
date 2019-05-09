#include"yRasterization.h"
#include<vector>
#define PI 3.14159265359f


using namespace std;


unsigned char img[600*600*3];


void draw(const char *name) {
	std::ofstream fout(name);
	fout << "P3\n" << Width << " " << Height << "\n255\n";

	for (int i = Height-1; i >=0 ; i--) {
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

void setPixel(int x, int y,unsigned char color) {
	unsigned char *p = img + (y*Width + x) * 3;
	p[0] = p[1] = p[2] = color;
}


void fillMesh(std::vector<Vector>& list) {
	for (int i = 0; i < list.size(); i+=3) {
		fillTriangle1(list[i], list[i+1], list[i+2]);
	}
}

void wareFrame(std::vector<Vector>& list) {
	for (int i = 0; i < list.size(); i +=3) {
		lineclip(list[i], list[i + 1]);
		lineclip(list[i], list[i + 2]);
		lineclip(list[i+1], list[i + 2]);
	}
}


int main() {
	memset(img, 255, sizeof(img));
	//float cx = Width * 0.5f - 0.5f, cy = Height * 0.5f - 0.5f;
	vertex a,b,c,d;
	a = vertex(100, 100); c = vertex(300, 400);
	b = vertex(400, 100); d = vertex(200, 400);
	
	//bresenham(lmin, lmax);
	//bresenham(lmin, rmin);
	//bresenham(lmax, rmax);
	//bresenham(rmin, rmax);

	//
	std::vector<Vector> list{ a,c,d,a,c,b };
	fillMesh(list);
	//fillTriangle2(a, b, c);
	//wareFrame(list);
	draw("color1.ppm");

	
}
















#include"yRasterization.h"



Vector VecCross(const Vector& a, const Vector& b) {
	return Vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);

}


//---------------------------------
//bresenham 直线画法			   |
//---------------------------------

void bresenham(const vertex &begin, const vertex& end) {

	int x1, y1, x2, y2;
	x1 = begin.x, y1 = begin.y, x2 = end.x, y2 = end.y;
	bresenham(x1, y1, x2, y2);
}



void bresenham(int x1,int y1, int x2,int y2) {
	int dx = abs(x2 - x1), dy = abs(y2 - y1);
	int sx, sy;
	sx = x2 > x1 ? 1 : -1;
	sy = y2 > y1 ? 1 : -1;
	int flag = 0;
	//判断x，y以谁为主轴
	if (dy > dx) {
		std::swap(dy, dx);
		flag = 1;
	}
	int p = 2 * dy - dx;
	for (int i = 0; i < dx; i++) {
		setPixel(x1, y1);
		if (p >= 0) {
			if (!flag)y1 += sy;
			else x1 += sx;
			p -= 2 * dx;
		}
		//主轴增加
		if (flag)y1 += sy;
		else x1 += sx;

		p += 2 * dy;
	}
}

//---------------------------------
//扫描线              光栅化三角形|
//---------------------------------
void fillTriangle1(vertex v0, vertex v1, vertex v2) {			//scan line
	//sort on the y axis
	if (v0.y > v1.y)std::swap(v0, v1);
	if (v0.y > v2.y)std::swap(v0, v2);
	if (v1.y > v2.y) std::swap(v1, v2);

	//scan flat-top triangel
	for (int y = v0.y; y <= v1.y; y++) {
		float x1 = (float)(y - v0.y)*(v1.x - v0.x) / (v1.y - v0.y + 1) + v0.x;
		float x2 = (float)(y - v0.y)*(v2.x - v0.x) / (v2.y - v0.y + 1) + v0.x;
		
		//if (x1 > x2)std::swap(x1, x2);
		vertex b(x1, y), e(x2, y);
		lineclip(b, e);
		/*for(int x = x1;x<=x2;x++)
			setPixel(x, y);*/
		
	}

	//scan flat triangle
	for (int y = v1.y; y <= v2.y; y++) {
		float x1 = (float)(y - v0.y)*(v2.x - v0.x) / (v2.y - v0.y + 1) + v0.x;
		float x2 = (float)(y - v1.y)*(v2.x - v1.x) / (v2.y - v1.y + 1) + v1.x;
	
		//if (x1 > x2)std::swap(x1, x2);
		/*for (int x = x1; x <= x2; x++)
			setPixel(x, y);*/
		vertex b(x1, y), e(x2, y);
		lineclip(b,e);
		
	}


}
//---------------------------------
//包围盒 判断重心坐标 光栅化三角形|
//---------------------------------

void fillTriangle2(const vertex& v0, const vertex& v1, const vertex& v2) {
	//create BoundingBox
	vertex lmin(Width-1, Height-1),rmax(0,0);
	lmin.x = std::min(v0.x, std::min(v1.x, v2.x));
	lmin.y = std::min(v0.y, std::min(v1.y, v2.y));
	rmax.x = std::max(v0.x, std::max(v1.x, v2.x));
	rmax.y = std::max(v0.y, std::max(v1.y, v2.y));

	lmin.x = std::max(0.0f, lmin.x);		lmin.y = std::max(0.0f, lmin.y);
	rmax.x = std::min(Width*1.0f, rmax.x);	rmax.y = std::min(Height*1.0f, rmax.y);

	for (int y = rmax.y; y >= lmin.y; y--) {
		for (int x = lmin.x; x <= rmax.x; x++) {
			Vector p = barycentric(v0, v1, v2, vertex(x, y));
			if (p.x < 0 || p.y < 0 || p.z < 0)
				continue;
			setPixel(x, y);
		}
	}

}

Vector barycentric(const vertex& v0, const vertex &v1, const vertex& v2,const vertex& v) {
	Vector bary = VecCross(Vector(v1.x-v0.x, v2.x - v0.x, v0.x - v.x), Vector(v1.y - v0.y, v2.y - v0.y, v0.y - v.y));
	if (bary.z < 1)return Vector(-1, 1, 1);
	float a = bary.x / bary.z, b = bary.y / bary.z;
	return Vector(1.0f - a - b, a, b);
}



//---------------------------------
//       直线裁剪				  |
//---------------------------------

int getcode(const vertex &v) {
	int p = 0;
	if (v.x < Lx) p |= L;
	if (v.y < By) p |= B;
	if (v.y > Ty) p |= T;
	if (v.x > Rx) p |= R;
	return p;

}

void clip(int& code, vertex& v, const vertex &u) {
	if (code&L) {
		v.y = v.y + (u.y - v.y)*(Lx - v.x) / (u.x - v.x);
		v.x = Lx;
	}
	else if (code&R) {
		v.y = v.y + (u.y - v.y)*(Rx - v.x) / (u.x - v.x);
		v.x = Rx;
	}
	else if (code&T) {
		v.x = v.x + (u.x - v.x)*(Ty - v.y) / (u.y - v.y);
		v.y = Ty;
	}
	else if (code&B) {
		v.x = v.x + (u.x - v.x)*(By - v.y) / (u.y - v.y);
		v.y = By;
	}
	code = getcode(v);

}

void lineclip(vertex& a, vertex& b) {

	//	int x0 = a.x, int y0 = a.y, int x1 = b.x, int y1=b.y;
	int p1 = 0, p2 = 0;
	p1 = getcode(a); p2 = getcode(b);
	if (p1&p2) return;
	//if (!p1 && !p2)

	while (p1 || p2) {			
		clip(p1, a, b);
		clip(p2, b, a);
	}
		
	
	bresenham(a, b);
}



//---------------------------------
//       矩阵相关函数			  |
//---------------------------------

Matrix ::Matrix(float x11, float x12, float x13, float x14,
	float x21, float x22, float x23, float x24,
	float x31, float x32, float x33, float x34,
	float x41, float x42, float x43, float x44) 
{
	value[0][0] = x11, value[0][1] = x12, value[0][2] = x13, value[0][3] = x14;
	value[1][0] = x21, value[1][1] = x22, value[1][2] = x23, value[1][3] = x24;
	value[2][0] = x11, value[2][1] = x22, value[2][2] = x33, value[2][3] = x34;
	value[3][0] = x11, value[3][1] = x22, value[3][2] = x43, value[3][3] = x44;
}

Matrix& operator *= (Matrix &lsh,Matrix &rsh) {
	lsh = lsh * rsh;
	return lsh;
	
}

Matrix operator * (const Matrix &lrh, const Matrix &rsh) {
	float x00 = lrh.value[0][0] * rsh(0, 0) + lrh.value[0][1] * rsh(1, 0) + lrh.value[0][2] * rsh(2, 0) + lrh.value[0][3] * rsh(3, 0);
	float x01 = lrh.value[0][0] * rsh(0, 1) + lrh.value[0][1] * rsh(1, 1) + lrh.value[0][2] * rsh(2, 1) + lrh.value[0][3] * rsh(3, 1);
	float x02 = lrh.value[0][0] * rsh(0, 2) + lrh.value[0][1] * rsh(1, 2) + lrh.value[0][2] * rsh(2, 2) + lrh.value[0][3] * rsh(3, 2);
	float x03 = lrh.value[0][0] * rsh(0, 3) + lrh.value[0][1] * rsh(1, 3) + lrh.value[0][2] * rsh(2, 3) + lrh.value[0][3] * rsh(3, 3);

	float x10 = lrh.value[1][0] * rsh(0, 0) + lrh.value[1][1] * rsh(1, 0) + lrh.value[1][2] * rsh(2, 0) + lrh.value[1][3] * rsh(3, 0);
	float x11 = lrh.value[1][0] * rsh(0, 1) + lrh.value[1][1] * rsh(1, 1) + lrh.value[1][2] * rsh(2, 1) + lrh.value[1][3] * rsh(3, 1);
	float x12 = lrh.value[1][0] * rsh(0, 2) + lrh.value[1][1] * rsh(1, 2) + lrh.value[1][2] * rsh(2, 2) + lrh.value[1][3] * rsh(3, 2);
	float x13 = lrh.value[1][0] * rsh(0, 3) + lrh.value[1][1] * rsh(1, 3) + lrh.value[1][2] * rsh(2, 3) + lrh.value[1][3] * rsh(3, 3);

	float x20 = lrh.value[2][0] * rsh(0, 0) + lrh.value[2][1] * rsh(1, 0) + lrh.value[2][2] * rsh(2, 0) + lrh.value[2][3] * rsh(3, 0);
	float x21 = lrh.value[2][0] * rsh(0, 1) + lrh.value[2][1] * rsh(1, 1) + lrh.value[2][2] * rsh(2, 1) + lrh.value[2][3] * rsh(3, 1);
	float x22 = lrh.value[2][0] * rsh(0, 2) + lrh.value[2][1] * rsh(1, 2) + lrh.value[2][2] * rsh(2, 2) + lrh.value[2][3] * rsh(3, 2);
	float x23 = lrh.value[2][0] * rsh(0, 3) + lrh.value[2][1] * rsh(1, 3) + lrh.value[2][2] * rsh(2, 3) + lrh.value[2][3] * rsh(3, 3);

	float x30 = lrh.value[3][0] * rsh(0, 0) + lrh.value[3][1] * rsh(1, 0) + lrh.value[3][2] * rsh(2, 0) + lrh.value[3][3] * rsh(3, 0);
	float x31 = lrh.value[3][0] * rsh(0, 1) + lrh.value[3][1] * rsh(1, 1) + lrh.value[3][2] * rsh(2, 1) + lrh.value[3][3] * rsh(3, 1);
	float x32 = lrh.value[3][0] * rsh(0, 2) + lrh.value[3][1] * rsh(1, 2) + lrh.value[3][2] * rsh(2, 2) + lrh.value[3][3] * rsh(3, 2);
	float x33 = lrh.value[3][0] * rsh(0, 3) + lrh.value[3][1] * rsh(1, 3) + lrh.value[3][2] * rsh(2, 3) + lrh.value[3][3] * rsh(3, 3);

	return Matrix(x00, x01, x02, x03,
		x10, x11, x12, x13,
		x20, x21, x22, x23,
		x30, x31, x32, x33);

}



Vector4 operator * (const Vector4 &lrh, const Matrix &rsh) {
	
	float x = lrh._x * rsh(0, 0) + lrh._y * rsh(1, 0) + lrh._z * rsh(2, 0) + lrh._w * rsh(3, 0);
	float y = lrh._x * rsh(0, 1) + lrh._y * rsh(1, 1) + lrh._z  * rsh(2, 1) + lrh._w * rsh(3, 1);
	float z = lrh._x * rsh(0, 2) + lrh._y * rsh(1, 2) + lrh._z  * rsh(2, 2) + lrh._w * rsh(3, 2);
	float w = lrh._x * rsh(0, 3) + lrh._y * rsh(1, 3) + lrh._z  * rsh(2, 3) + lrh._w * rsh(3, 3);
	return Vector4(x, y, z, w);
}


Vector4& operator *= (vertex4 &lsh, Matrix &rsh) {
	lsh = lsh * rsh;
	return lsh;

}





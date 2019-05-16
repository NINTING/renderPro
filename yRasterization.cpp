#include"yRasterization.h"



void VecCross(Vector4 *out,const Vector4& a, const Vector4& b) {
	 (*out)._w = 1.0f;
	 (*out)._x = a._y*b._z - a._z*b._y;
	 (*out)._y = a._z*b._x - a._x*b._z;
	 (*out)._z = a._x*b._y - a._y*b._x;
}


//---------------------------------
//bresenham 直线画法			   |
//---------------------------------

void bresenham(const vertex4 &begin, const vertex4& end) {

	int x1, y1, x2, y2;
	x1 = begin._x, y1 = begin._y, x2 = end._x, y2 = end._y;
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
		setPixel(x1, y1,Black);
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
void fillTriangle1(Vector4 v0, Vector4 v1, Vector4 v2) {			//scan line
	//sort on the y axis
	if (v0._y > v1._y)std::swap(v0, v1);
	if (v0._y > v2._y)std::swap(v0, v2);
	if (v1._y > v2._y) std::swap(v1, v2);

	//scan flat-top triangel
	for (int y = v0._y; y <= v1._y; y++) {
		float x1 = (float)(y - v0._y)*(v1._x - v0._x) / (v1._y - v0._y + 1) + v0._x;
		float x2 = (float)(y - v0._y)*(v2._x - v0._x) / (v2._y - v0._y + 1) + v0._x;

		if (x1 > x2)std::swap(x1, x2);
		Vector4 b(x1, y,0,1), e(x2, y,0,1);
		lineclip(b, e);
		for(int x = b._x;x<=e._x;x++)
			setPixel(x, y,Black);

	}

	//scan flat triangle
	for (int y = v1._y; y <= v2._y; y++) {
		float x1 = (float)(y - v0._y)*(v2._x - v0._x) / (v2._y - v0._y + 1) + v0._x;
		float x2 = (float)(y - v1._y)*(v2._x - v1._x) / (v2._y - v1._y + 1) + v1._x;

		//if (x1 > x2)std::swap(x1, x2);
		/*for (int x = x1; x <= x2; x++)
			setPixel(x, y);*/
		Vector4 b(x1, y,0,1), e(x2, y,0,1);
		lineclip(b,e);

	}


}
////---------------------------------
////包围盒 判断重心坐标 光栅化三角形|
////---------------------------------
//
//void fillTriangle2(const Vector4& v0, const Vector4& v1, const Vector4& v2) {
//	//create BoundingBox
//	Vector4 lmin(Width-1, Height-1,0,0),rmax(0,0,0,0);
//	lmin._x = std::min(v0._x, std::min(v1._x, v2._x));
//	lmin._y = std::min(v0._y, std::min(v1._y, v2._y));
//	rmax._x = std::max(v0._x, std::max(v1._x, v2._x));
//	rmax._y = std::max(v0._y, std::max(v1._y, v2._y));
//
//	lmin._x = std::max(0.0f, lmin._x);		lmin._y = std::max(0.0f, lmin._y);
//	rmax._x = std::min(Width*1.0f, rmax._x);	rmax._y = std::min(Height*1.0f, rmax._y);
//
//	for (int y = rmax._y; y >= lmin._y; y--) {
//		for (int x = lmin._x; x <= rmax._x; x++) {
//			Vector4 p = barycentric(v0, v1, v2, Vector4(x, y));
//			if (p._x < 0 || p._y < 0 || p.z < 0)
//				continue;
//			setPixel(x, y);
//		}
//	}
//
//}
//
//Vector4 barycentric(const Vector4& v0, const Vector4 &v1, const vertex4& v2,const vertex4& v) {
//	Vector4 bary = VecCross(Vector4(v1._x-v0._x, v2._x - v0._x, v0._x - v._x), Vector4(v1._y - v0._y, v2._y - v0._y, v0._y - v._y));
//	if (bary.z < 1)return Vector4(-1, 1, 1);
//	float a = bary._x / bary.z, b = bary._y / bary.z;
//	return Vector4(1.0f - a - b, a, b);
//}
//


//---------------------------------
//       2D直线裁剪				  |
//---------------------------------

int getcode(const vertex4 &v) {
	int p = 0;
	if (v._x < Lx) p |= L;
	if (v._y < By) p |= B;
	if (v._y > Ty) p |= T;
	if (v._x > Rx) p |= R;
	return p;

}

void clip(int& code, vertex4& v, const vertex4 &u) {
	if (code&L) {
		v._y = v._y + (u._y - v._y)*(Lx - v._x) / (u._x - v._x);
		v._x = Lx;
	}
	else if (code&R) {
		v._y = v._y + (u._y - v._y)*(Rx - v._x) / (u._x - v._x);
		v._x = Rx;
	}
	else if (code&T) {
		v._x = v._x + (u._x - v._x)*(Ty - v._y) / (u._y - v._y);
		v._y = Ty;
	}
	else if (code&B) {
		v._x = v._x + (u._x - v._x)*(By - v._y) / (u._y - v._y);
		v._y = By;
	}
	code = getcode(v);

}

void lineclip(vertex4 &a, vertex4 &b) {

	//	int x0 = a._x, int y0 = a._y, int x1 = b._x, int y1=b._y;
	int p1 = 0, p2 = 0;
	p1 = getcode(a); p2 = getcode(b);
	if (p1&p2) return;
	//if (!p1 && !p2)

	while (p1 || p2) {
		clip(p1, a, b);
		clip(p2, b, a);
	}



}

//---------------------------------
//       向量相关函数			  |
//---------------------------------

Vector4 operator - (const Vector4& lsh, const Vector4& rsh) {
	return Vector4(lsh._x - rsh._x, lsh._y - rsh._y, lsh._z - rsh._z,0);
	
}

Vector4& operator -= (Vector4& lsh, const Vector4& rsh) {
	lsh = lsh - rsh;
	return lsh;
}

Vector4 operator + (const Vector4& lsh, const Vector4& rsh) {
	return Vector4(lsh._x + rsh._x, lsh._y + rsh._y, lsh._z + rsh._z, 0);


}

Vector4& operator += (Vector4& lsh, const Vector4& rsh) {
	lsh = lsh + rsh;
	return lsh;

}

void VecNormalize(Vector4 & out,const Vector4& in) {
	float x = in._x, y = in._y, z = in._z;
	
	float k = sqrt(x * x + y * y + z * z);

	out._x = in._x / k, out._y = in._y / k, out._z = in._z / k;
}

float operator * (const Vector4& lsh, const Vector4& rsh) {
	return lsh._x*rsh._x + lsh._y*rsh._y + lsh._z*rsh._z;

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
	value[2][0] = x31, value[2][1] = x32, value[2][2] = x33, value[2][3] = x34;
	value[3][0] = x41, value[3][1] = x42, value[3][2] = x43, value[3][3] = x44;
}



void MatrxIdentity(Matrix &out) {
	out(0, 0) = 1.0f;
	out(1, 1) = 1.0f;
	out(2, 2) = 1.0f;
	out(3, 3) = 1.0f;

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

Matrix getIdentity() {
	return Matrix(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f );

}

void RotationAxisX(Matrix&out,float theta){

	out(0,0) = 1.0f;
	out(1,1) = cos(theta);
	out(1,2) = sin(theta);
	out(2,1) = -sin(theta);
	out(2,2) = cos(theta);
	out(3,3) = 1.0f;
}
void RotationAxisY(Matrix& out,float theta){

	out(1,1) = 1.0f;
	out(0,0) = cos(theta);
	out(0,2) = -sin(theta);
	out(2,0) = sin(theta);
	out(2,2) = cos(theta);
	out(3,3) = 1.0f;
}
void RotationAxisZ(Matrix& out,float theta){


	out(0,0) = cos(theta);
	out(0,1) = sin(theta);
	out(2,0) = -sin(theta);
	out(1,1) = cos(theta);
	out(2,2) = 1.0f;
	out(3,3) = 1.0f;
}



void RotationAxis(Matrix&out,const Vector4& axis,float theta){
    float cosT = cosf(theta);
    float sinT = cosf(theta);
    float oneSubc = 1.0f-cosT;

    float x = axis._x;
    float y = axis._y;
    float z = axis._z;

    out(0,0) = x*x*oneSubc+cosT; out(0,1) = x*y*oneSubc+x*sinT; out(0,2) = x*z*oneSubc-y*sinT;
	out(0, 0) = x * y*oneSubc - z*sinT; out(0, 1) = y * y*oneSubc + cosT; out(0, 2) = y * z*oneSubc - x * sinT;
	out(0, 0) = x * z*oneSubc + y*sinT; out(0, 1) = z * y*oneSubc - x * sinT; out(0, 2) = z * z*oneSubc + cosT;
	out(3, 3) = 1;
}



void viewMatrix(Matrix&out,const Vector4& at,const Vector4& view,const Vector4& up) {
	
	Vector4 Z = view - at;
	Vector4 X,Y;
	VecNormalize(Z,Z);
	VecCross(&X, up, Z);
	VecNormalize(X, X);
	VecCross(&Y, Z, X);
	VecNormalize(Y, Y);
	
	out(0, 0) = X._x;
	out(1, 0) = X._y;
	out(2, 0) = X._z;
	out(3, 0) = -(at * X);

	out(0, 1) = Y._x;
	out(1, 1) = Y._y;
	out(2, 1) = Y._z;
	out(3, 1) = -(at * Y);

	out(0, 2) = Z._x;
	out(1, 2) = Z._y;
	out(2, 2) = Z._z;
	out(3, 2) = -(at * Z);

	out(0, 3) = 0.0f;
	out(1, 3) = 0.0f;
	out(2, 3) = 0.0f;
	out(3, 3) = 1.0f;

}




//---------------------------------
//      背面消隐				  |
//---------------------------------


//正面为逆时针方向构成
bool backCull(const vertex4& a,const vertex4& b,const vertex4& c) {
	Vector4 u = b - a;
	Vector4 v = a - c;
	
	Vector4 zero(0, 0, 0, 0);
	Vector4 normal;
	VecCross(&normal, v, u);
	Vector4 view = a - zero;
	
	float cosT = view * normal;
	if (cosT <= 0) 
		return false;
	return true;
}

//---------------------------------
//       线性插值				  |
//---------------------------------

void slerp(float* s, float* t, vertex4 v0, vertex4 v1, vertex4 v2) {


}
void Zslerp(const vertexArr *va1, const vertexArr *va2, vertexArr *out, float t) {


}

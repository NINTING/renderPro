#include"yRasterization.h"

char img[(Width + 10)*(Height + 10)*3];

Vector4 Zero(0, 0, 0, 1);


void VecCross(Vector4 *out,const Vector4& a, const Vector4& b) {
	 (*out)._w = 0.0f;
	 (*out)._x = a._y*b._z - a._z*b._y;
	 (*out)._y = a._z*b._x - a._x*b._z;
	 (*out)._z = a._x*b._y - a._y*b._x;
}

//---------------------------------
//			��������			   |
//---------------------------------

int Presentfvf;		//��ǰ�Ķ����ʽ
Texture& PresentTex() {
	static Texture* PresentTex = new Texture();
	return *PresentTex;
}

Texture& Texture::operator = (const Texture& rsh) {
	width = rsh.width; height = rsh.height;
	tex = rsh.tex;
	return *this;
}

void releasePreTex() {
	delete[] PresentTex().tex;
	delete &PresentTex();
}


	//��ǰ������
void getTexPixel(const Texture &Texobj ,Color* c,float x,float y) {
	if (x <= 0.0f)x = 0.0f;
	if (x > 1.0f)x = 1.0f;
	if (y <= 0.0f)y = 0.0f;
	if (y > 1.0f)y = 1.0f;

	x *= Texobj.width - 1;
	y *= Texobj.height - 1;
	int lx = floorf(x), rx = ceilf(x);
	int ty = floorf(y), by = ceilf(y);



	float lt = (x - lx) * (y - ty), rt = (rx - x) * (y - ty);
	float lb = (x - lx) * (by - y), rb = (rx - x) * (by - y);
	float* p;
	Color c1, c2, c3, c4;
	p = Texobj.tex + (lx *Texobj.width  + ty)*3;
	c1._r *= p[0], c1._g *= p[1], c1._b *= p[2];

	 p = Texobj.tex + (rx * Texobj.width + ty) * 3;
	c2._r *= p[0], c2._g *= p[1], c2._b *= p[2];
	
	 p = Texobj.tex + (lx * Texobj.width + by) * 3;
	c3._r *= p[0], c3._g *= p[1], c3._b *= p[2];
	
	 p = Texobj.tex + (rx * Texobj.width + by) * 3;
	c4._r *= p[0], c4._g *= p[1], c4._b *= p[2];
	*c = c1 + c2 + c3 + c4;
	return;
}



void init_Texture(Texture * texobj, int w, int h) {
	texobj->width = w;	
	texobj->height = h;
	texobj->tex = new float[w * h * 3];
	for (int i = 0; i < h; i++) {
		float* p = texobj->tex + i * w * 3;
		for (int j = 0; j < w; j++, p += 3) {
			if ((i / 32 + j / 32) & 1)
				p[0] =1, p[1] = 1, p[2] = 1;
			else
				p[0] = (float)0x3f/255.0f, p[1] = (float)0xbc/ 255.0f, p[2] = (float)0xef/255.0f;
			
		}
	}
}

void RsetFVF(int FvF) {
	Presentfvf = 0;
	Presentfvf |= FvF;
}
void RsetTex(Texture *tex) {
	delete [] PresentTex().tex;
	
	PresentTex() = *tex;
	PresentTex().tex;
}
//---------------------------------
//			color ����			   |
//---------------------------------
Color operator +(const Color& a, const Color&b) {
	return Color(a._r + b._r, a._g + b._g, a._b + b._b);
}
Color& operator +=(Color& a, const Color&b) {
	a = a + b;
	return a;
}

Color operator -(const Color& a, const Color& b) {
	return Color(a._r - b._r, a._g - b._g, a._b - b._b);
}
Color& operator -=(Color& a, const Color& b) {
	a = a - b;
	return a;
}

Color operator *(const Color& a, float t) {
	return Color(a._r*t, a._g *t, a._b*t);
}
Color& operator *=(Color& a, float t) {
	return a = a * t;
	
}

Color operator *(float t, const Color& a) {
	return a * t;
}
Color operator *(const Color& a, const Color& b) {
	return Color(a._r * b._r, a._g * b._g, a._b * b._b);
}
Color& operator *=(Color& a, const Color b) {
	return a = a * b;
}

Color operator /(const Color& a, float t) {
	return Color(a._r / t, a._g / t, a._b / t);

}
Color& operator /=(Color& a, float t) {
	return a = a / t;
}


//---------------------------------
//bresenham ֱ�߻���			   |
//---------------------------------

void bresenham(const VertexAtrr &begin, const VertexAtrr& end) {

	int x1, y1, x2, y2;
	x1 = begin._v._x, y1 = begin._v._y, x2 = end._v._x,y2 = end._v._y;
	bresenham(x1, y1, x2, y2);
}



void bresenham(int x1,int y1, int x2,int y2) {
	int dx = abs(x2 - x1), dy = abs(y2 - y1);
	int sx, sy;
	sx = x2 > x1 ? 1 : -1;
	sy = y2 > y1 ? 1 : -1;
	int flag = 0;
	//�ж�x��y��˭Ϊ����
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
		//��������
		if (flag)y1 += sy;
		else x1 += sx;

		p += 2 * dy;
	}
}

//---------------------------------
//ɨ����        ��դ��������	  |
//---------------------------------

void wireFrame(std::vector<Triangle> list) {
	for (int i = 0; i < list.size(); i++) {
		lineclip(list[i].v0, list[i].v1);
		bresenham(list[i].v0, list[i].v1);

		lineclip(list[i].v0, list[i].v2);
		bresenham(list[i].v0, list[i].v2);

		lineclip(list[i].v1, list[i].v2);
		bresenham(list[i].v1, list[i].v2);
	}
}


void fillMesh(std::vector<VertexAtrr>& list) {
	size_t size = list.size();
	for (int i = 0; i < size; i += 3) {
		fillTriangle1(list[i], list[i + 1], list[i + 2]);
	}
}

void fillMesh(std::vector<Triangle>& list) {
	size_t size = list.size();
	if (Presentfvf == FVFwireframe)
		wireFrame(list);
	else
	for (int i = 0; i < size; i++) {
		fillTriangle1(list[i]);
	}
}






void checkColor(Color& c) {
	float maxc = std::max(c._r, std::max(c._b, c._g));
	if (maxc > 1.0f)
		c /= maxc;
	
}

void setPixel(int x, int y, Color color) {
	char *p = img + (y*Width + x)*3;
	
	checkColor(color);

	p[0] = (char)(color._r*255+0.5) ,p[1] = (char)(color._g*255+0.5) , p[2] = (char)(color._b*255+0.5);

}
void setPixel(const VertexAtrr& va) {
	int y = va._v._y, x = va._v._x;
	char *p = img + (y*Width + x) * 3;
	p[0] = va._c._r, p[1] = va._c._g, p[2] = va._c._b;
}

void fillTriangle1(Triangle& tri) {
	fillTriangle1(tri.v0, tri.v1, tri.v2);
}
void fillTriangle1(VertexAtrr v0, VertexAtrr v1, VertexAtrr v2) {			//scan line															//sort on the y axis
	if (v0._v._y > v1._v._y)std::swap(v0, v1);
	if (v0._v._y > v2._v._y)std::swap(v0, v2);
	if (v1._v._y > v2._v._y) std::swap(v1, v2);


	
	//scan flat-top triangel
	for (int y = ceilf(v0._v._y); y <= v2._v._y; y++) {
		//���²�ֵ
		float d1,d2;
		
		//���Ҷ˵�Ͳ�ֵ
		VertexAtrr lv, rv;
		d2 = (float)(y - v0._v._y) / (float)(v2._v._y - v0._v._y+1);
		vertexInterp(&rv, v0, v2, d2);

		if (y < v1._v._y) {
			d1 = (float)(y - v0._v._y) / (float)(v1._v._y - v0._v._y +1);
			vertexInterp(&lv, v0, v1, d1);
		}
		else {
			d1 = (float)(y - v1._v._y) / (float)(v2._v._y - v1._v._y + 1);
			vertexInterp(&lv, v1, v2, d1);
		}
		
		lv._v._y = rv._v._y = y;
		if (lv._v._x > rv._v._x)std::swap(lv, rv);
		
		lineclip(lv, rv);
		
		float overSubx = 1 / (rv._v._x - lv._v._x);
		if (rv._v._x - lv._v._x == 0) {
			overSubx = 0.0f;
			continue;
		}
		float zstep = (rv.rhw - lv.rhw) *overSubx;

		float ustep = 0;
		float vstep = 0;
		//float zbstep = (rv._v._z - lv._v._z) *overSubx;

		
		Color Coverz, cStep;		//��Ƚ�����ɫ������
		if (Presentfvf & FVFcolor) {
			//��ɫ��ֵ
			cStep = (rv._c - lv._c) * overSubx;
			/*cStep._b = (rv._c._b - lv._c._b)*overSubx;
			cStep._g = (rv._c._g - lv._c._g)*overSubx;
			cStep._r = (rv._c._r - lv._c._r)*overSubx;*/
			Coverz = lv._c;
		}
		else if (Presentfvf&FVFtexture)
		{
			ustep = (rv._tu - lv._tu)*overSubx;
			vstep = (rv._tv - lv._tv)*overSubx;

			//getTexPixel()
		}	
		Color Ldoverz(lv.ld),Lsoverz(lv.ls);
		Color LdStep((rv.ld - lv.ld)*overSubx), LsStep ((rv.ls-lv.ls)*overSubx);


		float overz,uoverz,voverz;
		int x;
		float z;
		for (x = lv._v._x+0.5f,z = lv._v._z,overz = lv.rhw,uoverz = lv._tu,voverz = lv._tv;
			x <= rv._v._x;
			x++,overz+=zstep,Ldoverz+=LdStep,Lsoverz+=LsStep) {
			Color c(1,1,1);
			float zz = 1.0f / overz;
			
			if (Presentfvf & FVFcolor) {
				c._r = Coverz._r * zz, c._g = Coverz._b * zz, c._b = Coverz._b * zz;
				Coverz += cStep;
			}
			if (Presentfvf & FVFtexture) {
				float u = uoverz * zz;
				float v = voverz * zz;
				getTexPixel(PresentTex(), &c, u,v);
				uoverz += ustep, voverz += vstep;
			}
			//ground��ɫ
			//ambient����==diffuse����==����������ɫ
			Color dc = Ldoverz * zz;
			c = c * (Ldoverz * zz)+ 0.6 * c * (Lsoverz*zz);
			
			if (overz > getZbuffer(x, y)) {
				setPixel(x, y, c);
				writeZbuffer(x,y,overz);
			}
		}

	}

}

void vertexInterp(VertexAtrr *out,const VertexAtrr &v0, const VertexAtrr &v1,float t) {

	out->_v._x = interp(v0._v._x, v1._v._x, t);
	out->_v._y = interp(v0._v._y, v1._v._y, t);
	out->_v._z = interp(v0._v._z, v1._v._z, t);
	float BoneOverZ = v0.rhw;
	float ToneOverZ = v1.rhw;

	
	out->rhw   = interp(BoneOverZ, ToneOverZ, t);
	out->_c    = ColorInterp(v0._c,v1._c,t);
	out->_tu   = interp(v0._tu, v1._tu, t);
	out->_tv   = interp(v0._tv, v1._tv, t);
	out->ld    = ColorInterp(v0.ld, v1.ld, t);
	out->ls    = ColorInterp(v0.ls, v1.ls, t);
}

////---------------------------------
////��Χ�� �ж��������� ��դ��������|
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
//       2Dֱ�߲ü�				  |
//---------------------------------

int getcode(const Vertex4 &v) {
	int p = 0;
	if (v._x < Lx) p |= L;
	if (v._y < By) p |= B;
	if (v._y > Ty) p |= T;
	if (v._x > Rx) p |= R;
	return p;

}

void clip(int& code, VertexAtrr& va, const VertexAtrr & ua) {
	if (code&L) {
		float t = Lx - va._v._x / (ua._v._x - va._v._x);
	//	va._v._y = va._v._y + (ua._v._y - va._v._y)*(Lx - va._v._x) / (ua._v._x - va._v._x);
	//	va._v._x = Lx;
		vertexInterp(&va, va, ua,t);
	}
	else if (code&R) {
		float t = Rx - va._v._x / (ua._v._x - va._v._x);
		/*va._v._y = va._v._y + (ua._v._y - va._v._y)*(Rx - va._v._x) / (ua._v._x - va._v._x);
		va._v._x = Rx;*/
		vertexInterp(&va, va, ua, t);
	}
	else if (code&T) {
		float t = Ty - va._v._y / (ua._v._y - va._v._y);
	/*	va._v._x = va._v._x + (ua._v._x - va._v._x)*(Ty - va._v._y) / (ua._v._y - va._v._y);
		va._v._y = Ty;*/
		vertexInterp(&va, va, ua, t);
	}
	else if (code&B) {
		float t = By - va._v._y / (ua._v._y - va._v._y);
		/*va._v._x = va._v._x + (ua._v._x - va._v._x)*(By - va._v._y) / (ua._v._y - va._v._y);
		va._v._y = By;*/
		vertexInterp(&va, va, ua, t);
	}
	code = getcode(va._v);

}

void lineclip(VertexAtrr &a, VertexAtrr &b) {

	//	int x0 = a._x, int y0 = a._y, int x1 = b._x, int y1=b._y;
	int p1 = 0, p2 = 0;
	p1 = getcode(a._v); p2 = getcode(b._v);
	if (p1&p2) return;
	//if (!p1 && !p2)

	while (p1 || p2) {
		clip(p1, a, b);
		clip(p2, b, a);
	}

}

//---------------------------------
//       ������غ���			  |
//---------------------------------
void AttrMulRhw(VertexAtrr& out, const VertexAtrr& in, float rhw) {
	out._c = in._c * rhw;
	out._tu = in._tu * rhw;
	out._tv = in._tv * rhw;
	out.ld *= rhw; out.ls *= rhw;
}

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
	
	float k =1.0/ sqrt(x * x + y * y + z * z);

	out._x = in._x * k, out._y = in._y * k, out._z = in._z * k;
}

float operator * (const Vector4& lsh, const Vector4& rsh) {
	return lsh._x*rsh._x + lsh._y*rsh._y + lsh._z*rsh._z;

}

Vector4 operator - (const Vector4& rsh) {
	return Vector4(-rsh._x, -rsh._y, -rsh._z, -rsh._w);
}


//---------------------------------
//       ������غ���			  |
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



void MatrixIdentity(Matrix &out) {
	MatrixZero(out);
	out(0, 0) = 1.0f;
	out(1, 1) = 1.0f;
	out(2, 2) = 1.0f;
	out(3, 3) = 1.0f;

}
void MatrixZero(Matrix& out) {
	out(0, 0) = 0.0f;
	out(0, 1) = 0.0f;
	out(0, 2) = 0.0f;
	out(0, 3) = 0.0f;
	
	out(1, 0) = 0.0f;
	out(1, 1) = 0.0f;
	out(1, 2) = 0.0f;
	out(1, 3) = 0.0f;
	
	out(2, 0) = 0.0f;
	out(2, 1) = 0.0f;
	out(2, 2) = 0.0f;
	out(2, 3) = 0.0f;

	out(3, 0) = 0.0f;
	out(3, 1) = 0.0f;
	out(3, 2) = 0.0f;
	out(3, 3) = 0.0f;
}

Matrix& operator *= (Matrix &lsh,const Matrix &rsh) {
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


Vector4& operator *= (Vector4 &lrh, const Matrix &rsh) {
	lrh = lrh * rsh;
	return lrh;

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

void MatrixTranslation(Matrix& out, float x, float y, float z) {
	out(0, 0) = 1.0f;
	out(1, 1) = 1.0f;
	out(2, 2) = 1.0f;
	out(3, 3) = 1.0f;
	out(3, 0) = x;
	out(3, 1) = y;
	out(3, 2) = z;
}


//---------------------------------
//      ��������				  |
//---------------------------------


void getNormal( Triangle * tri) {
	Vector4 u = tri->v1._v - tri->v0._v;
	Vector4 v = tri->v2._v - tri->v0._v;
	VecCross(&tri->normal, u, v);
	VecNormalize(tri->normal, tri->normal);
}

void getNormal(Vector4 *out, const Vertex4& a, const Vertex4& b, const Vertex4& c) {
	Vector4 u = b - a;
	Vector4 v = c - a;

	VecCross(out, u, v);
}



//����Ϊ��ʱ�뷽�򹹳�
bool backCull(const Vertex4& a,const Vertex4& b,const Vertex4& c) {
	Vector4 normal;
	getNormal(&normal, a, b, c);
	Vector4 zero(0, 0, 0, 0);
	Vector4 view =  zero - a;
	float cosT = view * normal;
	if (cosT >= 0) 
		return false;
	return true;
}

bool backCull(const Triangle& tri) {
	Vector4 zero(0, 0, 0, 0);
	Vector4 view = zero - tri.v0._v;
	float cosT = view * tri.normal;
	if (cosT >= 0)
		return false;
	return true;

}

//---------------------------------
//       ���Բ�ֵ				  |
//---------------------------------

float interp(float a,float b,float t) {
	return (b - a)*t + a;

}
Color ColorInterp(const Color& c0, const Color& c1, float t) {
	Color out;
	out._b = interp(c0._b, c1._b, t);
	out._r = interp(c0._r, c1._r, t);
	out._g = interp(c0._g ,c1._g, t);
	return out;
}



//---------------------------------
//      �ļ�����				  |
//---------------------------------

void draw(const char *name) {
	std::ofstream fout(name,std::ofstream::out|std::ofstream::binary);
	fout << "P6"<<std::endl << Width << " " << Height << std::endl<<"255"<<std::endl;
	

	for (int i = 0; i <= Height; i++) {
		char *p = img + (i * Width) * 3;
//		for (int j = 0; j < Width; j++) {

			fout.write(p, Width*3);
	//	}
	}
	fout.close();

}



void getTexture(const char* file) {
		
	std::ifstream fin(file, std::ifstream::in | std::ifstream::binary);
			
	fin.seekg(0, fin.end);//�ļ��ڲ�ָ�붨λ���ļ�ĩβ
	int length = fin.tellg();	//��ȡ�ڲ�ָ���ڵ�ǰλ�õ�ƫ����(�ļ��ĳ���)
	fin.seekg(0, fin.beg);	//���¶�λ���ļ���
	
	char *buffer = new char[length];
	fin.read(buffer, length);
	


	/*std::ofstream os("10.ppm", std::ofstream::out | std::ostream::binary);
	os.write(buffer, length);
	delete[] buffer;*/

}


//---------------------------------
//			  ���մ���			  |
//---------------------------------
Light initDirectionalLight(const Vector4& direction, const Color &color) {
	Light light;
	light.type = LightDirectional;
	
	light.direction = direction;

	light._Ambient = 0.2*color;
	light._Specular = color;
	light._Diffuse = color;
	return light;
}
void getSpecuC(Color* out, const Light& Light,const Matreial&M, const Vector4& view, const Vector4& normal) {
	Vector4 l = -Light.direction;
	Vector4 h = view + l;
	VecNormalize(h, h);
	float cost = std::max(h * normal, 0.0f);
	if (cost > 0)
		cost = powf(cost, M._power);
	*out = Light._Specular * cost  ;
}
//void getAmbientC(Color* out, const Matreial& M, const Light& Light) {
//	*out = Light._Ambient * M._Ambient;
//}

void getdiffuseC(Color* out, const Matreial& M, const Light& Light, const Vector4& normal) {
	float cost = std::max(-Light.direction * normal, 0.0f);
	*out = cost * Light._Diffuse;
}


Light* Lightarr() {
	static Light* Lightarr = new Light[10]();
	return Lightarr;
}
void setLight(int num, Light* light) {
	*(Lightarr() +num) = *light;
}
Light* getLight(int index) {
	return Lightarr() + index;
}
void releaseLight() {
	delete[] Lightarr();
}

void lightApplyMatrix(Light* light, const Matrix&M) 
{	
	light->position *= M;
	light->direction *= M;
}

void VertexShader(VertexAtrr &v,int lightNum) {
	v.ld = Black;
	v.ls = Black;
	Vector4 view = -v._v;
	VecNormalize(view, view);
	for (int i = 0;i<lightNum ; i++) {
		Color ls,ld;
		getSpecuC(&ls, *getLight(i),currentMatreial(),view,v.normal);
		getdiffuseC(&ld, currentMatreial(), *getLight(i), v.normal); 
		if (getLight(i)->type & (LightSpot | LightPoint)) {
			//˥������
		}
		v.ls += ls;
		v.ld += ld;
	}
	v.ld += getLight(0)->_Ambient;
}

//---------------------------------
//			Z-buffer			  |
//---------------------------------
float *Z_buffer;

void setZbuffer(int w, int h) {
	Z_buffer = new float[w*h]();
	//std::fill(Z_buffer, Z_buffer + w * h, 1.0f);
}

void releaseZbuffer() {
	delete[] Z_buffer;
}

void writeZbuffer(int x, int y, float z) {
	*(Z_buffer + y*Width+x) = z;
}
float getZbuffer(int x, int y) {
	return *(Z_buffer + y*Width+x);
}
//---------------------------------
//			��������			  |
//---------------------------------

void setTriNormal(Triangle& tri, const Vector4& Normal) {
	tri.v0.normal = Normal;
	tri.v1.normal = Normal;
	tri.v2.normal = Normal;
}



//---------------------------------
//			��������			  |
//---------------------------------

Matreial& currentMatreial() {
	static Matreial *CM = new Matreial();
	return *CM;
}

void setMatreial(Matreial *Mtr) {
	currentMatreial() = *Mtr;
}

void releaseCM() {
	delete &currentMatreial();
}


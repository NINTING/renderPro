#include"yRasterization.h"

class Device {
public:
	int screenW, screenH;
	

	Matrix tranMatrix;							//Õ∂”∞æÿ’Û
	Matrix ProjectMatr, ViewMatr, WorldMatr;	//world*view*pro
	Texture *texture;
	//int texw, texh;
	Light** lightArr;		//π‚’’
	int LCount;
	Matreial* presentMat;	//≤ƒ÷ 
	float* Z_buffer;		//zª∫≥Â
	char* frameBuffer;		//÷°ª∫¥Ê
	
	void setTexture(Texture* tex);
	void setframeBuffer(char *img);
	void setMatreial(Matreial* pm);
	void setLight(int index, Light* light);
	void setMatrix(int flag, Matrix* m);
	void getTranMatrix();
	void release();

};

void initDevice(Device*dev,int w,int h);


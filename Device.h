#include"yRasterization.h"

class Device {
public:
	int screenW, screenH;
	

	Matrix tranMatrix;							//ͶӰ����
	Matrix ProjectMatr, ViewMatr, WorldMatr;	//world*view*pro
	Texture *texture;
	//int texw, texh;
	Light** lightArr;		//����
	int LCount;
	Matreial* presentMat;	//����
	float* Z_buffer;		//z����
	char* frameBuffer;		//֡����
	
	void setTexture(Texture* tex);
	void setframeBuffer(char *img);
	void setMatreial(Matreial* pm);
	void setLight(int index, Light* light);
	void setMatrix(int flag, Matrix* m);
	void getTranMatrix();
	void release();

};

void initDevice(Device*dev,int w,int h);


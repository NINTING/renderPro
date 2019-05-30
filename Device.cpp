#include"Device.h"


void Device::setTexture(Texture* tex) {
	//texture->width = w; tex->height = h;
	texture = tex;
}
void Device:: setframeBuffer(char* img) {
	frameBuffer = img;
}
void Device:: setMatreial(Matreial* pm) {
	presentMat = pm;
}
void Device::setLight(int index, Light* light) {
	*(lightArr + index) = light;
}
void Device::setMatrix(int flag, Matrix* m) {
	if (flag & TS_WORLD)
		WorldMatr = *m;
	if (flag & TS_VIEW)
		ViewMatr = *m;
	if (flag & TS_PROJECTION)
		ProjectMatr = *m;
}
void Device:: getTranMatrix() {
	tranMatrix = WorldMatr * ViewMatr;
}
void Device ::release() {
	delete[] frameBuffer;
	delete[] Z_buffer;
	delete[] texture->tex;
	delete texture;
	delete lightArr;
	delete presentMat;	
	LCount = 0;
}



void initDevice(Device* dev, int w, int h) {
	dev->screenW = w;
	dev->screenH = h;
	dev->Z_buffer = new float[w * h];
	
}

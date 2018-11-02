#ifndef _DISPARITY_IO_H
#define _DISPARITY_IO_H

// value to use to represent unknown disparity
#define UNKNOWN_DISP 1e9

class CDispIO
{
public:
	// return whether disparity is valid
	template <class T>
	static bool InValid(T d){ return d < 0 || d >= UNKNOWN_DISP; };

// 	// read a flow file into 2-band image
// 	template <class T>
// 	static int ReadFlowFile(T* U, T* V, int* w, int* h, const char* filename);
// 
// 	// write a 2-band image into flow file 
// 	template <class T>
// 	static int WriteFlowFile(T* U, T* V, int w, int h, const char* filename);
// 
// 	// read a KITTI flow file into 2-band image
// 	template <class T>
// 	static int ReadKittiFlowFile(T* U, T* V, int* w, int* h, const char* filename);
// 
// 	// write a 2-band image into KITTI flow file 
// 	template <class T>
// 	static int WriteKittiFlowFile(T* U, T* V, int w, int h, const char* filename);

	// render the disparity to a 4-band BGRA color image
	template <class T>
	static double DispToColor(unsigned char* fillPix, T* U, T* V, int w, int h, float range = -1);

	template <class T>
	static void ShowDisp(const char* winname, T* U, T* V, int w, int h, float range = -1, int waittime = 1);
	template <class T>
	static void SaveDispAsImage(const char* imgName, T* U, T* V, int w, int h, float range = -1);

	template <class T1, class T2>
	static float CalcDispError(T1* u1, T1* v1, T2* u2, T2* v2, int w, int h);

private:

};

template <class T>
double CDispIO::DispToColor(unsigned char* fillPix, T* U, T* V, int w, int h, float range /*= -1*/)
{

// BGR color table
static unsigned char colorTbl[256][3] = {
{143,  0,  0},{147,  0,  0},{151,  0,  0},{155,  0,  0},{159,  0,  0},{163,  0,  0},{167,  0,  0},{171,  0,  0},{175,  0,  0},{179,  0,  0},{183,  0,  0},{187,  0,  0},{191,  0,  0},{195,  0,  0},{199,  0,  0},{203,  0,  0},{206,  0,  0},{210,  0,  0},{214,  0,  0},{218,  0,  0},{222,  0,  0},{226,  0,  0},{230,  0,  0},{234,  0,  0},{238,  0,  0},{242,  0,  0},{246,  0,  0},{250,  0,  0},{254,  0,  0},{255,  3,  0},{255,  7,  0},{255, 11,  0},
{255, 14,  0},{255, 18,  0},{255, 22,  0},{255, 26,  0},{255, 30,  0},{255, 34,  0},{255, 38,  0},{255, 42,  0},{255, 46,  0},{255, 50,  0},{255, 54,  0},{255, 58,  0},{255, 62,  0},{255, 66,  0},{255, 70,  0},{255, 74,  0},{255, 77,  0},{255, 81,  0},{255, 85,  0},{255, 89,  0},{255, 93,  0},{255, 97,  0},{255,101,  0},{255,105,  0},{255,109,  0},{255,113,  0},{255,117,  0},{255,121,  0},{255,125,  0},{255,129,  0},{255,133,  0},{255,137,  0},
{255,140,  0},{255,144,  0},{255,148,  0},{255,152,  0},{255,156,  0},{255,160,  0},{255,164,  0},{255,168,  0},{255,172,  0},{255,176,  0},{255,180,  0},{255,184,  0},{255,188,  0},{255,192,  0},{255,196,  0},{255,200,  0},{255,203,  0},{255,207,  0},{255,211,  0},{255,215,  0},{255,219,  0},{255,223,  0},{255,227,  0},{255,231,  0},{255,235,  0},{255,239,  0},{255,243,  0},{255,247,  0},{255,251,  0},{255,255,  0},{251,255,  4},{247,255,  8},
{244,255, 11},{240,255, 15},{236,255, 19},{232,255, 23},{228,255, 27},{224,255, 31},{220,255, 35},{216,255, 39},{212,255, 43},{208,255, 47},{204,255, 51},{200,255, 55},{196,255, 59},{192,255, 63},{188,255, 67},{184,255, 71},{181,255, 74},{177,255, 78},{173,255, 82},{169,255, 86},{165,255, 90},{161,255, 94},{157,255, 98},{153,255,102},{149,255,106},{145,255,110},{141,255,114},{137,255,118},{133,255,122},{129,255,126},{125,255,130},{122,255,133},
{118,255,137},{114,255,141},{110,255,145},{106,255,149},{102,255,153},{ 98,255,157},{ 94,255,161},{ 90,255,165},{ 86,255,169},{ 82,255,173},{ 78,255,177},{ 74,255,181},{ 70,255,185},{ 66,255,189},{ 62,255,193},{ 59,255,196},{ 55,255,200},{ 51,255,204},{ 47,255,208},{ 43,255,212},{ 39,255,216},{ 35,255,220},{ 31,255,224},{ 27,255,228},{ 23,255,232},{ 19,255,236},{ 15,255,240},{ 11,255,244},{  7,255,248},{  3,255,252},{  0,254,255},{  0,251,255},
{  0,247,255},{  0,243,255},{  0,239,255},{  0,235,255},{  0,231,255},{  0,227,255},{  0,223,255},{  0,219,255},{  0,215,255},{  0,211,255},{  0,207,255},{  0,203,255},{  0,199,255},{  0,195,255},{  0,191,255},{  0,188,255},{  0,184,255},{  0,180,255},{  0,176,255},{  0,172,255},{  0,168,255},{  0,164,255},{  0,160,255},{  0,156,255},{  0,152,255},{  0,148,255},{  0,144,255},{  0,140,255},{  0,136,255},{  0,132,255},{  0,128,255},{  0,125,255},
{  0,121,255},{  0,117,255},{  0,113,255},{  0,109,255},{  0,105,255},{  0,101,255},{  0, 97,255},{  0, 93,255},{  0, 89,255},{  0, 85,255},{  0, 81,255},{  0, 77,255},{  0, 73,255},{  0, 69,255},{  0, 65,255},{  0, 62,255},{  0, 58,255},{  0, 54,255},{  0, 50,255},{  0, 46,255},{  0, 42,255},{  0, 38,255},{  0, 34,255},{  0, 30,255},{  0, 26,255},{  0, 22,255},{  0, 18,255},{  0, 14,255},{  0, 10,255},{  0,  6,255},{  0,  2,255},{  0,  0,254},
{  0,  0,250},{  0,  0,246},{  0,  0,242},{  0,  0,238},{  0,  0,234},{  0,  0,230},{  0,  0,226},{  0,  0,222},{  0,  0,218},{  0,  0,214},{  0,  0,210},{  0,  0,206},{  0,  0,202},{  0,  0,198},{  0,  0,194},{  0,  0,191},{  0,  0,187},{  0,  0,183},{  0,  0,179},{  0,  0,175},{  0,  0,171},{  0,  0,167},{  0,  0,163},{  0,  0,159},{  0,  0,155},{  0,  0,151},{  0,  0,147},{  0,  0,143},{  0,  0,139},{  0,  0,135},{  0,  0,131},{  0,  0,128}
};

	// determine disparity range:
	double maxD = FLT_MIN;

	if (range > 0) {
		maxD = range;
	}else{	// obtain the disparity range according to the max disp.
		for (int i = 0; i < h; i++){
			for (int j = 0; j < w; j++){
				double u = U[i*w + j];
				if (!InValid(u))
					maxD = __max(maxD, u);
			}
		}
		if (maxD == 0) // if disp == 0 everywhere
			maxD = 1;
	}

	//printf("max motion: %.2f  motion range: u = [%.2f,%.2f];  v = [%.2f,%.2f]\n",
	//	maxrad, minu, maxu, minv, maxv);

	memset(fillPix, 0, 4 * w*h);
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int idx = i*w + j;
			double d = U[idx];
			if (!InValid(d)){
				unsigned char colorIdx = (d / maxD) * 255 + 0.5;
				memcpy(fillPix + idx * 4, colorTbl[colorIdx], 3);
			}
			fillPix[idx * 4 + 3] = 0xff; // alpha channel, only for alignment
		}
	}

	return maxD;
}

template <class T>
void CDispIO::ShowDisp(const char* winname, T* U, T* V, int w, int h, float range /*= -1*/, int waittime /*= 1*/)
{
	cv::Mat img(h, w, CV_8UC4);
	float maxDisp = DispToColor(img.data, U, V, w, h, range);

#if 1
	// get corner color
	int x = 10, y = 20;
	unsigned char color[4];
	unsigned char* pSrc = img.data + y*img.step + x * 4;
	color[0] = 255 - pSrc[0];
	color[1] = 255 - pSrc[1];
	color[2] = 255 - pSrc[2];
	char info[256];
	sprintf(info, "max: %.1f", maxDisp);
	cv::putText(img, info, cvPoint(x, y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cvScalar(color[0], color[1], color[2]));
#endif

	cv::imshow(winname, img);
	cv::waitKey(waittime);
}

template <class T>
void CDispIO::SaveDispAsImage(const char* imgName, T* U, T* V, int w, int h, float range /*= -1*/)
{
	cv::Mat img(h, w, CV_8UC4);
	float maxDisp = DispToColor(img.data, U, V, w, h, range);

#if 1
	// get corner color
	int x = 10, y = 20;
	unsigned char color[3];
	unsigned char* pSrc = img.data + y*img.step + x * 4;
	color[0] = 255 - pSrc[0];
	color[1] = 255 - pSrc[1];
	color[2] = 255 - pSrc[2];
	char info[256];
	sprintf(info, "max: %.1f", maxDisp);
	cv::putText(img, info, cvPoint(x, y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cvScalar(color[0], color[1], color[2]));
#endif

	cv::imwrite(imgName, img);
}

template <class T1, class T2>
float CDispIO::CalcDispError(T1* u1, T1* v1, T2* u2, T2* v2, int w, int h)
{
	return 0;
}

#endif
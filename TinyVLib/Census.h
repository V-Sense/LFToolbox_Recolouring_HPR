
#ifndef _CENSUS_H_
#define _CENSUS_H_

#define CENSUS_EPI (1e-03) //the threshold for the three-valued census transform
//#define CENSUS_EPI 0

// the number of bits that differ
int HammingDistance(unsigned char* sig1, unsigned char* sig2, int byteLen)
{
#if 1 // fast version
	int *s1 = (int*)sig1, *s2 = (int*)sig2;
	int iter = byteLen / 4;
	int diffCnt = 0;
	for (int i = 0; i < iter; i++){
		int n = (*s1)^(*s2); //XOR
		n = (n & 0x55555555) + ((n >> 1) & 0x55555555);
		n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
		n = (n & 0x0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f);
		n = (n & 0x00ff00ff) + ((n >> 8) & 0x00ff00ff);
		n = (n & 0x0000ffff) + ((n >> 16) & 0x0000ffff);
		diffCnt += n;
		s1++; s2++;
	}
	return diffCnt;

#else // slow version
	int diffCnt = 0;
	for (int i = 0; i < byteLen; i++){
		unsigned char bitDiff = sig1[i] ^ sig2[i];
		unsigned char d0 = (bitDiff >> 0) & 1;
		unsigned char d1 = (bitDiff >> 1) & 1;
		unsigned char d2 = (bitDiff >> 2) & 1;
		unsigned char d3 = (bitDiff >> 3) & 1;
		unsigned char d4 = (bitDiff >> 4) & 1;
		unsigned char d5 = (bitDiff >> 5) & 1;
		unsigned char d6 = (bitDiff >> 6) & 1;
		unsigned char d7 = (bitDiff >> 7) & 1;
		diffCnt += (d0 + d1 + d2 + d3 + d4 + d5 + d6 + d7);
	}
	return diffCnt;
#endif

}

/************************************************************************/
/*					         	LDP                                     */
/************************************************************************/
class CLdp
{
public:
	CLdp(){};
	~CLdp(){};
	void Transform(const FImage& img);
	inline UCImage* GetFeatures(){ return &m_feature; };
	void Test();
private:
	UCImage m_feature;
};

void CLdp::Transform(const FImage& img)
{
	int w = img.width();
	int h = img.height();
	int ch = img.nchannels();
	int featureCh = ch;
	if (!m_feature.matchDimension(w, h, featureCh)){
		m_feature.allocate(w, h, featureCh);
	}
	m_feature.setValue(0);

	float kernels[8][9] =
	{ { -3, -3, 5, -3, 0, 5, -3, -3, 5 },
	{ -3, 5, 5, -3, 0, 5, -3, -3, -3 },
	{ 5, 5, 5, -3, 0, -3, -3, -3, -3 },
	{ 5, 5, -3, 5, 0, -3, -3, -3, -3 },
	{ 5, -3, -3, 5, 0, -3, 5, -3, -3 },
	{ -3, -3, -3, 5, 0, -3, 5, 5, -3 },
	{ -3, -3, -3, -3, 0, -3, 5, 5, 5 },
	{ -3, -3, -3, -3, 0, 5, -3, 5, 5 } };

	FImage outConv(w, h, ch);
	for (int k = 0; k < 8; k++){
		ImageProcessing::filtering(img.pData, outConv.pData, w, h, ch, kernels[k], 1);
		//outConv.imshow("conv", 0);
		for (int i = 0; i < h; i++){
			for (int j = 0; j < w; j++){
				unsigned char* pOut = m_feature.pData + (i*w + j)*featureCh;
				for (int c = 0; c < ch; c++){
					unsigned char code = 0;
					if (outConv[(i*w + j)*ch + c] > 0){
						code = 1;
					}
					pOut[c] |= (code << k);
				}
			}
		}
	}
}

void CLdp::Test()
{

}

/************************************************************************/
/*					     Census Transform                               */
/************************************************************************/

class CCensus
{
public:
	CCensus(){};
	~CCensus(){};

	void Transform(const FImage& img);
	inline UCImage* GetFeatures(){ return &m_feature; };
	void Test();
private:
	UCImage m_feature;
};

void CCensus::Transform(const FImage& img)
{
//#define TWO_STATE

	int is5x5 = 0;	// default is 3x3

	int candiCnt = 8;
	int offset3x3[8][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 }, { -1, -1 }, { -1, 1 }, { 1, -1 }, { 1, 1 } };
	int offset5x5[24][2] = {
		{ -2, -2 }, { -2, -1 }, { -2, 0 }, { -2, 1 }, { -2, 2 },
		{ -1, -2 }, { -1, -1 }, { -1, 0 }, { -1, 1 }, { -1, 2 },
		{ 0, -2 }, { 0, -1 }, { 0, 1 }, { 0, 2 },
		{ 1, -2 }, { 1, -1 }, { 1, 0 }, { 1, 1 }, { 1, 2 },
		{ 2, -2 }, { 2, -1 }, { 2, 0 }, { 2, 1 }, { 2, 2 } };
	int (*offset)[2] = offset3x3;
	if (is5x5){
		candiCnt = 24;
		offset = offset5x5;
	}

	int w = img.width();
	int h = img.height();
	int ch = img.nchannels();

#ifdef TWO_STATE
	int bytePerItem = candiCnt / 8;
#else
	int bytePerItem = candiCnt * 2 / 8;
#endif

	int featureCh = bytePerItem * ch;
	if (!m_feature.matchDimension(w, h, featureCh)){
		m_feature.allocate(w, h, featureCh);
	}
	m_feature.setValue(0);

	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			unsigned char* pOut = m_feature.pData + (i*w + j)*featureCh;
			for (int c = 0; c < ch; c++){
				pOut += (bytePerItem * c);
#if 0
				float a = img[(i*w + j)*ch + c];
#else
				float a = img[(i*w + j)*ch + c];
				int validCnt = 1;
				for (int k = 0; k < candiCnt; k++){
					int ii = i + offset[k][0];
					int jj = j + offset[k][1];
					if (ii >= 0 && ii < h && jj >= 0 && jj < w){
						a += img[(ii*w + jj)*ch + c];
						validCnt++;
					}
				}
				a /= validCnt;
#endif
				unsigned char sig; //signature
				for (int k = 0; k < candiCnt; k++){
					float b = a;
					int ii = i + offset[k][0];
					int jj = j + offset[k][1];
					if (ii >= 0 && ii < h && jj >= 0 && jj < w){
						b = img[(ii*w + jj)*ch + c];
					}

					float diff = b - a;
					sig = 0x1;
#ifdef TWO_STATE
					if (diff < 0){
						sig = 0x0;
					}
					unsigned char* pChar = pOut + (int)(k / 8);
					int shiftCnt = (k % 8) * 1;
					(*pChar) |= (sig << shiftCnt);
#else
					if (diff > CENSUS_EPI){
						sig = 0x3;
					}
					else if (diff < -CENSUS_EPI){
						sig = 0x0;
					}
					unsigned char* pChar = pOut + (int)(k / 4);
					int shiftCnt = (k % 4) * 2;
					(*pChar) |= (sig << shiftCnt);
#endif
				}
			}
		}
	}
}

void CCensus::Test()
{
	FImage im;
	im.imread("d:/a.png");
	Transform(im);
	
	FImage tmpIm;
	m_feature.collapse(tmpIm, collapse_average);
	tmpIm.imshow("Census", 0);
}

/************************************************************************/
/*					 Scale-Robust Census Transform                      */
/************************************************************************/

class CSCensus
{
public:
	CSCensus();
	~CSCensus();

	void Init(float minScale, float maxScale, int scaleNum);
	void Transform(FImage& img);
	inline UCImage* GetFeatures(){ return m_features; };
	inline unsigned char* GetSig(int scaleIdx, int row, int col, int* sigLen);

	void Reset();
	void Test();
private:
	void ConstructCandiPts();
	UCImage* m_features;
	FImage m_outerOffsets;
	FImage m_innerOffsets;
	float* m_scales;
	int m_scaleNum;
};

CSCensus::CSCensus()
{
	m_features = NULL;
	m_scales = NULL;
	m_scaleNum = 0;
}

CSCensus::~CSCensus()
{
	Reset();
}

void CSCensus::Init(float minScale, float maxScale, int scaleNum)
{
	Reset();

	m_scaleNum = scaleNum;
	m_scales = new float[m_scaleNum];
	float step = (maxScale - minScale) / (m_scaleNum - 1);
	for (int i = 0; i < m_scaleNum; i++){
		m_scales[i] = minScale + i*step;
		//printf("%f\n", m_scales[i]);
	}
	ConstructCandiPts();
}

void CSCensus::Transform(FImage& img)
{
	int w = img.width();
	int h = img.height();
	int ch = img.nchannels();
	int featureCh = 6 * ch;
	
	// malloc space
	if (m_features){
		delete[] m_features;
		m_features = NULL;
	}
	m_features = new UCImage[m_scaleNum];
	for (int i = 0; i < m_scaleNum; i++){
		m_features[i].allocate(w, h, featureCh);
		m_features[i].setValue(0);
	}

	//
	float v[3];
	for (int s = 0; s < m_scaleNum; s++){
		UCImage* imf = &m_features[s];
		float* pi = m_innerOffsets.rowPtr(s);
		float* po = m_outerOffsets.rowPtr(s);
		for (int i = 0; i < h; i++){
			for (int j = 0; j < w; j++){
				// average of inner ring
				float innerAver[3] = { 0, 0, 0 };
				for (int k = 0; k < 12; k++){
					float dx = pi[2 * k];
					float dy = pi[2 * k + 1];
					ImageProcessing::BilinearInterpolate(img.pData, w, h, ch, j + dx, i + dy, v);
					for (int c = 0; c < ch; c++){
						innerAver[c] += v[c];
					}
				}
				for (int c = 0; c < ch; c++){
					innerAver[c] /= 12.;
				}

				// compare with outer ring
				unsigned char* pOut = imf->pData + (i*w + j)*featureCh;
				unsigned char sig; //signature
				for (int k = 0; k < 24; k++){
					float dx = po[2 * k];
					float dy = po[2 * k + 1];
					ImageProcessing::BilinearInterpolate(img.pData, w, h, ch, j + dx, i + dy, v);
					for (int c = 0; c < ch; c++){
						float diff = v[c] - innerAver[c];
						sig = 0x1;
						if (diff > CENSUS_EPI){
							sig = 0x3;
						}else if (diff < -CENSUS_EPI){
							sig = 0x0;
						}
						unsigned char* pChar = pOut + (int)(k / 4);
						int shiftCnt = (k % 4) * 2;
						(*pChar) |= (sig << shiftCnt);
					}
				}
			}
		}
	}
}

inline unsigned char* CSCensus::GetSig(int scaleIdx, int row, int col, int* sigLen)
{
	*sigLen = m_features[scaleIdx].nchannels();
	return m_features[scaleIdx].pixPtr(row, col);
}

void CSCensus::Reset()
{
	if (m_features){
		delete[] m_features;
		m_features = NULL;
	}
	if (m_scales){
		delete[] m_scales;
		m_scales = NULL;
	}
}

void CSCensus::Test()
{
	CTimer t;

	t.tic();
	Init(0.5, 2, 2);
	t.toc("SCensus init: ");

	FImage tmp(200, 200, 3);
	tmp.setValue(0);
	int x = 100, y = 100;
	float iv[3] = { 0, 0, 1 };
	float ov[3] = { 0, 1, 0 };
	for (int i = 0; i < 12; i++){
		float* pi = m_innerOffsets.rowPtr(0);
		float* po = m_outerOffsets.rowPtr(0);
		float nx = pi[2 * i] + x;
		float ny = pi[2 * i + 1] + y;
		tmp.setPixel(nx, ny, iv);
		float ox1 = po[4 * i] + x;
		float oy1 = po[4 * i + 1] + y;
		float ox2 = po[4 * i + 2] + x;
		float oy2 = po[4 * i + 3] + y;
		tmp.setPixel(ox1, oy1, ov);
		tmp.setPixel(ox2, oy2, ov);
	}
	tmp.imshow("SCensus", 0);

	FImage img;
	img.imread("d:/a.png");
	img.desaturate();

	t.tic();
	Transform(img);
	t.toc("SCensus transform: ");
}

// input is m_scales and scaleNum
// output is m_innerOffsets and m_outerOffsets
void CSCensus::ConstructCandiPts()
{
	// 360C / 12
	float candiAngles[12] = { 0, 30 * PI / 180, 60 * PI / 180, 90 * PI / 180, 120 * PI / 180, 150 * PI / 180, PI,
		210 * PI / 180, 240 * PI / 180, 270 * PI / 180, 300 * PI / 180, 330 * PI / 180 };

	m_innerOffsets.allocate(12*2, m_scaleNum);
	m_outerOffsets.allocate(24*2, m_scaleNum);

	int iidx, oidx; // inner and outer index
	for (int i = 0; i < m_scaleNum; i++){
		float s = m_scales[i];
		iidx = 0; oidx = 0;
		float* pi = m_innerOffsets.rowPtr(i);
		float* po = m_outerOffsets.rowPtr(i);
		for (int k = 0; k < 12; k++){
			float a = candiAngles[k];
			pi[iidx++] = cos(a) * s / 4.; // dx
			pi[iidx++] = sin(a) * s / 4.; // dy
			po[oidx++] = cos(a) * s; // dx
			po[oidx++] = sin(a) * s; // dy
			po[oidx++] = cos(a) * 2 * s; // dx
			po[oidx++] = sin(a) * 2 * s; // dy
		}
	}
}

#endif
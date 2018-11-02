
#include "Image.h"
#include "OpticFlowIO.h"

//extern "C" {
//	__declspec(dllexport) int match_number(char* fileName);
//	__declspec(dllexport) int match_density(char* fileName, int w, int h);
//	__declspec(dllexport) int match_accuracy(char* fileName, char* gtName, int w, int h);
//}

/* read matches, stored as x1 y1 x2 y2 per line (other values on the same is not taking into account */
void ReadMatches(const char *filename, FImage& outMat)
{
	float* tmp = new float[4 * 1000000]; // max number of match pair
	FILE *fid = fopen(filename, "r");
	int nmatch = 0;
	float x1, x2, y1, y2;
	while (!feof(fid) && fscanf(fid, "%f %f %f %f%*[^\n]", &x1, &y1, &x2, &y2) == 4){
		tmp[4 * nmatch] = x1;
		tmp[4 * nmatch + 1] = y1;
		tmp[4 * nmatch + 2] = x2;
		tmp[4 * nmatch + 3] = y2;
		nmatch++;
	}
	outMat.allocate(4, nmatch);
	memcpy(outMat.pData, tmp, nmatch * 4 * sizeof(float));
	fclose(fid);
	delete[] tmp;
}

void Match2Flow(FImage& inMat, FImage& ou, FImage& ov, int w, int h)
{
	if (!ou.matchDimension(w, h, 1)){
		ou.allocate(w, h, 1);
	}
	if (!ov.matchDimension(w, h, 1)){
		ov.allocate(w, h, 1);
	}
	ou.setValue(UNKNOWN_FLOW);
	ov.setValue(UNKNOWN_FLOW);
	int cnt = inMat.height();
	for (int i = 0; i < cnt; i++){
		float* p = inMat.rowPtr(i);
		int x = p[0];
		int y = p[1];
		int u = p[2] - p[0];
		int v = p[3] - p[1];

		for (int ii = -5; ii <= 4; ii++){
			for (int jj = -5; jj <= 4; jj++){
				int ni = ImageProcessing::EnforceRange(y + ii, h);
				int nj = ImageProcessing::EnforceRange(x + jj, w);
				ou[ni*w + nj] = u;
				ov[ni*w + nj] = v;

			}
		}
		//ou[y*w + x] = u;
		//ov[y*w + x] = v;
	}
}

float MatchPrecision(FImage& inMat, FImage& gtu, FImage& gtv)
{
	int validCnt = 0;
	int w = gtu.width();
	int h = gtu.height();
	int cnt = inMat.height();
	for (int i = 0; i < cnt; i++){
		float* p = inMat.rowPtr(i);
		int x = p[0];
		int y = p[1];
		if (x < 0 || x >= w || y < 0 || y >= h){
			continue;
		}
		int u = p[2] - p[0];
		int v = p[3] - p[1];

		float gu = gtu[y*w + x];
		float gv = gtv[y*w + x];

		float diff = sqrt((u - gu)*(u - gu) + (v - gv)*(v - gv));
		if (diff <= 5){
			validCnt++;
		}
	}
	return (float)validCnt / cnt;
}

int match_number(char* fileName)
{
	FImage mat;
	ReadMatches(fileName, mat);
	return mat.height();
}

int match_density(char* fileName, int w, int h)
{
	FImage mat;
	ReadMatches(fileName, mat);

	FImage mu, mv;
	Match2Flow(mat, mu, mv, w, h);

	int validCnt = 0;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			float u = mu[i*w + j];
			float v = mv[i*w + j];
			if (!OpticFlowIO::unknown_flow(u, v)){
				validCnt++;
			}
		}
	}
	float density = (float)validCnt / (w*h);

	return density * 10000;
}

int match_accuracy(char* fileName, char* gtName, int w, int h)
{
	FImage mat;
	ReadMatches(fileName, mat);
	FImage gtu(w, h), gtv(w, h);
	int gtw, gth;
	OpticFlowIO::ReadFlowFile(gtu.pData, gtv.pData, &gtw, &gth, gtName);
	if (gtw != w || gth != h){
		printf("Resolution is UnMatching !\n");
	}
	float accu = MatchPrecision(mat, gtu, gtv);
	return accu * 10000;
}

#ifndef _SUPER_PIXEL_H
#define _SUPER_PIXEL_H

#include <stdio.h>
#include <assert.h>
#include <stack>
#include <cmath>

#include "Image.h"
#include "SLIC.h"
#include "Heap.h"
#include "Util.h"
#include "OpticFlow.h"
#include "dt.h"
//#include "Clustering.h"

class CSuperPixelGraph
{
public:
	CSuperPixelGraph();
	~CSuperPixelGraph();

	int Init(FImage& img, const int K, const double compactness); // SLIC

	int Init(char* fileName, int width, int height);
	int Init(int* klabels, int width, int height);

	inline int NodeCnt(){ return m_nodeCnt; };
	inline int NodeSize(int nodeIdx){ return m_nodeSizes[nodeIdx]; };
	inline int* NodeItems(int nodeIdx) { return m_nodeItems[nodeIdx]; };
	inline int NeighborCnt(int nodeIdx){ return m_neighborCnts[nodeIdx]; };
	inline int* NeighborIndex(int nodeIdx){ return m_neighbors[nodeIdx]; };
	inline float* NeighborDis(int nodeIdx){ return m_nbDistances[nodeIdx]; };
	inline int* NodePos(int nodeIdx){ return m_centerPositions + nodeIdx * 2; };
	inline int* KLabels(){ return m_labels; };
	inline int* BorderFlags(){ return m_borderFlags; };
	inline int MaxNeighborCnt(){ return m_maxNeighBor; };

	int ComputeNeighborDis(FImage& img, FImage& costMap);
	int ComputeNeighborDis(FImage& disMap);
	void ComputeNeighborDis();
	
	int FindNN(int nodeIdx, int* validFlag, int* nnIdx, float* nnWeight, int maxCnt); // find nearest neighbor of some super pixel
	int FindNN_GD(int nodeIdx, int* validFlag, int* nnIdx, float* nnWeight, int maxCnt); // geodesic distance

	void AddSuperPixelMask(FImage& image);
	void ShowSuperPixel(FImage& image, char* winName);
	void ShowNN(char* winName, FImage& image, int nodeIdx, int* validFlag);
	void ShowNNClustering(FImage& image, int nodeIdx, FImage& seedsFlow);
	void ShowNNFlow(FImage& image, int nodeIdx, FImage& seedsFlow);

	static int ArrangeLabels(int*kLabMat, int width, int height);

public:
	void ComputeAverage(FImage& featureImg);
	void GetAverageItem(FImage& featureImg, int nodeIdx, float* c1, float* c2, float* c3);
	float GetWeight(int spIdx1, int spIdx2);
	float GetDistance(int spIdx1, int spIdx2);

	float GetColorDis(int spIdx1, int spIdx2);
	float GetSpaceDis(int spIdx1, int spIdx2);

	void ClearAllNeighbors();
public:
	void Reset();
	void AddNeighbor(int spIdx1, int spIdx2, float dis);
	void AddNeighbor0(int srcIdx, int dstIdx, float dis);
	void SetNeighborDis(int spIdx1, int spIdx2, float dis);
	void SetNeighborDis0(int spIdx1, int spIdx2, float dis);

	int* m_labels;
	int m_width, m_height;

	int m_nodeCnt;
	int** m_neighbors;
	float** m_nbDistances;
	int** m_nodeItems;
	int* m_nodeItemStore;
	int* m_neighborCnts;
	int* m_centerPositions;
	int* m_nodeSizes;
	int* m_borderFlags;

	float* m_averBGR; // save the average pixel value of every super pixel
	float* m_averLAB; 

	int m_maxNeighBor; // max neighbor number of each SuperPixel
};

int CSuperPixelGraph::Init(FImage& img, const int K, const double compactness)
{
	int* kLabels;
	SLIC slic;
	slic.OverSegmentaion(img, kLabels, K, compactness);
	Init(kLabels, img.width(), img.height());
	ComputeAverage(img);
	delete[] kLabels;
	return m_nodeCnt;
}

int CSuperPixelGraph::Init(char* fileName, int width, int height)
{
	FILE* fid = fopen(fileName, "rb");
	if (fid == NULL){
		fprintf(stderr, "open '%s' for CSuperPixel error!\n");
		return -1;
	}
	int nodeCnt = -1;
	int* data = new int[width*height];
	if (fread(data, sizeof(int), width*height, fid) == width*height){
		nodeCnt = Init(data, width, height);
	}
	delete[] data;
	fclose(fid);
	return nodeCnt;
}

int CSuperPixelGraph::Init(int* klabels, int width, int height)
{
	Reset();

	// re-arrange labels
	ArrangeLabels(klabels, width, height);

	// save label
	m_labels = new int[width*height];
	memcpy(m_labels, klabels, width*height*sizeof(int));

	//
	m_width = width;
	m_height = height;

	// check input labels
	int w = m_width;
	int h = m_height;
	m_nodeCnt = -1;
	for(int i=0; i<w*h; i++){
		if(m_labels[i]<0 || m_labels[i]>=1000*1000){
			printf("BAD superpixel label input: %d !\n", m_labels[i]);
			return -1;
		}
		if(m_labels[i] > m_nodeCnt){
			m_nodeCnt = m_labels[i];
		}
	}
	m_nodeCnt++;

	// construct the neighbor systems
	m_neighborCnts = new int[m_nodeCnt];
	memset(m_neighborCnts, 0, m_nodeCnt*sizeof(int));

	m_neighbors = (int**)malloc(m_nodeCnt*sizeof(int*));
	m_neighbors[0] = (int*)malloc(m_nodeCnt*m_maxNeighBor*sizeof(int));
	for(int i=1; i<m_nodeCnt; i++){
		m_neighbors[i] = m_neighbors[i-1] + m_maxNeighBor;
	}

	m_nbDistances = (float**)malloc(m_nodeCnt*sizeof(float*));
	m_nbDistances[0] = (float*)malloc(m_nodeCnt*m_maxNeighBor*sizeof(float));
	for (int i = 1; i < m_nodeCnt; i++){
		m_nbDistances[i] = m_nbDistances[i - 1] + m_maxNeighBor;
	}

	m_nodeItems = (int**)malloc(m_nodeCnt*sizeof(int*));
	m_nodeItemStore = (int*)malloc(w*h*2*sizeof(int));

	int nbOffset[4][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 } };
	//int nbOffset[8][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 }, { -1, -1 }, { -1, 1 }, { 1, -1 }, { 1, 1 } };
	m_centerPositions = new int[m_nodeCnt*2];

	m_borderFlags = new int[m_width*m_height];
	for (int i = 0; i < m_width*m_height; i++){
		m_borderFlags[i] = -1;
	}
	m_nodeSizes = new int[m_nodeCnt];
	memset(m_nodeSizes, 0, m_nodeCnt*sizeof(int));

	int* xvec = new int[w*h]; // simple stack
	int* yvec = new int[w*h];

	int* itemPos = m_nodeItemStore;
	int sumX, sumY;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int seedLabel = m_labels[i*w + j];
			if (seedLabel >= 0){
				//printf("%d ", seedLabel);
				sumX = j; sumY = i;
				m_labels[i*w + j] = -1; // flag

				xvec[0] = j;
				yvec[0] = i;

				m_nodeItems[seedLabel] = itemPos;
				itemPos[0] = j;
				itemPos[1] = i;
				itemPos += 2;

				int count = 1;
				for (int c = 0; c < count; c++){
					for (int n = 0; n < 4; n++){
						int x = xvec[c] + nbOffset[n][0];
						int y = yvec[c] + nbOffset[n][1];
						if ((x >= 0 && x < w) && (y >= 0 && y < h)){
							int nbIdx = y*w + x;
							if (m_labels[nbIdx] >= 0 && klabels[nbIdx] == seedLabel){
								xvec[count] = x;
								yvec[count] = y;
								count++;

								m_labels[nbIdx] = -1; // flag

								itemPos[0] = x;
								itemPos[1] = y;
								itemPos += 2;

								sumX += x;
								sumY += y; // for center calc;
							}
							if (klabels[nbIdx] != seedLabel) { // get neighbor
								int nbNum = m_neighborCnts[seedLabel];

								// exist?
								bool exist = false;
								for (int a = 0; a < nbNum; a++){
									if (m_neighbors[seedLabel][a] == klabels[nbIdx]){
										exist = true;
										break;
									}
								}
								if (!exist){	// added if not exist
									if (nbNum >= m_maxNeighBor){ // safety check
										printf("SuperPixel Construction Error: too many neighbors of some superpixels !\n");
										continue;
									}
									m_neighbors[seedLabel][nbNum] = klabels[nbIdx];
									m_neighborCnts[seedLabel]++;
								}

								// set border flag
								m_borderFlags[y*w + x] = seedLabel;
							}
						}
					}
				}
				m_nodeSizes[seedLabel] = count;
				m_centerPositions[seedLabel * 2 + 0] = (int)((float)sumX / count + 0.5);//(left+right)/2;
				m_centerPositions[seedLabel * 2 + 1] = (int)((float)sumY / count + 0.5);//(top+bottom)/2;
			}
		}
	}
	memcpy(m_labels, klabels, width*height*sizeof(int));

	delete[] xvec;
	delete[] yvec;

#if 1
	// check
	int cnt = NodeCnt();
	int sum = 0;
	for (int i = 0; i < cnt; i++){
		sum += NodeSize(i);
	}
	printf("%d / %d\n", sum, m_width*m_height);
	assert(sum == m_width*m_height);
#endif
	return m_nodeCnt;
}

CSuperPixelGraph::~CSuperPixelGraph()
{
	Reset();
}

CSuperPixelGraph::CSuperPixelGraph()
{
	memset(this, 0, sizeof(*this));
	m_maxNeighBor = 64;
}

void CSuperPixelGraph::Reset()
{
	if(m_labels){
		delete[] m_labels;
		m_labels = NULL;
	}
	if(m_neighbors){
		free(m_neighbors[0]);
		free(m_neighbors);
		m_neighbors = NULL;
	}
	if (m_nbDistances){
		free(m_nbDistances[0]);
		free(m_nbDistances);
		m_nbDistances = NULL;
	}
	if (m_nodeItems){
		free(m_nodeItemStore);
		free(m_nodeItems);
		m_nodeItems = NULL;
	}
	if(m_neighborCnts){
		delete[] m_neighborCnts;
		m_neighborCnts = NULL;
	}
	if(m_centerPositions){
		delete[] m_centerPositions;
		m_centerPositions = NULL;
	}
	if (m_nodeSizes){
		delete[] m_nodeSizes;
		m_nodeSizes = NULL;
	}
	if (m_borderFlags){
		delete[] m_borderFlags;
		m_borderFlags = NULL;
	}
	if (m_averBGR){
		delete[] m_averBGR;
		m_averBGR = NULL;
	}
	if (m_averLAB){
		delete[] m_averLAB;
		m_averLAB = NULL;
	}
}

void CSuperPixelGraph::AddNeighbor(int spIdx1, int spIdx2, float dis)
{
	AddNeighbor0(spIdx1, spIdx2, dis);
	AddNeighbor0(spIdx2, spIdx1, dis);
}

void CSuperPixelGraph::AddNeighbor0(int srcIdx, int dstIdx, float dis)
{
	int nbCnt = NeighborCnt(srcIdx);
	int* nbIdx = NeighborIndex(srcIdx);
	float* nbDis = NeighborDis(srcIdx);

	nbIdx[nbCnt] = dstIdx;
	nbDis[nbCnt] = dis;
	m_neighborCnts[srcIdx]++;
}

void CSuperPixelGraph::SetNeighborDis(int spIdx1, int spIdx2, float dis)
{
	SetNeighborDis0(spIdx1, spIdx2, dis);
	SetNeighborDis0(spIdx2, spIdx1, dis);
}

void CSuperPixelGraph::SetNeighborDis0(int spIdx1, int spIdx2, float dis)
{
	int nbCnt = NeighborCnt(spIdx1);
	int* nbIdx = NeighborIndex(spIdx1);
	float* nbDis = NeighborDis(spIdx1);
	for (int i = 0; i < nbCnt; i++){ // find the neighbor
		if (nbIdx[i] == spIdx2){
			nbDis[i] = dis;
			return;
		}
	}
}

void CSuperPixelGraph::ComputeAverage(FImage& featureImg)
{
	if (m_averBGR){
		delete[] m_averBGR;
		m_averBGR = NULL;
	}
	if (m_averLAB){
		delete[] m_averLAB;
		m_averLAB = NULL;
	}
	m_averBGR = (float*)malloc(m_nodeCnt * 3 * sizeof(float));
	m_averLAB = (float*)malloc(m_nodeCnt * 3 * sizeof(float));
	
	float bgr[3];
	float lab[3];
	for (int i = 0; i < m_nodeCnt; i++){
		GetAverageItem(featureImg, i, bgr, bgr + 1, bgr + 2);
		m_averBGR[3 * i + 0] = bgr[0];
		m_averBGR[3 * i + 1] = bgr[1];
		m_averBGR[3 * i + 2] = bgr[2];
		ImageProcessing::BGR2Lab(bgr, lab, 1, 1);
		m_averLAB[3 * i + 0] = lab[0];
		m_averLAB[3 * i + 1] = lab[1];
		m_averLAB[3 * i + 2] = lab[2];
	}
}

// Fast Approximate MBD (Minimum Barrier Distance) Transform
void MinimumBarrier_DT(FImage& img, float* seedsX, float* seedsY, int cnt, int* labels, float* dmap, float factor, int iter = 4)
{
#define GRAY

	int w = img.width();
	int h = img.height();
	int ch = img.nchannels();

	FImage g(w, h, ch);
	g.copyData(img);
	//img.smoothing(g, 2.0);
	//img.BoxFilter(g, 2);
	//g.imshow("g", 0);

	FImage U(w, h, ch), L(w, h, ch);

	//
	U.copyData(g);
	L.copyData(g);
	memset(labels, 0xFF, w*h*sizeof(int)); // -1
	memset(dmap, 0x7F, w*h*sizeof(float)); // MAX_FLT
	//FImage tmp(w, h);
	for (int i = 0; i < cnt; i++){
		int x = seedsX[i] + 0.5;
		int y = seedsY[i] + 0.5;
		int cIdx = y*w + x;
		dmap[cIdx] = 0;
		labels[cIdx] = i;
		//
		//tmp[cIdx] = 1;
	}
	//tmp.imshow("center", 0);

	for (int n = 0; n < iter; n++)
	{
		int startX = 0, endX = w;
		int startY = 0, endY = h;
		int step = 1;
		int ox[2], oy[2]; //offset
		ox[0] = 0; oy[0] = -1;
		ox[1] = -1; oy[1] = 0;
		if (n % 2 == 1){
			startX = w - 1;	endX = -1;
			startY = h - 1;	endY = -1;
			step = -1;
			ox[0] = 0; oy[0] = 1;
			ox[1] = 1; oy[1] = 0;
		}

		for (int i = startY; i != endY; i += step){
			for (int j = startX; j != endX; j += step){
				int idx = i*w + j;
				for (int k = 0; k <= 1; k++){
					int candix = j + ox[k];
					int candiy = i + oy[k];
					if (candix >= 0 && candix < w && candiy >= 0 && candiy < h){
						int canIdx = candiy*w + candix;
						int sIdx = labels[canIdx];
						if (sIdx >= 0){
#ifdef GRAY
							float maxCost = max(U[canIdx], g[idx]);
							float minCost = min(L[canIdx], g[idx]);
							float cDis = maxCost - minCost;
#else
							float maxCost[3], minCost[3], colorDis[3];
							maxCost[0] = __max(U[3 * canIdx + 0], g[3 * idx + 0]);
							maxCost[1] = __max(U[3 * canIdx + 1], g[3 * idx + 1]);
							maxCost[2] = __max(U[3 * canIdx + 2], g[3 * idx + 2]);
							minCost[0] = __min(L[3 * canIdx + 0], g[3 * idx + 0]);
							minCost[1] = __min(L[3 * canIdx + 1], g[3 * idx + 1]);
							minCost[2] = __min(L[3 * canIdx + 2], g[3 * idx + 2]);
							colorDis[0] = maxCost[0] - minCost[0];
							colorDis[1] = maxCost[1] - minCost[1];
							colorDis[2] = maxCost[2] - minCost[2];
							float cDis = __max(__max(colorDis[0], colorDis[1]), colorDis[2]);
							//float cDis = colorDis[0] + colorDis[1] + colorDis[2];
#endif
							float sDis2 = (seedsX[sIdx] - j)*(seedsX[sIdx] - j) + (seedsY[sIdx] - i)*(seedsY[sIdx] - i);
							float dis = cDis*cDis + factor*sDis2;

							//float sDis = len[canIdx] + 1;
							//float dis = cDis*cDis + m*m*((sDis*sDis) / spSize);

							if (dis < dmap[idx]){
								labels[idx] = sIdx;
								dmap[idx] = dis;
#ifdef GRAY
								U[idx] = maxCost;
								L[idx] = minCost;
#else
								memcpy(U.pData + 3 * idx, maxCost, 3 * sizeof(float));
								memcpy(L.pData + 3 * idx, minCost, 3 * sizeof(float));
#endif
							}
						}
					}
				}
			}
		}
	}
}

// a border between two int (neighbor) is represented as a long
static inline int64_t Key(int i, int j)
{
	// swap: always i<j
	if (i > j) { int t = i; i = j; j = t; };
        return (int64_t)i + ((int64_t)j << 32);
}
static inline int KeyFirst(int64_t i) { return int(i); }
static inline int KeySecond(int64_t i) { return int(i >> 32); }

int CSuperPixelGraph::ComputeNeighborDis(FImage& img, FImage& costMap)
{
	int w = costMap.width();
	int h = costMap.height();
	int sCnt = m_nodeCnt;

	//
	FImage dis(w, h);
	IntImage labels(w, h);

#if 0
	float* sx = new float[sCnt];
	float* sy = new float[sCnt];
	for (int i = 0; i < sCnt; i++) {
		sx[i] = m_centerPositions[2 * i];
		sy[i] = m_centerPositions[2 * i + 1];
	}
	MinimumBarrier_DT(img, sx, sy, sCnt, labels.pData, dis.pData, 0, 6);// 0.0001);
	delete[] sx;
	delete[] sy;
#else
	memset(dis.pData, 0x7F, w*h*sizeof(float)); // init approx. to FLT_MAX
	memset(labels.pData, 0xFF, w*h*sizeof(int)); // init to -1;
	for (int i = 0; i < sCnt; i++) {
		int x = m_centerPositions[2 * i];
		int y = m_centerPositions[2 * i + 1];
		int p = y*w + x;
		if (labels[p] >= 0){ // the same seeds
			//printf("ERROR in WDT.\n"); 
			// search in four directions for new center
			int off[4][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 } };
			float sc = 1;
			while (1){
				for (int k = 0; k < 4; k++){
					int tx = x + off[k][0]*sc;
					int ty = y + off[k][1]*sc;
					if (tx < 0 || tx >= w || ty < 0 || ty >= h) continue;
					if (labels[ty*w + tx] >= 0)continue;
					p = ty*w + tx;
					goto DONE;
				}
				sc += 1;
			}
		}
		DONE:
		labels[p] = i;
		dis[p] = costMap[p];
	}
	// some components will have the same label, 
	// and some labels will be merged, so disappeared
	GDT(costMap.pData, dis.pData, labels.pData, w, h);
#endif

#if 0
	costMap.imshow("cost Map");
	dis.imagesc("dis");
	labels.imagesc("labels", 0);
#endif

	CHeap H(w*h, true);
	for (int i = 1; i < h; i++){
		for (int j = 1; j < w; j++){
			int idx = i*w + j;

			int l0 = labels[idx];
			int l1 = labels[idx - 1];
			int l2 = labels[idx - w];

			if (l0 != l1){
                                int64_t k = Key(l0, l1);
				double v = dis[idx] + dis[idx - 1];
				H.Push(&v, (double*)&k, 1);
			}

			if (l0 != l2){
                                int64_t k = Key(l0, l2);
				double v = dis[idx] + dis[idx - w];
				H.Push(&v, (double*)&k, 1);
			}
		}
	}

#if 0
	// clear Neighbor system
	for (int i = 0; i < m_nodeCnt; i++){
		m_neighborCnts[i] = 0;
	}
#endif

	// init the neighbor distances
	for (int i = 0; i < sCnt;i++){
		int nbCnt = NeighborCnt(i);
		float* nbDis = NeighborDis(i);
		for (int i = 0; i < nbCnt; i++){
			nbDis[i] = 1e7;
		}
	}

	while (H.Size()){
                int64_t newK, currK;

		float minDis = H.Top((double*)&currK);
		int l0 = KeyFirst(currK);
		int l1 = KeySecond(currK);

		while (H.Size()){
			float dis = H.Top((double*)&newK);
			if (newK != currK){
				break;
			}
			H.Pop();

			if (dis < minDis){
				minDis = dis;
			}
		}

		//printf("%d, %d, %f\n", l0, l1, minDis);
		// AddNeighbor(l0, l1, minDis);
		SetNeighborDis(l0, l1, minDis);
		currK = newK;
	}

	return sCnt;
}

void CSuperPixelGraph::ComputeNeighborDis()
{
	for (int i = 0; i < m_nodeCnt; i++){
		int srcIdx = i;
		int nbCnt = NeighborCnt(srcIdx);
		int* nbIdx = NeighborIndex(srcIdx);
		float* nbDis = NeighborDis(srcIdx);
		for (int k = 0; k < nbCnt; k++){
			nbDis[k] = GetDistance(srcIdx, nbIdx[k]);
		}
	}
}

int CSuperPixelGraph::ComputeNeighborDis(FImage& disMap)
{
	int w = disMap.width();
	int h = disMap.height();
	int sCnt = m_nodeCnt;

	float* dis = disMap.data();

	CHeap H(w*h, true);
	for (int i = 1; i < h; i++){
		for (int j = 1; j < w; j++){
			int idx = i*w + j;

			int l0 = m_labels[idx];
			int l1 = m_labels[idx - 1];
			int l2 = m_labels[idx - w];

			if (l0 != l1){
                                int64_t k = Key(l0, l1);
				double v = dis[idx] + dis[idx - 1];
				H.Push(&v, (double*)&k, 1);
			}

			if (l0 != l2){
                                int64_t k = Key(l0, l2);
				double v = dis[idx] + dis[idx - w];
				H.Push(&v, (double*)&k, 1);
			}
		}
	}

	while (H.Size()){
                int64_t newK, currK;

		float minDis = H.Top((double*)&currK);
		int l0 = KeyFirst(currK);
		int l1 = KeySecond(currK);

		while (H.Size()){
			float dis = H.Top((double*)&newK);
			if (newK != currK){
				break;
			}
			H.Pop();

			if (dis < minDis){
				minDis = dis;
			}
		}

		//printf("%d, %d, %f\n", l0, l1, minDis);
		// AddNeighbor(l0, l1, minDis);
		SetNeighborDis(l0, l1, minDis);
		currK = newK;
	}

	return sCnt;
}

int CSuperPixelGraph::FindNN(int nodeIdx, int* validFlag, int* nnIdx, float* nnWeight, int maxCnt)
{
	CHeap H(m_nodeCnt);
	int* nodeVec = new int[m_nodeCnt];
	int* nodeFlag = new int[m_nodeCnt];
	memset(nodeFlag, 0, sizeof(int)*m_nodeCnt);

	int validCnt = validFlag[nodeIdx] ? 1 : 0;

	// in queue
	int currIdx = 0;
	nodeVec[currIdx++] = nodeIdx;
	H.Push(nodeIdx, 1); // max weight 1
	nodeFlag[nodeIdx] = 1;

	for (int i = 0; i < currIdx; i++){
		int nbNum = NeighborCnt(nodeVec[i]);
		int* nbIdx = NeighborIndex(nodeVec[i]);
		for (int k = 0; k < nbNum; k++){
			int nb = nbIdx[k];
			if (nodeFlag[nb] == 0){
				H.Push(nb, GetWeight(nodeIdx, nb));
				nodeVec[currIdx] = nb;
				currIdx++;
				if (validFlag[nb]){
					validCnt++;
				}
				if (validCnt >= maxCnt){
					goto NEXT;
				}
				nodeFlag[nb] = 1;
			}
		}
	}

NEXT:
	// out queue
	int nnCnt = 0;
	for (int i = 0; i < m_nodeCnt; i++){
		double wt;
		int idx = H.Pop(&wt);
		if (validFlag[idx]){
			nnIdx[nnCnt] = idx;
			nnWeight[nnCnt] = wt;
			nnCnt++;
			if (nnCnt >= maxCnt){
				break;
			}
		}
	}

	delete[] nodeVec;
	delete[] nodeFlag;
	return nnCnt;
}

int CSuperPixelGraph::FindNN_GD(int nodeIdx, int* validFlag, int* nnIdx, float* nnWeight, int maxCnt)
{
	CHeap H(m_nodeCnt, true); // min-heap

	float* currDis = new float[m_nodeCnt];
	memset(currDis, 0x7F, sizeof(float)*m_nodeCnt); // max float

	float* currColorDis = new float[m_nodeCnt];
	float* currSpaceDis = new float[m_nodeCnt];
	memset(currColorDis, 0x7F, sizeof(float)*m_nodeCnt); // max float
	memset(currSpaceDis, 0x7F, sizeof(float)*m_nodeCnt); // max float

	int nnCnt = 0;

	H.Push(nodeIdx, 0); // min distance
	currDis[nodeIdx] = 0;
	currColorDis[nodeIdx] = 0;
	currSpaceDis[nodeIdx] = 0;

	while(H.Size()){
		double dis;
		int idx = H.Pop(&dis);

		if (dis > currDis[idx]){
			continue;
		}

		if (validFlag[idx]){
			nnIdx[nnCnt] = idx;
			nnWeight[nnCnt] = dis;
			nnCnt++;
			if (nnCnt >= maxCnt){
				break;
			}
			if (nnCnt >= 10 && dis > 1.0){
				//break;
			}
		}

		int nbNum = NeighborCnt(idx);
		int* nbIdx = NeighborIndex(idx);
		float* nbDis = NeighborDis(idx);
		for (int k = 0; k < nbNum; k++){
			int nb = nbIdx[k];

			float colorDis = currColorDis[idx] + GetColorDis(idx, nb);
			float spaceDis = currSpaceDis[idx] + GetSpaceDis(idx, nb);
			float sigma1 = 200, sigma2 = 50;
			float newDis = 1 - exp(-(spaceDis*spaceDis / (2 * sigma1*sigma1)) - (colorDis*colorDis / (2 * sigma2*sigma2)));

			if (newDis < currDis[nb]){
				H.Push(nb, newDis);
				currDis[nb] = newDis;
				currColorDis[nb] = colorDis;
				currSpaceDis[nb] = spaceDis;
			}
		}
	}

// 	float sigma = 6;
// 	for (int i = 0; i < nnCnt; i++){
// 		nnWeight[i] = exp(-nnWeight[i] / (2 * sigma*sigma));
// 	}

	delete[] currDis;
	delete[] currColorDis;
	delete[] currSpaceDis;
	return nnCnt;
}

void CSuperPixelGraph::AddSuperPixelMask(FImage& image)
{
	if (!image.matchDimension(m_width, m_height, image.nchannels())){
		return;
	}
	int w = m_width;
	int h = m_height;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			if (m_borderFlags[i*w + j] >= 0){
				float* ptr = image.pixPtr(i, j);
				for (int c = 0; c < image.nchannels(); c++){
					ptr[c] = 0;
				}
			}
		}
	}

}

void CSuperPixelGraph::ShowSuperPixel(FImage& image, char* winName)
{
	int ch = image.nchannels();
	if (!image.matchDimension(m_width, m_height, ch)){
		return;
	}

	ComputeAverage(image);

	FImage showimg(image);
	for (int nodeIdx = 0; nodeIdx < m_nodeCnt; nodeIdx++)
	{
		float r, g, b;
		b = m_averBGR[3 * nodeIdx + 0];
		g = m_averBGR[3 * nodeIdx + 1];
		r = m_averBGR[3 * nodeIdx + 2];

		int cnt = NodeSize(nodeIdx);
		int* pt = NodeItems(nodeIdx);
		for (int i = 0; i < cnt; i++){
			int idx = m_width*pt[2 * i + 1] + pt[2 * i];
			showimg[idx * ch + 2] = r;
			showimg[idx * ch + 1] = g;
			showimg[idx * ch + 0] = b;
		}
	}

	showimg.imshow(winName);
}

void CSuperPixelGraph::ShowNN(char* winName, FImage& image, int nodeIdx, int* validFlag)
{
#define COLOR_SHOW
	int ch = image.nchannels();
	if (!image.matchDimension(m_width, m_height, ch)){
		return;
	}

	int maxNN = 100;
	int* nnIdx = new int[maxNN];
	float* nnWeight = new float[maxNN];

	int ptCnt = FindNN_GD(nodeIdx, validFlag, nnIdx, nnWeight, maxNN);

	float maxWeight = -1;
	for (int i = 0; i < ptCnt; i++){
		printf("%f, ", nnWeight[i]);
		if (nnWeight[i] > maxWeight){
			maxWeight = nnWeight[i];
		}
	}

	CColorTable colorTbl;
	unsigned char* pColor = NULL;

#ifdef COLOR_SHOW
	FImage showimg(image);
#else
	FImage showimg(m_width, m_height, 1);
#endif
	for (int i = 0; i < ptCnt; i++)
	{
		int idx = nnIdx[i];
		int colorIdx = (nnWeight[i] / maxWeight) * 255;
		pColor = colorTbl[colorIdx];

		int cnt = NodeSize(idx);
		int* pt = NodeItems(idx);
		for (int k = 0; k < cnt; k++){
			int idx = m_width*pt[2 * k + 1] + pt[2 * k];
			#ifdef COLOR_SHOW
			showimg[idx * ch + 0] = (float)pColor[0] / 255.;
			showimg[idx * ch + 1] = (float)pColor[1] / 255.;
			showimg[idx * ch + 2] = (float)pColor[2] / 255.;
			#else
			showimg[idx] = colorIdx / 255.;
			#endif
		}
	}
	showimg.imshow(winName);

	delete[] nnIdx;
	delete[] nnWeight;
}

void GenerateGaussPDF(Vector<float>& points, int ptCnt, float* weights, FImage& outPDF)
{
	int w = 300, h = 300;
	int hw = w / 2, hh = h / 2; // half
	outPDF.allocate(w, h, 1);

	int m = 5;
	float sigma = 3;
	for (int k = 0; k < ptCnt; k++){
		float x = points[2 * k];
		float y = points[2 * k + 1];
		float cx = x + hw, cy = y + hh;
		for (int dx = -m * sigma; dx <= m * sigma; dx++){
			for (int dy = -m * sigma; dy <= m * sigma; dy++){
				int candix = cx + dx;
				int candiy = cy + dy;
				if (candix >= 0 && candix < w && candiy >= 0 && candiy < h){
					float dis = dx*dx + dy*dy;
					float gaussWt = exp(-dis / (2 * sigma*sigma));
					//outPDF[candiy*w + candix] += gaussWt*weights[k];
					outPDF[candiy*w + candix] += gaussWt;
				}
			}
		}
	}
}

void CSuperPixelGraph::ShowNNClustering(FImage& image, int nodeIdx, FImage& seedsFlow)
{
	int ch = image.nchannels();
	if (!image.matchDimension(m_width, m_height, ch)){
		return;
	}

	int maxNN = 300;
	int* nnIdx = new int[maxNN];
	float* nnWeight = new float[maxNN];

	int* validFlag = new int[m_nodeCnt];
	memset(validFlag, 0, m_nodeCnt*sizeof(int));
	for (int i = 0; i < m_nodeCnt; i++){
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		if (!OpticFlowIO::unknown_flow(u, v)){
			validFlag[i] = 1;
		}
	}

	int ptCnt = FindNN_GD(nodeIdx, validFlag, nnIdx, nnWeight, maxNN);

#if 0
	{
		FImage points(2, ptCnt);
		for (int i = 0; i < ptCnt; i++){
			int idx = nnIdx[i];
			points[2 * i] = seedsFlow[2 * idx];
			points[2 * i + 1] = seedsFlow[2 * idx + 1];
		}
		Vector<int> clusterLabel(ptCnt);
		WeightedClustering(points, nnWeight, ptCnt, clusterLabel.data());

		int sw = 300, sh = 300; // show width and height
		int hw = sw / 2, hh = sh / 2; // half width and height
		FImage showFlowImg(sw, sh, 3);
		float colors[3];
		for (int i = 0; i < ptCnt; i++){
			int idx = nnIdx[i];
			int cu = seedsFlow[2 * idx];
			int cv = seedsFlow[2 * idx + 1];

			float* pDst = showFlowImg.pData + ((cv + hh)*sw + cu + hw) * 3;

			// generate random color according to the clustering label
			srand(clusterLabel[i]);
			colors[0] = (rand() % 6) * 50 / 255.;
			colors[1] = (rand() % 6) * 50 / 255.;
			colors[2] = (rand() % 6) * 50 / 255.;

			memcpy(pDst, colors, 3 * sizeof(float));
		}
		showFlowImg.imshow("cluster 1");
	}
#endif

#if 0
	//
	FImage showFlowImg;
	Vector<float> flowPts(2 * ptCnt);
	for (int i = 0; i < ptCnt; i++){
		int idx = nnIdx[i];
		flowPts[2 * i] = seedsFlow[2 * idx];
		flowPts[2 * i + 1] = seedsFlow[2 * idx + 1];
	}
	GenerateGaussPDF(flowPts, ptCnt, nnWeight, showFlowImg);
	//showFlowImg.imagesc("pdf");

	int tw = showFlowImg.width();
	int th = showFlowImg.height();
	ptCnt = 0;
	FImage points(2, tw*th);
	float* density = new float[tw*th];
	for (int i = 0; i < th; i += 1){
		for (int j = 0; j < tw; j += 1){
			density[ptCnt] = showFlowImg.pData[i*tw + j];
			if (density[ptCnt] > 1e-3){
				points[2 * ptCnt] = j - tw / 2;
				points[2 * ptCnt + 1] = i - th / 2;
				ptCnt++;
			}
		}
	}
	Vector<int> clusterLabel(ptCnt);
	Clustering(points, density, ptCnt, clusterLabel.data());

	//
	int sw = 300, sh = 300; // show width and height
	int hw = sw / 2, hh = sh / 2; // half width and height
	int showCh = 3;
	showFlowImg.allocate(sw, sh, showCh);
	float colors[3];
	showFlowImg.setValue(0); // background
	for (int i = 0; i < ptCnt; i++){
		int cu = points[2 * i];
		int cv = points[2 * i + 1];
		float* pDst = showFlowImg.pData + ((cv + hh)*sw + cu + hw)*showCh;
		// generate random color according to the clustering label
		srand(clusterLabel[i]);
		colors[0] = (rand() % 6) * 50 / 255.;
		colors[1] = (rand() % 6) * 50 / 255.;
		colors[2] = (rand() % 6) * 50 / 255.;
		memcpy(pDst, colors, showCh * sizeof(float));
	}
	showFlowImg.imshow("cluster 2");
	delete[] density;
#endif

	delete[] nnIdx;
	delete[] nnWeight;
	delete[] validFlag;
}

void CSuperPixelGraph::ShowNNFlow(FImage& image, int nodeIdx, FImage& seedsFlow)
{
	int ch = image.nchannels();
	if (!image.matchDimension(m_width, m_height, ch)){
		return;
	}

	int maxNN = 300;
	int* nnIdx = new int[maxNN];
	float* nnWeight = new float[maxNN];

	int* validFlag = new int[m_nodeCnt];
	memset(validFlag, 0, m_nodeCnt*sizeof(int));
	for (int i = 0; i < m_nodeCnt; i++){
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		if (!OpticFlowIO::unknown_flow(u, v)){
			validFlag[i] = 1;
		}
	}

	int ptCnt = FindNN_GD(nodeIdx, validFlag, nnIdx, nnWeight, maxNN);

	int sw = 300, sh = 300; // show width and height
	int hw = sw / 2, hh = sh / 2; // half width and height
	FImage showFlowImg(sw, sh);

	for (int i = 0; i < sw; i++){
		showFlowImg.pData[hh*sw + i] = 0.5;
	}
	for (int i = 0; i < sh; i++){
		showFlowImg.pData[i*sw + hw] = 0.5;
	}
	for (int i = 0; i < ptCnt; i++)
	{
		int idx = nnIdx[i];
		int cu = seedsFlow[2 * idx];
		int cv = seedsFlow[2 * idx + 1];
		if (!OpticFlowIO::unknown_flow(cu, cv)){
			showFlowImg.pData[(cv + hh)*sw + cu + hw] = 1;
		}
	}
	showFlowImg.imshow("stat 1");

	{
		FImage showFlowImg;
		Vector<float> flowPts(2 * ptCnt);
		for (int i = 0; i < ptCnt; i++){
			int idx = nnIdx[i];
			flowPts[2 * i] = seedsFlow[2 * idx];
			flowPts[2 * i + 1] = seedsFlow[2 * idx + 1];
		}
		GenerateGaussPDF(flowPts, ptCnt, nnWeight, showFlowImg);
		showFlowImg.imagesc("stat 2");
	}

	delete[] nnIdx;
	delete[] nnWeight;
	delete[] validFlag;
}

int CSuperPixelGraph::ArrangeLabels(int*kLabMat, int width, int height)
{
	int w = width;
	int h = height;
	int* tmpLabels = new int[w*h];
	int labelIdx = 0; // re-arrange the klables (start from 0)
	int nbOffset[4][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 } };
	//int nbOffset[8][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 }, { -1, -1 }, { -1, 1 }, { 1, -1 }, { 1, 1 } };

	std::stack<int> sti, stj;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int seedLabel = kLabMat[i*w + j];
			if (seedLabel >= 0){
				sti.push(i);
				stj.push(j);
				kLabMat[i*w + j] = -1; // flag
				tmpLabels[i*w + j] = labelIdx;
				while (!sti.empty()){  // search connected region
					int ti = sti.top(); sti.pop();
					int tj = stj.top(); stj.pop();
					for (int k = 0; k < 4; k++){
						int tii = ti + nbOffset[k][0];
						int tjj = tj + nbOffset[k][1];
						if (tii >= 0 && tii < h && tjj >= 0 && tjj < w){
							if (kLabMat[tii*w + tjj] == seedLabel){
								sti.push(tii);
								stj.push(tjj);
								kLabMat[tii*w + tjj] = -1; // flag
								tmpLabels[tii*w + tjj] = labelIdx;
							}
						}
					}
				}
				labelIdx++;
			}
		}
	}
	memcpy(kLabMat, tmpLabels, w*h*sizeof(int));
	delete[] tmpLabels;

	return labelIdx;
}

void CSuperPixelGraph::GetAverageItem(FImage& featureImg, int nodeIdx, float* c1, float* c2, float* c3)
{
	int cnt = NodeSize(nodeIdx);
	int* pt = NodeItems(nodeIdx);
	float averC1, averC2, averC3;
	averC1 = averC2 = averC3 = 0;
	for (int i = 0; i < cnt; i++){
		int idx = m_width*pt[2 * i + 1] + pt[2 * i];
		averC1 += featureImg[idx * 3 + 0];
		averC2 += featureImg[idx * 3 + 1];
		averC3 += featureImg[idx * 3 + 2];
	}
	*c1 = averC1 / cnt;
	*c2 = averC2 / cnt;
	*c3 = averC3 / cnt;
}

float CSuperPixelGraph::GetWeight(int spIdx1, int spIdx2)
{
	float spaceDis = GetSpaceDis(spIdx1, spIdx2);
	float colorDis = GetColorDis(spIdx1, spIdx2);
	//
	float sigma1 = 30, sigma2 = 5;
	//float sigma1 = 30, sigma2 = 0.1;
	return exp(-(spaceDis*spaceDis / (2 * sigma1*sigma1)) - (colorDis*colorDis / (2 * sigma2*sigma2)));
}

float CSuperPixelGraph::GetDistance(int spIdx1, int spIdx2)
{
	return 1 - GetWeight(spIdx1, spIdx2);
}

float CSuperPixelGraph::GetColorDis(int spIdx1, int spIdx2)
{
	float x1 = m_centerPositions[2 * spIdx1];
	float y1 = m_centerPositions[2 * spIdx1 + 1];
	float x2 = m_centerPositions[2 * spIdx2];
	float y2 = m_centerPositions[2 * spIdx2 + 1];
	//
	float aver1[3], aver2[3];
#if 1
	aver1[0] = m_averLAB[3 * spIdx1 + 0];
	aver1[1] = m_averLAB[3 * spIdx1 + 1];
	aver1[2] = m_averLAB[3 * spIdx1 + 2];
	//
	aver2[0] = m_averLAB[3 * spIdx2 + 0];
	aver2[1] = m_averLAB[3 * spIdx2 + 1];
	aver2[2] = m_averLAB[3 * spIdx2 + 2];
#else
	aver1[0] = m_averBGR[3 * spIdx1 + 0];
	aver1[1] = m_averBGR[3 * spIdx1 + 1];
	aver1[2] = m_averBGR[3 * spIdx1 + 2];
	//
	aver2[0] = m_averBGR[3 * spIdx2 + 0];
	aver2[1] = m_averBGR[3 * spIdx2 + 1];
	aver2[2] = m_averBGR[3 * spIdx2 + 2];
#endif
	float colorDis = (aver1[0] - aver2[0])*(aver1[0] - aver2[0])
		+ (aver1[1] - aver2[1])*(aver1[1] - aver2[1])
		+ (aver1[2] - aver2[2])*(aver1[2] - aver2[2]);
	
	return sqrt(colorDis);
}

float CSuperPixelGraph::GetSpaceDis(int spIdx1, int spIdx2)
{
	float x1 = m_centerPositions[2 * spIdx1];
	float y1 = m_centerPositions[2 * spIdx1 + 1];
	float x2 = m_centerPositions[2 * spIdx2];
	float y2 = m_centerPositions[2 * spIdx2 + 1];
	float spaceDis = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
	return sqrt(spaceDis);
}

void CSuperPixelGraph::ClearAllNeighbors()
{
	for (int i = 0; i < m_nodeCnt; i++){
		m_neighborCnts[i] = 0;
	}
}

void ConstructNewGraph(FImage& costMap, IntImage& seeds, CSuperPixelGraph& outGraph)
{
	int w = costMap.width();
	int h = costMap.height();
	int sCnt = seeds.height();

	//
	FImage dis(w, h);
	IntImage labels(w, h);

#if 0
	float* sx = new float[sCnt];
	float* sy = new float[sCnt];
	for (int i = 0; i < sCnt; i++) {
		sx[i] = seeds[2 * i];
		sy[i] = seeds[2 * i + 1];
	}
	MinimumBarrier_DT(costMap, sx, sy, sCnt, labels.pData, dis.pData, 0, 6);// 0.0001);
	delete[] sx;
	delete[] sy;
#else
	memset(dis.pData, 0x7F, w*h*sizeof(float)); // init approx. to FLT_MAX
	memset(labels.pData, 0xFF, w*h*sizeof(int)); // init to -1;
	for (int i = 0; i < sCnt; i++) {
		int p = seeds[2 * i] + seeds[2 * i + 1] * w;
		labels[p] = i;
		dis[p] = costMap[p];
	}
	GDT(costMap.pData, dis.pData, labels.pData, w, h);
#endif

	sCnt = CSuperPixelGraph::ArrangeLabels(labels.pData, w, h);

	// check
	for (int i = 0; i < w*h; i++){
		if (labels.pData[i] < 0 || labels.pData[i] >= sCnt){
			printf("WDT error !\n");
		}
	}
#if 0
	costMap.imshow("cost Map");
	dis.imagesc("dis");
	labels.imagesc("labels", 0);
#endif

	outGraph.Init(labels.pData, w, h);
	outGraph.ClearAllNeighbors();

	CHeap H(w*h, true);
	for (int i = 1; i < h; i++){
		for (int j = 1; j < w; j++){
			int idx = i*w + j;

			int l0 = labels[idx];
			int l1 = labels[idx - 1];
			int l2 = labels[idx - w];

			if (l0 != l1){
                                int64_t k = Key(l0, l1);
				double v = dis[idx] + dis[idx - 1];
				H.Push(&v, (double*)&k, 1);
			}

			if (l0 != l2){
                                int64_t k = Key(l0, l2);
				double v = dis[idx] + dis[idx - w];
				H.Push(&v, (double*)&k, 1);
			}
		}
	}

	while (H.Size()){
                int64_t newK, currK;

		float minDis = H.Top((double*)&currK);
		int l0 = KeyFirst(currK);
		int l1 = KeySecond(currK);

		while (H.Size()){
			float dis = H.Top((double*)&newK);
			if (newK != currK){
				break;
			}
			H.Pop();

			if (dis < minDis){
				minDis = dis;
			}
		}

		//printf("%d, %d, %f\n", l0, l1, minDis);
		outGraph.AddNeighbor(l0, l1, minDis);
		currK = newK;
	}
}

#endif //_SUPER_PIXEL_H

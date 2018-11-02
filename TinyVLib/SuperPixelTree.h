#ifndef _SUPER_PIXEL_TREE_H
#define _SUPER_PIXEL_TREE_H

#include "Util.h"
#include "Image.h"
#include "SuperPixelGraph.h"
#include "Graph.h"

class CSuperPixelTree
{
public:
	CSuperPixelTree(){};
	~CSuperPixelTree(){};

	void Init(FImage& featureImg, FImage& img, int spSize);
	void Init(FImage& featureImg, int* kLabels, int width, int height);
	void Init(FImage& featureImg, char* fileName, int width, int height);

	inline int NodeCnt(){ return m_graph.NodeCnt(); };
	inline int NodeSize(int nodeIdx){ return m_graph.NodeSize(nodeIdx); };
	inline int* NodeItems(int nodeIdx){ return m_graph.NodeItems(nodeIdx); };
	inline int* NodePos(int nodeIdx){ return m_graph.NodePos(nodeIdx); };
	
	inline CSuperPixelGraph* GetGraph(){ return &m_graph; };
	inline CMst* GetTree(){ return &m_tree; };

	void ShowTree(const FImage& img, char* winName, bool norm = false);
	void ShowByOrder(FImage& img, char* winName);

private:
	void ConstructTree(FImage& featureImg);
	float GetSuperPixelDis(FImage& featureImg, int spIdx1, int spIdx2);
	CSuperPixelGraph m_graph;
	CMst m_tree;
	// save the distance to the closest seed (superpixel center) for each pixel
	FImage m_dist;
};

void CSuperPixelTree::Init(FImage& featureImg, FImage& img, int spSize)
{
	int w = img.width();
	int h = img.height();
	m_graph.Init(img, w*h / spSize, 20);
	ConstructTree(featureImg);
}

void CSuperPixelTree::Init(FImage& featureImg, int* kLabels, int width, int height)
{
	m_graph.Init(kLabels, width, height);
	ConstructTree(featureImg);
}

void CSuperPixelTree::Init(FImage& featureImg, char* fileName, int width, int height)
{
	m_graph.Init(fileName, width, height);
	ConstructTree(featureImg);
}


void CSuperPixelTree::ConstructTree(FImage& featureImg)
{
	// construct tree from graph

	// construct the superpixel edges from neighborhood
	int numV = m_graph.NodeCnt();
	int numE = 0;
	int* vFlag = new int[numV]; // is the vertex visited?
	memset(vFlag, 0, sizeof(int)*numV);

	int maxEdgeCnt = numV * 16;
	int* edges = new int[maxEdgeCnt * 2];
	float* edgeLens = new float[maxEdgeCnt];

	for (int i = 0; i < numV; i++){
		int num = m_graph.NeighborCnt(i);
		int* idx = m_graph.NeighborIndex(i);
		for (int j = 0; j < num; j++){
			int candiIndex = idx[j];
			if (vFlag[candiIndex]){
				continue;
			}
			edges[numE * 2] = i;
			edges[numE * 2 + 1] = candiIndex;
			edgeLens[numE] = GetSuperPixelDis(featureImg, i, candiIndex);
			numE++;
		}
		vFlag[i] = 1;
	}
	m_tree.Init(numV, numE, edges, edgeLens);

	delete[] vFlag;
	delete[] edges;
	delete[] edgeLens;
}

void CSuperPixelTree::ShowTree(const FImage& img, char* winName, bool norm /* = false*/)
{
	CColorTable colorTbl;

	int w = img.width();
	int h = img.height();

	FImage tmpImg(img);
	m_graph.AddSuperPixelMask(tmpImg);
	cv::Mat cvImg = ImageIO::CvmatFromPixels(tmpImg.pData, w, h, 3);

	int numV = NodeCnt();

	// get the range of edge length
	float minV = FLT_MAX;
	float maxV = FLT_MIN;
	for (int i = 0; i < numV; i++){
		float paLen = m_tree.ParentDistance(i);
		if (paLen < minV){
			minV = paLen;
		}
		if (paLen > maxV){
			maxV = paLen;
		}
	}
	
	unsigned char whiteColor[3] = { 255, 255, 255 };
	unsigned char* color = whiteColor;
	for (int i = 0; i < numV; i++){
		int parent = m_tree.Parent(i);
		int* pos = NodePos(i);
		int* parentPos = NodePos(parent);
		float paLen = m_tree.ParentDistance(i);
		if (norm){
			color = colorTbl[(int)((paLen - minV) / (maxV - minV) * 255)];
		}
		cv::line(cvImg, cvPoint(pos[0], pos[1]), cvPoint(parentPos[0], parentPos[1]), cvScalar(color[0], color[1], color[2]));
	}

	// draw red centers
	for (int i = 0; i < numV; i++){
		int* pos = NodePos(i);
		int x = pos[0];
		int y = pos[1];
		cvImg.data[y*cvImg.step + 3 * x + 0] = 0;
		cvImg.data[y*cvImg.step + 3 * x + 1] = 0;
		cvImg.data[y*cvImg.step + 3 * x + 2] = 255;
	}
	cv::imshow(winName, cvImg);
	cv::waitKey(0);
}

void CSuperPixelTree::ShowByOrder(FImage& img, char* winName)
{
	FImage tmpImg(img);
	tmpImg.setValue(1);
	
	int w = img.width();
	int h = img.height();

	int numV = NodeCnt();
	int* kLabels = m_graph.KLabels();
	for (int i = 0; i < numV; i++){
		int idx = m_tree.NodeOrder(i);

		int cnt = NodeSize(idx);
		int* pt = NodeItems(idx);
		for (int k = 0; k < cnt; k++){
			int x = pt[2 * k];
			int y = pt[2 * k + 1];
			float* p = img.pixPtr(y, x);
			tmpImg.setPixel(y, x, p);
		}

		tmpImg.imshow(winName, false);
	}
}

#include "Util.h"
float CSuperPixelTree::GetSuperPixelDis(FImage& costMap, int spIdx1, int spIdx2)
{
#if 0
	int* borders = m_graph.BorderFlags();
	int* kLabels = m_graph.KLabels();

	int* rect1 = m_graph.NodeRect(spIdx1);
	int* rect2 = m_graph.NodeRect(spIdx2);
	int left = min(rect1[0], rect2[0]);
	int right = max(rect1[1], rect2[1]);
	int top = min(rect1[2], rect2[2]);
	int bottom = max(rect1[3], rect2[3]);

	int w = costMap.width();
	int h = costMap.height();
	float maxCost = FLT_MIN;
//	int mx, my;

// 	FImage tmp(w,h,3);
// 	float color[3] = { 0, 0, 1 };
// 	tmp.setValue(1);
// 	m_graph.AddSuperPixelMask(tmp);

	for (int i = top; i <= bottom; i++){
		for (int j = left; j <= right; j++){
			int idx = i*w + j;
			if (kLabels[idx] == spIdx1 || kLabels[idx] == spIdx2){
				float cost = costMap[idx];
				if (cost > maxCost){
					maxCost = cost;
					//mx = j; my = i;
				}
			}
		}
	}
	
// 	printf("%d, %d, %f\n", spIdx1, spIdx2, maxCost);
// 	tmp.setPixel(my, mx, color);
// 	tmp.imshow("tmp", 0);

	return maxCost;
#elif 0
	int* borders = m_graph.BorderFlags();
	int* kLabels = m_graph.KLabels();

	int* rect1 = m_graph.NodeRect(spIdx1);
	int* rect2 = m_graph.NodeRect(spIdx2);
	int left = min(rect1[0], rect2[0]);
	int right = max(rect1[1], rect2[1]);
	int top = min(rect1[2], rect2[2]);
	int bottom = max(rect1[3], rect2[3]);
	
	int w = costMap.width();
	int h = costMap.height();
	float maxCost = FLT_MIN;
	
	// debug show
// 	FImage tmp(w,h,3);
// 	float color[3] = { 0, 0, 1 };
// 	tmp.setValue(1);
// 	m_graph.AddSuperPixelMask(tmp);

	for (int i = top; i <= bottom; i++){
		for (int j = left; j <= right; j++){
			int idx = i*w + j;
			if (borders[idx] >= 0){
				if ((kLabels[idx] == spIdx1 && borders[idx] == spIdx2)
					|| (kLabels[idx] == spIdx2 && borders[idx] == spIdx1)){
					//tmp.setPixel(i, j, color);
					float cost = costMap[idx];
					if (cost > maxCost){
						maxCost = cost;
					}
				}
			}
		}
	}
// 	printf("%d, %d\n", spIdx1, spIdx2);
// 	tmp.imshow("tmp", 0);

	return maxCost;
#elif 1
	float r1, g1, b1, r2, g2, b2;
	//m_graph.GetAverRGB(costMap, spIdx1, &r1, &g1, &b1);
	//m_graph.GetAverRGB(costMap, spIdx2, &r2, &g2, &b2);
	float dis = sqrt((r1 - r2)*(r1 - r2) + (g1 - g2)*(g1 - g2) + (b1 - b2)*(b1 - b2));
	return dis;
#elif 0
	int w = costMap.width();
	CLine line;
	int* pos1 = m_graph.NodePos(spIdx1);
	int* pos2 = m_graph.NodePos(spIdx2);
	int lineLen = line.Draw(pos1[0], pos1[1], pos2[0], pos2[1]);
	int* x = line.GetDrawMaskX();
	int* y = line.GetDrawMaskY();
	float dis = 0;
	for (int i = 0; i < lineLen; i++){
		dis += costMap[y[i] * w + x[i]];
	}
	return dis;
#elif 0
	int* pos1 = m_graph.NodePos(spIdx1);
	int* pos2 = m_graph.NodePos(spIdx2);
	float* pix1 = costMap.pixPtr(pos1[1], pos1[0]);
	float* pix2 = costMap.pixPtr(pos2[1], pos2[0]);
	float dis = sqrt(pow(pix1[0] - pix2[0], 2)
		+ pow(pix1[1] - pix2[1], 2)
		+ pow(pix1[2] - pix2[2], 2));
	return dis;
#elif 0
	int w = costMap.width();
	int h = costMap.height();

	int nbOffset[4][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 } };
	//int nbOffset[8][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 }, { -1, -1 }, { -1, 1 }, { 1, -1 }, { 1, 1 } };
	int* rect1 = m_graph.NodeRect(spIdx1);
	int* rect2 = m_graph.NodeRect(spIdx2);
	int left = min(rect1[0], rect2[0]);
	int right = max(rect1[1], rect2[1]);
	int top = min(rect1[2], rect2[2]);
	int bottom = max(rect1[3], rect2[3]);
	int* kLabels = m_graph.KLabels();

	// construct the edges in the two superpixels
	int cropW = right - left + 1;
	int cropH = bottom - top + 1;
	int numV = cropW*cropH;
	int numE = 0;
	int* vFlag = new int[numV]; // is the vertex visited?
	memset(vFlag, 0, sizeof(int)*numV);
	int maxEdgeCnt = numV * 4;
	int* edges = new int[maxEdgeCnt * 2];
	float* edgeLens = new float[maxEdgeCnt];

	for (int i = 0; i < cropH; i++){
		for (int j = 0; j < cropW; j++){
			int idx = i*cropW + j;
			for (int k = 0; k < 4; k++){
				int nbRow = i + nbOffset[k][0];
				int nbCol = j + nbOffset[k][1];
				if (nbRow >= 0 && nbRow < cropH && nbCol >= 0 && nbCol < cropW){
					int nbIdx = nbRow*cropW + nbCol;
					if (vFlag[nbIdx]){
						continue;
					}

					// get edgelen
					float* pix1 = costMap.pixPtr(i + top, j + left);
					float dis = *pix1;

					edges[numE * 2] = idx;
					edges[numE * 2 + 1] = nbIdx;
					edgeLens[numE] = dis;
					numE++;

					edges[numE * 2] = nbIdx;
					edges[numE * 2 + 1] = idx;
					edgeLens[numE] = dis;
					numE++;
				}
			}
			vFlag[idx] = 1;
		}
	}

	int* pos1 = m_graph.NodePos(spIdx1);
	int* pos2 = m_graph.NodePos(spIdx2);
	int srcIdx = (pos1[1] - top)*cropW + pos1[0] - left;
	int dstIdx = (pos2[1] - top)*cropW + pos2[0] - left;
	
	CDijkstra dij;
	dij.Init(numV, numE, edges, edgeLens);
	dij.PathFinding(srcIdx);
	float dis = dij.PathDistance(dstIdx);

	delete[] vFlag;
	delete[] edges;
	delete[] edgeLens;

	return dis;
#endif
}

#endif
// [9/19/2015 Yinlin.Hu]

#ifndef _GRAPH_H
#define _GRAPH_H

#include <queue>
#include "Heap.h"
#include "Util.h"

/************************************************************************/
/*                         Minimum Spanning Tree                        */
/************************************************************************/
class CMst
{
public:
	CMst();
	~CMst();
	
	template <class T1>
	int Init(Image<T1>& image);

	int Init(int numV, int numE, int* edges, float* edgeLens);
	void Clean();

	inline int Root() { return 0; }
	inline int NodeCnt(){ return m_numV; }
	inline int Parent(int vIdx){ return m_parent[vIdx]; }
	inline float ParentDistance(int vIdx){ return m_parentDis[vIdx]; }
	inline int* ChildIndex(int vIdx){ return m_children[vIdx]; }
	inline int ChildNum(int vIdx){ return m_numChild[vIdx]; }
	inline int NodeOrder(int order){ return m_nodeOrder[order]; }

	void DistanceTransform(FImage& img, IntImage& seeds, FImage& outDis, IntImage& outLabel);

	int TreeFilter(float* data, int* validFlag, int len, float sigma);

	int FindNN(int vIdx, int* nnIdx, float* nnDis, int maxCnt);
	void ShowNN(char* winName, FImage& image, int vIdx);
	void ShowNN_UI(char* winName, FImage& image);

	void Test();

	void ShowMST(FImage& img);

private:
	int FindSet(int x);
	void Kruskal();
	void SortIncrease(int* sortIdx, unsigned short* data, int len);
	void BuildTree();

	int m_numV, m_numE;
	int* m_parent;
	int* m_edge;
	float* m_edgeLen;
	int* m_edgeOrder;
	int* m_nodeOrder;
	int* m_numChild;
	int* m_numConnect;
	int** m_connected;
	float** m_connectLength;
	float* m_parentDis;	// distance between this node and its parent
	int** m_children;
	int m_maxNeighBor;
};

void CMst::Test()
{
	// the graph from Figure 23.4 of <Introduction to Algorithms> 3e
	int numV = 9;
	int numE = 14;
	int edges[28] = { 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 0, 7, 1, 7, 8, 8, 6, 8, 2, 3, 5, 2, 5 };
	float edgeLens[14] = { 0.04, 0.08, 0.07, 0.09, 0.10, 0.02, 0.01, 0.08, 0.11, 0.07, 0.06, 0.02, 0.14, 0.04 };
	Init(numV, numE, edges, edgeLens);

	for (int i = 0; i < numV; i++){
		printf("%d, %f, %d", Parent(i), ParentDistance(i), NodeOrder(i));
		printf("< ");
		for (int j = 0; j < ChildNum(i); j++){
			printf("%d ", ChildIndex(i)[j]);
		}
		printf(">\n");
	}
}

CMst::CMst()
{
	memset(this, 0, sizeof(*this));
	m_maxNeighBor = 32;	// max neighbor count of each node
}

CMst::~CMst()
{
	Clean();
}

int CMst::Init(int numV, int numE, int* edges, float* edgeLens)
{
	Clean();

	m_numV = numV;	// vertex number
	m_numE = numE;	// edge number

	m_edge = new int[m_numE * 2];
	m_edgeLen = new float[m_numE];
	m_edgeOrder = new int[m_numE];
	memcpy(m_edge, edges, m_numE * 2 * sizeof(int));
	memcpy(m_edgeLen, edgeLens, numE * sizeof(float));

	// normalize the edgeLens
	float maxInEdgeLen = FLT_MIN;
	for (int i = 0; i < m_numE; i++){
		if (edgeLens[i] < 0){
			edgeLens[i] = 0;
		}
		if (edgeLens[i] > maxInEdgeLen){
			maxInEdgeLen = edgeLens[i];
		}
	}
	// quick sort
	unsigned short* tmpLen = new unsigned short[m_numE];
	for (int i = 0; i < m_numE; i++){
		tmpLen[i] = m_edgeLen[i] / maxInEdgeLen * 65535;
	}
	SortIncrease(m_edgeOrder, tmpLen, m_numE);
	delete[] tmpLen;

	m_parent = new int[m_numV];
	// assign every node's parent to itself.
	for (int i = 0; i < m_numV; i++){
		m_parent[i] = i;
	}

	m_connected = (int**)malloc(m_numV*sizeof(int*));
	m_connectLength = (float**)malloc(m_numV*sizeof(float*));
	m_children = (int**)malloc(m_numV*sizeof(int*));

	m_connected[0] = (int*)malloc(m_numV*m_maxNeighBor*sizeof(int));
	m_connectLength[0] = (float*)malloc(m_numV*m_maxNeighBor*sizeof(float));
	m_children[0] = (int*)malloc(m_numV*m_maxNeighBor*sizeof(int));
	for (int i = 1; i < numV; i++){
		m_connected[i] = m_connected[i - 1] + m_maxNeighBor;
		m_connectLength[i] = m_connectLength[i - 1] + m_maxNeighBor;
		m_children[i] = m_children[i - 1] + m_maxNeighBor;
	}

	m_numChild = new int[m_numV];
	m_numConnect = new int[m_numV];
	// there's no connection.
	memset(m_numConnect, 0, m_numV*sizeof(int));

	m_nodeOrder = new int[m_numV];
	m_parentDis = new float[m_numV];

	Kruskal();
	BuildTree();

	// debug show
	// 	for(int i=0; i<m_numV; i++){
	// 		printf("%d\t%d\t%d\t%d\n", m_parent[i], m_numChild[i], m_weights[i], m_nodeOrder[i]);
	// 	}

	return 0;
}

template <class T1>
int CMst::Init(Image<T1>& image)
{
	int w = image.width();
	int h = image.height();
	int ch = image.nchannels();
	T1* imgData = image.pData;

	int numV = w*h;
	int numE = (w - 1)*h + w*(h - 1); // 4-neighbor

	int* edges = new int[numE * 2];
	float* edgeLens = new float[numE];

	// edge construction
	int edgeIdx = 0;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w - 1; j++){
			int nodeIdx = i*w + j;
			int e1 = nodeIdx;
			int e2 = nodeIdx + 1;
			edges[2 * edgeIdx] = e1;
			edges[2 * edgeIdx + 1] = e2;
			float dis = 0;
			for (int k = 0; k < ch; k++) {
				dis = __max(dis, abs(imgData[e1*ch + k] - imgData[e2*ch + k]));
			}
			edgeLens[edgeIdx] = dis + 0.001;
			edgeIdx++;
		}
	}
	for (int i = 0; i < h - 1; i++){
		for (int j = 0; j < w; j++){
			int nodeIdx = i*w + j;
			int e1 = nodeIdx;
			int e2 = nodeIdx + w;
			edges[2 * edgeIdx] = e1;
			edges[2 * edgeIdx + 1] = e2;
			float dis = 0;
			for (int k = 0; k < ch; k++) {
				dis = __max(dis, abs(imgData[e1*ch + k] - imgData[e2*ch + k]));
			}
			edgeLens[edgeIdx] = dis + 0.001;
			edgeIdx++;
		}
	}

	Init(numV, numE, edges, edgeLens);

	delete[] edges;
	delete[] edgeLens;

	return 0;
}

void CMst::Clean()
{
	if (m_parent){
		delete[] m_parent;
		m_parent = NULL;
	}
	if (m_edge){
		delete[] m_edge;
		m_edge = NULL;
	}
	if (m_edgeLen){
		delete[] m_edgeLen;
		m_edgeLen = NULL;
	}
	if (m_edgeOrder){
		delete[] m_edgeOrder;
		m_edgeOrder = NULL;
	}
	if (m_nodeOrder){
		delete[] m_nodeOrder;
		m_nodeOrder = NULL;
	}
	if (m_numChild){
		delete[] m_numChild;
		m_numChild = NULL;
	}
	if (m_numConnect){
		delete[] m_numConnect;
		m_numConnect = NULL;
	}
	if (m_connected){
		free(m_connected[0]);
		free(m_connected);
		m_connected = NULL;
	}
	if (m_connectLength){
		free(m_connectLength[0]);
		free(m_connectLength);
		m_connectLength = NULL;
	}
	if (m_parentDis){
		delete[] m_parentDis;
		m_parentDis = NULL;
	}
	if (m_children){
		free(m_children[0]);
		free(m_children);
		m_children = NULL;
	}
}

void CMst::BuildTree()
{
	for (int i = 0; i < m_numV; i++){
		m_numChild[i] = 0;
		m_parent[i] = -1;
	}
	m_parent[0] = 0;// the root of the tree
	m_nodeOrder[0] = 0;// count from root
	m_parentDis[0] = 0;// there's no edge between the root and it's parent

	std::queue<int> Q;
	int parent;
	int len = 1;
	Q.push(0);
	while (!Q.empty())
	{
		parent = Q.front(); Q.pop();
		int numC = m_numConnect[parent];
		for (int i = 0; i < numC; i++){
			int nb = m_connected[parent][i];	// get NeighBor
			if (-1 == m_parent[nb]){
				Q.push(nb);
				m_parent[nb] = parent;
				m_parentDis[nb] = m_connectLength[parent][i];
				m_nodeOrder[len++] = nb;

				m_children[parent][m_numChild[parent]++] = nb;
			}
		}
	}
}

void CMst::SortIncrease(int* sortIdx, unsigned short* data, int len)
{
	const int binNum = 65536;
	int tmp[binNum + 1];
	int* histo = tmp + 1;	//histo[-1] for special use
	int hittedCnt[binNum];
	memset(histo, 0, binNum*sizeof(int));
	memset(hittedCnt, 0, binNum*sizeof(int));

	// get histogram
	for (int i = 0; i < len; i++){
		histo[data[i]]++;
	}
	// get accumulated histogram
	int accu = histo[0];
	for (int i = 1; i < binNum; i++){
		histo[i] += accu;
		accu = histo[i];
	}
	// sort according to the accumulated histogram
	histo[-1] = 0;
	for (int i = 0; i < len; i++){
		int v = data[i];
		int rank = histo[v - 1] + hittedCnt[v];
		sortIdx[rank] = i;
		hittedCnt[v]++;
	}

	// debug show
	// 	for(int i=0; i<len; i++){
	// 		printf("%d ", sortIdx[i]);
	// 	}
	// 	printf("\n");
}

// disjoint set operation
// return the representative of the set containing x
int CMst::FindSet(int x)
{
	// find the object with self parent
	int p = x;
	while (p != m_parent[p]){
		p = m_parent[p];
	}

	m_parent[x] = p; // update for quick ref
	return p;
}

// the kruskal algorithm
void CMst::Kruskal()
{
	for (int i = 0; i < m_numE; i++){
		int id = m_edgeOrder[i];

		int* edge = m_edge + id * 2;
		int u = edge[0];
		int v = edge[1];
		int pu = FindSet(u);
		int pv = FindSet(v);
		if (pu != pv){ // Union
			int idx = m_numConnect[u];
			m_connected[u][idx] = v;
			m_connectLength[u][idx] = m_edgeLen[id];
			m_numConnect[u]++;

			idx = m_numConnect[v];
			m_connected[v][idx] = u;
			m_connectLength[v][idx] = m_edgeLen[id];
			m_numConnect[v]++;

			m_parent[pu] = m_parent[pv];

			// safe check
			if (m_numConnect[u] > m_maxNeighBor
				|| m_numConnect[v] > m_maxNeighBor){
				printf("MST Error: some nodes have too many neighbors !\n");
				break;
			}
		}
	}
}

void CMst::DistanceTransform(FImage& img, IntImage& seeds, FImage& outDis, IntImage& outLabel)
{
	int w = img.width();
	int h = img.height();
	int ptCnt = seeds.height();
	int numV = w*h;

	outDis.allocate(w, h, 1);
	memset(outDis.pData, 0x7F, w*h*sizeof(float)); // max float
	outLabel.allocate(w, h);
	memset(outLabel.pData, 0xFF, w*h*sizeof(int)); // -1

	for (int i = 0; i < ptCnt; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		outDis[y*w + x] = 0;
		outLabel[y*w + x] = i;
	}

	// update parents from leaf to root
	for (int pos = numV - 1; pos >= 0; pos--){ // from the last node
		int idx = NodeOrder(pos);
		int p = Parent(idx);
		float pd = ParentDistance(idx);
		float nd = outDis[idx] + pd;
		if (nd < outDis[p]){
			outDis[p] = nd;
			outLabel[p] = outLabel[idx];
		}
	}

	// update child from root to leaf				 
	for (int pos = 0; pos < numV; pos++){ // from the first node
		int idx = NodeOrder(pos);
		int p = Parent(idx);
		float pd = ParentDistance(idx);
		float nd = outDis[p] + pd;
		if (nd < outDis[idx]){
			outDis[idx] = nd;
			outLabel[idx] = outLabel[p];
		}
	}
}

int CMst::TreeFilter(float* data, int* validFlag, int len, float sigma)
{
	if (len != m_numV){
		return -1;
	}

	float* tmpCost = new float[m_numV];
	float* weightTbl = new float[m_numV];
	for (int i = 0; i < m_numV; i++){
		weightTbl[i] = exp(-ParentDistance(i) / sigma);
	}

	// first step: from leaf to root
	for (int pos = m_numV - 1; pos >= 0; pos--){ // from last node
		int idx = NodeOrder(pos);

		int childNum = ChildNum(idx);
		int* childIdx = ChildIndex(idx);
		tmpCost[idx] = data[idx];
		for (int j = 0; j < childNum; j++){
			int idChild = childIdx[j];
			float simi = weightTbl[idChild];
			tmpCost[idx] += (simi*tmpCost[idChild]);
		}
	}

	// second step: from root to leaf
	data[0] = tmpCost[0];// the first is root node				 
	for (int pos = 1; pos < m_numV; pos++){ // from second node
		int idx = NodeOrder(pos);

		int parent = Parent(idx);
		float simi = weightTbl[idx];

		data[idx] = tmpCost[idx] + simi*(data[parent] - simi*tmpCost[idx]);
	}

	delete[] weightTbl;
	delete[] tmpCost;

	return 0;
}

int CMst::FindNN(int vIdx, int* nnIdx, float* nnDis, int maxCnt)
{
	CHeap H(m_numV, true); // min-heap
	int nnCnt = 0;

	int* vFlag = new int[m_numV];
	memset(vFlag, 0, m_numV*sizeof(int));

	H.Push(vIdx, 0); // min distance
	vFlag[vIdx] = 1;

	while (H.Size()){
		double dis;
		int idx = H.Pop(&dis);

		nnIdx[nnCnt] = idx;
		nnDis[nnCnt] = dis;
		nnCnt++;
		if (nnCnt >= maxCnt){
			break;
		}

		int p = Parent(idx);
		if (!vFlag[p]){
			H.Push(p, dis + ParentDistance(idx));
			vFlag[p] = 1;
		}
		int cn = ChildNum(idx);
		int* cIdx = ChildIndex(idx);
		for (int k = 0; k < cn; k++){
			int cd = cIdx[k];
			if (!vFlag[cd]){
				H.Push(cd, dis + ParentDistance(cd));
				vFlag[cd] = 1;
			}
		}
	}

	delete[] vFlag;
	return nnCnt;
}

void CMst::ShowNN(char* winName, FImage& image, int vIdx)
{
#define COLOR_SHOW

	int w = image.width();
	int h = image.height();
	int ch = image.nchannels();

	int maxNN = 200*200;
	int* nnIdx = new int[maxNN];
	float* nnDis = new float[maxNN];

// 	CTimer t;
// 	for (int i = 0; i < m_numV; i++){
// 		if (i % 100 == 0)
// 		{
// 			FindNN(i, nnIdx, nnDis, maxNN);
// 		}
// 	}
// 	t.toc("FindNN: ");

	int ptCnt = FindNN(vIdx, nnIdx, nnDis, maxNN);

	float maxDis = -1;
	for (int i = 0; i < ptCnt; i++){
		if (nnDis[i] > maxDis){
			maxDis = nnDis[i];
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
		int colorIdx = (nnDis[i] / maxDis) * 255;
		pColor = colorTbl[colorIdx];

#ifdef COLOR_SHOW
		showimg[idx * ch + 0] = (float)pColor[0] / 255.;
		showimg[idx * ch + 1] = (float)pColor[1] / 255.;
		showimg[idx * ch + 2] = (float)pColor[2] / 255.;
#else
		showimg[idx] = colorIdx / 255.;
#endif
	}
	showimg.imshow(winName);

	delete[] nnIdx;
	delete[] nnDis;
}

FImage g_mst_image;
char g_mst_win_name[256];
static void MST_OnMouse(int event, int x, int y, int flags, void* param)
{
	CMst* mst = (CMst*)param;
	int w = g_mst_image.width();
	switch (event)
	{
	case CV_EVENT_LBUTTONDOWN:
		int vIdx = y*w + x;
		mst->ShowNN(g_mst_win_name, g_mst_image, vIdx);
		break;
	}
}

void CMst::ShowNN_UI(char* winName, FImage& image)
{
	strcpy(g_mst_win_name, winName);
	g_mst_image.copy(image);
	g_mst_image.imshow(winName);
	cvSetMouseCallback(winName, MST_OnMouse, this);
	cvWaitKey(0);
}

void CMst::ShowMST(FImage& img)
{
	int w = img.width();
	int h = img.height();
	float ratio = 2;

	FImage showImg;
	img.imresize(showImg, ratio, INTER_NN);

	cv::Mat imgMat = ImageIO::CvmatFromPixels(showImg.pData, showImg.width(), showImg.height(), showImg.nchannels());

	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int idx = i*w + j;
			int parent = Parent(idx);
			int px = parent % w;
			int py = parent / w;
//			cv::line(imgMat, CvPoint(ratio*j, ratio*i), CvPoint(ratio*px, ratio*py),
//				cv::Scalar(255, 255, 255), 1);
		}
	}

	cv::imshow("MST", imgMat);

}


/************************************************************************/
/*                              Graph                                   */
/************************************************************************/
class CGraph
{
public:
	CGraph();
	~CGraph();

	// node id is from 0 to numV-1
	int Init(int numV, int numE, int* edges, float* edgeLens);

	void Reset();

	inline int NodeCnt(){ return m_numV; };
	inline int NeighborCnt(int nodeIdx){ return m_neighborCnts[nodeIdx]; };
	inline int* NeighborIndex(int nodeIdx){ return m_neighborIDs[nodeIdx]; };
	inline float* NeighborDis(int nodeIdx){ return m_neighborDis[nodeIdx]; };
	inline int MaxNeighborCnt(){ return m_maxNeighBor; };

private:
	void _AddNeighbor(int spIdx1, int spIdx2, float dis);
	void _AddNeighbor0(int srcIdx, int dstIdx, float dis);

	int m_numV, m_numE;
	int** m_neighborIDs;
	float** m_neighborDis;
	int* m_neighborCnts;
	int m_maxNeighBor; // max neighbor number of each SuperPixel
};

CGraph::CGraph()
{
	memset(this, 0, sizeof(*this));
	m_maxNeighBor = 64;
}

CGraph::~CGraph()
{
	Reset();
}

int CGraph::Init(int numV, int numE, int* edges, float* edgeLens)
{
	Reset();

	m_numV = numV;
	m_numE = numE;

	// malloc space for the neighbor system
	m_neighborCnts = new int[m_numV];
	memset(m_neighborCnts, 0, m_numV*sizeof(int));

	m_neighborIDs = (int**)malloc(m_numV*sizeof(int*));
	m_neighborIDs[0] = (int*)malloc(m_numV*m_maxNeighBor*sizeof(int));
	for (int i = 1; i < m_numV; i++){
		m_neighborIDs[i] = m_neighborIDs[i - 1] + m_maxNeighBor;
	}

	m_neighborDis = (float**)malloc(m_numV*sizeof(float*));
	m_neighborDis[0] = (float*)malloc(m_numV*m_maxNeighBor*sizeof(float));
	for (int i = 1; i < m_numV; i++){
		m_neighborDis[i] = m_neighborDis[i - 1] + m_maxNeighBor;
	}

	// construct the neighbor system from the edges
	for (int i = 0; i < m_numE; i++){
		int e0 = edges[2 * i];
		int e1 = edges[2 * i + 1];
		float dis = edgeLens[i];
		_AddNeighbor(e0, e1, dis);
	}

	return m_numV;
}

void CGraph::Reset()
{
	if (m_neighborIDs){
		free(m_neighborIDs[0]);
		free(m_neighborIDs);
		m_neighborIDs = NULL;
	}
	if (m_neighborDis){
		free(m_neighborDis[0]);
		free(m_neighborDis);
		m_neighborDis = NULL;
	}
	if (m_neighborCnts){
		delete[] m_neighborCnts;
		m_neighborCnts = NULL;
	}
}

void CGraph::_AddNeighbor(int spIdx1, int spIdx2, float dis)
{
	_AddNeighbor0(spIdx1, spIdx2, dis);
	_AddNeighbor0(spIdx2, spIdx1, dis);
}

void CGraph::_AddNeighbor0(int srcIdx, int dstIdx, float dis)
{
	int nbCnt = m_neighborCnts[srcIdx];
	int* nbIDs = m_neighborIDs[srcIdx];
	float* nbDis = m_neighborDis[srcIdx];

	if (nbCnt >= m_maxNeighBor){
		printf("Too many neighbors for node %d in CGraph!\n", srcIdx);
		return;
	}

	// check if the neighbor exists
	for (int i = 0; i < nbCnt; i++){
		if (nbIDs[i] == dstIdx){
			printf("Duplicate edges <%d - %d> in CGraph!\n", srcIdx, dstIdx);
			return;
		}
	}

	nbIDs[nbCnt] = dstIdx;
	nbDis[nbCnt] = dis;
	m_neighborCnts[srcIdx]++;
}

/************************************************************************/
/*  Single Source Shortest Path: Dijkstra                               */
/************************************************************************/
class CDijkstra
{
public:
	CDijkstra();
	~CDijkstra();

	void Init(int numV, int numE, int* edges, float* edgeLens);
	void Clean();

	// find the path to dstIdx from srcIdx
	// float PathFinding(int srcIdx, int dstIdx);

	// find the paths to all nodes from scrIdx
	void PathFinding(int srcIdx);
	inline float PathDistance(int dstIdx){ return m_dist[dstIdx]; }
	inline int PathPrevNode(int idx){ return m_prev[idx]; }

	void Test();
private:
	int m_numV;
	int* m_inS;
	int* m_prev;
	float* m_dist;
	int* m_numConnect;
	int** m_connected;
	float** m_connectLength;
	int m_maxNeighBor;
};

CDijkstra::CDijkstra()
{
	memset(this, 0, sizeof(*this));
	m_maxNeighBor = 32; // max neighbor count of each node
}

CDijkstra::~CDijkstra()
{
	Clean();
}

void CDijkstra::Init(int numV, int numE, int* edges, float* edgeLens)
{
	Clean();

	m_numV = numV;

	m_inS = new int[m_numV];
	m_prev = new int[m_numV];
	m_dist = new float[m_numV];
	m_numConnect = new int[m_numV];

	m_connected = (int**)malloc(m_numV*sizeof(int*));
	m_connectLength = (float**)malloc(m_numV*sizeof(float*));
	m_connected[0] = (int*)malloc(m_numV*m_maxNeighBor*sizeof(int));
	m_connectLength[0] = (float*)malloc(m_numV*m_maxNeighBor*sizeof(float));
	for (int i = 1; i < m_numV; i++){
		m_connected[i] = m_connected[i - 1] + m_maxNeighBor;
		m_connectLength[i] = m_connectLength[i - 1] + m_maxNeighBor;
	}

	// construct the inner neighbor system from edge informations
	memset(m_numConnect, 0, sizeof(int)*m_numV);
	for (int i = 0; i < numE; i++){
		int e0 = edges[2 * i + 0];
		int e1 = edges[2 * i + 1];
		int num0 = m_numConnect[e0];
		if (num0 == m_maxNeighBor){
			fprintf(stderr, "Too many edges for dijkstra.");
			exit(1);
		}
		m_connected[e0][num0] = e1;
		m_connectLength[e0][num0] = edgeLens[i];
		m_numConnect[e0]++;
	}

	//debug show
#if 0
	for (int i = 0; i < m_numV; i++){
		int nbNum = m_numConnect[i];
		int* nbIdx = m_connected[i];
		printf("%d > ", i);
		for (int j = 0; j < nbNum; j++){
			printf("%d(%.1f) ", nbIdx[j], m_connectLength[i][j]);
		}
		printf("\n");
	}
#endif
}

void CDijkstra::Clean()
{
	if (m_inS){
		delete[] m_inS;
		m_inS = NULL;
	}
	if (m_prev){
		delete[] m_prev;
		m_prev = NULL;
	}
	if (m_dist){
		delete[] m_dist;
		m_dist = NULL;
	}
	if (m_numConnect){
		delete[] m_numConnect;
		m_numConnect = NULL;
	}
	if (m_connected){
		free(m_connected[0]);
		free(m_connected);
		m_connected = NULL;
	}
	if(m_connectLength){
		free(m_connectLength[0]);
		free(m_connectLength);
		m_connectLength = NULL;
	}
}

void CDijkstra::PathFinding(int srcIdx)
{
	// init: there's no node in S
	for (int i = 0; i < m_numV; i++){
		m_inS[i] = 0;
		m_prev[i] = -1;
		m_dist[i] = FLT_MAX;
	}

	// add srcIdx into S and update data
	m_inS[srcIdx] = 1;
	m_dist[srcIdx] = 0;
	int nbNum = m_numConnect[srcIdx];
	int* nbIdx = m_connected[srcIdx];
	for (int i = 0; i < nbNum; i++){
		m_prev[nbIdx[i]] = srcIdx;
		m_dist[nbIdx[i]] = m_connectLength[srcIdx][i];
	}

	// add other node (m_numV-1 nodes) into S
	for (int n = 0; n < m_numV - 1; n++){
		// find the node with minimum dist not in S
		float minDist = FLT_MAX;
		int minIdx = -1;
		for (int i = 0; i < m_numV; i++){
			if (!m_inS[i] && m_dist[i] < minDist){
				minDist = m_dist[i];
				minIdx = i;
			}
		}
		// add this node to S and update its neighbors
		m_inS[minIdx] = 1;
		nbNum = m_numConnect[minIdx];
		nbIdx = m_connected[minIdx];
		for (int i = 0; i < nbNum; i++){
			if (!m_inS[nbIdx[i]]){
				float newDis = m_dist[minIdx] 
					+ m_connectLength[minIdx][i];
				if (newDis < m_dist[nbIdx[i]]){
					m_dist[nbIdx[i]] = newDis;
					m_prev[nbIdx[i]] = minIdx;
				}
			}
		}
	}
}

void CDijkstra::Test()
{
	// the graph from Figure 24.6 of <Introduction to Algorithms> 3e
	int numV = 5;
	int numE = 10;
	int edges[20] = { 0, 1, 2, 0, 0, 3, 3, 0, 3, 1, 1, 4, 4, 1, 2, 3, 4, 2, 3, 4 };
	float edgeLens[10] = { 1, 10, 2, 3, 9, 4, 6, 5, 7, 2 };
	Init(numV, numE, edges, edgeLens);

	int startIdx = 2;
	PathFinding(startIdx);
	for (int i = 0; i < numV; i++){
		int node = i;
		while (node != startIdx){
			printf("%d ", node);
			node = m_prev[node];
		}
		printf("%d\n", startIdx);
	}
}

#endif

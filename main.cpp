
#include <queue>

#include "SuperPixelTree.h"
#include "Util.h"
#include "SLIC.h"
#include "TreeFilter.h"
#include "OpticFlowIO.h"
#include "Census.h"
#include "ImageFeature.h"
#include "ImagePyramid.h"
#include "Fitting.h"
//#include "matcher.h"

#define MAX_DISPLACEMENT 400

//#define DEBUG_SHOW

int STEP = 3;
int MAXR = 4;
int LEVELS = 5;
int ITER_N = 6;
float CHECK_TH = 2 * STEP;
int CHECK_LEVELS = 2;
int BORDER_WIDTH = 0;

int FB_CHECK = 1;

int* gKLabels = NULL;
int* gKLabels2 = NULL;
int gw = 0;
int gh = 0;

#if 0

#include "epic_aux.h"
/* Compute the neighboring matrix between seeds as well as the closest seed for each pixel */
void distance_transform_and_graph(const int_image* seeds, const float_image* cost, dt_params_t* dt_params,
	int_image* labels, float_image* dmap, csr_matrix* ngh, int n_thread);
int find_nn_graph_arr(csr_matrix* graph, int seed, int nmax, int* best, float* dist);

#endif 

template<class T>
void ToImage(Image<T>& outImg, T* pData, int width, int height)
{
	if (!outImg.matchDimension(width, height, 1)){
		outImg.allocate(width, height);
	}
	memcpy(outImg.pData, pData, width*height*sizeof(T));
}

#if 0
// the graph should be sysmetric
void Graph2Tree(csr_matrix* graph, CMst* tree)
{
	assert(graph->nr == graph->nc);
	int numV = graph->nr;
	int numE = 0;
	int maxEdgeCnt = numV * 32;
	int* edges = new int[maxEdgeCnt*2];
	float* edgeLens = new float[maxEdgeCnt];

	int* indptr = graph->indptr;
	for (int i = 0; i < numV; i++){
		for (int nb = indptr[i]; nb < indptr[i + 1]; nb++){
			edges[2 * numE] = i;
			edges[2 * numE + 1] = graph->indices[nb];
			edgeLens[numE] = graph->data[nb];
			numE++;
		}
	}

	tree->Init(numV, numE, edges, edgeLens);
	delete[] edges;
	delete[] edgeLens;
}

void Segmentation(IntImage& seeds, FImage& costMap, IntImage& kLabels, IntImage& nnf, FImage& dis)
{
	int w = costMap.width();
	int h = costMap.height();
	int numV = seeds.height();

	// reset the tree using geodestic distance
	int_image interSeeds; empty_image(&interSeeds, int, 2, numV);
	for (int i = 0; i < numV; i++){
		//if (ptx[i] == 718)
		//	printf("%d: %d, %d\n", i, ptx[i], pty[i]);
		interSeeds.pixels[2 * i] = seeds[2 * i];
		interSeeds.pixels[2 * i + 1] = seeds[2 * i + 1];
	}
	float_image fcost; empty_image(&fcost, float, w, h);
	memcpy(fcost.pixels, costMap.pData, w*h*sizeof(float));
	float_image dmap = { 0 };
	int_image labels; empty_image(&labels, int, w, h);
	// compute distance transform and build graph
	csr_matrix ngh = { 0 };
	distance_transform_and_graph(&interSeeds, &fcost, NULL, &labels, &dmap, &ngh, 1);

	if (!kLabels.matchDimension(w, h, 1)){
		kLabels.allocate(w, h);
	}
	memcpy(kLabels.pData, labels.pixels, w*h*sizeof(int));

	int nbCnt = 20;
	if (!nnf.matchDimension(nbCnt, numV, 1)){
		nnf.allocate(nbCnt, numV, 1);
	}
	if (!dis.matchDimension(nbCnt, numV, 1)){
		dis.allocate(nbCnt, numV, 1);
	}
	// compute nearest neighbors using the graph
	for (int i = 0; i < numV; i++)
		find_nn_graph_arr(&ngh, i, nbCnt, nnf.pData + i*nbCnt, dis.pData + i*nbCnt);

	free(interSeeds.pixels);
	free(fcost.pixels);
	free(labels.pixels);
	free(dmap.pixels);
	free(ngh.data);
	free(ngh.indices);
	free(ngh.indptr);
}

void ConstructTree(int* ptx, int* pty, int len, FImage& costMap, CMst* tree, IntImage& kLabels, IntImage& nnf, FImage& dis)
{
	int w = costMap.width();
	int h = costMap.height();
	int numV = len;

	// reset the tree using geodestic distance
	int_image seeds; empty_image(&seeds, int, 2, numV);
	for (int i = 0; i < numV; i++){
		//if (ptx[i] == 718)
		//	printf("%d: %d, %d\n", i, ptx[i], pty[i]);
		seeds.pixels[2 * i] = ptx[i];
		seeds.pixels[2 * i + 1] = pty[i];
	}
	float_image fcost; empty_image(&fcost, float, w, h);
	memcpy(fcost.pixels, costMap.pData, w*h*sizeof(float));
	float_image dmap = { 0 };
	int_image labels; empty_image(&labels, int, w, h);
	// compute distance transform and build graph
	csr_matrix ngh = { 0 };
	distance_transform_and_graph(&seeds, &fcost, NULL, &labels, &dmap, &ngh, 1);

	int* indptr = ngh.indptr;
	int maxNeighborCnt = -1;
	int noNbCnt = 0;
	for (int i = 0; i < numV; i++){
		//printf("%d neighbors: ", i);
		int nbCnt = indptr[i + 1] - indptr[i];
		if (nbCnt > maxNeighborCnt){
			maxNeighborCnt = nbCnt;
		}
		if (nbCnt == 0){
			//printf(": %d ", i);
			noNbCnt++;
		}
		for (int nb = indptr[i]; nb < indptr[i + 1]; nb++){
			//printf("%d(%.3f) ", ngh.indices[nb], ngh.data[nb]);
		}
		//printf("\n");
	}
	if (noNbCnt > 0){
		printf("Tree Error !!! \n");
	}
	printf("%d, %d\n", maxNeighborCnt, noNbCnt);

#if 0
	FImage hcost, hdmap;
	IntImage hlabels;
	ToImage(hcost, fcost.pixels, fcost.tx, fcost.ty);
	ToImage(hlabels, labels.pixels, labels.tx, labels.ty);
	ToImage(hdmap, dmap.pixels, dmap.tx, dmap.ty);
	hcost.imshow("cost");
	hdmap.imshow("dmap");
	hlabels.imshow("labels", 0);
#endif
	Graph2Tree(&ngh, tree);

	if (!kLabels.matchDimension(w, h, 1)){
		kLabels.allocate(w, h);
	}
	int nbCnt = 500;
	if (!nnf.matchDimension(nbCnt, numV, 1)){
		nnf.allocate(nbCnt, numV, 1);
	}
	if (!dis.matchDimension(nbCnt, numV, 1)){
		dis.allocate(nbCnt, numV, 1);
	}
	CTimer t;
	// compute nearest neighbors using the graph
	for (int i = 0; i < numV; i++)
		find_nn_graph_arr(&ngh, i, nbCnt, nnf.pData + i*nbCnt, dis.pData+i*nbCnt);
	t.toc("nnf: ");

	memcpy(kLabels.pData, labels.pixels, w*h*sizeof(int));

	free(seeds.pixels);
	free(fcost.pixels);
	free(labels.pixels);
	free(dmap.pixels);
	free(ngh.data);
	free(ngh.indices);
	free(ngh.indptr);
}

int VisoMatch(FImage& img1, FImage& img2, FImage& outMat)
{
	Matcher::parameters param;
	param.nms_n = 10;   // non - max - suppression: min.distance between maxima(in pixels)
	param.nms_tau = 50;  // non - max - suppression: interest point peakiness threshold
	param.match_binsize = 50;  // matching bin width / height(affects efficiency only)
	param.match_radius = 200; // matching radius(du / dv in pixels)
	param.match_disp_tolerance = 1;   // du tolerance for stereo matches(in pixels)
	param.outlier_disp_tolerance = 5;   // outlier removal : disparity tolerance(in pixels)
	param.outlier_flow_tolerance = 5;   // outlier removal : flow tolerance(in pixels)
	param.multi_stage = 1; // 0 = disabled, 1 = multistage matching(denser and faster)
	param.half_resolution = 0; // 0 = disabled, 1 = match at half resolution, refine at full resolution
	param.refinement = 1;   // refinement(0 = none, 1 = pixel, 2 = subpixel)
	Matcher *M = new Matcher(param);

	UCImage im1, im2;
	im1.copy(img1);
	im2.copy(img2);
	im1.desaturate();
	im2.desaturate();
	int dim[3];
	dim[0] = im1.width();
	dim[1] = im1.height();
	dim[2] = im1.width();

	M->pushBack(im1.pData, dim, 0);
	M->pushBack(im2.pData, dim, 0);
	M->matchFeatures(0);

	std::vector<Matcher::p_match> matches = M->getMatches();
	delete M;

	int nmatch = matches.size();
	outMat.allocate(4, nmatch);
	for (int i = 0; i < nmatch; i++){
		outMat.pData[4 * i + 0] = matches[i].u1p;
		outMat.pData[4 * i + 1] = matches[i].v1p;
		outMat.pData[4 * i + 2] = matches[i].u1c;
		outMat.pData[4 * i + 3] = matches[i].v1c;
	}
	return nmatch;
}
#endif

/* read matches, stored as x1 y1 x2 y2 per line (other values on the same is not taking into account */
void ReadMatches(const char *filename, FImage& outMat)
{
	float* tmp = new float[4 * 100000]; // max number of match pair
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

void WriteMatches(const char *filename, FImage& inMat)
{
	int len = inMat.height();
	FILE *fid = fopen(filename, "w");
	for (int i = 0; i < len; i++){
		float x1 = inMat[4 * i + 0];
		float y1 = inMat[4 * i + 1];
		float x2 = inMat[4 * i + 2];
		float y2 = inMat[4 * i + 3];
		fprintf(fid, "%.0f %.0f %.0f %.0f\n", x1, y1, x2, y2);
		//fprintf(fid, "%.3f %.3f %.3f %.3f 1 100\n", x1, y1, x2, y2);
		//fprintf(fid, "%d %d %d %d\n", (int)x1, (int)y1, (int)x2, (int)y2);
	}
	fclose(fid);
}

void WriteMatches(const char *filename, IntImage& seeds, float* u, float* v)
{
	int len = seeds.height();
	FILE *fid = fopen(filename, "w");
	for (int i = 0; i < len; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		fprintf(fid, "%d %d %.3f %.3f 1 100\n", x, y, x+u[i], y+v[i]);
	}
	fclose(fid);
}

void DrawMatches(FImage& outImg, FImage& inMat)
{
	int w = outImg.width();
	int h = outImg.height();
	int ch = outImg.nchannels();
	cv::Mat tmpImg = ImageIO::CvmatFromPixels(outImg.pData, w, h, ch);
	int nbMatch = inMat.rows();
	for (int i = 0; i < nbMatch; i++){
		float* mat = inMat.rowPtr(i);
		cv::line(tmpImg, cvPoint(mat[0], mat[1]), cvPoint(mat[2], mat[3]), cvScalar(0, 0, 255));
		cv::line(tmpImg, cvPoint(mat[0], mat[1]), cvPoint(mat[0], mat[1]), cvScalar(255, 255, 255));
	}
	cv::imshow("matches", tmpImg);
}


int borderWidth = 0;
float MatchCost(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f, int x1, int y1, int x2, int y2)
{
#if 1
	int w = im1f->width();
	int h = im1f->height();
	int ch = im1f->nchannels();
	float totalDiff;

// 	int b = 2;
// 	if (x1 < b || x1 >= w - b || y1 < b || y1 >= h - b)
// 		return FLT_MAX;
// 	if (x2 < b || x2 >= w - b || y2 < b || y2 >= h - b)
// 		return FLT_MAX;

	// so fast
	x1 = ImageProcessing::EnforceRange(x1, w);
	x2 = ImageProcessing::EnforceRange(x2, w);
	y1 = ImageProcessing::EnforceRange(y1, h);
	y2 = ImageProcessing::EnforceRange(y2, h);

	unsigned char* p1 = im1f->pixPtr(y1, x1);
	unsigned char* p2 = im2f->pixPtr(y2, x2);

	// Hamming
// 	totalDiff = HammingDistance(p1, p2, ch);
// 	return totalDiff;

#ifdef WITH_SSE
	totalDiff = 0;
	__m128i *r1 = (__m128i*)p1, *r2 = (__m128i*)p2, r3;
	assert(ch % 16 == 0);
	int iter = ch / 16;
	for (int i = 0; i < iter; i++){
		r3 = _mm_sad_epu8(*r1, *r2);
		totalDiff += (r3.m128i_u16[0] + r3.m128i_u16[4]);
		r1++; r2++;
	}
	return totalDiff;
#else
	totalDiff = 0;
	for (int i = 0; i < ch; i++){
		totalDiff += abs(p1[i] - p2[i]);
	}
	return totalDiff;
#endif

#elif 0
	int w = im1f->width();
	int h = im1f->height();
	int ch = im1f->nchannels();

	int r = 4;  // patch radius
	float totalDiff = 0;
	for (int i = -r; i < r; i+=1){
		int r1 = ImageProcessing::EnforceRange(y1 + i, h);
		int r2 = ImageProcessing::EnforceRange(y2 + i, h);
		for (int j = -r; j < r; j+=1){
			int c1 = ImageProcessing::EnforceRange(x1 + j, w);
			int c2 = ImageProcessing::EnforceRange(x2 + j, w);
			unsigned char* sig1 = im1f->pixPtr(r1, c1);
			unsigned char* sig2 = im2f->pixPtr(r2, c2);
			totalDiff += HammingDistance(sig1, sig2, ch);
		}
	}
	return totalDiff;

#elif 0 // RGB 8x8
	int w = img1.width();
	int h = img1.height();
	int ch = img1.nchannels();

	int r = 4;  // patch radius
	float totalDiff = 0;
	for (int i = -r; i < r; i+=1){
		int r1 = ImageProcessing::EnforceRange(y1 + i, h);
		int r2 = ImageProcessing::EnforceRange(y2 + i, h);
		for (int j = -r; j < r; j+=1){
			int c1 = ImageProcessing::EnforceRange(x1 + j, w);
			int c2 = ImageProcessing::EnforceRange(x2 + j, w);
			float* p1 = img1.pixPtr(r1, c1);
			float* p2 = img2.pixPtr(r2, c2);
			float diff = 0;
			for (int k = 0; k < ch; k++){
				diff += pow(p1[k] - p2[k], 2);
			}
			totalDiff += sqrt(diff);
		}
	}
	return totalDiff;
#elif 0
	int w = img1.width();
	int h = img1.height();
	int ch = img1.nchannels();

	float minDiff = FLT_MAX;

	int r = 8; // patch radius
	float scaleFactor[5] = { 1.0, 1.25, 1.5, 1.75, 2 };
	float interV[3];
	for (int s = 0; s < 5; s++){
		float totalDiff = 0;
		for (int i = -r; i <= r; i += 1){
			int r1 = ImageProcessing::EnforceRange(y1 + i, h);
			float r2 = ImageProcessing::EnforceRange(y2 + i*scaleFactor[s], h);
			for (int j = -r; j <= r; j += 1){
				int c1 = ImageProcessing::EnforceRange(x1 + j, w);
				float c2 = ImageProcessing::EnforceRange(x2 + j*scaleFactor[s], w);
				ImageProcessing::BilinearInterpolate(img2.pData, w, h, ch, c2, r2, interV);
				float* p1 = img1.pixPtr(r1, c1);
				float diff = 0;
				for (int k = 0; k < ch; k++){
					diff += pow(p1[k] - interV[k], 2);
				}
				totalDiff += sqrt(diff);
			}
		}
		if (totalDiff < minDiff){
			minDiff = totalDiff;
		}
	}
	return minDiff;
#endif
}

void match2flow(int* kLabels, int w, int h, FImage& seedsFlow, FImage& u, FImage& v)
{
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int id = kLabels[i*w + j];
			u[i*w + j] = seedsFlow[2 * id];
			v[i*w + j] = seedsFlow[2 * id + 1];
		}
	}
}

void match2flow(IntImage& seeds, FImage& seedsFlow, int w, int h, FImage& u, FImage& v)
{
	u.setValue(UNKNOWN_FLOW);
	v.setValue(UNKNOWN_FLOW);
	int seedCnt = seeds.height();
	int r = STEP / 2;
	for (int k = 0; k < seedCnt; k++){
		int x = seeds[2 * k];
		int y = seeds[2 * k + 1];
		for (int i = x - r; i <= x + r; i++){
			for (int j = y - r; j <= y + r; j++){
				if (i >= 0 && i < w && j >= 0 && j < h){
					u[j*w + i] = seedsFlow[2 * k];
					v[j*w + i] = seedsFlow[2 * k + 1];
				}
			}
		}
	}
}

void flow2match(FImage& u, FImage& v, FImage& outMat)
{
	int w = u.width();
	int h = u.height();
	int matcnt = 0;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int idx = i*w + j;
			if (!OpticFlowIO::unknown_flow(u[idx], v[idx])){
				matcnt++;
			}
		}
	}
	outMat.allocate(4, matcnt, 1);
	int cntIdx = 0;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int idx = i*w + j;
			if (!OpticFlowIO::unknown_flow(u[idx], v[idx])){
				outMat[4 * cntIdx + 0] = j;
				outMat[4 * cntIdx + 1] = i;
				outMat[4 * cntIdx + 2] = j + u[idx];
				outMat[4 * cntIdx + 3] = i + v[idx];
				cntIdx++;
			}
		}
	}
}

void ShowSparseFlow(char* winName, const FImage& img, IntImage& seeds, FImage& seedsFlow)
{
	int w = img.width();
	int h = img.height();
	int ch = img.nchannels();
	int len = seeds.height();
	cv::Mat tmpImg = ImageIO::CvmatFromPixels(img.pData, w, h, ch);
	for (int i = 0; i < len; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		if (!OpticFlowIO::unknown_flow(u,v)){
			cv::line(tmpImg, cvPoint(x, y), cvPoint(x + u, y + v), cvScalar(0, 0, 255));
		}
	}
	// draw white seeds
	for (int i = 0; i < len; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		//cv::line(tmpImg, cvPoint(x, y), cvPoint(x, y), cvScalar(255, 255, 255));
	}
	cv::imshow(winName, tmpImg);
}

void ShowSeedsFlow(char* winName, int w, int h, IntImage& seeds, FImage& seedsFlow)
{
	FImage u(w, h), v(w, h);
	match2flow(seeds, seedsFlow, w, h, u, v);
	OpticFlowIO::ShowFlow(winName, u.pData, v.pData, w, h, 82.81);
}

void SetSuperPixel(CSuperPixelGraph& graph, int spIdx, float* v, FImage& outImg)
{
	int cnt = graph.NodeSize(spIdx);
	int* pt = graph.NodeItems(spIdx);
	for (int i = 0; i < cnt; i++){
		outImg.setPixel(pt[2 * i + 1], pt[2 * i], v);
	}
}

void FilteringOnSegmentation(CMst* tree, int* ptx, int* pty, int len,
	FImage& img1, FImage& img2, UCImage* im1f, UCImage*im2f,
	int* isSubRoots, float* bestU, float* bestV)
{
	CTimer t;

	int numV = tree->NodeCnt();

	// construct sub tree
	int subEdgeNum;
	int* subEdges = new int[numV * 2];
	float* subEdgeLens = new float[numV * 4];
	int* sub2RawTable = new int[numV];
	int* raw2SubTable = new int[numV];

	int candiW = 2 * MAX_DISPLACEMENT + 1;
	int candiH = 2 * MAX_DISPLACEMENT + 1;
	IntImage candiFlowHisto(candiW, candiH);
	int oy = candiH / 2;
	int ox = candiW / 2;

	for (int i = 0; i < numV; i++){
		if (!isSubRoots[i])
			continue;

		int subIdx;
		int rootId = i;
		int parent;
		subIdx = 0;
		subEdgeNum = 0;
		
		sub2RawTable[subIdx] = rootId;
		raw2SubTable[rootId] = subIdx;
		subIdx++;

		candiFlowHisto.setValue(0);

		std::queue<int> Q;
		Q.push(rootId);
		while (!Q.empty()){
			parent = Q.front(); Q.pop();

			// collect candidate flow
			if (!OpticFlowIO::unknown_flow(bestU[parent], bestV[parent])){
				candiFlowHisto.pData[(oy + (int)bestV[parent])*candiW + ox + (int)bestU[parent]] ++;
			}

			int childNum = tree->ChildNum(parent);
			int* childIdx = tree->ChildIndex(parent);
			for (int c = 0; c < childNum; c++){
				int child = childIdx[c];
				if (!isSubRoots[child]){
					sub2RawTable[subIdx] = child;
					raw2SubTable[child] = subIdx;
					subIdx++;

					subEdges[2 * subEdgeNum] = raw2SubTable[parent];
					subEdges[2 * subEdgeNum + 1] = raw2SubTable[child];
					subEdgeLens[subEdgeNum] = tree->ParentDistance(child);
					subEdgeNum++;

					Q.push(child);
				}
			}
		}
		int subNumV = subIdx;
		CMst subTree;
		subTree.Init(subNumV, subEdgeNum, subEdges, subEdgeLens);
		CTreeFilter subTreeFilter;
		subTreeFilter.Init(&subTree);
		subTreeFilter.SetSigma(0.1);
		
		//candiFlowHisto.GaussianSmoothing(0.5, 2);
		//candiFlowHisto.imagesc("candi flow", 0);

		int candiCnt = 0;
		for (int ii = 0; ii < candiH; ii++){
			for (int jj = 0; jj < candiW; jj++){
				if (candiFlowHisto[ii*candiW + jj] > 1){
					candiCnt++;
				}
			}
		}
		printf("candidate count: %d\n", candiCnt);

		FImage candiFlow(2, candiCnt);
		int tmpIdx = 0;
		for (int ii = 0; ii < candiH; ii++){
			for (int jj = 0; jj < candiW; jj++){
				if (candiFlowHisto[ii*candiW + jj] > 1){
					candiFlow[2 * tmpIdx + 0] = jj - ox;
					candiFlow[2 * tmpIdx + 1] = ii - oy;
					tmpIdx++;
				}
			}
		}

		FImage costV(subNumV, candiCnt);
		for (int i = 0; i < candiCnt; i++){
			float u = candiFlow[i * 2 + 0];
			float v = candiFlow[i * 2 + 1];
			for (int j = 0; j < subNumV; j++){
				float cost;
#if 1
				if (OpticFlowIO::unknown_flow(bestU[sub2RawTable[j]], bestV[sub2RawTable[j]])){
					cost = 0;
				}
				else
#endif
				{
					int rawIdx = sub2RawTable[j];
					cost = MatchCost(img1, img2, im1f, im2f, 
						ptx[rawIdx], pty[rawIdx], ptx[rawIdx] + u, pty[rawIdx] + v);
				}
				costV[i*subNumV + j] = cost;
			}
		}

		for (int i = 0; i < candiCnt; i++){
			subTreeFilter.Filter(costV.rowPtr(i), subNumV);
		}

		for (int i = 0; i < subNumV; i++){
			float minCost = FLT_MAX;
			int minCandi = -1;
			for (int j = 0; j < candiCnt; j++){
				float cost = costV[j*subNumV + i];
				if (cost < minCost){
					minCost = cost;
					minCandi = j;
				}
			}
			bestU[sub2RawTable[i]] = candiFlow[minCandi * 2];
			bestV[sub2RawTable[i]] = candiFlow[minCandi * 2 + 1];
		}	

		//ShowSuperPixelFlow("tf", *graph, img1, bestU, bestV, numV);
	}
	t.toc("seg filtering: ");

	//ShowSuperPixelFlow("tf", *graph, img1, bestU, bestV, numV);

	delete[] subEdges;
	delete[] subEdgeLens;
	delete[] sub2RawTable;
	delete[] raw2SubTable;
}

void TreeSegmentation(CMst* tree, int* isSubRoots, float th, int minTreeNodeNum)
{
	th = 1.0;
	minTreeNodeNum = 50;

	CTimer t;

	// merge from bottom to up
	int numV = tree->NodeCnt();

	memset(isSubRoots, 0, sizeof(int)*numV);

	isSubRoots[0] = 1;
	for (int i = 0; i < numV; i++){
		if (tree->ParentDistance(i) >= th){
			isSubRoots[i] = 1;
		}
	}

	// is the tree too small?
	for (int i = numV - 1; i >= 0; i--){ // bottom to up
		if (!isSubRoots[i]){
			continue;
		}
		int rootId = i;
		int parent;
		int treeNodeNum = 0;
		std::queue<int> Q;
		Q.push(rootId);
		while (!Q.empty()){
			parent = Q.front(); Q.pop();
			treeNodeNum++;
			if (treeNodeNum >= minTreeNodeNum){
				break;
			}
			int childNum = tree->ChildNum(parent);
			int* childIdx = tree->ChildIndex(parent);
			for (int c = 0; c < childNum; c++){
				int child = childIdx[c];
				if (!isSubRoots[child]){
					Q.push(child);
				}
			}
		}
		if (treeNodeNum < minTreeNodeNum){
			isSubRoots[i] = 0;
		}
	}

	t.toc("tree seg: ");
}

#include <unordered_map>

// a flow is represented as a int
static inline int key(short u, short v) { return u + (v << 16); }
static inline short first(int i) { return short(i); }
static inline short second(int i) { return short(i >> 16); }
void NNF_Filtering(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f,
	IntImage& seeds, FImage& seedsFlow, IntImage& nnf, FImage& dis)
{
	CTimer t;
	double sigma = 0.05;

	int numV = seeds.height();
	int maxNbCnt = nnf.width();

	int* candiCnt = new int[numV];
	FImage candiU(maxNbCnt, numV);
	FImage candiV(maxNbCnt, numV);
	FImage costVol(maxNbCnt, numV);

	int* unstable = new int[numV];
	int cnt = 0;
	memset(unstable, 0, sizeof(int)*numV);
	for (int i = 0; i < numV; i++){
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		if (OpticFlowIO::unknown_flow(u, v)){
			unstable[i] = 1;
			cnt++;
		}
	}
	printf("unstable count: %d\n", cnt);

	typedef std::unordered_map<int, float> umap;
	umap* costMap = new umap[numV];
	int totalCandiCnt = 0;
	for (int i = 0; i < numV; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		int idx = 0;
		for (int j = 0; j < maxNbCnt; j++){
			int nbIdx = nnf[i*maxNbCnt + j];
			float d = dis[i*maxNbCnt + j];
			if (nbIdx < 0){// || d > 1.0){
				break;
			}
			short u = seedsFlow[nbIdx * 2 + 0];
			short v = seedsFlow[nbIdx * 2 + 1];
			int k = key(u, v);
			umap::const_iterator iter = costMap[i].find(k);
			if (iter == costMap[i].end()){ // not found
				float cost = 0;
				if (!unstable[i]){
					cost = MatchCost(img1, img2, im1f, im2f, x, y, x + u, y + v);
				}
				costMap[i].insert(umap::value_type(k, cost));

				candiU[i*maxNbCnt + idx] = u;
				candiV[i*maxNbCnt + idx] = v;
				idx++;
			}
		}
		candiCnt[i] = idx;
		totalCandiCnt += idx;
		//printf("%d ", idx);
	}
	//printf("%d\n", totalCandiCnt);
	for (int i = 0; i < numV; i++){
		int cnt = candiCnt[i];
		for (int j = 0; j < cnt; j++){
			short u = candiU[i*maxNbCnt + j];
			short v = candiV[i*maxNbCnt + j];
			int k = key(u, v);
			float wCost = 0;
			int validnbCnt = 0;
			for (int n = 0; n < maxNbCnt; n++){
				int nbIdx = nnf[i*maxNbCnt + n];
				float d = dis[i*maxNbCnt + n];
				if (nbIdx < 0)
					break;
				if (unstable[nbIdx])
					continue;
				int x = seeds[2 * nbIdx];
				int y = seeds[2 * nbIdx + 1];
				float cost;
				umap::const_iterator iter = costMap[nbIdx].find(k);
				if (iter == costMap[nbIdx].end()){
					cost = MatchCost(img1, img2, im1f, im2f, x, y, x + u, y + v);
					costMap[nbIdx].insert(umap::value_type(k, cost));
				}else{
					cost = iter->second;
				}
				wCost += (exp(-d / sigma) * cost);
				validnbCnt++;
				if (validnbCnt > 50){
					break;
				}
			}
			costVol[i*maxNbCnt + j] = wCost;
		}
	
	}

	// get the minimum cost
	for (int i = 0; i < numV; i++){
		float minCost = FLT_MAX;
		short bestu = 0;
		short bestv = 0;
		int cnt = candiCnt[i];
		for (int j = 0; j < cnt; j++){
			float cost = costVol[i*maxNbCnt + j];
			if (cost < minCost){
				minCost = cost;
				bestu = candiU[i*maxNbCnt + j];
				bestv = candiV[i*maxNbCnt + j];
			}
		}
		seedsFlow[2 * i + 0] = bestu;
		seedsFlow[2 * i + 1] = bestv;
	}

	delete[] candiCnt;
	delete[] unstable;
	delete[] costMap;

	t.toc("nnf filtering: ");
}

// smooth seeds flow
void NonLocalFiltering(CMst* tree, IntImage& seeds, FImage& seedsFlow)
{
	CTreeFilter tf;
	tf.Init(tree);
	tf.SetSigma(0.1);

	int numV = seeds.height();
	float* outU = new float[numV];
	float* outV = new float[numV];

	// prepare data
	for (int i = 0; i < numV; i++){
		outU[i] = seedsFlow[2 * i];
		outV[i] = seedsFlow[2 * i + 1];
		if (OpticFlowIO::unknown_flow(outU[i], outV[i])){
			outU[i] = 0;
			outV[i] = 0;
		}
	}

	tf.Filter(outU, numV);
	tf.Filter(outV, numV);

	// write back
	for (int i = 0; i < numV; i++){
		seedsFlow[2 * i] = outU[i];
		seedsFlow[2 * i + 1] = outV[i];
	}

	delete[] outU;
	delete[] outV;
}

void TreeFilter(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f,
	CMst* tree, IntImage& seeds, FImage& seedsFlow, FImage& candiFlow)
{
	CTreeFilter tf;
	tf.Init(tree);
	tf.SetSigma(0.1);

	int numV = seeds.height();
	int candiCnt = candiFlow.height();

	FImage costV(numV, candiCnt);
	for (int i = 0; i < candiCnt; i++){
		float u = candiFlow[i * 2 + 0];
		float v = candiFlow[i * 2 + 1];
		for (int j = 0; j < numV; j++){
			int x = seeds[j * 2];
			int y = seeds[j * 2 + 1];
			float cost ;
#if 1
			if (OpticFlowIO::unknown_flow(seedsFlow[2 * j], seedsFlow[2 * j + 1])){
				cost = 0;
			}else
#endif
			{
				cost = MatchCost(img1, img2, im1f, im2f, x, y, x + u, y + v);
			}
			costV[i*numV + j] = cost;
		}
	}

	for (int i = 0; i < candiCnt; i++){
		tf.Filter(costV.rowPtr(i), numV);
	}

	for (int i = 0; i < numV; i++){
		float minCost = FLT_MAX;
		int minCandi = -1;
		for (int j = 0; j < candiCnt; j++){
			float cost = costV[j*numV + i];
			if (cost < minCost){
				minCost = cost;
				minCandi = j;
			}
		}
		seedsFlow[2 * i] = candiFlow[minCandi * 2];
		seedsFlow[2 * i + 1] = candiFlow[minCandi * 2 + 1];
	}
}

// #define DEPRIVE_DUP_SEEDS
// a good initialization is already stored in bestU & bestV
void Propogate(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* pyd1f, UCImage* pyd2f, int level,
	int maxR, int iterCnt, IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow, float* bestCosts, 
	IntImage& seedGuessStore, bool guessStore = false)
{
	int nLevels = pyd1.nlevels();
	float ratio = pyd1.ratio();

	FImage im1 = pyd1[level];
	FImage im2 = pyd2[level];
	UCImage* im1f = pyd1f + level;
	UCImage* im2f = pyd2f + level;
	IntImage* seeds = pydSeeds + level;
	FImage* seedsFlow = pydSeedsFlow + level;

	int w = im1.width();
	int h = im1.height();
	int ptNum = seeds->height();

#ifdef DEPRIVE_DUP_SEEDS
	int* dupSeeds = new int[w*h]; //handle duplicated seeds
#endif

	int* seedGuessNumber = new int[ptNum];
	memset(seedGuessNumber, 0, ptNum*sizeof(int));

	int maxNb = neighbors.width();
	int* vFlags = new int[ptNum];

	// init cost
	for (int i = 0; i < ptNum; i++){
		int x = seeds->pData[2 * i];
		int y = seeds->pData[2 * i + 1];
		float u = seedsFlow->pData[2 * i];
		float v = seedsFlow->pData[2 * i + 1];
		bestCosts[i] = MatchCost(im1, im2, im1f, im2f, x, y, x + u, y + v);
	}

	for (int iter = 0; iter < iterCnt; iter++)
	{
		//ShowSuperPixelFlow(spt, img1, bestU, bestV, ptNum);
		//ShowSparseFlow("propogate", im1, seeds, bestU, bestV);
		//cvWaitKey();

		memset(vFlags, 0, sizeof(int)*ptNum);

#ifdef DEPRIVE_DUP_SEEDS
		memset(dupSeeds, 0xFF, sizeof(int)*w*h);
#endif

		int startPos = 0, endPos = ptNum, step = 1;
		if (iter % 2 == 1){
			startPos = ptNum - 1; endPos = -1; step = -1;
		}
		for (int pos = startPos; pos != endPos; pos += step){
			int idx = pos;

			int x = seeds->pData[2 * idx];
			int y = seeds->pData[2 * idx + 1];

#ifdef DEPRIVE_DUP_SEEDS
			if (dupSeeds[y*w + x] >= 0){
				int dupIdx = dupSeeds[y*w + x];
				seedsFlow->pData[2 * idx] = seedsFlow->pData[2 * dupIdx];
				seedsFlow->pData[2 * idx + 1] = seedsFlow->pData[2 * dupIdx + 1];
				bestCosts[idx] = bestCosts[dupIdx];
				vFlags[idx] = 1;
				continue;
			}
#endif
			int* nbIdx = neighbors.rowPtr(idx);

			// Propagation: Improve current guess by trying instead correspondences from neighbors
			for (int i = 0; i < maxNb; i++){
				if (nbIdx[i] < 0){
					break;
				}
				if (!vFlags[nbIdx[i]]){ // unvisited yet
					continue;
				}
				float tu = seedsFlow->pData[2 * nbIdx[i]];
				float tv = seedsFlow->pData[2 * nbIdx[i] + 1];
				float tc = MatchCost(im1, im2, im1f, im2f, x, y, x + tu, y + tv);
				if (tc < bestCosts[idx]){
					bestCosts[idx] = tc;
					seedsFlow->pData[2 * idx] = tu;
					seedsFlow->pData[2 * idx + 1] = tv;
				}

				// save guess
				if (guessStore){
					int* seedGuess = seedGuessStore.rowPtr(idx);
					int guessIdx = seedGuessNumber[idx];
					seedGuess[3 * guessIdx] = tu;
					seedGuess[3 * guessIdx + 1] = tv;
					seedGuess[3 * guessIdx + 2] = tc;
					seedGuessNumber[idx]++;
				}
			}

			// Random search: Improve current guess by searching in boxes
			// of exponentially decreasing size around the current best guess.
			for (int mag = maxR; mag >= 1; mag /= 2) {
				/* Sampling window */
				float tu = seedsFlow->pData[2 * idx] + rand() % (2 * mag + 1) - mag;
				float tv = seedsFlow->pData[2 * idx + 1] + rand() % (2 * mag + 1) - mag;
				float tc = MatchCost(im1, im2, im1f, im2f, x, y, x + tu, y + tv);
				if (tc < bestCosts[idx]){
					bestCosts[idx] = tc;
					seedsFlow->pData[2 * idx] = tu;
					seedsFlow->pData[2 * idx + 1] = tv;
				}

				// save guess
				if (guessStore){
					int* seedGuess = seedGuessStore.rowPtr(idx);
					int guessIdx = seedGuessNumber[idx];
					seedGuess[3 * guessIdx] = tu;
					seedGuess[3 * guessIdx + 1] = tv;
					seedGuess[3 * guessIdx + 2] = tc;
					seedGuessNumber[idx]++;
				}
			}
			vFlags[idx] = 1;
#ifdef DEPRIVE_DUP_SEEDS
			dupSeeds[y*w + x] = idx;
#endif
			//ShowSuperPixelFlow(spt, img1, bestU, bestV, ptNum);
		}
		//printf("iter %d: %f [s]\n", iter, t.toc());
	}

#ifdef DEPRIVE_DUP_SEEDS
	delete[] dupSeeds;
#endif

	// remove outliers
#if 0
	int b = 0;
	for (int i = 0; i < ptNum; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		int x2 = x + u;
		int y2 = y + v;
		if (x < b || x >= w - b || y < b || y >= h - b
			|| x2 < b || x2 >= w - b || y2 < b || y2 >= h - b
			|| sqrt(u*u + v*v)>MAX_DISPLACEMENT){
			seedsFlow[2 * i] = UNKNOWN_FLOW;
			seedsFlow[2 * i + 1] = UNKNOWN_FLOW;
		}
	}
#endif

	delete[] vFlags;
	delete[] seedGuessNumber;
}
inline bool IsFlowValid(float u, float v)
{
	bool unValid = OpticFlowIO::unknown_flow(u, v) 
		|| abs(u) > MAX_DISPLACEMENT || abs(v) > MAX_DISPLACEMENT;
	return !unValid;
}
#if 0
void SeedMatch(FImage& im1, FImage& im2, FImage& costMap1, FImage& costMap2,
	UCImage* im1f, UCImage* im2f, int* ptx, int* pty, int& len,
	float* bestU, float* bestV)
{
	int w = im1.width();
	int h = im1.height();

	// check the invalid points
	int* validFlag = new int[len];
	memset(validFlag, 0, sizeof(int)*len);
	int validNum = 0;
	for (int i = 0; i < len; i++){
		int npx = ptx[i] + bestU[i];
		int npy = pty[i] + bestV[i];
		if (ptx[i] < 0 || ptx[i] >= w || pty[i] < 0 || pty[i] >= h
			|| npx<0 || npx >= w || npy < 0 || npy >= h
			|| abs(bestU[i])>MAX_DISPLACEMENT || abs(bestV[i])>MAX_DISPLACEMENT){
			continue;
		}
		validFlag[i] = 1;
		validNum++;
	}
	IntImage seeds(2, validNum);
	FImage seedsFlow(2, validNum);
	int* validPtx = new int[validNum];
	int* validPty = new int[validNum];
	float* bu = new float[validNum];
	float* bv = new float[validNum];
	int idx = 0;
	for (int i = 0; i < len; i++){
		if (validFlag[i]){
			validPtx[idx] = ptx[i];
			validPty[idx] = pty[i];
			seeds[2 * idx] = ptx[i];
			seeds[2 * idx + 1] = pty[i];
			bu[idx] = bestU[i];
			bv[idx] = bestV[i];
			seedsFlow[2 * idx] = bestU[i];
			seedsFlow[2 * idx + 1] = bestV[i];
			idx++;
		}
	}
	assert(idx == validNum);

	//
	CMst tree;
	
	IntImage kLabels, nnf;
	FImage dis;
	ConstructTree(validPtx, validPty, validNum, costMap1, &tree, kLabels, nnf, dis);

	IntImage neighbors(32, validNum);
	neighbors.setValue(-1);
	for (int i = 0; i < validNum; i++){
		int* nbIdx = nnf.rowPtr(i);
		for (int j = 0; j < neighbors.width(); j++){
			if (nbIdx[j] < 0 || j>=32)
				break;
			neighbors[i*neighbors.width() + j] = nbIdx[j];
		}
	}

	int* isSubRoots = new int[validNum];

#ifdef DEBUG_SHOW
	ShowSuperPixelFlow("reverse", kLabels.pData, w, h, bu, bv);
#endif

	Propogate(im1, im2, im1f, im2f, 10, seeds, neighbors, bu, bv);

//	NNF_Filtering(im1, im2, im1f, im2f, seeds, seedsFlow, nnf, dis);
// 	for (int i = 0; i < validNum; i++){
// 		bu[i] = seedsFlow[2 * i];
// 		bv[i] = seedsFlow[2 * i + 1];
// 	}

	//TreeSegmentation(&tree, isSubRoots, 0.04, 50);
	
	//FilteringOnSegmentation(&tree, validPtx, validPty, validNum, im1, im2, im1f, im2f, isSubRoots, bu, bv);

#ifdef DEBUG_SHOW
	ShowSuperPixelFlow("reverse filtering", kLabels.pData, w, h, bu, bv);
#endif

	// copy back
	idx = 0;
	for (int i = 0; i < len; i++){
		if (validFlag[i]){
			bestU[i] = bu[idx];
			bestV[i] = bv[idx];
			idx++;
		}else{
			bestU[i] = UNKNOWN_FLOW;
			bestV[i] = UNKNOWN_FLOW;
		}
	}

	delete[] isSubRoots;
	delete[] validPtx;
	delete[] validPty;
	delete[] bu;
	delete[] bv;
	delete[] validFlag;
}
#endif

void PyramidRandomSearch(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f,
	IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow)
{
	int nLevels = pyd1.nlevels();
	float ratio = pyd1.ratio();

	FImage rawImg1 = pyd1[0];
	FImage rawImg2 = pyd2[0];
	srand(0);

	int w = rawImg1.width();
	int h = rawImg1.height();
	int numV = pydSeeds[0].height();

	float* bestCosts = new float[numV];

	// random Initialization on coarsest level
	// int initR = MAX_DISPLACEMENT * pow(ratio, nLevels - 1);
	int initR = max(pyd1[nLevels - 1].width(), pyd1[nLevels - 1].height());
	for (int i = 0; i < numV; i++){
		pydSeedsFlow[nLevels - 1][2 * i] = rand() % (2 * initR + 1) - initR;
		pydSeedsFlow[nLevels - 1][2 * i + 1] = rand() % (2 * initR + 1) - initR;
	}

	// added libviso2 Initialization
#if 0
	t.tic();
	FImage mat;
	VisoMatch(rawImg1, rawImg2, mat);
	t.toc("viso: ");
	DrawMatches(rawImg1, mat);
	for (int i = 0; i < mat.height(); i++){
		int x = mat[i * 4 + 0];
		int y = mat[i * 4 + 1];
		int labIdx = gKLabels[y*w + x];
		bestU[labIdx] = mat[i * 4 + 2] - x;
		bestV[labIdx] = mat[i * 4 + 3] - y;
		bestCosts[labIdx] = MatchCost(rawImg1, rawImg2, im1f, im2f, x, y, x + bestU[labIdx], y + bestV[labIdx]);
	}
	//ShowSuperPixelFlow("init + viso2", graph, im1, bestU, bestV, numV);
	//cvWaitKey();
#endif

	int* searchR = new int[nLevels];
	int* iterCnts = new int[nLevels];
	searchR[nLevels - 1] = initR;
	iterCnts[nLevels - 1] = 8;
	for (int i = 0; i < nLevels - 1; i++){
		searchR[i] = MAXR;
		iterCnts[i] = ITER_N;
	}

	IntImage seedGuessStore(300, numV); // 100 max number of guesses. format:(u,v,cost)
	seedGuessStore.setValue(-1);

	for (int l = nLevels - 1; l >= 0; l--){ // coarse-to-fine
		//show the patch flow
// 		FImage t1(pyd1[l]), t2(pyd2[l]);
// 		t1.imresize(pow(1. / ratio, l), INTER_NN);
// 		t2.imresize(pow(1. / ratio, l), INTER_NN);
// 		t1.imshow("pyd1");
// 		t2.imshow("pyd2");
// 		FImage tmpFlow(pydSeedsFlow[l]);
// 		tmpFlow.Multiplywith(pow(1. / ratio, l));
// 		ShowSuperPixelFlow("pyd", gKLabels, w, h, tmpFlow);
// 		cvWaitKey(0);

		if (l == 0){
			Propogate(pyd1, pyd2, im1f, im2f, l, searchR[l], iterCnts[l], pydSeeds, neighbors, 
				pydSeedsFlow, bestCosts, seedGuessStore, true);
		}else{
			Propogate(pyd1, pyd2, im1f, im2f, l, searchR[l], iterCnts[l], pydSeeds, neighbors,
				pydSeedsFlow, bestCosts, seedGuessStore, false);
			pydSeedsFlow[l - 1].copyData(pydSeedsFlow[l]);
			pydSeedsFlow[l - 1].Multiplywith(1./ratio);
		}
	}

#if 0
	// show guess store
	int tmpw = 200, tmph = 200;
	IntImage tmpGuessImg(tmpw, tmph);

	for (int vIdx = 0; vIdx < numV; vIdx++){
		FImage tmpImg(pyd1[0]);

		float bestu = pydSeedsFlow[0][2 * vIdx];
		float bestv = pydSeedsFlow[0][2 * vIdx + 1];
		tmpGuessImg.setValue(0);
		int* guessIdx = seedGuessStore.rowPtr(vIdx);
		int idx = 0;
		while (1){
			int tu = guessIdx[3 * idx];
			int tv = guessIdx[3 * idx + 1];
			int tc = guessIdx[3 * idx + 2];
			if (tc >= 0){
				int x = tu - bestu + tmpw / 2;
				int y = tv - bestv + tmph / 2;
				tmpGuessImg[y*tmpw + x] = tc;
				tmpGuessImg[(tmph / 2)*tmpw + (tmph / 2)] = 0;
				idx++;
			}else{
				break;
			}
		}

		// set superpixel
		float tmpVal[3] = { 0, 0, 1 };
		for (int ii = 0; ii < h; ii++){
			for (int jj = 0; jj < w; jj++){
				int id = gKLabels[ii*w + jj];
				if (id == vIdx){
					tmpImg.setPixel(ii, jj, tmpVal);
				}
			}
		}

		printf("%d count!\n", idx);
		tmpImg.imshow("tmpImg");
		tmpGuessImg.imagesc("guess", 0);
	}
#endif

#if 0
	Vector<float> coe(5);
	Vector<float> tx(1024), ty(1024), tz(1024);

	for (int i = 0; i < numV; i++){
		float bestu = pydSeedsFlow[0][2 * i];
		float bestv = pydSeedsFlow[0][2 * i + 1];
		int* guessIdx = seedGuessStore.rowPtr(i);
		int idx = 0;
		while (1){
			int tc = guessIdx[3 * idx + 2];
			if (tc >= 0){
				idx++;
			}else{
				break;
			}
		}

		// paraboloid fitting
		
		int fitCnt = idx;
		if (fitCnt >= 5){
			for (int ii = 0; ii < fitCnt; ii++){
				tx[ii] = guessIdx[3 * ii];
				ty[ii] = guessIdx[3 * ii + 1];
				tz[ii] = guessIdx[3 * ii + 2];
			}
			coe = ParaboloidFitting(tx, ty, tz, fitCnt);
			float rho1 = coe[0];
			float rho2 = coe[1];
			pydSeedsFlow[0][2 * i] = UNKNOWN_FLOW;
			pydSeedsFlow[0][2 * i + 1] = UNKNOWN_FLOW;
			if (rho1 > 0 && rho2 > 0){
				float tu = -coe[2] / (2 * rho1);
				float tv = -coe[3] / (2 * rho2);
				float dis = sqrt((tu - bestu)*(tu - bestu) + (tv - bestv)*(tv - bestv));
				//printf("<x,y>:<%f,%f>\t<sx,sy>:<%f,%f>\tdis: %f\n", bestu, bestv, tu, tv, dis);
				if (dis < 1)
				{
 					pydSeedsFlow[0][2 * i] = bestu;
 					pydSeedsFlow[0][2 * i + 1] = bestv;
					//pydSeedsFlow[0][2 * i] = tu;
					//pydSeedsFlow[0][2 * i + 1] = tv;
				}
			}
		}
	}
#endif

	//FImage matchImg(pyd1[0]);
	//ShowSparseFlow("match", matchImg, pydSeeds[0], pydSeedsFlow[0]);

#if 0
	FImage outlierImg(w, h);
	memset(outlierImg.pData, 0, w*h*sizeof(float));

	for (int i = 0; i < numV; i++){
		float x = pydSeeds[0][2 * i];
		float y = pydSeeds[0][2 * i + 1];
		float bestu = pydSeedsFlow[0][2 * i];
		float bestv = pydSeedsFlow[0][2 * i + 1];
		float tx = ImageProcessing::EnforceRange(x + bestu, w);
		float ty = ImageProcessing::EnforceRange(y + bestv, h);
		outlierImg[ty*w + tx] = 1;
	}
	outlierImg.imshow("outlier", 0);
#endif

#if 0
	// remove outliers ( if more pixels point to the same position )
	IntImage targetFlag(w, h);
	memset(targetFlag.pData, 0xFF, w*h*sizeof(int)); //-1
	for (int seedIdx = 0; seedIdx < numV; seedIdx++){
		float x = pydSeeds[0][2 * seedIdx];
		float y = pydSeeds[0][2 * seedIdx + 1];
		float bestu = pydSeedsFlow[0][2 * seedIdx];
		float bestv = pydSeedsFlow[0][2 * seedIdx + 1];
		int tx = (int)(x + bestu + 0.5);
		int ty = (int)(y + bestv + 0.5);
		float currCost = bestCosts[seedIdx];
		if (tx < 0 || tx >= w || ty < 0 || ty >= h){ // out of range
			pydSeedsFlow[0][2 * seedIdx] = UNKNOWN_FLOW;
			pydSeedsFlow[0][2 * seedIdx + 1] = UNKNOWN_FLOW;
			continue;
		}

		int r = 3;// STEP / 2; // radius
		bool valid = true;
		for (int i = ty - r; i <= ty + r; i++){
			for (int j = tx - r; j <= tx + r; j++){
				if (j < 0 || j >= w || i < 0 || i >= h){
					continue;
				}
				int nbIdx = targetFlag.pData[i*w + j];
				if (nbIdx >= 0){
					if (bestCosts[nbIdx] < currCost)
					{
						valid = false;
						goto CHECK_DONE;
					}
				}
			}
		}
	CHECK_DONE:
		if (valid){
			targetFlag.pData[ty*w + tx] = seedIdx;
		}else{
			pydSeedsFlow[0][2 * seedIdx] = UNKNOWN_FLOW;
			pydSeedsFlow[0][2 * seedIdx + 1] = UNKNOWN_FLOW;
		}
	}

#endif

	delete[] searchR;
	delete[] bestCosts;
}

void OnePass(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, FImage& costMap1, FImage& costMap2,
	IntImage& seeds, IntImage& neighbors, FImage* pydSeedsFlow, IntImage& nnf, FImage& dis, bool filtering = true)
{
	FImage rawImg1 = pyd1[0];
	FImage rawImg2 = pyd2[0];

	int nLevels = pyd1.nlevels();
	float ratio = pyd1.ratio();

	int numV = seeds.height();

	IntImage* pydSeeds = new IntImage[nLevels];
	for (int i = 0; i < nLevels; i++){
		pydSeeds[i].allocate(2, numV);
		int sw = pyd1[i].width();
		int sh = pyd1[i].height();
		for (int n = 0; n < numV; n++){
			pydSeeds[i][2 * n] = ImageProcessing::EnforceRange(seeds[2 * n] * pow(ratio, i), sw);
			pydSeeds[i][2 * n + 1] = ImageProcessing::EnforceRange(seeds[2 * n + 1] * pow(ratio, i), sh);
		}
	}

	PyramidRandomSearch(pyd1, pyd2, im1f, im2f, pydSeeds, neighbors, pydSeedsFlow);

	// scale
	int b = borderWidth;
	for (int i = 0; i < nLevels; i++){
// 		int sw = pyd1[i].width();
// 		int sh = pyd1[i].height();
// 		for (int n = 0; n < numV; n++){
// 			int sx = pydSeeds[i][2 * n];
// 			int sy = pydSeeds[i][2 * n + 1];
// 			if (sx < b || sx >= sw-b || sy < b || sy >= sh-b){
// 				pydSeedsFlow[i][2 * n] = 0;
// 				pydSeedsFlow[i][2 * n + 1] = 0;
// 			}
// 		}
		pydSeedsFlow[i].Multiplywith(pow(1./ratio, i));
	}

// 	if (filtering){
// 		int numV = seeds.height();
// 		float* candiU = new float[numV];
// 		float* candiV = new float[numV];
// 		int idx = 0;
// 		typedef std::unordered_map<int, float> umap;
// 		umap flowMap;
// 		for (int i = 0; i < numV; i++){
// 			short u = seedsFlow[2 * i];
// 			short v = seedsFlow[2 * i + 1];
// 			if (OpticFlowIO::unknown_flow(u, v))
// 				continue;
// 			int k = key(u, v);
// 			umap::const_iterator iter = flowMap.find(k);
// 			if (iter == flowMap.end()){
// 				candiU[idx] = u;
// 				candiV[idx] = v;
// 				idx++;
// 				flowMap.insert(umap::value_type(k, 1));
// 			}
// 		}
// 		FImage candiFlow(2, idx);
// 		for (int i = 0; i < idx; i++){
// 			candiFlow[2 * i] = candiU[i];
// 			candiFlow[2 * i + 1] = candiV[i];
// 		}
// 		delete[] candiU;
// 		delete[] candiV;
// 		TreeFilter(rawImg1, rawImg1, im1f, im2f, tree, seeds, seedsFlow, candiFlow);
// 
// 		//NNF_Filtering(rawImg1, rawImg2, im1f, im2f, seeds, seedsFlow, nnf, dis);
// 	}
	delete[] pydSeeds;
}
/*
void InitMatch(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, FImage& costMap1, FImage& costMap2,
	IntImage& seeds, IntImage& neighbors, float* bestU, float* bestV, int* unstable)
{
	CTimer t;

	FImage rawImg1 = pyd1[0];
	FImage rawImg2 = pyd2[0];

	int w = rawImg1.width();
	int h = rawImg1.height();
	int numV = seeds.height();

	PyramidRandomSearch(pyd1, pyd2, im1f, im2f, seeds, neighbors, bestU, bestV);

#ifdef DEBUG_SHOW
	//ShowSparseFlow("1", rawImg1, ptx, pty, len, bestU, bestV);
	ShowSuperPixelFlow("1", gKLabels, w, h, bestU, bestV);
#endif

	// reverse check
	float* bestU2 = new float[numV];
	float* bestV2 = new float[numV];
	int* ptx2 = new int[numV];
	int* pty2 = new int[numV];
	for (int i = 0; i < numV; i++){
		int x2 = seeds[2 * i] + bestU[i];
		int y2 = seeds[2 * i + 1] + bestV[i];
		ptx2[i] = x2;
		pty2[i] = y2;
		bestU2[i] = -bestU[i];
		bestV2[i] = -bestV[i];
	}

// 	FImage tmpImg(rawImg2);
// 	float color[3] = { 0, 0, 1 };
// 	for (int i = 0; i < numV; i++){
// 		//tmpImg.setPixel(pty[i], ptx[i], v);
// 		ptx2[i] = ImageProcessing::EnforceRange(ptx2[i], w);
// 		pty2[i] = ImageProcessing::EnforceRange(pty2[i], h);
// 		tmpImg.setPixel(pty2[i], ptx2[i], color);
// 	}
// 	tmpImg.imshow("tmp", 0);

	//ShowSuperPixelFlow("1", graph, im1, bestU, bestV, numV);

#if 1
	FImage u(w, h);
	FImage v(w, h);
	SeedMatch(rawImg2, rawImg1, costMap2, costMap1, im2f + 0, im1f + 0, ptx2, pty2, numV, bestU2, bestV2);
	for (int i = 0; i < numV; i++){
		float du = abs(bestU[i] + bestU2[i]);
		float dv = abs(bestV[i] + bestV2[i]);
		float diff = sqrt(du*du + dv*dv);
		if (diff > 5){
			unstable[i] = 1;
			bestU[i] = UNKNOWN_FLOW;
			bestV[i] = UNKNOWN_FLOW;
		}
	}
	// deprive outlier
	for (int i = 0; i < numV; i++){
		int* nbIdx = neighbors.rowPtr(i);
		int validNbCnt = 0;
		for (int j = 0; j < neighbors.width(); j++){
			int nb = nbIdx[j];
			if (nb < 0)
				break;
			if (!OpticFlowIO::unknown_flow(bestU[nb], bestV[nb])){
				validNbCnt++;
			}
		}
		if (validNbCnt < 3){
			unstable[i] = 1;
			bestU[i] = UNKNOWN_FLOW;
			bestV[i] = UNKNOWN_FLOW;
		}
	}
#endif

#ifdef DEBUG_SHOW
	ShowSuperPixelFlow("afterCheck", gKLabels, w, h, bestU, bestV);
#endif
	delete[] bestU2;
	delete[] bestV2;
	delete[] ptx2;
	delete[] pty2;
}
*/

#if 0
void ReConstructTree(CSuperPixelTree& spt, FImage& costMap, IntImage& nnf, FImage& dis)
{
	CSuperPixelGraph* graph = spt.GetGraph();
	CMst* tree = spt.GetTree();
	int numV = graph->NodeCnt();

	int* ptx = new int[numV];
	int* pty = new int[numV];

	for (int i = 0; i < numV; i++){
		int* p = graph->NodePos(i);
		ptx[i] = p[0];
		pty[i] = p[1];
	}

	IntImage kLabels;
	ConstructTree(ptx, pty, numV, costMap, tree, kLabels, nnf, dis);

	delete[] ptx;
	delete[] pty;
}

IntImage gNNF;
FImage gDis;
void OnMouse(int Event, int x, int y, int flags, void* param)
{
	switch (Event){
	case CV_EVENT_MOUSEMOVE:
		FImage* tmp = (FImage*) param;
		FImage newTmp(*tmp);
		newTmp.desaturate();

		int idx = gKLabels2[y*gw + x];
		printf("%d: ", idx);
		int nnfw = gNNF.width();
		for (int i = 0; i < 5; i++){
			int nbIdx = gNNF.pData[idx*nnfw + i];
			for (int k = 0; k < gw*gh; k++){
				if (gKLabels2[k] == nbIdx){
					newTmp[k] = 1;
				}
			}
			printf("%d: %.2f ", gNNF.pData[idx*nnfw + i], gDis.pData[idx*nnfw + i]);
			//printf("%.2f ", gDis.pData[idx*nnfw + i]);
		}
		printf("\n");
		newTmp.imshow("nnf");
		break;
	}
}

void WeightedAverage(IntImage& seeds, FImage& seedsFlow, IntImage& nnf, FImage& dis)
{
	int numV = seeds.height();
	int maxNbCnt = nnf.width();
	FImage tmpFlow(seedsFlow);
	for (int i = 0; i < numV; i++){
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		float sumu = 0, sumv = 0;
		float wt = 0;
		for (int j = 0; j <= 8; j++){
			int nb = nnf[i*maxNbCnt + j];
			float d = dis[i*maxNbCnt + j];
			if (nb < 0)
				break;
			float nu = seedsFlow[2 * nb];
			float nv = seedsFlow[2 * nb + 1];
			sumu += d*nu;
			sumv += d*nv;
			wt += d;
		}
		tmpFlow[2 * i] = sumu / wt;
		tmpFlow[2 * i + 1] = sumv / wt;
	}
	seedsFlow.copyData(tmpFlow);
}

void HistogramFilter(FImage& seedsFlow)
{
	int candiW = 2 * MAX_DISPLACEMENT + 1;
	int candiH = 2 * MAX_DISPLACEMENT + 1;
	IntImage candiFlowHisto(candiW, candiH);
	candiFlowHisto.setValue(0);
	int oy = candiH / 2;
	int ox = candiW / 2;
	int numV = seedsFlow.height();

	for (int i = 0; i < numV; i++){
		int u = seedsFlow[2 * i];
		int v = seedsFlow[2 * i + 1];
		if (!OpticFlowIO::unknown_flow(u, v)){
			candiFlowHisto.pData[(oy + v)*candiW + ox + u] = 255;
		}
	}
	candiFlowHisto.imshow("histo");
}

#endif

void CrossCheck(IntImage& seeds, FImage& seedsFlow, FImage& seedsFlow2, 
	IntImage& kLabel2, int* valid, float th)
{
	int w = kLabel2.width();
	int h = kLabel2.height();
	int numV = seeds.height();
	for (int i = 0; i < numV; i++){
		valid[i] = 1;
	}

	// cross check (1st step)
	int b = borderWidth;
	for (int i = 0; i < numV; i++){
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		int x2 = x + u;
		int y2 = y + v;
		if (x < b || x >= w - b || y < b || y >= h - b
			|| x2 < b || x2 >= w - b || y2 < b || y2 >= h - b
			|| sqrt(u*u + v*v)>MAX_DISPLACEMENT){
			valid[i] = 0;
			continue;
		}
#if 1
		int idx2 = kLabel2[y2*w + x2];
		float u2 = seedsFlow2[2 * idx2];
		float v2 = seedsFlow2[2 * idx2 + 1];
		float diff = sqrt((u + u2)*(u + u2) + (v + v2)*(v + v2));
		if (diff > th){
			valid[i] = 0;
		}
#endif
	}
}

void ShowResponseMap(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f, int x, int y, float scale)
{
	int w = img1.width();
	int h = img1.height();
	FImage map(w, h);
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			map[i*w + j] = MatchCost(img1, img2, im1f, im2f, x, y, j, i);
		}
	}
	map.Multiplywith(-1);
	//map.imwrite("d:/x.png");
	map.imresize(scale, INTER_NN);
	map.imagesc("response map", 0);
}

void SubSample(FImage& u, FImage& v)
{
	int w = u.width();
	int h = u.height();
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			double du = u.data()[i*w + j];
			double dv = v.data()[i*w + j];
			if (!OpticFlowIO::unknown_flow(du, dv)){
				for (int ii = -2; ii <= 2; ii++){
					for (int jj = -2; jj <= 2; jj++){
						if (ii == 0 && jj == 0)
							continue;
						int ni = ImageProcessing::EnforceRange(i + ii, h);
						int nj = ImageProcessing::EnforceRange(j + jj, w);
						u.data()[ni*w + nj] = UNKNOWN_FLOW;
						v.data()[ni*w + nj] = UNKNOWN_FLOW;
					}
				}
			}
		}
	}
}

void CPM_Flow(FImage& img1, FImage& img2, FImage& costMap1, FImage& costMap2, IntImage& kLabMat, 
	FImage& u, FImage& v, FImage& matches)
{
	CTimer t;

	int w = img1.width();
	int h = img1.height();

	float ratio = 0.5;
	int nLevels = LEVELS;
	//int nLevels = 1;
	FImagePyramid pyd1, pyd2;

    t.tic();
	pyd1.ConstructPyramidLevels(img1, ratio, nLevels);
	pyd2.ConstructPyramidLevels(img2, ratio, nLevels);
	
	t.toc("construct pyramid: ");

	UCImage* im1f = new UCImage[nLevels];
	UCImage* im2f = new UCImage[nLevels];
	for (int i = 0; i < nLevels; i++){
// 		ImageFeature::imCensus(pyd1[i], im1f[i], 8);
// 		ImageFeature::imCensus(pyd2[i], im2f[i], 8);
		ImageFeature::imSIFT(pyd1[i], im1f[i], 2, 1, true, 8);
		ImageFeature::imSIFT(pyd2[i], im2f[i], 2, 1, true, 8);
	}
	t.toc("get feature: ");

#if 0
	int ptx = 495;
	int pty = 82;
	int level = 0;
	int inx = ptx*pow(ratio, level);
	int iny = pty*pow(ratio, level);
	float scale = 1./pow(ratio, level);
	ShowResponseMap(pyd1[level], pyd2[level], im1f+level, im2f+level, inx, iny, scale);

	FImage crop1;
	img1.crop(crop1, ptx - 80, pty - 80, 160, 160);
	crop1.imshow("crop");

	cv::Mat cvImg1 = ImageIO::CvmatFromPixels(img1.pData, img1.width(), img1.height(), img1.nchannels());
	int matchRadius = 8;
	cv::rectangle(cvImg1, cvRect(ptx - matchRadius / 2, pty - matchRadius / 2, matchRadius, matchRadius), cvScalar(255, 255, 255));
	matchRadius = 16;
	cv::rectangle(cvImg1, cvRect(ptx - matchRadius / 2, pty - matchRadius / 2, matchRadius, matchRadius), cvScalar(255, 255, 255));
	matchRadius = 32;
	cv::rectangle(cvImg1, cvRect(ptx - matchRadius / 2, pty - matchRadius / 2, matchRadius, matchRadius), cvScalar(255, 255, 255));
	matchRadius = 64;
	cv::rectangle(cvImg1, cvRect(ptx - matchRadius / 2, pty - matchRadius / 2, matchRadius, matchRadius), cvScalar(255, 255, 255));
	cv::imshow("xx", cvImg1);
#endif
	
#if 0
	CSuperPixelGraph tmpGraph;
	tmpGraph.Init(img1, 200, 10);
	FImage tmp(img1);
	tmpGraph.AddSuperPixelMask(tmp);
	tmp.imshow("xx",0);

	float val[3] = {0,0,1};
	SetSuperPixel(tmpGraph,20,val,tmp);
	tmp.imshow("xx",0);
#endif

#if 0
	t.tic();
	FImage lab1, lab2;
	img1.ToLab(lab1);
	img2.ToLab(lab2);
	//lab1.imshow("lab1");
	//lab2.imshow("lab2", 0);
	t.toc("to lab: ");

	t.tic();
	CSuperPixelGraph graph1, graph2;
	graph1.Init(img1, w*h / (STEP*STEP), 20);
	graph2.Init(img2, w*h / (STEP*STEP), 20);
	t.toc("super pixel graph: ");

	graph1.ShowSuperPixel(img1);
	graph2.ShowSuperPixel(img2);

	int numV = graph1.NodeCnt();
	int numV2 = graph2.NodeCnt();

	FImage* pydSeedsFlow = new FImage[nLevels];
	FImage* pydSeedsFlow2 = new FImage[nLevels];
	for (int i = 0; i < nLevels; i++){
		pydSeedsFlow[i].allocate(2, numV);
		pydSeedsFlow2[i].allocate(2, numV2);
	}
	FImage seedsFlow(2, numV);

	IntImage seeds(2, numV), seeds2(2,numV2);
	IntImage neighbors(graph1.MaxNeighborCnt(), numV);
	IntImage neighbors2(graph2.MaxNeighborCnt(), numV2);
	neighbors.setValue(-1);
	neighbors2.setValue(-1);
	for (int i = 0; i < numV; i++){
		int* p = graph1.NodePos(i);
		seeds[2 * i] = p[0];
		seeds[2 * i + 1] = p[1];
		int nbNums = graph1.NeighborCnt(i);
		int* nbIdx = graph1.NeighborIndex(i);
		for (int j = 0; j < nbNums; j++){
			neighbors[i*neighbors.width() + j] = nbIdx[j];
		}
	}
	for (int i = 0; i < numV2; i++){
		int* p = graph2.NodePos(i);
		seeds2[2 * i] = p[0];
		seeds2[2 * i + 1] = p[1];
		int nbNums = graph2.NeighborCnt(i);
		int* nbIdx = graph2.NeighborIndex(i);
		for (int j = 0; j < nbNums; j++){
			neighbors2[i*neighbors.width() + j] = nbIdx[j];
		}
	}
	IntImage nnf, nnf2;
	FImage dis, dis2;

	IntImage kLabels(w, h), kLabels2(w, h);
	memcpy(kLabels.pData, graph1.KLabels(), w*h*sizeof(int));
	memcpy(kLabels2.pData, graph2.KLabels(), w*h*sizeof(int));
	//kLabels.imshow("kLabels", 0);
	gKLabels = kLabels.pData;
	gKLabels2 = kLabels2.pData;
#else
	int step = STEP;
	int gridw = w / step;
	int gridh = h / step;
	int xoffset = (w - (gridw - 1)*step) / 2;
	int yoffset = (h - (gridh - 1)*step) / 2;
	int numV = gridw * gridh;
	int numV2 = numV;

	FImage* pydSeedsFlow = new FImage[nLevels];
	FImage* pydSeedsFlow2 = new FImage[nLevels];
	for (int i = 0; i < nLevels; i++){
		pydSeedsFlow[i].allocate(2, numV);
		pydSeedsFlow2[i].allocate(2, numV2);
	}

	IntImage seeds(2, numV);
	IntImage neighbors(12, numV);
	neighbors.setValue(-1);
	int nbOffset[8][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 }, { -1, -1 }, { -1, 1 }, { 1, -1 }, { 1, 1 } };
	for (int i = 0; i < numV; i++){
		int gridX = i % gridw;
		int gridY = i / gridw;
		seeds[2 * i] = gridX * step + xoffset;
		seeds[2 * i + 1] = gridY * step + yoffset;
		int nbIdx = 0;
		for (int j = 0; j < 4; j++){
			int nbGridX = gridX + nbOffset[j][0];
			int nbGridY = gridY + nbOffset[j][1];
			if (nbGridX < 0 || nbGridX >= gridw || nbGridY < 0 || nbGridY >= gridh)
				continue;
			neighbors[i*neighbors.width() + nbIdx] = nbGridY*gridw + nbGridX;
			nbIdx++;
		}
	}
	IntImage seeds2(seeds);
	IntImage neighbors2(neighbors);
	FImage seedsFlow(2, numV), seedsFlow2(2, numV2);
	
	IntImage nnf, nnf2;
	FImage dis, dis2;

	IntImage kLabels(w,h), kLabels2;
	//Segmentation(seeds, costMap1, kLabels, nnf, dis);
	//Segmentation(seeds2, costMap2, kLabels2, nnf2, dis2);
	for (int i = 0; i < numV; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		int r = step / 2;
		for (int ii = -r; ii <= r; ii++){
			for (int jj = -r; jj <= r; jj++){
				int xx = ImageProcessing::EnforceRange(x + ii, w);
				int yy = ImageProcessing::EnforceRange(y + jj, h);
				kLabels[yy*w + xx] = i;
			}
		}
	}
	kLabels2.copy(kLabels);
	//kLabels.imshow("kLabels", 0);
	gKLabels = kLabels.pData;
	gKLabels2 = kLabels2.pData;

	//t.toc("generate seeds: ");
#endif

#ifdef DEBUG_SHOW
	FImage seedShow(w, h);
	seedShow.setValue(1);
	for (int i = 0; i < numV; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		float val = 0;
		seedShow.setPixel(y, x, &val);
	}
	seedShow.imshow("seeds");

#endif
	//printf("seeds num: %d\n", numV);
	//printf("seeds num2: %d\n", numV2);

	//t.tic();
	//IntImage nnf, nnf2;
	//FImage dis, dis2;
	//ReConstructTree(spt1, costMap1, nnf, dis);
	//ReConstructTree(spt2, costMap2, nnf2, dis2);
	//t.toc("reconstruct tree: ");

#ifdef DEBUG_SHOW
	//spt1.ShowTree(img1, "tree1", true);
	//spt2.ShowTree(img2, "tree2", true);
// 	FImage tmp(img2);
// 	gNNF.copyData(nnf2);
// 	gDis.copyData(dis2);
// 	graph2->AddSuperPixelMask(tmp);
// 	tmp.imshow("nnf");
// 	cvSetMouseCallback("nnf", OnMouse, &tmp);
// 	cvWaitKey();
#endif

	t.tic();
	OnePass(pyd1, pyd2, im1f, im2f, costMap1, costMap2, seeds, neighbors, pydSeedsFlow, nnf, dis, false);
	t.toc("forward matching: ");
	//WeightedAverage(seeds, seedsFlow, nnf, dis);
	OnePass(pyd2, pyd1, im2f, im1f, costMap2, costMap1, seeds2, neighbors2, pydSeedsFlow2, nnf2, dis2, false);
	t.toc("backward matching: ");
	//WeightedAverage(seeds2, seedsFlow2, nnf2, dis2);

#ifdef DEBUG_SHOW
// 	for (int i = nLevels - 1; i >= 0; i--){
// 		ShowSuperPixelFlow("1", gKLabels, w, h, pydSeedsFlow[i]);
// 		ShowSuperPixelFlow("2", gKLabels2, w, h, pydSeedsFlow2[i]);
// 		cvWaitKey();
// 	}
#endif

	// cross check
	IntImage labs2(w, h);
	memcpy(labs2.pData, gKLabels2, w*h*sizeof(int));
	IntImage validFlag;
	validFlag.allocate(numV, nLevels);
	float baseTh = CHECK_TH;
	float* th = new float[nLevels];
	for (int i = 0; i < nLevels; i++){
		th[i] = baseTh;// *pow(1. / ratio, i);
	}
	for (int i = 0; i < nLevels; i++){
		CrossCheck(seeds, pydSeedsFlow[i], pydSeedsFlow2[i], labs2, validFlag.rowPtr(i), th[i]);
	}

#ifdef DEBUG_SHOW
// 	for (int l = 0; l < nLevels; l++){
// 		for (int i = 0; i < numV; i++){
// 			if (!validFlag[l*numV + i]){
// 				pydSeedsFlow[l][2 * i] = UNKNOWN_FLOW;
// 				pydSeedsFlow[l][2 * i + 1] = UNKNOWN_FLOW;
// 			}
// 		}
// 	}
	for (int i = nLevels - 1; i >= 0; i--){
		ShowSeedsFlow("1", w, h, seeds, pydSeedsFlow[i]);
		ShowSeedsFlow("2", w, h, seeds2, pydSeedsFlow2[i]);
		//FImage tmp(img1);
		//ShowSparseFlow("1", tmp, seeds, pydSeedsFlow[i]);
		//ShowSuperPixelFlow("1", gKLabels, w, h, pydSeedsFlow[i]);
		//ShowSuperPixelFlow("2", gKLabels2, w, h, pydSeedsFlow2[i]);
		cvWaitKey();
	}
#endif

	seedsFlow.copyData(pydSeedsFlow[0]);

	if (FB_CHECK){
		int* errorCnt = new int[numV];
		int* unstable = new int[numV];
		memset(errorCnt, 0, sizeof(int)*numV);
		memset(unstable, 0, sizeof(int)*numV);
		for (int i = 0; i < numV; i++){
			for (int l = 0; l < nLevels; l++){
				if (!validFlag[l*numV + i]){
					errorCnt[i]++;
				}
			}
			//printf("%d ", validCnt[i]);
		}

		int b = BORDER_WIDTH;
		int lvlCnt = min(nLevels, CHECK_LEVELS);
		for (int i = 0; i < numV; i++){

			for (int j = 0; j < lvlCnt; j++){
				if (!validFlag[j * numV + i]){
					unstable[i] = 1;
					break;
				}
			}

			int x = seeds[2 * i];
			int y = seeds[2 * i + 1];
			if (x < b || x >= w - b || y < b || y >= h - b){
				unstable[i] = 1;
			}
		}

		for (int i = 0; i < numV; i++){
			if (unstable[i]){
				seedsFlow[2 * i] = UNKNOWN_FLOW;
				seedsFlow[2 * i + 1] = UNKNOWN_FLOW;
			}
		}
		delete[] errorCnt;
		delete[] unstable;
	}

	delete[] th;
	
	// cross check
// 	IntImage labs2(w, h);
// 	memcpy(labs2.pData, gKLabels2, w*h*sizeof(int));
// 	int* valid = new int[numV];
// 	CrossCheck(seeds, pydSeedsFlow[s], pydSeedsFlow2[s], labs2, valid, th);
// 	seedsFlow.copyData(pydSeedsFlow[s]);
// 	for (int i = 0; i < numV; i++){
// 		if(!valid[i]){
// 			seedsFlow[2*i] = UNKNOWN_FLOW;
// 			seedsFlow[2*i+1] = UNKNOWN_FLOW;
// 		}
// 	}
// 	delete[] valid;

	// deprive outlier
#if 0
	for (int i = 0; i < numV; i++){
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		if (OpticFlowIO::unknown_flow(u,v))
			continue;
		int* nbIdx = neighbors.rowPtr(i);
		int validNbCnt = 0;
		for (int j = 0; j < neighbors.width(); j++){
			int nb = nbIdx[j];
			if (nb < 0)
				break;
			u = seedsFlow[2 * nb];
			v = seedsFlow[2 * nb + 1];
			if (!OpticFlowIO::unknown_flow(u,v)){
				validNbCnt++;
			}
		}
		if (validNbCnt == 0){
			seedsFlow[2 * i] = UNKNOWN_FLOW;
			seedsFlow[2 * i + 1] = UNKNOWN_FLOW;
		}
	}
#endif

	//HistogramFilter(seedsFlow);

#if 0
	float* candiU = new float[numV];
	float* candiV = new float[numV];
	int idx = 0;
	typedef std::unordered_map<int, float> umap;
	umap flowMap;
	for (int i = 0; i < numV; i++){
		short u = seedsFlow[2 * i];
		short v = seedsFlow[2 * i + 1];
		if (OpticFlowIO::unknown_flow(u, v))
			continue;
		int k = key(u, v);
		umap::const_iterator iter = flowMap.find(k);
		if (iter == flowMap.end()){
			candiU[idx] = u;
			candiV[idx] = v;
			idx++;
			flowMap.insert(umap::value_type(k, 1));
		}
	}
	FImage candiFlow(2, idx);
	for (int i = 0; i < idx; i++){
		candiFlow[2 * i] = candiU[i];
		candiFlow[2 * i + 1] = candiV[i];
	}
	delete[] candiU;
	delete[] candiV;
	printf("candi flow count: %d\n", idx);
	TreeFilter(img1, img2, im1f, im2f, spt1.GetTree(), seeds, seedsFlow, candiFlow);
#endif
	//NNF_Filtering(img1, img2, im1f, im2f, seeds, seedsFlow, nnf, dis);

	//t.toc("cross check: ");

// 	for (int i = 0; i < numV; i++){
// 		seedsFlow[2 * i] = bestU[i];
// 		seedsFlow[2 * i + 1] = bestV[i];
// 	}
// 	//ShowSuperPixelFlow("before", gKLabels, w, h, bestU, bestV);
// 	NNF_Filtering(img1, img2, im1f, im2f, seeds, seedsFlow, nnf, dis);
	//ShowSuperPixelFlow("after", gKLabels, w, h, bestU, bestV);
	
	//TreeSegmentation(tree, isSubRoots, 0.04, 50);

	// show tree segmentation
#if 0
	FImage tmpImg(w, h, 3);
	float color[3] = { 0, 1, 1 };
	srand(0);
	for (int i = 0; i < numV; i++){
		if (!isSubRoots[i])
			continue;
		int rootId = i;
		int parent;
		std::queue<int> Q;
		Q.push(rootId);
		while (!Q.empty()){
			parent = Q.front(); Q.pop();
			SetSuperPixel(*graph1, parent, color, tmpImg);
			int childNum = tree->ChildNum(parent);
			int* childIdx = tree->ChildIndex(parent);
			for (int c = 0; c < childNum; c++){
				int child = childIdx[c];
				if (!isSubRoots[child]){
					Q.push(child);
				}
			}
		}
		color[0] = rand() % 80 / 100. + 0.2;
		color[1] = rand() % 80 / 100. + 0.2;
		color[2] = rand() % 80 / 100. + 0.2;

	}
	tmpImg.imshow("seg");
#endif

	//FilteringOnSegmentation(tree, ptx, pty, numV, img1, img2, im1f, im2f, isSubRoots, bestU, bestV);

// 	for (int i = 0; i < numV; i++){
// 		int x2 = seeds[2 * i] + bestU[i];
// 		int y2 = seeds[2 * i + 1] + bestV[i];
// 		if (x2 < 0 || x2 >= w || y2 < 0 || y2 >= h){
// 			bestU[i] = UNKNOWN_FLOW;
// 			bestV[i] = UNKNOWN_FLOW;
// 		}
// 	}
	
	//NonLocalFiltering(spt1.GetTree(), seeds, seedsFlow);

#if 0
	t.tic();
	spt1.Init(costMap1, img1, SP_SIZE);
	int* spLabels = spt1.GetGraph()->KLabels();
	t.toc("super pixel: ");
	IntImage seedMap(w, h);
	memset(seedMap.pData, 0xff, w*h*sizeof(int));
	for (int i = 0; i < numV; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		if (!OpticFlowIO::unknown_flow(u, v)){
			seedMap[y*w + x] = i;
		}
	}
	//seedMap.imshow("seedmap", 0);
	int spCnt = spt1.NodeCnt();
	IntImage spSeeds(2, spCnt);
	FImage spFlow(2, spCnt);
	for (int i = 0; i < spCnt; i++){
		spSeeds[2 * i] = -1;
		spSeeds[2 * i + 1] = -1;
		spFlow[2 * i] = UNKNOWN_FLOW;
		spFlow[2 * i + 1] = UNKNOWN_FLOW;

		int *p = spt1.NodePos(i);
		int x = p[0], y = p[1];
		float s = 0;
		float d = 0;
		float u = 0;
		float v = 0;
		int* rect = spt1.GetGraph()->NodeRect(i);
		int left = rect[0];
		int right = rect[1];
		int top = rect[2];
		int bottom = rect[3];
		for (int ii = top; ii <= bottom; ii++){
			for (int jj = left; jj <= right; jj++){
				if (spLabels[ii*w + jj] == i){
					int seedIdx = seedMap[ii*w + jj];
					if (seedIdx >= 0){
						int sx = seeds[2 * seedIdx];
						int sy = seeds[2 * seedIdx + 1];
						float su = seedsFlow[2 * seedIdx];
						float sv = seedsFlow[2 * seedIdx + 1];
						d = 0;// sqrt((x - sx)*(x - sx) + (y - sy)*(y - sy));
						d = exp(-d);
						u += d*su;
						v += d*sv;
						s += d;
					}
				}
			}
		}

		if (s > 0){
			u /= s;
			v /= s;
			float x2 = x + u;
			float y2 = y + v;
			if (x2 >= 0 && x2 < w && y2 >= 0 && y2 < h){
				spSeeds[2 * i] = x;
				spSeeds[2 * i + 1] = y;
				spFlow[2 * i] = u;
				spFlow[2 * i + 1] = v;
			}

		}
	}
	t.toc("sp flow: ");

	// save match
	match2flow(spLabels, w, h, spFlow, u, v);
	FImage tmpMatch(4, spCnt);
	tmpMatch.setValue(-1);
	int validMatCnt = 0;
	for (int i = 0; i < spCnt; i++){
		int x = spSeeds[2 * i];
		int y = spSeeds[2 * i + 1];
		float u = spFlow[2 * i];
		float v = spFlow[2 * i + 1];
		float x2 = x + u;
		float y2 = y + v;
		if (!OpticFlowIO::unknown_flow(u, v)){
			tmpMatch[4 * i + 0] = x;
			tmpMatch[4 * i + 1] = y;
			tmpMatch[4 * i + 2] = x2;
			tmpMatch[4 * i + 3] = y2;
			validMatCnt++;
		}
	}
	if (!matches.matchDimension(4, validMatCnt, 1)){
		matches.allocate(4, validMatCnt, 1);
	}
	int tmpIdx = 0;
	for (int i = 0; i < spCnt; i++){
		if (tmpMatch[4 * i + 0] >= 0){
			memcpy(matches.rowPtr(tmpIdx), tmpMatch.rowPtr(i), sizeof(float) * 4);
			tmpIdx++;
		}
	}
#else
	//match2flow(seeds, seedsFlow, w, h, u, v);
	match2flow(gKLabels, w, h, seedsFlow, u, v);
	FImage tmpMatch(4, numV);
	tmpMatch.setValue(-1);
	int validMatCnt = 0;
	for (int i = 0; i < numV; i++){
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		float x2 = x + u;
		float y2 = y + v;
		if (!OpticFlowIO::unknown_flow(u, v)){
			tmpMatch[4 * i + 0] = x;
			tmpMatch[4 * i + 1] = y;
			tmpMatch[4 * i + 2] = x2;
			tmpMatch[4 * i + 3] = y2;
			validMatCnt++;
		}
	}
	if (!matches.matchDimension(4, validMatCnt, 1)){
		matches.allocate(4, validMatCnt, 1);
	}
	int tmpIdx = 0;
	for (int i = 0; i < numV; i++){
		if (tmpMatch[4 * i + 0] >= 0){
			memcpy(matches.rowPtr(tmpIdx), tmpMatch.rowPtr(i), sizeof(int) * 4);
			tmpIdx++;
		}
	}
#endif

#if 0
	// subsample
	if (STEP == 1){
		SubSample(u, v);
		flow2match(u, v, matches);
	}
#endif

	//t.toc("generate matches: ");
	//ShowSparseFlow("final match", img1, seeds, bestU, bestV);

	//printf("%f\n", tmp);

	delete[] im1f;
	delete[] im2f;
	delete[] pydSeedsFlow;
	delete[] pydSeedsFlow2;
}

void FullySearch(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f, FImage& u, FImage& v)
{
	CSuperPixelTree spt;
	spt.Init(img1, img1, 100);

	CMst* tree = spt.GetTree();

	float color[3] = { 0, 0, 1 };

	int numV = spt.NodeCnt();
	int R = MAX_DISPLACEMENT; // search window radius

	CTimer t;

	float* nodeU = new float[numV];
	float* nodeV = new float[numV];
	for (int n = 0; n < numV; n++){
		int* pos = spt.NodePos(n);
		int x = pos[0];
		int y = pos[1];

		float minCost = FLT_MAX;
		float uBest, vBest;
		for (int i = -R; i < R; i+=1){
			for (int j = -R; j < R; j+=1){
				float cost = MatchCost(img1, img2, im1f, im2f, x, y, x + i, y + j);
				if (cost < minCost){
					uBest = i;
					vBest = j;
					minCost = cost;
				}
			}
		}
		nodeU[n] = uBest;
		nodeV[n] = vBest;
	}

	t.toc("fullly search: ");

	//ShowSuperPixelFlow("fullly search", *spt.GetGraph(), img1, nodeU, nodeV, numV);
	delete[] nodeU;
	delete[] nodeV;
}

void ReadEdges(char* fileName, FImage& edge)
{
	int w = edge.width();
	int h = edge.height();
	FILE* fp = fopen(fileName, "rb");
	int size = fread(edge.pData, sizeof(float), w*h, fp);
	assert(size == w*h);
	fclose(fp);
}

// read the K Labels of over-segmentation
void ReadKlabels(const char *fileName, IntImage& kLabMat)
{
	int w = kLabMat.width();
	int h = kLabMat.height();
	FILE* fp = fopen(fileName, "rb");
	int size = fread(kLabMat.pData, sizeof(int), w*h, fp);
	assert(size == w*h);
	fclose(fp);
}

// 	img1.imread("d:/1.png");
// 	img2.imread("d:/2.png");
// 	img1.imread("F:/Database_OF/Middlebury/eval-data/Backyard/frame10.png");
// 	img2.imread("F:/Database_OF/Middlebury/eval-data/Backyard/frame11.png");
// 	img1.imread("F:/Database_OF/KITTI/testing/image_0/000114_10.png");
// 	img2.imread("F:/Database_OF/KITTI/testing/image_0/000114_11.png");
// 	img1.imread("F:/Database_OF/MPI-Sintel/test/final/market_4/frame_0025.png");
// 	img2.imread("F:/Database_OF/MPI-Sintel/test/final/market_4/frame_0026.png");

// char* img1Name = "D:/aa.png";
// char* img2Name = "D:/bb.png";

// char* img1Name = "F:/Database_OF/Middlebury/eval-data/Teddy/frame10.png";
// char* img2Name = "F:/Database_OF/Middlebury/eval-data/Teddy/frame11.png";

// char* img1Name = "F:/Database_OF/MPI-Sintel/training/clean/ambush_5/frame_0002.png";
// char* img2Name = "F:/Database_OF/MPI-Sintel/training/clean/ambush_5/frame_0003.png";

// char* img1Name = "F:/Database_OF/MPI-Sintel/training/clean/ambush_6/frame_0005.png";
// char* img2Name = "F:/Database_OF/MPI-Sintel/training/clean/ambush_6/frame_0006.png";

// char* img1Name = "F:/Database_OF/MPI-Sintel/training/final/alley_1/frame_0014.png";
// char* img2Name = "F:/Database_OF/MPI-Sintel/training/final/alley_1/frame_0015.png";

// char* img1Name = "F:/Database_OF/MPI-Sintel/test/final/market_1/frame_0012.png";
// char* img2Name = "F:/Database_OF/MPI-Sintel/test/final/market_1/frame_0013.png";

char* img1Name = "F:/Database_OF/MPI-Sintel/training/final/mountain_1/frame_0045.png";
char* img2Name = "F:/Database_OF/MPI-Sintel/training/final/mountain_1/frame_0046.png";

// char* img1Name = "F:/Database_OF/MPI-Sintel/training/final/temple_3/frame_0005.png";
// char* img2Name = "F:/Database_OF/MPI-Sintel/training/final/temple_3/frame_0006.png";

// char* img1Name = "F:/Database_OF/KITTI2015/training/image_2/000019_10.png";
// char* img2Name = "F:/Database_OF/KITTI2015/training/image_2/000019_11.png";

int main(int argc, char** argv)
{
	CTimer t;
#if 1
	FImage img1, img2, img1g, img2g;

#ifdef DEBUG_SHOW
	img1.imread(img1Name);
	img2.imread(img2Name);
#else

	if (argc < 4){
        printf("USAGE: ./cpm img1Name img2Name outMatchName <c> <d> <r> <k> <n>\n");
		printf("Where:\n");
		printf("c: forward-backward check flag <default: 1>\n");
		printf("d: grid spacing <default: 3>\n");
		printf("r: search radius <default: 4>\n");
		printf("k: pyramid levels <default: 5>\n");
		printf("n: iteration times <default: 6>\n");
        return 1;
	}

	if (argc >= 5){
		FB_CHECK = atoi(argv[4]);
	}
	if (argc >= 6){
		STEP = atoi(argv[5]);
		CHECK_TH = STEP;
	}
	if (argc >= 7)
		MAXR = atoi(argv[6]);
	if (argc >= 8)
		LEVELS = atoi(argv[7]);
	if (argc >= 9)
		ITER_N = atoi(argv[8]);

	printf("------------\nCPM v1.1\n------------\n");
	printf("c: %d\t", FB_CHECK);
	printf("d: %d\t", STEP);
	printf("r: %d\t", MAXR);
	printf("k: %d\t", LEVELS);
	printf("n: %d\n", ITER_N);

	img1.imread(argv[1]);
	img2.imread(argv[2]);
	//img1EdgeName = argv[3];
	//img2EdgeName = argv[4];
#endif

	int w = img1.width();
	int h = img2.height();
#ifdef DEBUG_SHOW
	img1.imshow("img1");
	img2.imshow("img2");
#endif
	//img1.desaturate(img1g);
	//img2.desaturate(img2g);

	FImage lab1, lab2;
// 	img1.ToLab(lab1);
// 	img2.ToLab(lab2);

// 	img1g.imshow("gray1");
// 	img2g.imshow("gray2");

	FImage edgeR1(w, h, 1);
	FImage edgeR2(w, h, 1);
// 	ReadEdges(img1EdgeName, edgeR1);
// 	ReadEdges(img2EdgeName, edgeR2);
// 	edgeR1.Add(0.01f); // the raw 0.001f is too small
// 	edgeR2.Add(0.01f);
#ifdef DEBUG_SHOW
 	//edgeR1.imshow("edgeR1");
 	//edgeR2.imshow("edgeR2");
#endif
	IntImage kLabMat(w, h);
// 	ReadKlabels("d:/sp.bin", kLabMat);
// 	//kLabMat.imshow("kLabels", 0);
// 	CSuperPixelGraph::ArrangeLabels(kLabMat.pData, w, h);
	//kLabMat.imshow("kLabels1", 0);
	//t.toc("prepare data: ");

// 
#ifdef DEBUG_SHOW
// 	FImage tmp1f, tmp2f;
// 	im1f->collapse(tmp1f, collapse_max);
// 	im2f->collapse(tmp2f, collapse_max);
// 	tmp1f.imshow("f1");
// 	tmp2f.imshow("f2", 0);
#endif

	CTimer totalT;

    totalT.tic();
	FImage u(w, h), v(w, h);
	FImage matches;
	//FullySearch(img1, img2, im1f, im2f, u, v);
	//TreeFlow(img1, img2, &im1f, &im2f, edge1, kLabMat, u, v);
	CPM_Flow(img1, img2, edgeR1, edgeR2, kLabMat, u, v, matches);
	
	totalT.toc("total time: ");

#ifdef DEBUG_SHOW
	//OpticFlowIO::WriteFlowFile(u.pData, v.pData, w, h, "d:/in.flo");
	OpticFlowIO::ShowFlow("final", u.pData, v.pData, w, h);
	cvWaitKey(0);
#else
	char outName[256];
	strcpy(outName, argv[3]);
	strcat(outName, ".png");
	//OpticFlowIO::WriteFlowFile(mu.pData, mv.pData, w, h, "unrefined_flow.flo");
	OpticFlowIO::SaveFlowAsImage(outName, u.pData, v.pData, w, h);
	//strcpy(outName, argv[3]);
	//strcat(outName, ".flo");
	//OpticFlowIO::WriteFlowFile(u.pData, v.pData, w, h, outName);
	WriteMatches(argv[3], matches);
#endif
	
    return 1;

#endif

	///////////////////////////////
	FImage img;
	img.imread("d:/a.png");

	CSuperPixelTree spTree;
	t.tic();
	spTree.Init(img, img, 100);
	t.toc("spt: ");
	//spTree.ShowTree(img, "spt");

	CTreeFilter tf;
	tf.Init(spTree.GetTree());
	tf.SetSigma(0.1);
	int numV = spTree.NodeCnt();
	float* cost = new float[numV];
	for (int i = 0; i < numV; i++){
		cost[i] = i;
	}
	t.tic();
	for (int i = 0; i < 10000; i++){
		tf.Filter(cost, numV);
	}
	t.toc("tf: ");
	delete[] cost;
// 	UCImage ct;
// 	img.desaturate();
// 	img.imshow("img", false);
// 	t.tic();
// 	img.CensusTransform(ct);
// 	t.toc("ct: ");
// 	ct.imshow("ct");

	//printf("%d\n", hd);
// 	FImage matches;
// 	ReadMatches("d:/mat.txt", matches);
// 	DrawMatches(img, matches);
// 	img.imshow("match", 0);
    return 1;
}


#include "Image.h"
#include "OpticFlowIO.h"


#include <stdlib.h>
#include <string.h>


#include "epic.h"
extern "C" {
#include "image.h"
#include "io.h"
#include "variational.h"
}

/* read matches, stored as x1 y1 x2 y2 per line (other values on the same is not taking into account */
int ReadMatches(const char *filename, FImage& outMat)
{
	float* tmp = new float[4 * 100000]; // max number of match pair
	FILE *fid = fopen(filename, "r");
	int nmatch = 0;
	float x1, x2, y1, y2;
	while (!feof(fid) && fscanf(fid, "%f %f %f %f%*[^\n]", &x1, &y1, &x2, &y2) == 4)
	{
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
	return nmatch;
}

void WriteMatches(const char *filename, FImage& inMat)
{
	int len = inMat.height();
	FILE *fid = fopen(filename, "w");
	for (int i = 0; i < len; i++)
	{
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

void Match2Flow(FImage& inMat, FImage& ou, FImage& ov,  int w, int h)
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

		for (int jj = -5; jj <= 4; jj++){
			for (int ii = -5; ii <= 4; ii++){
				int nj = ImageProcessing::EnforceRange(y + jj, h);
				int ni = ImageProcessing::EnforceRange(x + ii, w);
				ou[nj*w + ni] = u;
				ov[nj*w + ni] = v;
			}
		}
		//ou[y*w + x] = u;
		//ov[y*w + x] = v;
	}
	for (int m=0; m<ou.nelements(); m++) {
		if(ou[m]>9999999)
			ou[m] = 0;
		if(ov[m]>9999999)
			ov[m] = 0;
	}
	
}

void FImage2image_t(FImage mu, FImage mv, image_t* wx, image_t* wy, int w, int h, int num_match) {
/*
    for (int j = 0; j < h; j++) {
	for (int i = 0; i < w; i++) {
	    wx->data[j*w+i] = mu.pData[j*w+i];
	    wy->data[j*w+i] = mv.pData[j*w+i];
	}
    }
*/
	//memcpy(wx->data, mu.pData, mu.nelements() * sizeof(float));
	//memcpy(wy->data, mv.pData, mv.nelements() * sizeof(float));
/*
	wx->data = mu.pData;
	wy->data = mv.pData;
*/
/*
	OpticFlowIO::WriteFlowFile(mu.pData, mv.pData, w, h, "unrefined_flow.flo");
	image_t** unrefine_flow = readFlowFile("unrefined_flow.flo");
	wx = unrefine_flow[0];
	wy = unrefine_flow[1];
	//writeFlowFile("output.flo", wx, wy);

*/
}


void image_t2FImage(image_t* wx, FImage mu) {
	mu.pData = wx->data;
}


int main(int argc, char** argv)
{
    if (argc != 6){
        printf("USAGE: ./Match2Flow img1 img2 matches.txt edges.tmp outputfile.flo\n");
        return 0;
    }



	///////////////////////////////////////////////////////
	//set parameters for refinement from EpicFlow
    	color_image_t *im1 = color_image_load(argv[1]);
    	color_image_t *im2 = color_image_load(argv[2]);


	//original code	
	//FImage mat;
	//int num_match = ReadMatches(argv[1], mat);
	float_image matches = read_matches(argv[3]);
	float_image edges = read_edges(argv[4], im1->width, im1->height);
    	//int w = atoi(argv[3]);
    	//int h = atoi(argv[4]);
	const char *outputfile = argv[5];


	epic_params_t epic_params;
	epic_params_default(&epic_params);
	variational_params_t flow_params;
    	variational_params_default(&flow_params);
        image_t *wx = image_new(im1->width, im1->height), *wy = image_new(im1->width, im1->height);

	//origianl code
	//slightly modefied original code
	/*
	cout<<"bp1!"<<endl;
	FImage mu, mv;
	Match2Flow(mat, mu, mv, w, h);
	*/


	///////////////////////////////////////////////////////

	//cout<<"elements:"<<mu.nelements()<<endl;
	cout<<"bp2!"<<endl;
	//FImage2image_t(mu, mv, wx, wy, w, h, num_match);

	//add interpolation from EpicFlow
	//OpticFlowIO::WriteFlowFile(mu.pData, mv.pData, w, h, "unrefined_flow.flo");
	//image_t **unrefine_flow = readFlowFile("unrefined_flow.flo");
	//wx = unrefine_flow[0];
	//wy = unrefine_flow[1];
	
	color_image_t *imlab = rgb_to_lab(im1);
    	epic(wx, wy, imlab, &matches, &edges, &epic_params, 1);

	//add refinement from EpicFlow
	cout<<"bp2.5!"<<endl;
	variational(wx, wy, im1, im2, &flow_params);
	cout<<"bp3!"<<endl;
	writeFlowFile(outputfile, wx, wy);
	cout<<"bp4!"<<endl;

	/*
	FImage mu2, mv2;
	image_t2FImage(wx,mu2);
	image_t2FImage(wy,mv2);
	cout<<"bp2"<<endl;
	*/
	///////////////////////////////////////////////////////

	//original code
	/*
	float maxFlow = -1;
	bool needUpdate = true;
	while (1){
		if (needUpdate){
			maxFlow = OpticFlowIO::ShowFlow("mat flow", mu.pData, mv.pData, w, h, maxFlow);
			OpticFlowIO::WriteFlowFile(mu.pData, mv.pData, w, h, "flow.flo");
			OpticFlowIO::SaveFlowAsImage("flow.png", mu.pData, mv.pData, w, h);
			printf("Max Flow: %.2f\n", maxFlow);
			needUpdate = false;
		}
		int key = cv::waitKey(20);
		int lowCode = key & 0xFFFF;
		int highCode = key >> 16;
		//printf("%d,%d\n", highCode, lowCode);

		switch (lowCode){
		case 0: // non-char
			switch (highCode)
			{
			case 38: // up
				maxFlow = (int)maxFlow / 5 * 5 + 5;
				break;
			case 40: // down
				maxFlow = (int)maxFlow / 5 * 5 - 5;
				break;
			}
			needUpdate = true;
			break;
		case 27: //ESC
			break;
		}
	}
	*/

	return 0;
}

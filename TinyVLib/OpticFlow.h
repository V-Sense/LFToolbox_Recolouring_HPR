#ifndef _OPTIC_FLOW_H
#define _OPTIC_FLOW_H

#include "Image.h"
#include "OpticFlowIO.h"

void ConsistencyMap(FImage& u1, FImage& v1, FImage& u2, FImage& v2, FImage& outMap)
{
	int w = u1.width();
	int h = u1.height();

	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			float du = u1[i*w + j];
			float dv = v1[i*w + j];

			float rDu, rDv;
			rDu = ImageProcessing::BilinearInterpolate(u2.data(), w, h, j + du, i + dv);
			rDv = ImageProcessing::BilinearInterpolate(v2.data(), w, h, j + du, i + dv);

			float offset = sqrt((du + rDu)*(du + rDu) + (dv + rDv)*(dv + rDv));
			outMap[i*w + j] = offset;
		}
	}
}

void ConsistencyCheck(FImage& u1, FImage& v1, FImage& u2, FImage& v2, float th)
{
	int w = u1.width();
	int h = u1.height();

	FImage beliefMap(w, h, 1);
	ConsistencyMap(u1, v1, u2, v2, beliefMap);
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			float offset = beliefMap[i*w + j];
			if (offset >= th){
				u1[i*w + j] = UNKNOWN_FLOW;
				v1[i*w + j] = UNKNOWN_FLOW;
			}
		}
	}
}

#endif
#pragma once

#include "Image.h"
#include "NoiseModel.h"
#include "Vector.h"
#include <vector>

typedef float _FlowPrecision;

class OpticalFlow
{
public:
	static bool IsDisplay;
public:
	enum InterpolationMethod {Bilinear,Bicubic};
	static InterpolationMethod interpolation;
	enum NoiseModel {GMixture,Lap};
	OpticalFlow(void);
	~OpticalFlow(void);
	static GaussianMixture GMPara;
	static Vector<float> LapPara;
	static NoiseModel noiseModel;
public:
	static void getDxs(FImage& imdx,FImage& imdy,FImage& imdt,const FImage& im1,const FImage& im2);
	static void getDxxs(FImage& imdxx, FImage& imdyy, FImage& imdxy, FImage& imdxt, FImage& imdyt, const FImage& im1, const FImage& im2);
	static void SanityCheck(const FImage& imdx,const FImage& imdy,const FImage& imdt,float du,float dv);
	static void warpFL(FImage& warpIm2,const FImage& Im1,const FImage& Im2,const FImage& vx,const FImage& vy);
	static void warpFL(FImage& warpIm2,const FImage& Im1,const FImage& Im2,const FImage& flow);


	static void genConstFlow(FImage& flow,float value,int width,int height);
	static void genInImageMask(FImage& mask,const FImage& vx,const FImage& vy,int interval = 0);
	static void genInImageMask(FImage& mask,const FImage& flow,int interval =0 );
	static void SmoothFlowPDE(const FImage& Im1,const FImage& Im2, FImage& warpIm2,FImage& vx,FImage& vy,
														 float alpha,int nOuterFPIterations,int nInnerFPIterations,int nCGIterations);
	
	static void SmoothFlowSOR(const FImage& Im1,const FImage& Im2, FImage& warpIm2, FImage& vx, FImage& vy,
														 float alpha,int nOuterFPIterations,int nInnerFPIterations,int nSORIterations);
	static void OneLevelFlow(const FImage& Im1, const FImage& Im2, FImage& warpIm2, FImage& vx, FImage& vy,
		float alpha, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations);
	static void OneLevelFlow_Brox(const FImage& Im1, const FImage& Im2, FImage& warpIm2, FImage& vx, FImage& vy,
		float alpha, float gamma, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations);

	static void estGaussianMixture(const FImage& Im1,const FImage& Im2,GaussianMixture& para,float prior = 0.9);
	static void estLaplacianNoise(const FImage& Im1,const FImage& Im2,Vector<float>& para);
	static void Laplacian(FImage& output,const FImage& input,const FImage& weight);
	static void testLaplacian(int dim=3);

	// function of coarse to fine optical flow
	static void Coarse2FineFlow(FImage& vx,FImage& vy,FImage &warpI2,const FImage& Im1,const FImage& Im2,float alpha,float ratio,int minWidth,
															int nOuterFPIterations,int nInnerFPIterations,int nCGIterations);

	static void Coarse2FineFlowLevel(FImage& vx,FImage& vy,FImage &warpI2,const FImage& Im1,const FImage& Im2,float alpha,float ratio,int nLevels,
															int nOuterFPIterations,int nInnerFPIterations,int nCGIterations);

	// function to convert image to features
	static void im2feature(FImage& imfeature,const FImage& im);

	// function to load optical flow
	static bool LoadOpticalFlow(const char* filename,FImage& flow);

	static bool LoadOpticalFlow(ifstream& myfile,FImage& flow);

	static bool SaveOpticalFlow(const FImage& flow, const char* filename);

	static bool SaveOpticalFlow(const FImage& flow,ofstream& myfile);

	static bool showFlow(const FImage& vx,const char* filename);

	// function to assemble and dissemble flows
	static void AssembleFlow(const FImage& vx,const FImage& vy,FImage& flow)
	{
		if(!flow.matchDimension(vx.width(),vx.height(),2))
			flow.allocate(vx.width(),vx.height(),2);
		for(int i = 0;i<vx.npixels();i++)
		{
			flow.data()[i*2] = vx.data()[i];
			flow.data()[i*2+1] = vy.data()[i];
		}
	}
	static void DissembleFlow(const FImage& flow,FImage& vx,FImage& vy)
	{
		if(!vx.matchDimension(flow.width(),flow.height(),1))
			vx.allocate(flow.width(),flow.height());
		if(!vy.matchDimension(flow.width(),flow.height(),1))
			vy.allocate(flow.width(),flow.height());
		for(int i =0;i<vx.npixels();i++)
		{
			vx.data()[i] = flow.data()[i*2];
			vy.data()[i] = flow.data()[i*2+1];
		}
	}
	static void ComputeOpticalFlow(const FImage& Im1,const FImage& Im2,FImage& flow)
	{
		if(!Im1.matchDimension(Im2))
		{
			cout<<"The input images for optical flow have different dimensions!"<<endl;
			return;
		}
		if(!flow.matchDimension(Im1.width(),Im1.height(),2))
			flow.allocate(Im1.width(),Im1.height(),2);

		float alpha=0.01;
		float ratio=0.75;
		int minWidth=30;
		int nOuterFPIterations=15;
		int nInnerFPIterations=1;
		int nCGIterations=40;

		FImage vx,vy,warpI2;
		OpticalFlow::Coarse2FineFlow(vx,vy,warpI2,Im1,Im2,alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nCGIterations);
		AssembleFlow(vx,vy,flow);
	}
};

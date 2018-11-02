#include "OpticalFlow.h"
#include "ImageProcessing.h"
#include "ImagePyramid.h"
#include <cstdlib>
#include <iostream>

#include "OpticFlowIO.h"

using namespace std;

bool OpticalFlow::IsDisplay=false;


//OpticalFlow::InterpolationMethod OpticalFlow::interpolation = OpticalFlow::Bicubic;
OpticalFlow::InterpolationMethod OpticalFlow::interpolation = OpticalFlow::Bilinear;
OpticalFlow::NoiseModel OpticalFlow::noiseModel = OpticalFlow::Lap;
GaussianMixture OpticalFlow::GMPara;
Vector<float> OpticalFlow::LapPara;


OpticalFlow::OpticalFlow(void)
{
}

OpticalFlow::~OpticalFlow(void)
{
}

//--------------------------------------------------------------------------------------------------------
//  function to compute dx, dy and dt for motion estimation
//--------------------------------------------------------------------------------------------------------
void OpticalFlow::getDxs(FImage &imdx, FImage &imdy, FImage &imdt, const FImage &im1, const FImage &im2)
{
	//float gfilter[5]={0.01,0.09,0.8,0.09,0.01};
	float gfilter[5]={0.02,0.11,0.74,0.11,0.02};
	//float gfilter[5]={0,0,1,0,0};
	if(1)
	{
		//DImage foo,Im;
		//Im.Add(im1,im2);
		//Im.Multiplywith(0.5);
		////foo.imfilter_hv(Im,gfilter,2,gfilter,2);
		//Im.dx(imdx,true);
		//Im.dy(imdy,true);
		//imdt.Subtract(im2,im1);
		FImage Im1,Im2,Im;

		im1.imfilter_hv(Im1,gfilter,2,gfilter,2);
		im2.imfilter_hv(Im2,gfilter,2,gfilter,2);
		Im.copyData(Im1);
		Im.Multiplywith(0.4);
		Im.Add(Im2,0.6);
		//Im.Multiplywith(0.5);
		//Im1.copyData(im1);
		//Im2.copyData(im2);

		Im.dx(imdx,true);
		Im.dy(imdy,true);
		imdt.Subtract(Im2,Im1);
	}
	else
	{
		// Im1 and Im2 are the smoothed version of im1 and im2
		FImage Im1,Im2;

		im1.imfilter_hv(Im1,gfilter,2,gfilter,2);
		im2.imfilter_hv(Im2,gfilter,2,gfilter,2);

		//Im1.copyData(im1);
		//Im2.copyData(im2);

		Im2.dx(imdx,true);
		Im2.dy(imdy,true);
		imdt.Subtract(Im2,Im1);
	}


	imdx.setDerivative();
	imdy.setDerivative();
	imdt.setDerivative();
}

//--------------------------------------------------------------------------------------------------------
//  Yinlin.Hu: function to compute dxx, dyy, dxy, dxt and dyt for motion estimation
//--------------------------------------------------------------------------------------------------------
void OpticalFlow::getDxxs(FImage& imdxx, FImage& imdyy, FImage& imdxy, FImage& imdxt, FImage& imdyt, const FImage& im1, const FImage& im2)
{
	float gfilter[5] = { 0.02, 0.11, 0.74, 0.11, 0.02 };

	// Im1 and Im2 are the smoothed version of im1 and im2
	FImage Im1, Im2;
	FImage dx1, dy1, dx2, dy2;

	im1.imfilter_hv(Im1, gfilter, 2, gfilter, 2);
	im2.imfilter_hv(Im2, gfilter, 2, gfilter, 2);

	//Im1.copyData(im1);
	//Im2.copyData(im2);

	Im2.dxx(imdxx);
	Im2.dyy(imdyy);
	Im2.dxy(imdxy);

	Im1.dx(dx1, true);
	Im1.dy(dy1, true);
	Im2.dx(dx2, true);
	Im2.dy(dy2, true);
	imdxt.Subtract(dx2, dx1);
	imdyt.Subtract(dy2, dy1);

	imdxx.setDerivative();
	imdyy.setDerivative();
	imdxy.setDerivative();
	imdxt.setDerivative();
	imdyt.setDerivative();
}

//--------------------------------------------------------------------------------------------------------
// function to do sanity check: imdx*du+imdy*dy+imdt=0
//--------------------------------------------------------------------------------------------------------
void OpticalFlow::SanityCheck(const FImage &imdx, const FImage &imdy, const FImage &imdt, float du, float dv)
{
	if(imdx.matchDimension(imdy)==false || imdx.matchDimension(imdt)==false)
	{
		cout<<"The dimensions of the derivatives don't match!"<<endl;
		return;
	}
	const _FlowPrecision* pImDx,*pImDy,*pImDt;
	pImDx=imdx.data();
	pImDy=imdy.data();
	pImDt=imdt.data();
	float error=0;
	for(int i=0;i<imdx.height();i++)
		for(int j=0;j<imdx.width();j++)
			for(int k=0;k<imdx.nchannels();k++)
			{
				int offset=(i*imdx.width()+j)*imdx.nchannels()+k;
				float temp=pImDx[offset]*du+pImDy[offset]*dv+pImDt[offset];
				error+=fabs(temp);
			}
	error/=imdx.nelements();
	cout<<"The mean error of |dx*u+dy*v+dt| is "<<error<<endl;
}

//--------------------------------------------------------------------------------------------------------
// function to warp image based on the flow field
//--------------------------------------------------------------------------------------------------------
void OpticalFlow::warpFL(FImage &warpIm2, const FImage &Im1, const FImage &Im2, const FImage &vx, const FImage &vy)
{
	if(warpIm2.matchDimension(Im2)==false)
		warpIm2.allocate(Im2.width(),Im2.height(),Im2.nchannels());
	ImageProcessing::warpImage(warpIm2.data(),Im1.data(),Im2.data(),vx.data(),vy.data(),Im2.width(),Im2.height(),Im2.nchannels());
}

void OpticalFlow::warpFL(FImage &warpIm2, const FImage &Im1, const FImage &Im2, const FImage &Flow)
{
	if(warpIm2.matchDimension(Im2)==false)
		warpIm2.allocate(Im2.width(),Im2.height(),Im2.nchannels());
	ImageProcessing::warpImageFlow(warpIm2.data(),Im1.data(),Im2.data(),Flow.data(),Im2.width(),Im2.height(),Im2.nchannels());
}


//--------------------------------------------------------------------------------------------------------
// function to generate mask of the pixels that move inside the image boundary
//--------------------------------------------------------------------------------------------------------
void OpticalFlow::genInImageMask(FImage &mask, const FImage &vx, const FImage &vy,int interval)
{
	int imWidth,imHeight;
	imWidth=vx.width();
	imHeight=vx.height();
	if(mask.matchDimension(vx)==false)
		mask.allocate(imWidth,imHeight);
	const _FlowPrecision *pVx,*pVy;
	_FlowPrecision *pMask;
	pVx=vx.data();
	pVy=vy.data();
	mask.reset();
	pMask=mask.data();
	float x,y;
	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
		{
			int offset=i*imWidth+j;
			y=i+pVx[offset];
			x=j+pVy[offset];
			if(x<interval  || x>imWidth-1-interval || y<interval || y>imHeight-1-interval)
				continue;
			pMask[offset]=1;
		}
}

void OpticalFlow::genInImageMask(FImage &mask, const FImage &flow,int interval)
{
	int imWidth,imHeight;
	imWidth=flow.width();
	imHeight=flow.height();
	if(mask.matchDimension(flow.width(),flow.height(),1)==false)
		mask.allocate(imWidth,imHeight);
	else
		mask.reset();

	const _FlowPrecision *pFlow;
	_FlowPrecision *pMask;
	pFlow = flow.data();;
	pMask=mask.data();
	float x,y;
	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
		{
			int offset=i*imWidth+j;
			y=i+pFlow[offset*2+1];
			x=j+pFlow[offset*2];
			if(x<interval  || x>imWidth-1-interval || y<interval || y>imHeight-1-interval)
				continue;
			pMask[offset]=1;
		}
}

//--------------------------------------------------------------------------------------------------------
// function to compute optical flow field using two fixed point iterations
// Input arguments:
//     Im1, Im2:						frame 1 and frame 2
//	warpIm2:						the warped frame 2 according to the current flow field u and v
//	u,v:									the current flow field, NOTICE that they are also output arguments
//
//--------------------------------------------------------------------------------------------------------
void OpticalFlow::SmoothFlowSOR(const FImage &Im1, const FImage &Im2, FImage &warpIm2, FImage &u, FImage &v,
																    float alpha, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations)
{
	FImage mask,imdx,imdy,imdt;
	int imWidth,imHeight,nChannels,nPixels;
	imWidth=Im1.width();
	imHeight=Im1.height();
	nChannels=Im1.nchannels();
	nPixels=imWidth*imHeight;

	FImage du(imWidth,imHeight),dv(imWidth,imHeight);
	FImage uu(imWidth,imHeight),vv(imWidth,imHeight);
	FImage uux(imWidth,imHeight),uuy(imWidth,imHeight);
	FImage vvx(imWidth,imHeight),vvy(imWidth,imHeight);
	FImage Phi_1st(imWidth,imHeight);
	FImage Psi_1st(imWidth,imHeight,nChannels);

	FImage imdxdy,imdx2,imdy2,imdtdx,imdtdy;
	FImage ImDxDy,ImDx2,ImDy2,ImDtDx,ImDtDy;
	FImage foo1,foo2;

	float prob1,prob2,prob11,prob22;

	float varepsilon_phi=pow(0.001,2);
	float varepsilon_psi=pow(0.001,2);

	//--------------------------------------------------------------------------
	// the outer fixed point iteration
	//--------------------------------------------------------------------------
	for(int count=0;count<nOuterFPIterations;count++)
	{
		// compute the gradient
		getDxs(imdx,imdy,imdt,Im1,warpIm2);

		// generate the mask to set the weight of the pixels moving outside of the image boundary to be zero
		genInImageMask(mask,u,v);

		// set the derivative of the flow field to be zero
		du.reset();
		dv.reset();

		//--------------------------------------------------------------------------
		// the inner fixed point iteration
		//--------------------------------------------------------------------------
		for(int hh=0;hh<nInnerFPIterations;hh++)
		{
			// compute the derivatives of the current flow field
			if(hh==0)
			{
				uu.copyData(u);
				vv.copyData(v);
			}
			else
			{
				uu.Add(u,du);
				vv.Add(v,dv);
			}
			uu.dx(uux);
			uu.dy(uuy);
			vv.dx(vvx);
			vv.dy(vvy);

			// compute the weight of phi
			Phi_1st.reset();
			_FlowPrecision* phiData=Phi_1st.data();
			float temp;
			const _FlowPrecision *uuxData,*uuyData,*vvxData,*vvyData;
			uuxData=uux.data();
			uuyData=uuy.data();
			vvxData=vvx.data();
			vvyData=vvy.data();
			float power_alpha = 0.5;
			for(int i=0;i<nPixels;i++)
			{
				temp=uuxData[i]*uuxData[i]+uuyData[i]*uuyData[i]+vvxData[i]*vvxData[i]+vvyData[i]*vvyData[i];
				//phiData[i]=power_alpha*pow(temp+varepsilon_phi,power_alpha-1);
				phiData[i] = 0.5/sqrt(temp+varepsilon_phi);
				//phiData[i] = 1/(power_alpha+temp);
			}

			// compute the nonlinear term of psi
			Psi_1st.reset();
			_FlowPrecision* psiData=Psi_1st.data();
			const _FlowPrecision *imdxData,*imdyData,*imdtData;
			const _FlowPrecision *duData,*dvData;
			imdxData=imdx.data();
			imdyData=imdy.data();
			imdtData=imdt.data();
			duData=du.data();
			dvData=dv.data();

			float _a  = 10000, _b = 0.1;
			if(nChannels==1)
				for(int i=0;i<nPixels;i++)
				{
					temp=imdtData[i]+imdxData[i]*duData[i]+imdyData[i]*dvData[i];
					//if(temp*temp<0.04)
					// psiData[i]=1/(2*sqrt(temp*temp+varepsilon_psi));
					//psiData[i] = _a*_b/(1+_a*temp*temp);

					// the following code is for log Gaussian mixture probability model
					temp *= temp;
					switch(noiseModel)
					{
					case GMixture:
						prob1 = GMPara.Gaussian(temp,0,0)*GMPara.alpha[0];
						prob2 = GMPara.Gaussian(temp,1,0)*(1-GMPara.alpha[0]);
						prob11 = prob1/(2*GMPara.sigma_square[0]);
						prob22 = prob2/(2*GMPara.beta_square[0]);
						psiData[i] = (prob11+prob22)/(prob1+prob2);
						break;
					case Lap:
						if(LapPara[0]<1E-20)
							continue;
						//psiData[i]=1/(2*sqrt(temp+varepsilon_psi)*LapPara[0]);
                        psiData[i]=1/(2*sqrt(temp+varepsilon_psi));
						break;
					}
				}
			else
				for(int i=0;i<nPixels;i++)
					for(int k=0;k<nChannels;k++)
					{
						int offset=i*nChannels+k;
						temp=imdtData[offset]+imdxData[offset]*duData[i]+imdyData[offset]*dvData[i];
						//if(temp*temp<0.04)
						 // psiData[offset]=1/(2*sqrt(temp*temp+varepsilon_psi));
						//psiData[offset] =  _a*_b/(1+_a*temp*temp);
						temp *= temp;
						switch(noiseModel)
						{
						case GMixture:
							prob1 = GMPara.Gaussian(temp,0,k)*GMPara.alpha[k];
							prob2 = GMPara.Gaussian(temp,1,k)*(1-GMPara.alpha[k]);
							prob11 = prob1/(2*GMPara.sigma_square[k]);
							prob22 = prob2/(2*GMPara.beta_square[k]);
							psiData[offset] = (prob11+prob22)/(prob1+prob2);
							break;
						case Lap:
							if(LapPara[k]<1E-20)
								continue;
							//psiData[offset]=1/(2*sqrt(temp+varepsilon_psi)*LapPara[k]);
                            psiData[offset]=1/(2*sqrt(temp+varepsilon_psi));
							break;
						}
					}
			// prepare the components of the large linear system
			ImDxDy.Multiply(Psi_1st,imdx,imdy);
			ImDx2.Multiply(Psi_1st,imdx,imdx);
			ImDy2.Multiply(Psi_1st,imdy,imdy);
			ImDtDx.Multiply(Psi_1st,imdx,imdt);
			ImDtDy.Multiply(Psi_1st,imdy,imdt);

			ImDxDy.collapse(imdxdy);
			ImDx2.collapse(imdx2);
			ImDy2.collapse(imdy2);
			ImDtDx.collapse(imdtdx);
			ImDtDy.collapse(imdtdy);

			// laplacian filtering of the current flow field
		    Laplacian(foo1,u,Phi_1st);
			Laplacian(foo2,v,Phi_1st);

			for(int i=0;i<nPixels;i++)
			{
				imdtdx.data()[i] = -imdtdx.data()[i] + alpha*foo1.data()[i];
				imdtdy.data()[i] = -imdtdy.data()[i] + alpha*foo2.data()[i];
			}

			// here we start SOR

			// set omega
			float omega = 1.9;

			du.reset();
			dv.reset();

			for(int k = 0; k<nSORIterations; k++)
				for(int i = 0; i<imHeight; i++)
					for(int j = 0; j<imWidth; j++)
					{
						int offset = i * imWidth+j;
						float sigma1 = 0, sigma2 = 0, coeff = 0;
                        float _weight;


						if(j>0)
						{
                            _weight = phiData[offset-1];
							sigma1  += _weight*du.data()[offset-1];
							sigma2  += _weight*dv.data()[offset-1];
							coeff   += _weight;
						}
						if(j<imWidth-1)
						{
                            _weight = phiData[offset];
							sigma1 += _weight*du.data()[offset+1];
							sigma2 += _weight*dv.data()[offset+1];
							coeff   += _weight;
						}
						if(i>0)
						{
                            _weight = phiData[offset-imWidth];
							sigma1 += _weight*du.data()[offset-imWidth];
							sigma2 += _weight*dv.data()[offset-imWidth];
							coeff   += _weight;
						}
						if(i<imHeight-1)
						{
                            _weight = phiData[offset];
							sigma1  += _weight*du.data()[offset+imWidth];
							sigma2  += _weight*dv.data()[offset+imWidth];
							coeff   += _weight;
						}
						sigma1 *= -alpha;
						sigma2 *= -alpha;
						coeff *= alpha;
						 // compute du
						sigma1 += imdxdy.data()[offset]*dv.data()[offset];
						du.data()[offset] = (1-omega)*du.data()[offset] +
							omega/(imdx2.data()[offset] /*+ alpha*0.05*/ + coeff)*(imdtdx.data()[offset] - sigma1);
						// compute dv
						sigma2 += imdxdy.data()[offset]*du.data()[offset];
						dv.data()[offset] = (1-omega)*dv.data()[offset] +
							omega/(imdy2.data()[offset] /*+ alpha*0.05*/ + coeff)*(imdtdy.data()[offset] - sigma2);
					}
		}

		u.Add(du);
		v.Add(dv);
		if(interpolation == Bilinear)
			warpFL(warpIm2,Im1,Im2,u,v);
		else
		{
			Im2.warpImageBicubicRef(Im1,warpIm2,u,v);
			warpIm2.threshold();
		}

		//Im2.warpImageBicubicRef(Im1,warpIm2,BicubicCoeff,u,v);

		// estimate noise level
		switch(noiseModel)
		{
		case GMixture:
			estGaussianMixture(Im1,warpIm2,GMPara);
			break;
		case Lap:
			estLaplacianNoise(Im1,warpIm2,LapPara);
		}
	}

}

void OpticalFlow::OneLevelFlow(const FImage &Im1, const FImage &Im2, FImage &warpIm2, FImage &u, FImage &v,
	float alpha, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations)
{
	FImage mask, imdx, imdy, imdt;
	int imWidth, imHeight, nChannels, nPixels;
	imWidth = Im1.width();
	imHeight = Im1.height();
	nChannels = Im1.nchannels();
	nPixels = imWidth*imHeight;

	FImage du(imWidth, imHeight), dv(imWidth, imHeight);
	FImage uu(imWidth, imHeight), vv(imWidth, imHeight);
	FImage uux(imWidth, imHeight), uuy(imWidth, imHeight);
	FImage vvx(imWidth, imHeight), vvy(imWidth, imHeight);
	FImage psi_d(imWidth, imHeight, nChannels);
	FImage psi_s(imWidth, imHeight);

	FImage D, Du, Dv, Au, Av;
	FImage dD, dDu, dDv, dAu, dAv;
	FImage div_u, div_v;

	float prob1, prob2, prob11, prob22;

	float varepsilon_psi_d = pow(0.001, 2);
	float varepsilon_psi_s = pow(0.001, 2);

	//--------------------------------------------------------------------------
	// the outer fixed point iteration
	//--------------------------------------------------------------------------
	for (int count = 0; count < nOuterFPIterations; count++)
	{
		// compute the gradient
		getDxs(imdx, imdy, imdt, Im1, warpIm2);

		// generate the mask to set the weight of the pixels moving outside of the image boundary to be zero
		// genInImageMask(mask, u, v);

		// set the derivative of the flow field to be zero
		du.reset();
		dv.reset();

		//--------------------------------------------------------------------------
		// the inner fixed point iteration
		//--------------------------------------------------------------------------
		for (int hh = 0; hh < nInnerFPIterations; hh++)
		{
			// compute the derivatives of the current flow field
			if (hh == 0)
			{
				uu.copyData(u);
				vv.copyData(v);
			}
			else
			{
				uu.Add(u, du);
				vv.Add(v, dv);
			}
			uu.dx(uux);
			uu.dy(uuy);
			vv.dx(vvx);
			vv.dy(vvy);

			// compute the weight of phi
			psi_s.reset();
			_FlowPrecision* psi_sData = psi_s.data();
			float temp;
			const _FlowPrecision *uuxData, *uuyData, *vvxData, *vvyData;
			uuxData = uux.data();
			uuyData = uuy.data();
			vvxData = vvx.data();
			vvyData = vvy.data();
			float power_alpha = 0.5;
			for (int i = 0; i < nPixels; i++)
			{
				temp = uuxData[i] * uuxData[i] + uuyData[i] * uuyData[i] + vvxData[i] * vvxData[i] + vvyData[i] * vvyData[i];
				//phiData[i]=power_alpha*pow(temp+varepsilon_phi,power_alpha-1);
				psi_sData[i] = 0.5 / sqrt(temp + varepsilon_psi_s);
				//phiData[i] = 1/(power_alpha+temp);
			}

			// compute the nonlinear term of psi
			psi_d.reset();
			_FlowPrecision* psi_dData = psi_d.data();
			const _FlowPrecision *imdxData, *imdyData, *imdtData;
			const _FlowPrecision *duData, *dvData;
			imdxData = imdx.data();
			imdyData = imdy.data();
			imdtData = imdt.data();
			duData = du.data();
			dvData = dv.data();

			float _a = 10000, _b = 0.1;
			if (nChannels == 1)
				for (int i = 0; i < nPixels; i++)
				{
					temp = imdtData[i] + imdxData[i] * duData[i] + imdyData[i] * dvData[i];
					//if(temp*temp<0.04)
					// psiData[i]=1/(2*sqrt(temp*temp+varepsilon_psi));
					//psiData[i] = _a*_b/(1+_a*temp*temp);

					// the following code is for log Gaussian mixture probability model
					temp *= temp;
					switch (noiseModel)
					{
					case GMixture:
						prob1 = GMPara.Gaussian(temp, 0, 0)*GMPara.alpha[0];
						prob2 = GMPara.Gaussian(temp, 1, 0)*(1 - GMPara.alpha[0]);
						prob11 = prob1 / (2 * GMPara.sigma_square[0]);
						prob22 = prob2 / (2 * GMPara.beta_square[0]);
						psi_dData[i] = (prob11 + prob22) / (prob1 + prob2);
						break;
					case Lap:
						if (LapPara[0] < 1E-20)
							continue;
						//psiData[i]=1/(2*sqrt(temp+varepsilon_psi)*LapPara[0]);
						psi_dData[i] = 1 / (2 * sqrt(temp + varepsilon_psi_d));
						break;
					}
				}
			else
				for (int i = 0; i < nPixels; i++)
					for (int k = 0; k < nChannels; k++)
					{
						int offset = i*nChannels + k;
						temp = imdtData[offset] + imdxData[offset] * duData[i] + imdyData[offset] * dvData[i];
						//if(temp*temp<0.04)
						// psiData[offset]=1/(2*sqrt(temp*temp+varepsilon_psi));
						//psiData[offset] =  _a*_b/(1+_a*temp*temp);
						temp *= temp;
						switch (noiseModel)
						{
						case GMixture:
							prob1 = GMPara.Gaussian(temp, 0, k)*GMPara.alpha[k];
							prob2 = GMPara.Gaussian(temp, 1, k)*(1 - GMPara.alpha[k]);
							prob11 = prob1 / (2 * GMPara.sigma_square[k]);
							prob22 = prob2 / (2 * GMPara.beta_square[k]);
							psi_dData[offset] = (prob11 + prob22) / (prob1 + prob2);
							break;
						case Lap:
							if (LapPara[k] < 1E-20)
								continue;
							//psiData[offset]=1/(2*sqrt(temp+varepsilon_psi)*LapPara[k]);
							psi_dData[offset] = 1 / (2 * sqrt(temp + varepsilon_psi_d));
							break;
						}
					}
			// prepare the components of the large linear system
			dD.Multiply(psi_d, imdx, imdy);
			dDu.Multiply(psi_d, imdx, imdx);
			dDv.Multiply(psi_d, imdy, imdy);
			dAu.Multiply(psi_d, imdx, imdt);
			dAv.Multiply(psi_d, imdy, imdt);

			dD.collapse(D);
			dDu.collapse(Du);
			dDv.collapse(Dv);
			dAu.collapse(Au);
			dAv.collapse(Av);

			// laplacian filtering of the current flow field
			Laplacian(div_u, u, psi_s);
			Laplacian(div_v, v, psi_s);

			for (int i = 0; i < nPixels; i++)
			{
				Au.data()[i] = -Au.data()[i] + alpha*div_u.data()[i];
				Av.data()[i] = -Av.data()[i] + alpha*div_v.data()[i];
			}

			// here we start SOR

			// set omega
			float omega = 1.9;

			du.reset();
			dv.reset();

			for (int k = 0; k < nSORIterations; k++)
				for (int i = 0; i < imHeight; i++)
					for (int j = 0; j < imWidth; j++)
					{
						int offset = i*imWidth+j;

						float div_du = 0, div_dv = 0, div_d = 0;
						float _weight;

						if (j > 0){
							_weight = psi_sData[offset - 1];
							div_du += _weight*du.data()[offset - 1];
							div_dv += _weight*dv.data()[offset - 1];
							div_d += _weight;
						}
						if (j < imWidth - 1){
							_weight = psi_sData[offset];
							div_du += _weight*du.data()[offset + 1];
							div_dv += _weight*dv.data()[offset + 1];
							div_d += _weight;
						}
						if (i > 0){
							_weight = psi_sData[offset - imWidth];
							div_du += _weight*du.data()[offset - imWidth];
							div_dv += _weight*dv.data()[offset - imWidth];
							div_d += _weight;
						}
						if (i < imHeight - 1){
							_weight = psi_sData[offset];
							div_du += _weight*du.data()[offset + imWidth];
							div_dv += _weight*dv.data()[offset + imWidth];
							div_d += _weight;
						}

						// compute du
						du.data()[offset] =
							(1 - omega)*du.data()[offset] + omega *
							(Au.data()[offset] - D.data()[offset]*dv.data()[offset] + alpha*div_du) /
							(Du.data()[offset] + alpha*div_d);
						// compute dv
						dv.data()[offset] =
							(1 - omega)*dv.data()[offset] + omega *
							(Av.data()[offset] - D.data()[offset]*du.data()[offset] + alpha*div_dv) /
							(Dv.data()[offset] + alpha*div_d);
					}
		}

		u.Add(du);
		v.Add(dv);

		//
// 		u.MedianFiltering(1);
// 		v.MedianFiltering(1);

		if (interpolation == Bilinear)
			warpFL(warpIm2, Im1, Im2, u, v);
		else
		{
			Im2.warpImageBicubicRef(Im1, warpIm2, u, v);
			warpIm2.threshold();
		}

		//Im2.warpImageBicubicRef(Im1,warpIm2,BicubicCoeff,u,v);

		// estimate noise level
		switch (noiseModel)
		{
		case GMixture:
			estGaussianMixture(Im1, warpIm2, GMPara);
			break;
		case Lap:
			estLaplacianNoise(Im1, warpIm2, LapPara);
		}
	}

}


void OpticalFlow::OneLevelFlow_Brox(const FImage &Im1, const FImage &Im2, FImage &warpIm2, FImage &u, FImage &v,
	float alpha, float gamma, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations)
{
	FImage mask, imdx, imdy, imdt;
	FImage imdxx, imdyy, imdxy, imdxt, imdyt;

	int imWidth, imHeight, nChannels, nElements, nPixels;
	imWidth = Im1.width();
	imHeight = Im1.height();
	nChannels = Im1.nchannels();
	nElements = Im1.nelements();
	nPixels = imWidth*imHeight;

	FImage du(imWidth, imHeight), dv(imWidth, imHeight);
	FImage uu(imWidth, imHeight), vv(imWidth, imHeight);
	FImage uux(imWidth, imHeight), uuy(imWidth, imHeight);
	FImage vvx(imWidth, imHeight), vvy(imWidth, imHeight);
	FImage psi_d(imWidth, imHeight, nChannels);
	FImage psi_g(imWidth, imHeight, nChannels);
	FImage psi_s(imWidth, imHeight);

	FImage Au(imWidth, imHeight);
	FImage Av(imWidth, imHeight);
	FImage Du(imWidth, imHeight);
	FImage Dv(imWidth, imHeight);
	FImage D(imWidth, imHeight);

	FImage dAu, dAv, dDu, dDv, dD;
	FImage gAu, gAv, gDu, gDv, gD;
	FImage div_u, div_v;

	float prob1, prob2, prob11, prob22;

	float varepsilon_psi_d = pow(0.001, 2);
	float varepsilon_psi_g = pow(0.001, 2);
	float varepsilon_psi_s = pow(0.001, 2);

	//--------------------------------------------------------------------------
	// the outer fixed point iteration
	//--------------------------------------------------------------------------
	for (int count = 0; count < nOuterFPIterations; count++)
	{
		// compute the gradient
		getDxs(imdx, imdy, imdt, Im1, warpIm2);
		getDxxs(imdxx, imdyy, imdxy, imdxt, imdyt, Im1, warpIm2);

		// generate the mask to set the weight of the pixels moving outside of the image boundary to be zero
		// genInImageMask(mask, u, v);

		// set the derivative of the flow field to be zero
		du.reset();
		dv.reset();

		//--------------------------------------------------------------------------
		// the inner fixed point iteration
		//--------------------------------------------------------------------------
		for (int hh = 0; hh < nInnerFPIterations; hh++)
		{
			// compute the derivatives of the current flow field
			if (hh == 0)
			{
				uu.copyData(u);
				vv.copyData(v);
			}
			else
			{
				uu.Add(u, du);
				vv.Add(v, dv);
			}
			uu.dx(uux);
			uu.dy(uuy);
			vv.dx(vvx);
			vv.dy(vvy);

			// compute the weight of smoothness
			psi_s.reset();
			_FlowPrecision* psi_sData = psi_s.data();
			float temp;
			const _FlowPrecision *uuxData, *uuyData, *vvxData, *vvyData;
			uuxData = uux.data();
			uuyData = uuy.data();
			vvxData = vvx.data();
			vvyData = vvy.data();
			float power_alpha = 0.5;
			for (int i = 0; i < nPixels; i++)
			{
				temp = uuxData[i] * uuxData[i] + uuyData[i] * uuyData[i] + vvxData[i] * vvxData[i] + vvyData[i] * vvyData[i];
				//phiData[i]=power_alpha*pow(temp+varepsilon_phi,power_alpha-1);
				psi_sData[i] = 0.5 / sqrt(temp + varepsilon_psi_s);
				//phiData[i] = 1/(power_alpha+temp);
			}

			// compute the weight of gradient
			psi_g.reset();
			_FlowPrecision* psi_gData = psi_g.data();
			const _FlowPrecision *imdxxData, *imdyyData, *imdxyData, *imdxtData, *imdytData;
			const _FlowPrecision *duData, *dvData;
			imdxxData = imdxx.data();
			imdyyData = imdyy.data();
			imdxyData = imdxy.data();
			imdxtData = imdxt.data();
			imdytData = imdyt.data();
			duData = du.data();
			dvData = dv.data();
			float tempu, tempv;
			for (int i = 0; i < nPixels; i++){
				for (int k = 0; k < nChannels; k++){
					int offset = i*nChannels + k;
					tempu = imdxtData[offset] + imdxxData[offset] * duData[i] + imdxyData[offset] * dvData[i];
					tempv = imdytData[offset] + imdxyData[offset] * duData[i] + imdyyData[offset] * dvData[i];
					tempu *= tempu;
					tempv *= tempv;
					psi_gData[offset] = 1 / (2 * sqrt(tempu + tempv + varepsilon_psi_g));
				}
			}

			// compute the weight of data
			psi_d.reset();
			_FlowPrecision* psi_dData = psi_d.data();
			const _FlowPrecision *imdxData, *imdyData, *imdtData;
			imdxData = imdx.data();
			imdyData = imdy.data();
			imdtData = imdt.data();

			float _a = 10000, _b = 0.1;
			if (nChannels == 1)
				for (int i = 0; i < nPixels; i++)
				{
					temp = imdtData[i] + imdxData[i] * duData[i] + imdyData[i] * dvData[i];
					//if(temp*temp<0.04)
					// psiData[i]=1/(2*sqrt(temp*temp+varepsilon_psi));
					//psiData[i] = _a*_b/(1+_a*temp*temp);

					// the following code is for log Gaussian mixture probability model
					temp *= temp;
					switch (noiseModel)
					{
					case GMixture:
						prob1 = GMPara.Gaussian(temp, 0, 0)*GMPara.alpha[0];
						prob2 = GMPara.Gaussian(temp, 1, 0)*(1 - GMPara.alpha[0]);
						prob11 = prob1 / (2 * GMPara.sigma_square[0]);
						prob22 = prob2 / (2 * GMPara.beta_square[0]);
						psi_dData[i] = (prob11 + prob22) / (prob1 + prob2);
						break;
					case Lap:
						if (LapPara[0] < 1E-20)
							continue;
						//psiData[i]=1/(2*sqrt(temp+varepsilon_psi)*LapPara[0]);
						psi_dData[i] = 1 / (2 * sqrt(temp + varepsilon_psi_d));
						break;
					}
				}
			else
				for (int i = 0; i < nPixels; i++)
					for (int k = 0; k < nChannels; k++)
					{
						int offset = i*nChannels + k;
						temp = imdtData[offset] + imdxData[offset] * duData[i] + imdyData[offset] * dvData[i];
						//if(temp*temp<0.04)
						// psiData[offset]=1/(2*sqrt(temp*temp+varepsilon_psi));
						//psiData[offset] =  _a*_b/(1+_a*temp*temp);
						temp *= temp;
						switch (noiseModel)
						{
						case GMixture:
							prob1 = GMPara.Gaussian(temp, 0, k)*GMPara.alpha[k];
							prob2 = GMPara.Gaussian(temp, 1, k)*(1 - GMPara.alpha[k]);
							prob11 = prob1 / (2 * GMPara.sigma_square[k]);
							prob22 = prob2 / (2 * GMPara.beta_square[k]);
							psi_dData[offset] = (prob11 + prob22) / (prob1 + prob2);
							break;
						case Lap:
							if (LapPara[k] < 1E-20)
								continue;
							//psiData[offset]=1/(2*sqrt(temp+varepsilon_psi)*LapPara[k]);
							psi_dData[offset] = 1 / (2 * sqrt(temp + varepsilon_psi_d));
							break;
						}
					}
			// prepare the components of the large linear system

			gAu.allocate(Im1);
			gAv.allocate(Im1);
			gDu.allocate(Im1);
			gDv.allocate(Im1);
			gD.allocate(Im1);
			for (int i = 0; i < nElements; i++){
				gAu.pData[i] = psi_gData[i] * (imdxtData[i] * imdxxData[i] + imdytData[i] * imdxyData[i]);
				gAv.pData[i] = psi_gData[i] * (imdxtData[i] * imdxyData[i] + imdytData[i] * imdyyData[i]);
				gDu.pData[i] = psi_gData[i] * (imdxxData[i] * imdxxData[i] + imdxyData[i] * imdxyData[i]);
				gDv.pData[i] = psi_gData[i] * (imdyyData[i] * imdyyData[i] + imdxyData[i] * imdxyData[i]);
				gD.pData[i] = psi_gData[i] * (imdxxData[i] + imdyyData[i]) * imdxyData[i];
			}

			dD.Multiply(psi_d, imdx, imdy);
			dDu.Multiply(psi_d, imdx, imdx);
			dDv.Multiply(psi_d, imdy, imdy);
			dAu.Multiply(psi_d, imdx, imdt);
			dAv.Multiply(psi_d, imdy, imdt);

			dD.collapse();
			dDu.collapse();
			dDv.collapse();
			dAu.collapse();
			dAv.collapse();

			gD.collapse();
			gDu.collapse();
			gDv.collapse();
			gAu.collapse();
			gAv.collapse();

			// laplacian filtering of the current flow field
			Laplacian(div_u, u, psi_s);
			Laplacian(div_v, v, psi_s);

			for (int i = 0; i < nPixels; i++){
				Au.pData[i] = -dAu.pData[i] + alpha*div_u.pData[i] - gamma*gAu.pData[i];
				Av.pData[i] = -dAv.pData[i] + alpha*div_v.pData[i] - gamma*gAv.pData[i];
				Du.pData[i] = dDu.pData[i] + gamma*gDu.pData[i];
				Dv.pData[i] = dDv.pData[i] + gamma*gDv.pData[i];
				D.pData[i] = dD.pData[i] + gamma*gD.pData[i];
			}

			// here we start SOR

			// set omega
			float omega = 1.9;

			du.reset();
			dv.reset();

			for (int k = 0; k < nSORIterations; k++)
				for (int i = 0; i < imHeight; i++)
					for (int j = 0; j < imWidth; j++)
					{
						int offset = i*imWidth + j;

						float div_du = 0, div_dv = 0, div_d = 0;
						float _weight;

						if (j > 0){
							_weight = psi_sData[offset - 1];
							div_du += _weight*du.data()[offset - 1];
							div_dv += _weight*dv.data()[offset - 1];
							div_d += _weight;
						}
						if (j < imWidth - 1){
							_weight = psi_sData[offset];
							div_du += _weight*du.data()[offset + 1];
							div_dv += _weight*dv.data()[offset + 1];
							div_d += _weight;
						}
						if (i > 0){
							_weight = psi_sData[offset - imWidth];
							div_du += _weight*du.data()[offset - imWidth];
							div_dv += _weight*dv.data()[offset - imWidth];
							div_d += _weight;
						}
						if (i < imHeight - 1){
							_weight = psi_sData[offset];
							div_du += _weight*du.data()[offset + imWidth];
							div_dv += _weight*dv.data()[offset + imWidth];
							div_d += _weight;
						}

						// compute du
						du.data()[offset] =
							(1 - omega)*du.data()[offset] + omega *
							(Au.data()[offset] - D.data()[offset] * dv.data()[offset] + alpha*div_du) /
							(Du.data()[offset] + alpha*div_d);
						// compute dv
						dv.data()[offset] =
							(1 - omega)*dv.data()[offset] + omega *
							(Av.data()[offset] - D.data()[offset] * du.data()[offset] + alpha*div_dv) /
							(Dv.data()[offset] + alpha*div_d);
					}
		}

		u.Add(du);
		v.Add(dv);

		//
// 		u.MedianFiltering(1);
// 		v.MedianFiltering(1);

		if (interpolation == Bilinear)
			warpFL(warpIm2, Im1, Im2, u, v);
		else
		{
			Im2.warpImageBicubicRef(Im1, warpIm2, u, v);
			warpIm2.threshold();
		}

		//Im2.warpImageBicubicRef(Im1,warpIm2,BicubicCoeff,u,v);

		// estimate noise level
		switch (noiseModel)
		{
		case GMixture:
			estGaussianMixture(Im1, warpIm2, GMPara);
			break;
		case Lap:
			estLaplacianNoise(Im1, warpIm2, LapPara);
		}
	}

}

//--------------------------------------------------------------------------------------------------------
// function to compute optical flow field using two fixed point iterations
// Input arguments:
//     Im1, Im2:						frame 1 and frame 2
//	warpIm2:						the warped frame 2 according to the current flow field u and v
//	u,v:									the current flow field, NOTICE that they are also output arguments
//
//--------------------------------------------------------------------------------------------------------
void OpticalFlow::SmoothFlowPDE(const FImage &Im1, const FImage &Im2, FImage &warpIm2, FImage &u, FImage &v,
																    float alpha, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations)
{
	FImage mask, imdx, imdy, imdt;
	int imWidth, imHeight, nChannels, nPixels;
	imWidth = Im1.width();
	imHeight = Im1.height();
	nChannels = Im1.nchannels();
	nPixels = imWidth*imHeight;

	FImage du(imWidth, imHeight), dv(imWidth, imHeight);
	FImage uu(imWidth, imHeight), vv(imWidth, imHeight);
	FImage ux(imWidth, imHeight), uy(imWidth, imHeight);
	FImage vx(imWidth, imHeight), vy(imWidth, imHeight);
	FImage Phi_1st(imWidth, imHeight);
	FImage Psi_1st(imWidth, imHeight, nChannels);

	FImage imdxy, imdx2, imdy2, imdtdx, imdtdy;
	FImage ImDxy, ImDx2, ImDy2, ImDtDx, ImDtDy;
	FImage A11, A12, A22, b1, b2;
	FImage foo1, foo2;

	// variables for conjugate gradient
	FImage r1, r2, p1, p2, q1, q2;
	float* rou;
	rou = new float[nCGIterations];

	float varepsilon_phi = pow(0.001, 2);
	float varepsilon_psi = pow(0.001, 2);

	//--------------------------------------------------------------------------
	// the outer fixed point iteration
	//--------------------------------------------------------------------------
	for (int count = 0; count < nOuterFPIterations; count++)
	{
		// compute the gradient
		getDxs(imdx, imdy, imdt, Im1, warpIm2);

		// generate the mask to set the weight of the pxiels moving outside of the image boundary to be zero
		genInImageMask(mask, vx, vy);

		// set the derivative of the flow field to be zero
		du.reset();
		dv.reset();

		//--------------------------------------------------------------------------
		// the inner fixed point iteration
		//--------------------------------------------------------------------------
		for (int hh = 0; hh < nInnerFPIterations; hh++)
		{
			// compute the derivatives of the current flow field
			if (hh == 0)
			{
				uu.copyData(u);
				vv.copyData(v);
			}
			else
			{
				uu.Add(u, du);
				vv.Add(v, dv);
			}
			uu.dx(ux);
			uu.dy(uy);
			vv.dx(vx);
			vv.dy(vy);

			// compute the weight of phi
			Phi_1st.reset();
			float* phiData = Phi_1st.data();
			float temp;
			const float *uxData, *uyData, *vxData, *vyData;
			uxData = ux.data();
			uyData = uy.data();
			vxData = vx.data();
			vyData = vy.data();
			for (int i = 0; i < nPixels; i++)
			{
				temp = uxData[i] * uxData[i] + uyData[i] * uyData[i] + vxData[i] * vxData[i] + vyData[i] * vyData[i];
				phiData[i] = 1 / (2 * sqrt(temp + varepsilon_phi));
			}

			// compute the nonlinear term of psi
			float* psiData = Psi_1st.data();
			const float *imdxData, *imdyData, *imdtData;
			const float *duData, *dvData;
			imdxData = imdx.data();
			imdyData = imdy.data();
			imdtData = imdt.data();
			duData = du.data();
			dvData = dv.data();

			if (nChannels == 1)
			{
				for (int i = 0; i < nPixels; i++)
				{
					temp = imdtData[i] + imdxData[i] * duData[i] + imdyData[i] * dvData[i];
					psiData[i] = 1 / (2 * sqrt(temp*temp + varepsilon_psi));
				}
			}
			else
			{
				for (int i = 0; i < nPixels; i++)
					for (int k = 0; k < nChannels; k++)
					{
						int offset = i*nChannels + k;
						temp = imdtData[offset] + imdxData[offset] * duData[i] + imdyData[offset] * dvData[i];
						psiData[offset] = 1 / (2 * sqrt(temp*temp + varepsilon_psi));
					}
			}

			// prepare the components of the large linear system
			ImDxy.Multiply(Psi_1st, imdx, imdy);
			ImDx2.Multiply(Psi_1st, imdx, imdx);
			ImDy2.Multiply(Psi_1st, imdy, imdy);
			ImDtDx.Multiply(Psi_1st, imdx, imdt);
			ImDtDy.Multiply(Psi_1st, imdy, imdt);

			ImDxy.collapse(imdxy);
			ImDx2.collapse(imdx2);
			ImDy2.collapse(imdy2);
			ImDtDx.collapse(imdtdx);
			ImDtDy.collapse(imdtdy);

			// filtering
			imdx2.smoothing(A11, 3);
			imdxy.smoothing(A12, 3);
			imdy2.smoothing(A22, 3);

			// add epsilon to A11 and A22
			A11.Add(alpha*0.1);
			A22.Add(alpha*0.1);

			// form b
			imdtdx.smoothing(b1, 3);
			imdtdy.smoothing(b2, 3);
			// laplacian filtering of the current flow field
			Laplacian(foo1, u, Phi_1st);
			Laplacian(foo2, v, Phi_1st);
			float *b1Data, *b2Data;
			const float *foo1Data, *foo2Data;
			b1Data = b1.data();
			b2Data = b2.data();
			foo1Data = foo1.data();
			foo2Data = foo2.data();

			for (int i = 0; i < nPixels; i++)
			{
				b1Data[i] = -b1Data[i] + alpha*foo1Data[i];
				b2Data[i] = -b2Data[i] + alpha*foo2Data[i];
			}

			// for debug only, displaying the matrix coefficients
			//A11.imwrite("A11.bmp",ImageIO::normalized);
			//A12.imwrite("A12.bmp",ImageIO::normalized);
			//A22.imwrite("A22.bmp",ImageIO::normalized);
			//b1.imwrite("b1.bmp",ImageIO::normalized);
			//b2.imwrite("b2.bmp",ImageIO::normalized);

			//-----------------------------------------------------------------------
			// conjugate gradient algorithm
			//-----------------------------------------------------------------------
			r1.copyData(b1);
			r2.copyData(b2);
			du.reset();
			dv.reset();

			for (int k = 0; k < nCGIterations; k++)
			{
				rou[k] = r1.norm2() + r2.norm2();
				//cout<<rou[k]<<endl;
				if (rou[k] < 1E-10)
					break;
				if (k == 0)
				{
					p1.copyData(r1);
					p2.copyData(r2);
				}
				else
				{
					float ratio = rou[k] / rou[k - 1];
					p1.Add(r1, p1, ratio);
					p2.Add(r2, p2, ratio);
				}
				// go through the large linear system
				foo1.Multiply(A11, p1);
				foo2.Multiply(A12, p2);
				q1.Add(foo1, foo2);
				Laplacian(foo1, p1, Phi_1st);
				q1.Add(foo1, alpha);

				foo1.Multiply(A12, p1);
				foo2.Multiply(A22, p2);
				q2.Add(foo1, foo2);
				Laplacian(foo2, p2, Phi_1st);
				q2.Add(foo2, alpha);

				float beta;
				beta = rou[k] / (p1.innerproduct(q1) + p2.innerproduct(q2));

				du.Add(p1, beta);
				dv.Add(p2, beta);

				r1.Add(q1, -beta);
				r2.Add(q2, -beta);
			}
			//-----------------------------------------------------------------------
			// end of conjugate gradient algorithm
			//-----------------------------------------------------------------------
		}// end of inner fixed point iteration

		// the following procedure is merely for debugging
		//cout<<"du "<<du.norm2()<<" dv "<<dv.norm2()<<endl;
		// update the flow field
		u.Add(du, 1);
		v.Add(dv, 1);
		warpFL(warpIm2, Im1, Im2, u, v);
	}// end of outer fixed point iteration
	delete rou;
}

void OpticalFlow::estGaussianMixture(const FImage& Im1,const FImage& Im2,GaussianMixture& para,float prior)
{
	int nIterations = 3, nChannels = Im1.nchannels();
	FImage weight1(Im1),weight2(Im1);
	float *total1,*total2;
	total1 = new float[nChannels];
	total2 = new float[nChannels];
	for(int count = 0; count<nIterations; count++)
	{
		float temp;
		memset(total1,0,sizeof(float)*nChannels);
		memset(total2,0,sizeof(float)*nChannels);

		// E step
		for(int i = 0;i<weight1.npixels();i++)
			for(int k=0;k<nChannels;k++)
			{
				int offset = i*weight1.nchannels()+k;
				temp = Im1[offset]-Im2[offset];
				temp *= temp;
				weight1[offset] = para.Gaussian(temp,0,k)*para.alpha[k];
				weight2[offset] = para.Gaussian(temp,1,k)*(1-para.alpha[k]);
				temp = weight1[offset]+weight2[offset];
				weight1[offset]/=temp;
				weight2[offset]/=temp;
				total1[k] += weight1[offset];
				total2[k] += weight2[offset];
			}

		// M step
		para.reset();


		for(int i = 0;i<weight1.npixels();i++)
			for(int k =0;k<nChannels;k++)
			{
				int offset = i*weight1.nchannels()+k;
				temp = Im1[offset]-Im2[offset];
				temp *= temp;
				para.sigma[k]+= weight1[offset]*temp;
				para.beta[k] += weight2[offset]*temp;
			}

		for(int k =0;k<nChannels;k++)
		{
			para.alpha[k] = total1[k]/(total1[k]+total2[k])*(1-prior)+0.95*prior; // regularize alpha
			para.sigma[k] = sqrt(para.sigma[k]/total1[k]);
			para.beta[k]   = sqrt(para.beta[k]/total2[k])*(1-prior)+0.3*prior; // regularize beta
		}
		para.square();
		count = count;
	}
}

void OpticalFlow::estLaplacianNoise(const FImage& Im1,const FImage& Im2,Vector<float>& para)
{
	int nChannels = Im1.nchannels();
	if(para.dim()!=nChannels)
		para.allocate(nChannels);
	else
		para.reset();
	float temp;
	Vector<float> total(nChannels);
	for(int k = 0;k<nChannels;k++)
		total[k] = 0;

	for(int i =0;i<Im1.npixels();i++)
		for(int k = 0;k<nChannels;k++)
		{
			int offset = i*nChannels+k;
			temp= abs(Im1.data()[offset]-Im2.data()[offset]);
			if(temp>0 && temp<1000000)
			{
				para[k] += temp;
				total[k]++;
			}
		}
	for(int k = 0;k<nChannels;k++)
	{
		if(total[k]==0)
		{
			cout<<"All the pixels are invalid in estimation Laplacian noise!!!"<<endl;
			cout<<"Something severely wrong happened!!!"<<endl;
			para[k] = 0.001;
		}
		else
			para[k]/=total[k];
	}
}

void OpticalFlow::Laplacian(FImage &output, const FImage &input, const FImage& weight)
{
	if(output.matchDimension(input)==false)
		output.allocate(input);
	output.reset();

	if(input.matchDimension(weight)==false)
	{
		cout<<"Error in image dimension matching OpticalFlow::Laplacian()!"<<endl;
		return;
	}

	const _FlowPrecision *inputData=input.data(),*weightData=weight.data();
	int width=input.width(),height=input.height();
	FImage foo(width,height);
	_FlowPrecision *fooData=foo.data(),*outputData=output.data();


	// horizontal filtering
	for(int i=0;i<height;i++)
		for(int j=0;j<width-1;j++)
		{
			int offset=i*width+j;
			fooData[offset]=(inputData[offset+1]-inputData[offset])*weightData[offset];
		}
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)
		{
			int offset=i*width+j;
			if(j<width-1)
				outputData[offset]+=fooData[offset];
			if(j>0)
				outputData[offset]-=fooData[offset-1];
		}
	foo.reset();
	// vertical filtering
	for(int i=0;i<height-1;i++)
		for(int j=0;j<width;j++)
		{
			int offset=i*width+j;
			fooData[offset]=(inputData[offset+width]-inputData[offset])*weightData[offset];
		}
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)
		{
			int offset=i*width+j;
			if(i<height-1)
				outputData[offset]+=fooData[offset];
			if(i>0)
				outputData[offset]-=fooData[offset-width];
		}
}

void OpticalFlow::testLaplacian(int dim)
{
	// generate the random weight
	FImage weight(dim,dim);
	for(int i=0;i<dim;i++)
		for(int j=0;j<dim;j++)
			//weight.data()[i*dim+j]=(float)rand()/RAND_MAX+1;
			weight.data()[i*dim+j]=1;
	// go through the linear system;
	FImage sysMatrix(dim*dim,dim*dim);
	FImage u(dim,dim),du(dim,dim);
	for(int i=0;i<dim*dim;i++)
	{
		u.reset();
		u.data()[i]=1;
		Laplacian(du,u,weight);
		for(int j=0;j<dim*dim;j++)
			sysMatrix.data()[j*dim*dim+i]=du.data()[j];
	}
	// test whether the matrix is symmetric
	for(int i=0;i<dim*dim;i++)
	{
		for(int j=0;j<dim*dim;j++)
		{
			if(sysMatrix.data()[i*dim*dim+j]>=0)
				printf(" ");
			printf(" %1.0f ",sysMatrix.data()[i*dim*dim+j]);
		}
		printf("\n");
	}
}

//--------------------------------------------------------------------------------------
// function to perform coarse to fine optical flow estimation
//--------------------------------------------------------------------------------------
void OpticalFlow::Coarse2FineFlow(FImage &vx, FImage &vy, FImage &warpI2,const FImage &Im1, const FImage &Im2, float alpha, float ratio, int minWidth,
																	 int nOuterFPIterations, int nInnerFPIterations, int nCGIterations)
{
	// first build the pyramid of the two images
	FImagePyramid GPyramid1;
	FImagePyramid GPyramid2;
	if(IsDisplay)
		cout<<"Constructing pyramid...";
	GPyramid1.ConstructPyramid(Im1,ratio,minWidth);
	GPyramid2.ConstructPyramid(Im2,ratio,minWidth);
	if(IsDisplay)
		cout<<"done!"<<endl;

	// now iterate from the top level to the bottom
	FImage Image1,Image2,WarpImage2;
	//GaussianMixture GMPara(Im1.nchannels()+2);

	// initialize noise
	switch(noiseModel){
	case GMixture:
		GMPara.reset(Im1.nchannels()+2);
		break;
	case Lap:
		LapPara.allocate(Im1.nchannels()+2);
		for(int i = 0;i<LapPara.dim();i++)
			LapPara[i] = 0.02;
		break;
	}

	for(int k=GPyramid1.nlevels()-1;k>=0;k--)
	{
		if(IsDisplay)
			cout<<"Pyramid level "<<k;
		int width=GPyramid1[k].width();
		int height=GPyramid1[k].height();
// 		im2feature(Image1,GPyramid1.Image(k));
// 		im2feature(Image2,GPyramid2.Image(k));
// 		Image1.copy(GPyramid1.Image(k));
// 		Image2.copy(GPyramid2.Image(k));
		GPyramid1[k].desaturate(Image1);
		GPyramid2[k].desaturate(Image2);

		if (GPyramid1.nlevels() == 1){ // only one level, the vx,vy are initialized
			if (interpolation == Bilinear)
				warpFL(WarpImage2, Image1, Image2, vx, vy);
			else
				Image2.warpImageBicubicRef(Image1, WarpImage2, vx, vy);
		}else{ // multi level
			if (k == GPyramid1.nlevels() - 1) // if at the top level
			{
				vx.allocate(width, height);
				vy.allocate(width, height);
				//warpI2.copyData(Image2);
				WarpImage2.copyData(Image2);
			}
			else
			{
				vx.imresize(width, height);
				vx.Multiplywith(1 / ratio);
				vy.imresize(width, height);
				vy.Multiplywith(1 / ratio);
				//warpFL(warpI2,GPyramid1.Image(k),GPyramid2.Image(k),vx,vy);
				if (interpolation == Bilinear)
					warpFL(WarpImage2, Image1, Image2, vx, vy);
				else
					Image2.warpImageBicubicRef(Image1, WarpImage2, vx, vy);
			}
		}
// 		if (k == 0){
// 			OpticFlowIO::ShowFlow("k-1", vx.pData, vy.pData, width, height);
// 		}
		//SmoothFlowPDE(GPyramid1.Image(k),GPyramid2.Image(k),warpI2,vx,vy,alpha,nOuterFPIterations,nInnerFPIterations,nCGIterations);
		//SmoothFlowPDE(Image1,Image2,WarpImage2,vx,vy,alpha*pow((1/ratio),k),nOuterFPIterations,nInnerFPIterations,nCGIterations,GMPara);
		//SmoothFlowPDE(Image1,Image2,WarpImage2,vx,vy,alpha,nOuterFPIterations,nInnerFPIterations,nCGIterations);
		
		//SmoothFlowSOR(Image1, Image2, WarpImage2, vx, vy, alpha, nOuterFPIterations + k, nInnerFPIterations, nCGIterations + k * 3);
		//OneLevelFlow(Image1, Image2, WarpImage2, vx, vy, alpha, nOuterFPIterations + k, nInnerFPIterations, nCGIterations + k * 3);
		OneLevelFlow_Brox(Image1, Image2, WarpImage2, vx, vy, alpha, 1.5, nOuterFPIterations + k, nInnerFPIterations, nCGIterations + k * 3);
		
// 		if (k == 0){
// 			OpticFlowIO::ShowFlow("k", vx.pData, vy.pData, width, height, true);
// 		}

		//GMPara.display();
		if(IsDisplay)
			cout<<endl;
	}
	//warpFL(warpI2,Im1,Im2,vx,vy);
	Im2.warpImageBicubicRef(Im1,warpI2,vx,vy);
	warpI2.threshold();
}

void OpticalFlow::Coarse2FineFlowLevel(FImage &vx, FImage &vy, FImage &warpI2,const FImage &Im1, const FImage &Im2, float alpha, float ratio, int nLevels,
																	 int nOuterFPIterations, int nInnerFPIterations, int nCGIterations)
{
	// first build the pyramid of the two images
	FImagePyramid GPyramid1;
	FImagePyramid GPyramid2;
	FImagePyramid GFlow;
	FImage flow;
	AssembleFlow(vx,vy,flow);
	if(IsDisplay)
		cout<<"Constructing pyramid...";
	GPyramid1.ConstructPyramidLevels(Im1,ratio,nLevels);
	GPyramid2.ConstructPyramidLevels(Im2,ratio,nLevels);
	GFlow.ConstructPyramidLevels(flow,ratio,nLevels);
	flow= GFlow[nLevels-1];
	flow.Multiplywith(pow(ratio,nLevels-1));
	DissembleFlow(flow,vx,vy);

	if(IsDisplay)
		cout<<"done!"<<endl;

	// now iterate from the top level to the bottom
	FImage Image1,Image2,WarpImage2;

	// initialize noise
	switch(noiseModel){
	case GMixture:
		GMPara.reset(Im1.nchannels()+2);
		break;
	case Lap:
		LapPara.allocate(Im1.nchannels()+2);
		for(int i = 0;i<LapPara.dim();i++)
			LapPara[i] = 0.02;
		break;
	}


	for(int k=GPyramid1.nlevels()-1;k>=0;k--)
	{
		if(IsDisplay)
			cout<<"Pyramid level "<<k;
		int width=GPyramid1[k].width();
		int height=GPyramid1[k].height();
		im2feature(Image1,GPyramid1[k]);
		im2feature(Image2,GPyramid2[k]);

		if(k<GPyramid1.nlevels()-1) // if at the top level
		{
			vx.imresize(width,height);
			vx.Multiplywith(1/ratio);
			vy.imresize(width,height);
			vy.Multiplywith(1/ratio);
		}
		if(interpolation == Bilinear)
			warpFL(WarpImage2,Image1,Image2,vx,vy);
		else
			Image2.warpImageBicubicRef(Image1,WarpImage2,vx,vy);
		//SmoothFlowPDE(GPyramid1.Image(k),GPyramid2.Image(k),warpI2,vx,vy,alpha,nOuterFPIterations,nInnerFPIterations,nCGIterations);
		//SmoothFlowPDE(Image1,Image2,WarpImage2,vx,vy,alpha*pow((1/ratio),k),nOuterFPIterations,nInnerFPIterations,nCGIterations,GMPara);

		SmoothFlowPDE(Image1,Image2,WarpImage2,vx,vy,alpha,nOuterFPIterations,nInnerFPIterations,nCGIterations);
		//GMPara.display();
		if(IsDisplay)
			cout<<endl;
	}
	//warpFL(warpI2,Im1,Im2,vx,vy);
	Im2.warpImageBicubicRef(Im1,warpI2,vx,vy);
	warpI2.threshold();
}

//---------------------------------------------------------------------------------------
// function to convert image to feature image
//---------------------------------------------------------------------------------------
void OpticalFlow::im2feature(FImage &imfeature, const FImage &im)
{
	int width=im.width();
	int height=im.height();
	int nchannels=im.nchannels();
	if(nchannels==1)
	{
		imfeature.allocate(im.width(),im.height(),3);
		FImage imdx,imdy;
		im.dx(imdx,true);
		im.dy(imdy,true);
		_FlowPrecision* data=imfeature.data();
		for(int i=0;i<height;i++)
			for(int j=0;j<width;j++)
			{
				int offset=i*width+j;
				data[offset*3]=im.data()[offset];
				data[offset*3+1]=imdx.data()[offset];
				data[offset*3+2]=imdy.data()[offset];
			}
	}
	else if(nchannels==3)
	{
		FImage grayImage;
		im.desaturate(grayImage);

		imfeature.allocate(im.width(),im.height(),5);
		FImage imdx,imdy;
		grayImage.dx(imdx,true);
		grayImage.dy(imdy,true);
		_FlowPrecision* data=imfeature.data();
		for(int i=0;i<height;i++)
			for(int j=0;j<width;j++)
			{
				int offset=i*width+j;
				data[offset*5]=grayImage.data()[offset];
				data[offset*5+1]=imdx.data()[offset];
				data[offset*5+2]=imdy.data()[offset];
				data[offset*5+3]=im.data()[offset*3+1]-im.data()[offset*3];
				data[offset*5+4]=im.data()[offset*3+1]-im.data()[offset*3+2];
			}
	}
	else
		imfeature.copyData(im);
}

bool OpticalFlow::LoadOpticalFlow(const char* filename,FImage &flow)
{
	Image<unsigned short int> foo;
	if(foo.loadImage(filename) == false)
		return false;
	if(!flow.matchDimension(foo))
		flow.allocate(foo);
	for(int  i = 0;i<flow.npixels();i++)
	{
		flow.data()[i*2] = (float)foo.data()[i*2]/160-200;
		flow.data()[i*2+1] = (float)foo.data()[i*2+1]/160-200;
	}
	return true;
}

bool OpticalFlow::LoadOpticalFlow(ifstream& myfile,FImage& flow)
{
	Image<unsigned short int> foo;
	if(foo.loadImage(myfile) == false)
		return false;
	if(!flow.matchDimension(foo))
		flow.allocate(foo);
	for(int  i = 0;i<flow.npixels();i++)
	{
		flow.data()[i*2] = (float)foo.data()[i*2]/160-200;
		flow.data()[i*2+1] = (float)foo.data()[i*2+1]/160-200;
	}
	return true;
}

bool OpticalFlow::SaveOpticalFlow(const FImage& flow, const char* filename)
{
	Image<unsigned short int> foo;
	foo.allocate(flow);
	for(int i =0;i<flow.npixels();i++)
	{
		foo.data()[i*2] = (__min(__max(flow.data()[i*2],-200),200)+200)*160;
		foo.data()[i*2+1] = (__min(__max(flow.data()[i*2+1],-200),200)+200)*160;
	}
	return foo.saveImage(filename);
}

bool OpticalFlow::SaveOpticalFlow(const FImage& flow,ofstream& myfile)
{
	Image<unsigned short int> foo;
	foo.allocate(flow);
	for(int i =0;i<flow.npixels();i++)
	{
		foo.data()[i*2] = (__min(__max(flow.data()[i*2],-200),200)+200)*160;
		foo.data()[i*2+1] = (__min(__max(flow.data()[i*2+1],-200),200)+200)*160;
	}
	return foo.saveImage(myfile);
}

bool OpticalFlow::showFlow(const FImage& flow,const char* filename)
{
	if(flow.nchannels()!=1)
	{
		cout<<"The flow must be a single channel image!"<<endl;
		return false;
	}
	Image<unsigned char> foo;
	foo.allocate(flow.width(),flow.height());
	float Max = flow.maximum();
	float Min = flow.minimum();
	for(int i = 0;i<flow.npixels(); i++)
		foo[i] = (flow[i]-Min)/(Max-Min)*255;
	foo.imwrite(filename);
}

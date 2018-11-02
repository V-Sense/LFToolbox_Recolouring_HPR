// ʹ��Flann����������ƥ��.cpp : �������̨Ӧ�ó������ڵ㡣
//
#include "Util.h"
#include <opencv2/opencv.hpp>
#include <highgui/highgui.hpp>
#include <xfeatures2d/nonfree.hpp>
#include <vector>
using namespace cv;
using namespace xfeatures2d;
using namespace std;

int main(int argc, char* argv[])
{
// 	if (argc != 4){
// 		printf("USAGE: sift.exe image1.png image2.png out.txt\n");
// 		return 0;
// 	}

// 	Mat input1 = imread(argv[1], 1);
// 	Mat input2 = imread(argv[2], 1);
	Mat input1 = imread("d:/frame_0024.png", 1);
	Mat input2 = imread("d:/frame_0025.png", 1);
	if (input1.empty() || input2.empty())
	{
		cout << "������������ͼƬ" << endl;
		system("pause");
		return -1;
	}
	CTimer t;

	/************************************************************************/
	/*���������ȡ������*/
	/************************************************************************/
	vector<KeyPoint> kerpoints1;

	Ptr<Feature2D> sift;
	sift = SIFT::create();
	sift->detect(input1, kerpoints1);

// 	Mat output1;
// 	drawKeypoints(input1, kerpoints1, output1);
	vector<KeyPoint> kerpoints2;
	sift->detect(input2, kerpoints2);

// 	Mat output2;
// 	drawKeypoints(input2, kerpoints2, output2);
// 	imshow("��ȡ��������box.png", output1);
// 	imshow("��ȡ��������box_in_scene.png", output2);
// 	imwrite("��ȡ��������box.png", output1);
// 	imwrite("��ȡ��������box_in_scene.png", output2);
// 	cout << "box��ȡ����������Ϊ:" << kerpoints1.size() << endl;
// 	cout << "box_in_scene����������Ϊ:" << kerpoints2.size() << endl;
	/************************************************************************/
	/* �����������������ȡ */
	/************************************************************************/
	Mat description1;
	sift->compute(input1, kerpoints1, description1);
	Mat description2;
	sift->compute(input2, kerpoints2, description2);
	/************************************************************************/
	/* ����������������ٽ�ƥ�� */
	/************************************************************************/
	vector<DMatch> matches;
	FlannBasedMatcher matcher;
	Mat image_match;
	matcher.match(description1, description2, matches);
	if (matches.size() != kerpoints1.size()){
		printf("SIFT Exception!\n");
	}
	// save match
	printf("%d\n", matches.size());
// 	FILE *fid = fopen(argv[3], "w");
// 	for (int i = 0; i < matches.size(); i++){
// 		int srcIdx = matches[i].queryIdx;
// 		int dstIdx = matches[i].trainIdx;
// 		float x1 = kerpoints1[srcIdx].pt.x;
// 		float y1 = kerpoints1[srcIdx].pt.y;
// 		float x2 = kerpoints2[dstIdx].pt.x;
// 		float y2 = kerpoints2[dstIdx].pt.y;
// 		fprintf(fid, "%.3f %.3f %.3f %.3f 1 100\n", x1, y1, x2, y2);
// 	}
// 	fclose(fid);
	t.toc("sift: ");
	return 0;

	/************************************************************************/
	/* �������������������ֵ����Сֵ */
	/************************************************************************/
	double max_dist = 0, min_dist = 10000;
	for (int i = 0; i < description1.rows; i++)
	{
		if (matches.at(i).distance>max_dist)
		{
			max_dist = matches[i].distance;
		}
		if (matches[i].distance < min_dist)
		{
			min_dist = matches[i].distance;
		}
	}
	cout << "��С����Ϊ" << min_dist << endl;
	cout << "������Ϊ" << max_dist << endl;
	/************************************************************************/
	/* �õ�����С�ڶ�V����С�����ƥ�� */
	/************************************************************************/
	vector<DMatch> good_matches;
	for (int i = 0; i < matches.size(); i++)
	{
		if (matches[i].distance < 0.8*max_dist)
		{
			good_matches.push_back(matches[i]);
			cout << "��һ��ͼ�е�" << matches[i].queryIdx << "ƥ���˵ڶ���ͼ�е�" << matches[i].trainIdx << endl;
		}
	}
	drawMatches(input1, kerpoints1, input2, kerpoints2, good_matches, image_match);
	imshow("ƥ����ͼƬ", image_match);
	imwrite("ƥ����ͼƬ.png", image_match);
	cout << "ƥ�����������Ϊ:" << good_matches.size() << endl;
	waitKey(0);
	return 0;
}
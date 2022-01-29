#include <opencv2/opencv.hpp>
#include <iostream>
#include "fluids.h"
using namespace std;
using namespace cv;

Mat ShowScalar(int width, int height, float* scalar, float scale);

int main()
{
	Mat region_flag_img = imread("7.png", IMREAD_GRAYSCALE);
	imshow("region_flag_img", region_flag_img);

	int width = region_flag_img.size().width;
	int height = region_flag_img.size().height;

	float* u = new float[(width + 1) * height];
	float* v = new float[width * (height + 1)];
	float* u_temp = new float[(width + 1) * height];
	float* v_temp = new float[width * (height + 1)];
	float* p = new float[width * height];       // this can be omitted by reusing an idle array(e.g. u_temp)
	float* s1 = new float[width * height];
	float* s1_temp = new float[width * height];
	unsigned char* region_flag = new unsigned char[width * height];

	Fluids2D fluids;
	fluids.u = u;
	fluids.v = v;
	fluids.u_temp = u_temp;
	fluids.v_temp = v_temp;
	fluids.p = p;
	fluids.scalars.push_back(s1);
	fluids.scalars_temp.push_back(s1_temp);

	memset(u, 0, (width + 1) * height * sizeof(float));
	memset(v, 0, width * (height + 1) * sizeof(float));
	memset(u_temp, 0, (width + 1) * height * sizeof(float));
	memset(v_temp, 0, width * (height + 1) * sizeof(float));
	memset(p, 0, width * height * sizeof(float));
	memset(s1, 0, width * height * sizeof(float));
	memset(s1_temp, 0, width * height * sizeof(float));

	// TODO : init region_flag
	
	memcpy(region_flag, region_flag_img.data, height * width * sizeof(unsigned char));

	Mat init = ShowScalar(width, height, fluids.scalars.at(0), 20);
	imshow("s1", init);
	waitKey(0);

	//VideoWriter out;
	//out.open("D:\\research\\video\\my_stable_fluids" + to_string(rand() * 200) + ".avi", CV_FOURCC('X', 'V', 'I', 'D'), 60, init.size());

	for (int iter = 0; iter < 3000; iter++)
	{
		int x = 300, y = height / 2, r1 = 20, r2 = 20;
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				if (pow(i - x, 2) + pow(j - y, 2) <= r1 * r1)
				{
					s1[j * width + i] += 4;
				}
				if (pow(i - x, 2) + pow(j - y, 2) <= r2 * r2)
				{
					u[j * (width + 1) + i] = -8;
					v[j * width + i] = 0;
				}
			}
		}

		UpdateFluids(width, height, region_flag, &fluids);

		Mat s1_img = ShowScalar(width, height, fluids.scalars.at(0), 20);
		//out << s1_img;
		imshow("s1", s1_img);
		if (waitKey(1) == 'q')
			break;
	}
	//out.release();

	cout << "over" << endl;
	waitKey(0);

	delete[] u;
	delete[] v;
	delete[] u_temp;
	delete[] v_temp;
	delete[] p;
	delete[] s1;
	delete[] s1_temp;
	delete[] region_flag;

	return 0;
}

Mat ShowScalar(int width, int height, float* scalar, float scale)
{
	Mat img = Mat::zeros(Size(width, height), CV_8UC3);

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			float val = scalar[j * width + i] * scale;
			if (val > 255) val = 255;

			img.at<Vec3b>(Point(i, j)) = Vec3b(val, val, val);
		}
	}
	return img;
}
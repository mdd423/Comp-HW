
#include <iostream>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <cmath>

std::vector<float> forwardFirst(std::vector<float> xs, std::vector<float> ys, int size)
{
	std::vector<float> deriv;
	deriv.resize(size - 1);
	for (int iii{ 0 }; iii < size - 1; ++iii)
	{
		deriv[iii] = (ys[iii + 1] - ys[iii]) / (xs[iii + 1] - xs[iii]);
	}
	return deriv;
}

std::vector<float> centralFirst(std::vector<float> xs, std::vector<float> ys, int size)
{
	std::vector<float> deriv;
	deriv.resize(size - 2);
	for (int iii{ 0 }; iii < size - 2; ++iii)
	{
		deriv[iii] = (ys[iii + 2] - ys[iii]) / (xs[iii + 2] - xs[iii]);
	}
	return deriv;
}

std::vector<float> centralFirstFive(std::vector<float> xs, std::vector<float> ys, int size, float len)
{
	std::vector<float> deriv;
	deriv.resize(size - 4);
	for (int iii{ 0 }; iii < size - 4; ++iii)
	{
		deriv[iii] = (-ys[iii + 4] + 8 * ys[iii + 3] - 8 * ys[iii + 1] + ys[iii]) / (12 * len);
	}
	return deriv;
}

float getError(float numerical, float value)
{
	return abs((numerical - value) / value);
}


int main()
{
	int N{ 15 };
	float x{ 0.1 };
	float y = cos(x);
	float d = -sin(x);
	float h{ 10 };

	std::ofstream deriveFile;
	deriveFile.open("derivative_cos_10.txt");

	for (int iii{ 0 }; iii < N; ++iii)
	{
		float err_for{ 0.0 };
		float err_cen{ 0.0 };
		float err_cen5{ 0.0 };
		std::vector<float> xs_for{ x, x + h };
		std::vector<float> ys_for;
		ys_for.resize(2);
		for (int jjj{ 0 };  jjj < 2; ++jjj)
		{
			ys_for[jjj] = cos(xs_for[jjj]);
		}
		std::vector<float> xs_cen{ x - h, x, x + h };
		std::vector<float> ys_cen;
		ys_cen.resize(3);
		for (int jjj{ 0 }; jjj < 3; ++jjj)
		{
			ys_cen[jjj] = cos(xs_cen[jjj]);
		}
		std::vector<float> xs_cen5{ x - (2*h), x - h, x, x + h, x + (2*h) };
		std::vector<float> ys_cen5;
		ys_cen5.resize(5);
		for (int jjj{ 0 }; jjj < 5; ++jjj)
		{
			ys_cen5[jjj] = cos(xs_cen5[jjj]);
		}

		float num_for = (forwardFirst(xs_for, ys_for, 2))[0];
		float num_cen = (centralFirst(xs_cen, ys_cen, 3))[0];
		float num_cen5 = (centralFirstFive(xs_cen5, ys_cen5, 5, h))[0];

		err_for = getError(num_for, d);
		err_cen = getError(num_cen, d);
		err_cen5 = getError(num_cen5, d);

		deriveFile << iii << "\t" << h << "\t" << err_for << "\t" << err_cen << "\t" << err_cen5 << std::endl;
		h *= .5;
	}

	deriveFile.close();
	return 0;
}
�

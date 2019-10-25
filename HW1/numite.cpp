
#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

float simpsons(std::vector<float> xs, int size, float len)
{
	float integr{ 0. };
	integr += exp(-xs[0]);
	integr += exp(-xs[size-1]);
	for (int iii{ 1 }; iii < size; iii += 2)
	{
		integr += 4 * exp(-xs[iii]);
		// std::cout << iii << std::endl;
	}
	for (int iii{ 2 }; iii < size-1; iii += 2)
	{
		integr += 2 * exp(-xs[iii]);
		// std::cout << iii << std::endl;
	}
	return integr * len / 3.0;
}

int main()
{
	int times{ 22 };
	int N{ 2 };

	float anlyt_result = 1. - exp(-1.);

	std::ofstream intFile;
	intFile.open("integration.txt");
	for (int mmm{ 0 }; mmm < times; ++mmm)
	{
		float h = 1 / static_cast<float>(N);
		std::vector<float> x;
		x.resize(N+1);
		for (int iii{ 0 }; iii < N+1; ++iii) x[iii] = iii / static_cast<float>(N);

		float midpt{ 0 };
		for (int iii{ 0 }; iii < N - 1; ++iii) midpt += exp(-(x[iii] + (h / 2.0))) * h;

		float traps{ 0 };
		for (int iii{ 0 }; iii < N - 1; ++iii) traps += (exp(-x[iii]) + exp(-x[iii + 1])) / (2.0) * h;

		float simp{ 0 };
		simp = simpsons(x, N, h);

		intFile << mmm << "\t" << N << "\t" << h << "\t";
		intFile	<< abs((midpt - anlyt_result)/anlyt_result) << "\t" << abs((traps - anlyt_result)/anlyt_result) << "\t";
		intFile << abs((simp - anlyt_result)/anlyt_result) << "\t";
		intFile << midpt << "\t" << traps << "\t" << simp << "\t" << anlyt_result << std::endl;

		N *= 2;
	}
	intFile.close();
	return 0;
}

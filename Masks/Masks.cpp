#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
using namespace std;

const int lookup_size = 512;
float lookup[lookup_size][lookup_size] = {0};

/*Ustawianie lookup table dla funkcji gaussa*/
void set_lookup_table(float sampling_x, float sampling_y, float w0 = 7.0, float lambda = 0.632)
{
	//float z0 = M_PI * w0 * w0 / lambda;
	float x = 0, y = 0;

	for (int i = 0; i < lookup_size; i++)
	{
		x = (i - lookup_size / 2)*sampling_x;
		for (int j = 0; j < lookup_size; j++)
		{
			y = (j - lookup_size / 2)*sampling_y;
			lookup[i][j] = exp(-(x*x + y * y) / w0 / w0);
		}
	}
}

/*Zwraca wartoœæ funkcji cos^2 dla danego piksela*/
void cos2(int px_x, int px_y, float sampling_x, float sampling_y, float &value_px) //sampling in um//
{
	float x = px_x * sampling_x;
	float y = px_y * sampling_y;
	x = cos(x / 8 * 2 * M_PI) * cos(x / 8 * 2 * M_PI);
	y = cos(y / 8 * 2 * M_PI) * cos(y / 8 * 2 * M_PI);
	value_px = x * y;
}

/*Zwraca wartoœæ funkcji gaussowskiej dla danego piksela*/
void gauss(int px_x, int px_y, float sampling_x, float sampling_y, float &value_px) //sampling in um//
{
	float x = px_x * sampling_x;
	float y = px_y * sampling_y;
	int lookup_x = (int)(fmod(x, 8)/8 * 512);
	int lookup_y = (int)(fmod(y, 8)/8 * 512);
	value_px = lookup[lookup_x][lookup_y];
}

int main()
{
	cout << "START\n";

	float sampling_x = 1; //um
	float sampling_y = 1; //um
	float w0 = 7; //um
	float lambda = 0.632; //um
	const int mask_size = 128; // 32768;

	set_lookup_table(sampling_x, sampling_y, w0, lambda);
	cout << "Lookup table set!\n";

	float gauss_mask[mask_size][mask_size] = { 0 };
	float value_gauss = 0;
	float cos_mask[mask_size][mask_size] = { 0 };
	float value_cos = 0;

	cout << "Ready to create masks!\n";

	for (int i = 0; i < mask_size; i++)
	{
		for (int j = 0; j < mask_size; j++)
		{
			gauss(i, j, sampling_x, sampling_y, value_gauss);
			gauss_mask[i][j] = value_gauss;
		}
	}

	cout << "Gaussian mask ready!\n";

	for (int i = 0; i < mask_size; i++)
	{
		for (int j = 0; j < mask_size; j++)
		{
			cos2(i, j, sampling_x, sampling_y, value_cos);
			cos_mask[i][j] = value_cos;
		}
	}
	
	cout << "Cos^2 mask ready!\nSTOP";
}
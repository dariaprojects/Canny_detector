#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>


#pragma pack(1)

#define THR 3

typedef struct tagRGBTRIPLE
{
	char rgbBlue;
	char rgbGreen;
	char rgbRed;
} RGBTRIPLE;

struct f_info
{
	unsigned char signature[2];
	unsigned int sizefile;
	unsigned int reserved;
	unsigned int addr_offset;
};

struct pic_info
{
	unsigned int Size;
	unsigned int Width;
	unsigned int Height;
	unsigned short int  Planes;
	unsigned short int  BitCount;
	unsigned int Compression;
	unsigned int SizeImage;
	unsigned int XPelsPerMeter;
	unsigned int YPelsPerMeter;
	unsigned int ClrUsed;
	unsigned int ClrImportant;
};

void Filtr_Cob(float *First, float *Second, int height, int width) {

	int window[9] = { -1,0,1,
		-2,0,2,
		-1,0,1 };

	int window2[9] = {
		-1,-2,-1,
		0,0,0,
		1,2,1 };

	int Mx = 3, My = 3, r,g=9,winf2, winf, Pm, Qm;
	r = 9;
	winf = 0;
	winf2 = 0;
	while (r--)winf += window[r];
	if (winf == 0)winf = 1;

	while (g--)winf2 += window2[g];
	if (winf2 == 0)winf2 = 1;

	Pm = Mx / 2;
	Qm = My / 2;
	for (unsigned int i = 1; i < height - 1; i++) {
		for (unsigned int j = 1; j < width - 1; j++) {
			r = 0;
			g = 0;
			float SumI = 0;
			float SumI1 = 0;
			for (int l = -Qm; l <= Qm; l++) {
				for (int k = -Pm; k <= Pm; k++) {
					float I;
					I = First[(i + l)*width + (j + k)];
					SumI += I * window[r++];
					SumI1 += I * window2[g++];
				}
			}
			Second[i*width + j] = sqrt(SumI * SumI + SumI1 * SumI1);
		}
	}
}

void Filtr_G(float *First, float *Second, int height, int width) {

	int window[9] = { 1,1,1,
		1,2,1,
		1,1,1 };

	int Mx = 3, My = 3, r, winf, Pm, Qm;
	r = 9;
	winf = 0;
	while (r--)winf += window[r];
	if (winf == 0)winf = 1;

	Pm = Mx / 2;
	Qm = My / 2;
	for (unsigned int i = 1; i<height - 1; i++) {
		for (unsigned int j = 1; j<width - 1; j++) {
			r = 0;
			float SumI = 0;
			for (int l = -Qm; l <= Qm; l++) {
				for (int k = -Pm; k <= Pm; k++) {
					float I;
					I = First[(i + l)*width + (j + k)];
					SumI += I*window[r];
					r += 1;
				}
			}
			SumI = SumI / (winf);
			Second[i*width + j] = SumI;
		}
	}
}


int Check(float inf1, float inf2) {
	if (inf1 <= inf2)
		return 1;
	else
		return 0;
}

void NonMaxSupp(float *Matrix1, int height, int width, float *Matrix2) {
	for (int i = 1; i<height - 1; i++) {
		for (int j = 1; j<width - 1; j++) {
			int dx, dy;
			dx = (int)fabs(cos(Matrix1[i*width + j]));
			dy = (int)-fabs((sin(Matrix1[i*width + j])));
			if (Check(Matrix1[(i + dy)*width + (j + dx)], Matrix1[i*width + j]) == 1) {
				Matrix2[(i + dy)*width + (j + dx)] = 0;
			}
			if (Check(Matrix1[(i - dy)*width + (j - dx)], Matrix1[i*width + j]) == 1) {
				Matrix2[(i - dy)*width + (j - dx)] = 0;
			}
			Matrix2[i*width + j] = Matrix1[i*width + j];
		}
	}

}

void Double_p_f(float *Matrix1, int height, int width, float *Matrix2, float low_p, float high_p) {
	float down, up;
	down = low_p * 255;
	up = high_p * 255;
	for (int i = 1; i<height - 1; i++) {
		for (int j = 1; j<width - 1; j++) {
			if (Matrix1[i*width + j] >= up) {
				Matrix2[i*width + j] = 255;
			}
			else if (Matrix1[i*width + j] <= down) {
				Matrix2[i*width + j] = 0;
			}
			else {
				Matrix2[i*width + j] = 127;
			}
		}
	}
}



void Blob_analis(float *Matrix1, int height, int width, float *Matrix2, int high, int clear) {
	int MoveDir[2][8] = { { -1,-1,-1,0,0,1,1,1 },
						{ -1,0,1,-1,1,-1,0,1 } };
	for (int i = 0; i<height - 1; i++) {
		for (int j = 0; j<width - 1; j++) {
			if (Matrix1[i*width + j] == high) {
				Matrix2[i*width + j] = high;
				for (int k = 0; k<8; k++) {
					int dx = MoveDir[0][k];
					int dy = MoveDir[1][k];
					int im = i, jm = j;
					while (1) {
						im += dx;
						jm += dy;
						if (((im*width + jm)<0) || ((im*width + jm)>(height - 1)*(width - 1))) {
							break;
						}
						
						if ((Matrix1[im*width + jm] == clear) || (Matrix1[(i + dx)*width + (j + dy)] == high)) {
							float x = Matrix1[im*width + jm];
							break;
						}
						Matrix2[im*width + jm] = high;
					}
				}
			}
			else {
				Matrix2[i*width + j] = clear;
			}
		}
	}
}

int main(int argc, char *argv[])
{
	char infile[100] = "E:\\VSProjects\\Canny_method\\11.bmp";
	char outfile[100] = "new1.bmp";
	if (argc>1) strcpy_s(infile, argv[1]);
	if (argc>2) strcpy_s(outfile, argv[2]);

	std::ifstream in(infile, std::ios::in | std::ios::binary);
	if (!in) { std::cout << "File not found..."; exit(1); }
	struct f_info f_i;
	struct pic_info pic_i;
	in.read((char*)&f_i, sizeof(f_info));
	in.read((char*)&pic_i, sizeof(pic_info));
	int ln_str = (f_i.sizefile - 54) / pic_i.Height;
	float* image = new float[pic_i.Height * pic_i.Width];
	std::cout << "Width - " << pic_i.Width << ", Height - " << pic_i.Height << std::endl;

	char d;
	for (unsigned int i = 0; i < pic_i.Height; i++)
	{
		for (unsigned int j = 0; j < pic_i.Width; j++)
		{
			in.get(d);
			image[i * pic_i.Width + j] = (unsigned char)d * 0.114;
			in.get(d);
			image[i * pic_i.Width + j] += (unsigned char)d * 0.587;
			in.get(d);
			image[i * pic_i.Width + j] += (unsigned char)d * 0.299;
		}
		for (unsigned int k = 0; k < (ln_str - pic_i.Width * 3); k++)
		{
			char d;

			in.get(d);
		}
	}
	in.close();

	// float *image2 = new float[pic_i.Height*pic_i.Width];
	float *image3 = new float[pic_i.Height*pic_i.Width];
	float *image4 = new float[pic_i.Height*pic_i.Width];
	float *image5 = new float[pic_i.Height*pic_i.Width];
	float *image6 = new float[pic_i.Height*pic_i.Width];
	float *image7 = new float[pic_i.Height*pic_i.Width];

	
	Filtr_G(image, image3, pic_i.Height, pic_i.Width);
	Filtr_Cob(image3, image4, pic_i.Height, pic_i.Width);
	NonMaxSupp(image4, pic_i.Height, pic_i.Width, image5);
	Double_p_f(image5, pic_i.Height, pic_i.Width, image6, 0.3, 1);
	Blob_analis(image6, pic_i.Height, pic_i.Width, image7, 127, 0);


	std::ofstream out;
	out.open(outfile, std::ios::out | std::ios::binary);
	out.write((char*)&f_i, sizeof(f_info));
	out.write((char*)&pic_i, sizeof(pic_info));
	for (unsigned int i = 0; i<pic_i.Height; i++)
	{
		for (unsigned int j = 0; j<pic_i.Width; j++)
		{
			out.put(image7[i*pic_i.Width + j]);
			out.put(image7[i*pic_i.Width + j]);
			out.put(image7[i*pic_i.Width + j]);
		}
		for (unsigned int k = 0; k<(ln_str - pic_i.Width * 3); k++) out.put(0);
	}

	std::cout << "Press any key...\n";
	std::cout << "Width - " << pic_i.Width << ", Height - " << pic_i.Height << std::endl;

	delete image;
	delete image3;
	delete image4;
	delete image5;
	delete image6;
	delete image7;
	return 0;
}


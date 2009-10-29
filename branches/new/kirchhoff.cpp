/*
 * kirchhoff.cpp
 *
 *  Created on: 22-Oct-2009
 *      Author: dima
 */

#include <cmath>

using std::sqrt;

#include "kirchhoff.h"

const double PI=3.1415926535897932384626433832795028841971693993751;

void Kirchhoff2d::create(){

	if(!need_create)return;

	cout << " Kirchhoff2d::create()\n";

	// create fourier

	if(fourier != NULL){
		delete fourier;
		fourier = NULL;
	}

	fourier = new Fourier2d<double_complex>(comm, 2*height, 2*width);

	Matrix<double_complex>* image = new Matrix<double_complex>(2*height/comm_size, 2*width);

	// pre-compute
	double dis2 = distance*distance;
	double k = 2 * PI / lambda;

	int y_start = comm_rank * 2*height/comm_size;

	for (int y = 0; y < 2*height/comm_size; y++)
	{
		//              /\
		// we need:    /  \
		//
		double dy = (height - abs(y+y_start-height)) / resolution;

		for (int x = 0; x < 2*width; x++)
		{
		double dx = (width - abs(x-width)) / resolution;

			double r = sqrt( dx*dx + dy*dy + dis2);
			if( r == 0.0 ) continue;

			double kr = k*r;
			double CosKR = cos(kr);
			double SinKR = sin(kr);

			(*image)[y][x].re = CosKR/r/lambda;
			(*image)[y][x].im = SinKR/r/lambda;
		}// y
	}// x

	// and fourier it!
	// create image

	if(f_image != NULL){
		delete f_image;
		f_image = NULL;
	}

	f_image = fourier->exec(image);
	delete image;

	need_create = false;

	cout << " Kirchhoff2d::~create()\n";
}// create


Matrix<double_complex>* Kirchhoff2d::exec(const Matrix<double_complex>* arg){
	assert(arg->getHeight() == height/comm_size*2 && arg->getWidth() == width*2);

	create();

	//matrix must be already doubled!

	cout << " Kirchhoff2d::exec()\n";

	// four 1
	cout << " Kirchhoff2d::exec()::forward\n";
	Matrix<double_complex>* f_in = fourier->exec(arg);
	cout << " Kirchhoff2d::exec()::~forward\n";

	// multiply f_in*f_image->f_in
	cout << " Kirchhoff2d::exec()::mult\n";
	for(int i=0; i<2*height/comm_size; i++){
		for(int j=0; j<2*width; j++){
			double r1 = (*f_in)[i][j].re;
			double i1 = (*f_in)[i][j].im;
			double r2 = (*f_image)[i][j].re;
			double i2 = (*f_image)[i][j].im;

			double re = r1*r2 - i1*i2;
			double im = r1*i2 + r2*i1;
			(*f_in)[i][j].re = re;
			(*f_in)[i][j].im = im;
		}
	}
	cout << " Kirchhoff2d::exec()::~mult\n";

	fourier->setInverse(true);

	// fourier back
	cout << " Kirchhoff2d::exec()::backward\n";
	Matrix<double_complex>* out = fourier->exec(f_in);
	cout << " Kirchhoff2d::exec()::~backward\n";
	delete f_in;

	cout << " Kirchhoff2d::~exec()\n";

	return out;
}// exec!

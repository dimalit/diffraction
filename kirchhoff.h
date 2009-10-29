/*
 * kirchhoff.h
 *
 *  Created on: 21-Oct-2009
 *      Author: dima
 */

#ifndef KIRCHHOFF_H_
#define KIRCHHOFF_H_

#include "fourier.h"
#include "mpiclass.h"

struct double_complex{
	double re;
	double im;
};

class Kirchhoff2d: private MpiClass{
private:
	int height, width;
	double resolution, distance, lambda;

	Fourier2d<double_complex>* fourier;
	Matrix<double_complex>*	   f_image;

	bool need_create;

	// re-create fourier if comm or sizes were changed
	// re-create image   if lambda, distance or res were changed
	void create();
public:
	Kirchhoff2d():MpiClass(MPI_COMM_WORLD){
		lambda = 0; distance = 0;
		resolution = 0;
		width = 0;  height = 0;

//		comm = MPI_COMM_WORLD;

		fourier = NULL;
		f_image = NULL;
		need_create = true;
	}

	Kirchhoff2d(MPI_Comm cm, int h, int w, double res, double dst, double l):MpiClass(cm){

		height = h;
		width  = w;
		resolution = res;
		distance = dst;
		lambda = res;

		fourier = NULL;
		f_image = NULL;
		need_create = true;
	}

	void setLambda(double l)		{ lambda = l;       need_create = true;}
	void setDistance(double d)		{ distance = d;     need_create = true;}
	void setResolution(double res)  { resolution = res; need_create = true;}
	void setWidth(int w) 			{width = w;  		need_create = true;}
	void setHeight(int h) 			{height = h; 		need_create = true;}
	void setComm(MPI_Comm c){
		comm = c;
		MPI_Comm_size(comm, &comm_size);
		MPI_Comm_rank(comm, &comm_rank);

		need_create = true;
	}

	double getLambda()		const {return lambda;}
	double getDistance()	const { return distance; }
	double getResolution()	const { return resolution;}
	int getWidth()			const { return width; }
	int getHeight()			const { return height; }
	MPI_Comm getComm()		const {return comm;}


	Matrix<double_complex>* exec(const Matrix<double_complex>* arg);

	~Kirchhoff2d(){
		if(fourier != NULL)
			delete fourier;
		if(f_image != NULL)
			delete f_image;
	}
};

#endif /* KIRCHHOFF_H_ */

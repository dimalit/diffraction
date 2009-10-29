/*
 * fourier.h
 *
 *  Created on: 21-Oct-2009
 *      Author: dima
 */

#ifndef FOURIER_H_
#define FOURIER_H_

#include <iostream.h>

#include <mpi.h>
#include <mkl_cdft.h>

#include "mpiclass.h"
#include "matrix.h"

template<class T> class Fourier2d: private MpiClass{
private:
	int height, width;
	bool inverse;

	bool need_create;			// whether in exec we will re-create plan

// MKL-specific:
	DFTI_DESCRIPTOR_DM_HANDLE mkl_desc;

private:

	// prepare for FFT
	void create();

public:
	Fourier2d(MPI_Comm comm, int nrows, int ncols);

	void setInverse(bool inv){
		inverse = inv;
	}

	Matrix<T>* exec(const Matrix<T>* in);

	~Fourier2d(){
	}
};

template<class T> void Fourier2d<T>::create(){
	if(!need_create)
		return;

	// free if needed
	if(mkl_desc){
		DftiFreeDescriptorDM(&mkl_desc);
		mkl_desc = NULL;
	}

	MKL_LONG len[2] = {height, width};

	// select precision
	DFTI_CONFIG_VALUE prec;
	switch(sizeof(T)){
	case 8:
		prec = DFTI_SINGLE;
		break;
	case 16:
		prec = DFTI_DOUBLE;
		break;
	}

	DftiCreateDescriptorDM((MPI_Comm)MPI_Comm_c2f(comm), &mkl_desc, prec, DFTI_COMPLEX, 2, len);

	// assert right sizes
	MKL_LONG v,nx,start;
	DftiGetValueDM(mkl_desc, CDFT_LOCAL_SIZE, &v);
	assert(v == width*height/comm_size);
	DftiGetValueDM(mkl_desc, CDFT_LOCAL_NX, &nx);
	assert(nx*width == v);
	DftiGetValueDM(mkl_desc, CDFT_LOCAL_X_START, &start);
	assert(start == comm_rank*nx);

	DftiSetValueDM(mkl_desc,DFTI_PLACEMENT,DFTI_NOT_INPLACE);	// in->out
	DftiCommitDescriptorDM(mkl_desc);

	need_create = false;
}

template<class T> Fourier2d<T>::Fourier2d(MPI_Comm comm, int nrows, int ncols):MpiClass(comm){
	assert(nrows > 1 && ncols > 1);

	assert(nrows % comm_size == 0);

	height = nrows;
	width  = ncols;
	inverse = false;

	mkl_desc = NULL;
	need_create = true;
}

template<class T> Matrix<T>* Fourier2d<T>::exec(const Matrix<T>* in){

//	cout << ":" << endl;
//	cout << "height = "<< height << " comm_size = " << comm_size << " width = " << width << endl;
//	cout << "getHeight() = " << in->getHeight() << " getWidth() = " << in->getWidth() << endl << endl;

	assert(in->getHeight() == height / comm_size && in->getWidth() == width);

	create();

	Matrix<T>* out = new Matrix<T>(height / comm_size, width);

	// temp pointers
	void* a =  in->getContiguos();
	void* b = out->getContiguos();
	if(!inverse)
		DftiComputeForwardDM(mkl_desc, a, b);
	else
		DftiComputeBackwardDM(mkl_desc, a, b);

	return out;		// feel free to free it!
}

#endif /* FOURIER_H_ */

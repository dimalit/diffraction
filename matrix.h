/*
 * matrix.h
 *
 *  Created on: 21-Oct-2009
 *      Author: dima
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <memory.h>
#include <assert.h>

template<class T> class Matrix{

private:
	int width, height;
	T*  buffer;
	T** data;

public:
	Matrix(int nrows, int ncols);

	int getHeight() const {return height;}
	int getWidth()  const {return width;}

	T* operator[](int row) const{
		assert(row >= 0 && row < height);
		return data[row];
	}

	T* operator[](int row){
		assert(row >= 0 && row < height);
		return data[row];
	}

	T** getStructured() const{
		return data;
	}

	T* getContiguos() const{
		return buffer;
	}

	void clear(){
		memset(buffer, 0, sizeof(T)*height*width);
	}

	~Matrix(){
		delete[] data;
		delete[] buffer;
	}
};

template<class T> Matrix<T>::Matrix(int nrows, int ncols){
	assert(nrows > 1 && ncols > 1);

	height = nrows;
	width  = ncols;

	// prepare continuous buffer
	buffer = new T[width*height];

	// prepare structured
	data = new T*[height];
	for(int i=0; i<height; i++){
		data[i] = &buffer[i*width];
	}
}

#endif /* MATRIX_H_ */

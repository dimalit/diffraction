/*
 * main.cpp
 *
 *  Created on: 21-Oct-2009
 *      Author: dima
 */

#include <unistd.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

#include "kirchhoff.h"

int comm_rank;
int comm_size;

int test_fourier(int argc, char** argv){

	if(argc != 3){
		cout << "USAGE: diffraction_new size cycle_count\n";
		return -1;
	}

	int dim = atoi(argv[1]);
	int cnt = atoi(argv[2]);

	Fourier2d<double_complex> fourier(MPI_COMM_WORLD, dim, dim);

	Matrix<double_complex> in(dim, dim);
	in.clear();

	Matrix<double_complex>* out;
	for(int i=0; i<cnt; i++){
		out = fourier.exec(&in);
		delete out;
	}

	return 0;
}

int test_kirchhoff(int argc, char** argv){

	if(argc != 3){
		cout << "USAGE: diffraction_new size cycle_count\n";
		return -1;
	}

	int dim = atoi(argv[1]);
	int cnt = atoi(argv[2]);

	Kirchhoff2d kir(MPI_COMM_WORLD, dim, dim, 200, 0.006, 0.00063);

	Matrix<double_complex> in(2*dim/comm_size, 2*dim);
	in.clear();

	// create two points:
	if(comm_rank == 0){
		in[0][dim/2 - 2].re = 1;
		in[0][dim/2 + 2].re = 1;
	}

	Matrix<double_complex>* out;

	for(int i=0; i < cnt; i++){
		out = kir.exec(&in);
		delete out;
	}

	return 0;
}

int file_main(int argc, char** argv){
	//////////////// PARSE ARGS ////////////////////

	char *file_in, *file_out;
	int height, width;
	double res;
	double distance;
	double lambda;

	if( argc != 8 ){
		puts("USAGE: diffraction file_in file_out height width resolution distance lambda\n");
		return -1;
	}

	if(comm_size % 2 != 0){
		puts("-np must be divisible by two!\n");
		return -1;
	}

	int i=1;
	file_in  = argv[i++];
	file_out = argv[i++];
	height   = atoi(argv[i++]);
	width    = atoi(argv[i++]);
	res      = atof(argv[i++]);
	distance = atof(argv[i++]);
	lambda   = atof(argv[i++]);

	// prepare
	Kirchhoff2d algo(MPI_COMM_WORLD, height, width, res, distance, lambda);

	Matrix<double_complex> in(2*height / comm_size, 2*width);		// 2x matrix!
	in.clear();

	/////////////////// READ DATA //////////////////////

	// read and distribute
	if( comm_rank==0 ){

		FILE* fp = fopen(file_in, "rb");
		if(fp == NULL){
			cout << " Cannot open input file " << file_in << " - generating..."<< endl;
			fp = NULL;
		}

		// read full image and send only to comm_size/2 procs
		for(int i = 0; i < comm_size / 2; i++){
			for(int y = 0; y < in.getHeight(); y++){

					if(fp == NULL) continue;

					if(fread(in[y], sizeof(double_complex), width, fp) != width){
						cout << " Bad input file\n";
						fclose(fp);
						return -1;
					}// if read

			}// for y
			if(i!=0)
				MPI_Send(in.getContiguos(), in.getHeight()*in.getWidth() * 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}// for i

		// read for myself :)
		if(fp != NULL)
			fseek(fp, 0, 0);
		for(int y = 0; y < in.getHeight(); y++){

			if(fp == NULL) continue;

			if(fread(in[y], sizeof(double_complex), width, fp) != width){
				cout << " Bad input file\n";
				fclose(fp);
				return -1;
			}// if read

		}// for y

		if(fp != NULL){
			fclose(fp);
			cout << " loaded file " << file_in << endl;
		}

	}// if 0
	else if( comm_rank < comm_size / 2){		// receive
		MPI_Status st;
		MPI_Recv(in.getContiguos(), in.getHeight()*in.getWidth() * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);
		// TODO: return -1?
	}

	////////////// COMPUTE ////////////////////////
	cout << " started computing.\n";
	time_t tm = time(NULL);
	Matrix<double_complex>* out = algo.exec(&in);
	tm = time(NULL) - tm;

	cout << " Finished in " << tm << " s\n";

	/////////////// WRITE //////////////////////
	if(comm_rank == 0){

		FILE* fp = fopen(file_out, "wb");
		if(fp == NULL){
			cout << " Cannot open output file!\n";
			delete out;
			return -1;
		}

		for(int i=0; i < comm_size / 2; i++){

			MPI_Status st;
			if(i != 0)
				MPI_Recv(out->getContiguos(), out->getHeight()*out->getWidth()*2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &st);

			for(int y = 0; y < out->getHeight(); y++){
				if(fwrite((*out)[y], sizeof(double_complex), out->getWidth() / 2, fp) != out->getWidth() / 2){
					cout << " Cannot write output file!\n";
					delete out;
					fclose(fp);
					return -1;
				}
			}// for y
		}// for i

		fclose(fp);

		cout << " stored file " << file_out << endl;
	}
	else if( comm_rank < comm_size / 2){		// send
		MPI_Send(out->getContiguos(), out->getHeight()*out->getWidth()*2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	delete out;

	return 0;
}

int main(int argc, char** argv){

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

	int i = -1;
	while(i==comm_rank){
		sleep(5);
	}

//	int ret = file_main(argc, argv);
	int ret = test_kirchhoff(argc, argv);

	if(ret < 0)
		MPI_Abort(MPI_COMM_WORLD, ret);

	cout << comm_rank << " MPI_Finalize()\n";
	MPI_Finalize();

	return ret;
}

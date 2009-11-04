/*
 * my_fft.h
 *
 *  Created on: 02-Nov-2009
 *      Author: dima
 */

#ifndef MY_FFT_H_
#define MY_FFT_H_

#include <cmath>
#include <sys/time.h>
#include <time.h>
using namespace std;

double dtime(){
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int log2(int arg){
	// compute log
	int res = 0;
	for(; arg >> res != 1; res++);
	return res;
}

template<class T> struct complex{
	T re;
	T im;

	T mod2(){
		return re*re+im*im;
	}

	T mod(){
		return sqrt(mod2());
	}

	complex<T>& operator=(const complex& arg){
		this->re = arg.re;
		this->im = arg.im;
		return *this;
	}

	complex<T> operator+(const complex& arg){
		complex<T> res;
		res.re = this->re + arg.re;
		res.im = this->im + arg.im;

		return res;
	}

	complex<T> operator-(const complex& arg){
		complex<T> res;
		res.re = this->re - arg.re;
		res.im = this->im - arg.im;

		return res;
	}

	complex<T> operator*(const complex& arg){
		complex<T> res;
		res.re = this->re*arg.re - this->im*arg.im;
		res.im = this->re*arg.im + this->im*arg.re;
		//TODO: optimize?

		return res;
	}
};

const double PI=3.1415926535897932384626433832795028841971693993751;
const int BUF_SIZE = 16*1024;
//const int BUF_SIZE = 4;

template<class T> class my_fft{
private:
	int dim;
	int dim_log;
	T one;						// for +-1 in inverse
	complex<T> *exp_array;

	void create_exp_array(){

		// delete if needed
		if(exp_array != NULL){
			delete[] exp_array;
			exp_array = NULL;
		}

		// create
		exp_array = new complex<T>[dim/2];		// we don't need 2nd half

		// fill
		for(int i=0; i < dim/2; i++){
			exp_array[i].re = cos(2*PI*i/dim);
			exp_array[i].im = one*sin(2*PI*i/dim);
		}
	}

	// return i-th from 0..n
	complex<T> exp(int i, int n){
//		complex<T> res;
//		res.re = cos(2.0*PI*i/n);
//		res.im = one*sin(2.0*PI*i/n);
//		return res;

		return exp_array[i*(dim/n)];
	}

	int bit_reverse(int arg){
		int res = 0;
		for(int i=0; i<dim_log; i++){
			res <<= 1;
			res |= arg&1;
			arg >>= 1;
		}
		return res;
	}

	// 2-inversion
	void reorder(complex<T>* arg){
		for(int i=0; i < dim; i++){
			int j = bit_reverse(i);

			// swap only upwards
			if(j > i){
				complex<T> tmp = arg[i];
				arg[i] = arg[j];
				arg[j] = tmp;
			}
		}
	}

	// do one thick butterfly with 2*half_size elements
	void group_butterfly(complex<T>* arg, int half_size){

		int size = half_size * 2;

//		cout << dim/size << endl;
		for(int i=0; i<half_size; i++){

			complex<T> tmp = arg[i + half_size] * exp(i, size);
			arg[i + half_size] = arg[i] - tmp;
			arg[i]             = arg[i] + tmp;

//			cout << "EXP: " << exp(i, size).re << " " << exp(i, size).im << endl;

		}// for i
	}// goup_butterfly

	void remap(complex<T>* arg, int p){
		// p denotes num procs
		int b = dim / p;		// cache size

		// initially we have b blocks p each
		int ind1 = 0;		// source
		int ind2 = 0;		// destination
		for(int i=0; i < b; i++){
			for(int j=0; j < p; j++, ind1++){
				ind2 = b*j+i;		// ind1 = p*i+j;
				if(ind2 < ind1)
					continue;
				// else swap
				complex<T> tmp = arg[ind1];
				arg[ind1] = arg[ind2];
				arg[ind2] = tmp;
			}// for j

		}// for i

	}// b 2 c

public:
	my_fft(int dim, bool inverse=false){
		this->dim = dim;
		this->one = inverse ? +1.0 : -1.0;

		dim_log = log2(dim);

		exp_array = NULL;
		create_exp_array();
	}

	void setInverse(bool inv){
		double old_one = one;
		this->one = inv ? +1.0 : -1.0;
		if(one != old_one)
			create_exp_array();
	}

	void exec_inplace(complex<T>* arg){

		reorder(arg);

		int num_blocks = dim / BUF_SIZE;
		int num_steps  = log2(BUF_SIZE);
		int hs;
		for(int blk=0; blk<num_blocks; blk++){
			hs = 1;
			for(int step=0; step<num_steps; step++, hs*=2){

				int num_butts = BUF_SIZE / (hs * 2);

				for(int butt=0; butt<num_butts; butt++){
					group_butterfly(arg + blk*BUF_SIZE + (hs * 2)*butt, hs);
				}
			}
		}// for block

		// block to cyclic!
		remap(arg, num_blocks);

		// super-fast again ;):
		num_steps = dim_log - num_steps;
		for(int blk=0; blk<num_blocks; blk++){
			hs = 1;
			for(int step=0; step<num_steps; step++, hs*=2){

				int num_butts = BUF_SIZE / (hs * 2);

				for(int butt=0; butt<num_butts; butt++){
					group_butterfly(arg + blk*BUF_SIZE + (hs * 2)*butt, hs);
				}
			}
		}// for block

		// back to block
		remap(arg, BUF_SIZE);

/*
		// without block 2 cyclic
		for(int step=num_steps; step<dim_log; step++, hs*=2){
			int num_butts = dim / (hs * 2);

			for(int butt=0; butt<num_butts; butt++){
				group_butterfly(arg + (hs * 2)*butt, hs);
			}
		}// for step
*/

/*
		// without blocks
		int hs = 1;
		for(int step=0; step<dim_log; step++, hs*=2){
			int num_butts = dim / (hs * 2);

			for(int butt=0; butt<num_butts; butt++){
				group_butterfly(arg + (hs * 2)*butt, hs);
			}
		}// for step
*/
	}
};

#endif /* MY_FFT_H_ */

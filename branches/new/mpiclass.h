/*
 * mpiclass.h
 *
 *  Created on: 28-Oct-2009
 *      Author: dima
 */

#ifndef MPICLASS_H_
#define MPICLASS_H_

#include <mpi.h>
#include <iostream>
#include <iomanip>

class MpiClass {
private:
	static double start_wtime;
public:
	MPI_Comm comm;
	int comm_size;
	int comm_rank;

	class logstream{
	private:
		MpiClass& owner;
	public:
		logstream(MpiClass& ow):owner(ow){}

		template<class T> std::ostream& operator <<(const T& arg){
			std::cout << std::fixed << std::setprecision(3) << MPI_Wtime() - MpiClass::start_wtime << " " << owner.comm_rank << arg << std::flush;
			return std::cout;
		}
	};

	logstream cout;
public:
	MpiClass(MPI_Comm comm): cout(*this){
		// init start time in 1-st call
		if(start_wtime == -1){
			start_wtime = MPI_Wtime();
		}
		this->comm = comm;
		MPI_Comm_size(comm, &comm_size);
		MPI_Comm_rank(comm, &comm_rank);
	}
	virtual ~MpiClass(){}
};

#endif /* MPICLASS_H_ */

#include <stdio.h>
#include <mpi.h>
#include <omp.h>

int main(int argc, char** argv) {
  int rank, size, thread_num, num_threads;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  #pragma omp parallel private(thread_num, num_threads)
  {
    thread_num = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    
    printf("Hello from rank %d, thread %d of %d\n", rank, thread_num, num_threads);
  }
  
  MPI_Finalize();
  
  return 0;
}

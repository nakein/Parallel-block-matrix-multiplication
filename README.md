# Parallel-block-matrix-multiplication
Block matrix multiplication using Pthreads, OpenMP and MPI

Statement:
1. Solve sequentially and using the block technique NxN square matrix multiplication:
$$C = A*B$$

Test for different block sizes in powers of 2 and determine which block maximizes performance.

2. Use the optimal block size from Exercise 1 to solve, sequentially and in parallel, the following equation:
$$R = PromP*P$$
Where :

$$P = MaxD*(ABC) + MinA*(DCB)$$

R, P, A, B, C and D are square matrices of NxN
MaxD and MinA are the maximum and minimum value of the elements of matrices D and A, respectively.
matrices D and A, respectively.
PromR is the average of the P values obtained after solving the equation for P.

Delivery guidelines:
- Do not perform mathematical simplifications for the equations.
- Exercise 1 and 2 should be tested for values of N = 1024, 2048, and 4096.
- Exercise 2 should be solved in the Shared Memory model (Pthreads and Pthreads). (Pthreads and OpenMP) and in the Distributed Memory (MPI) model.
- In shared memory run for 4 and 8 threads.
- In distributed memory run for:
• 4 cores: using 2 machines, 2 processes on each machine.
• 8 cores: using 2 machines, 4 processes on each machine.
• 16 cores: using 2 machines, 8 processes on each machine.

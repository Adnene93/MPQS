# MPQS
MPQS : The quadratic sieve algorithm (QS) is an integer factorization algorithm, MPQS is the parallel version of it, this version is programmed in C using MPI. 
#Compilation
`gcc -o factorize factorize.c -lgmp -lm`

#execute
```mpirun -np 1 --hostfile=hostfile ./factorize n```
where `np` is the number of process and `n` is the semiprime to factorize

import sys
sys.path.append("../variant_effect_analysis")

from mpi4py.futures import MPIPoolExecutor


# both of the examples run with the following setup:
# salloc --partition=normal --mem=1G --ntasks=3
## the default loaded module gnu9/9.3.0 and openmpi4/4.0.4 does not work with following mpirun command. so load as below which works.
# module load gnu10 #gnu10/10.3.0-ya
# module load openmpi # openmpi/4.1.2-4a
# mpirun -np 3 python -m mpi4py.futures mpi_practice/mpi_working_example.py

# if __name__ == '__main__':
#     executor = MPIPoolExecutor()
#     iterable = ((2, n) for n in range(32))
#     for result in executor.starmap(pow, iterable):
#        print(result)
#     executor.shutdown()

def my_f(y):
    r = pow(2, y)
    print(r)
    # return r

if __name__ == '__main__':
    executor = MPIPoolExecutor()
    executor.map(my_f, list(range(32)))
    # for result in executor.map(pow, [2]*32, list(range(32))):
    #     print(result)

    executor.shutdown()


from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

a = 0.6
b = 0.5

if rank == 0:
    print('-' * 60)
    print('a = {} and b = {}'.format(a, b))
    print('-' * 60)
    ope = 'a - b'
    res = a - b

elif rank == 1:
    ope = 'a + b'
    res = a + b

elif rank == 2:
    ope = 'a * b'
    res = a * b

else:
    ope = None
    res = None

print('Task {} out of {}: {} = {}'.format(rank, size, ope, res))
comm.barrier()


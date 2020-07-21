from mpi4py import MPI
import sys

comm = MPI.COMM_WORLD
size = comm.Get_size()

if size != 2:
    print('This example requires exactly 2 ranks')
    sys.exit(1)

rank = comm.Get_rank()

def hello_world(n):
    return 'Hello World from {}!'.format(n)

data_out = hello_world(rank)
data_in = None

if rank == 0:
    comm.send(data_out, dest=1)
    data_in = comm.recv(source=1)

elif rank == 1:
    data_in = comm.recv(source=0)
    comm.send(data_out, dest=0)

print('Task {} out of {} has sent {} and received {}'.format(rank, size, data_out, data_in))
comm.barrier()


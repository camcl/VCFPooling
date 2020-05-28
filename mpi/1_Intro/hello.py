from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def hello_world():
    return 'Hello World!'

print('This is the message from task {} out of {}: {}'.format(rank, size, hello_world()))


from mpi4py import MPI
import sys
import numpy as np

n_numbers = 1000000

# Get my rank and the number of ranks
rank = MPI.COMM_WORLD.Get_rank()
n_ranks = MPI.COMM_WORLD.Get_size()

# Check that there are exactly two ranks
if n_ranks != 2:
    print("This example requires exactly two ranks")
    sys.exit(1)

# Call the other rank the neighbour
if rank == 0:
    neighbour = 1
else:
    neighbour = 0

# Generate numbers to send
send_message = np.arange(n_numbers, dtype=int) * (rank + 1)
buff = np.empty(n_numbers, dtype=int)

# Send the message to other rank
MPI.COMM_WORLD.Isend([send_message, MPI.INT], dest=neighbour, tag=0)

# Receive the message from the other rank
req = MPI.COMM_WORLD.Irecv([buff, MPI.INT], source=neighbour, tag=0)
req.Wait()  # makes sure this communication is finished.
print("Message received by rank {}: {}".format(rank, np.max(buff)))


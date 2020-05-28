from mpi4py import MPI
import sys

n_numbers = 1000

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
send_message = []
for i in range(n_numbers):
    j = i * (rank + 1)
    send_message.append(j)

# Send the message to other rank
MPI.COMM_WORLD.isend(send_message, dest=neighbour, tag=0)

# Receive the message from the other rank
req = MPI.COMM_WORLD.irecv(source=neighbour, tag=0)
recv_message = req.wait()  # makes sure this communication is finished.
print("Message received by rank {}: {}".format(rank, max(recv_message)))


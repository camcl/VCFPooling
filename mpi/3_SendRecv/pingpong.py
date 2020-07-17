from mpi4py import MPI
import sys

maxpts = int(1e06)

comm = MPI.COMM_WORLD

# Get my rank and the number of ranks
rank = comm.Get_rank()
n_ranks = comm.Get_size()

# Check that there are exactly two ranks
if n_ranks != 2:
    print("This example requires exactly two ranks")
    sys.exit(1)

score = 0

if rank == 0:
    print('Player 0 starts the game')
    player = 0
    target = 1
    comm.send(score, dest=target)  # this line is important for starting the game

elif rank == 1:
    player = 1
    target = 0

while score < maxpts:
    score = comm.recv(source=target) + 1  # recv MUST happen before sending "the ball" back
    comm.send(score, dest=target)

print("Player {} got {} points".format(player, score))
comm.barrier()

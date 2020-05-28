# PDC/PRACE Online Course: Writing Parallel Applications Using MPI

See page at: https://pdc-support.github.io/introduction-to-mpi/index.html

### Environment preparation on Rackham at UPPMAX
Connect to Rackham cluster. Following lines should be included in the sh file to be sbatched.

Make Python available
```bash
module load python/3.6.8
```

Make openmpi available
```bash
module load openmpi/4.0.2
```
Verification test for OpenMPI environement activation
```bash
mpirun -n 4 echo Hello World!

Hello World!
Hello World!
Hello World!
Hello World!
```

Prepare working directory, access to the 1000GP scripts among others

```bash
[camcl609@rackham2 ~]$ mkdir PDC-MPI
[camcl609@rackham2 ~]$ ln -s ~/1000Genomes/scripts/VCFPooling/python/ ./1000Gpy
[camcl609@rackham2 ~]$ ln -s ~/1000Genomes/scripts/VCFPooling/mpi/ ./1000Gmpi
```

Create virtual environment

Remember for transferring files and directories:
* Directory from local to Rackham
```bash
[LOCAL] $ scp -r ~/1000Genomes/scripts/VCFPooling/mpi/ camcl609@rackham.uppmax.uu.se:/home/camcl609/1000Gmpi
```
* Directory from Rackham to local
```bash
[LOCAL] $ scp -r camcl609@rackham.uppmax.uu.se:/home/camcl609/1000Gmpi ~/1000Genomes/scripts/VCFPooling/mpi/
```

MPI-bindings installation for Python via `pip`
```bash
$ sudo apt install libopenmpi-dev  # achieves prerequisite installation
$ pip3 install mpi4py
```
NB: this works on local only (`sudo`), where `pip3` is my default version for installing Python packages, `pip` itself is not installed on my system. On Rackham:
```bash
source ~/1000Genomes/venv3.6/bin/activate
module load openmpi/4.0.2  # this provides the prerequisites on Rackham
pip install mpi4py
```


### Template for submitting a job to SLURM that will be run on devel node (hence high priority)
```bash
#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p devel
#SBATCH --qos=short
#SBATCH -n 4
#SBATCH -t 00:15:00
#SBATCH -J beagle_impute_v1.3
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
source ~/1000Genomes/venv3.6/bin/activate
module load openmpi/4.0.2
pip install mpi4py
mpirun -n 4 python3 -u myscript.py <args>
```

Note the `pip`here, not `pip3`



## 1. Introduction to Parallel Computing
One (Open)MPI version is usually always installed on HPC clusters.

What does the `mpirun` command do?
* Each process is a copy of the script/piece of code and works as if it would be on this own,
* I need to know from the beginning how many tasks we will run
* MPI will not allow fro running more tasks/ranks than the number of core on the machine ('not enough slots available' error)
* All the processes can communicate together via the MPI communicator

_Exercise_: **Hello World!**
Trial for running on login node (!!!) on Rackham
```bash
(venv3.6) [camcl609@rackham2 1_Intro]$ mpirun -n 2 python3 hello.py
--------------------------------------------------------------------------
By default, for Open MPI 4.0 and later, infiniband ports on a device
are not used by default.  The intent is to use UCX for these devices.
You can override this policy by setting the btl_openib_allow_ib MCA parameter
to true.

  Local host:              rackham2
  Local adapter:           mlx4_0
  Local port:              2

--------------------------------------------------------------------------
--------------------------------------------------------------------------
WARNING: There was an error initializing an OpenFabrics device.

  Local host:   rackham2
  Local device: mlx4_0
--------------------------------------------------------------------------
This is the message from task 0 out of 2: Hello World!
This is the message from task 1 out of 2: Hello World!
```

With `sbatch`:
See __hello_batch.sh__

_Exercise_: **Using ranks**
```bash
 mpirun -n 4 python3 using_ranks.py
```
Q: If I want to decide upon a and b at run time, I should modify my python script such that it accepts a an b as command line arguments then? 
A: Yes, or it can read it from an input file, which is common practice in scientific code 
I'll just hop over to BR4 to see if there are any questions there.


## 2. Serial and Parallel regions
Crucial to identify before parallelizing: which steps interfere with each others (not independent then?) 
Only independent processes can be easily parallelized --> atomicity of the parallelization problem.

## 3. MPI_Send and MPI_Recv
Message TAG (`tag` kwarg in Python): The message tag is used to differentiate messages, in case rank A has sent multiple pieces of data to rank B. When rank B requests for a message with the correct tag, the data buffer will be overwritten by that message.

You can never now which ranks will be terminated before the others. The script has to be written taking that into account.

Q: Does the buffer have to have exactly the memory size of the data to be sent? Can I “oversize” my buffer?
A: oversize is ok. `count` determines how many elements should be sent/received starting from the buffer’s base address.

If communication volume becomes a bottle-neck (as with NumPy), use the buffer-like communication (upper-case methods for NumPy)

_Exercise_ **blocking.py**:
This code freezes!
Possible solution: use isend and ireceive for no-blocking? Yes, it works (see __nonblocking.py__)

## 5. Parallel Paradigms and Parallel Algorithms
_Exercise_ **pingpong.py**:
Consider what rank starts and the order of `send`/`recv` requests.

## 6. Non-blocking Communication
```python
recv_message = req.wait()
```
If one prints the received message before the `wait()` method, there is no guarantee this will be the right object.

Use non-blocking communication carefully for this reason! Always make sure you get the data you want.

_Exercise_ **nonblocking.py**:
Obviously there is a limit in the size of the list mpi4py can handle: `n_numbers = 1000` works, but `n_numbers = 10000` and above fail. Get 'EOFError: Ran out of input' or even segmentation fault from openmpi. 
Implementing a NumPy version with buffered message exchanges (see **nonblocknumpy.py**) solves the issue

## Procedures for message passing
Rank 0 will read the data for example, and dispatch it to all the other ranks that will process it, and then rank 0 gathers back the result


## Brainstorming
* Project for using MPI framework in the 1000GP research topic: compute the metrics on the final files (reduce operations might be interesting to bin by MAF?)
* Process more tasks than the number of cores available: read about __Embarrassingly Parallel Problems__,  Pools and Workers https://mpi4py.readthedocs.io/en/stable/mpi4py.futures.html
import multiprocessing as mp
import numpy as np

print('# of CPU -->', mp.cpu_count())

# Drawing without replacement
a = np.arange(1, 100)
n = len(a)//16
print('\n# of sets to draw: ', n)
set  = np.random.choice(a, size =(6, 16), replace=False)
print('\n6 sets of 16 samples: \n', set)
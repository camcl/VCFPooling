from src.VCFPooling.poolSNPs.pooler import Encoder
from src.VCFPooling.python.archived.cydecoder import *
import timeit

"""
Using cydecoder.pyx
Not very promising... cythonizatio makes a fuss in the code
"""

np.random.seed(123)

## Encoding
nb_blocks = 1200
var = np.random.binomial(1, 0.85, (16*nb_blocks, 2))
dse = Design(blocks=nb_blocks)
dme = dse.matrix
enc = Encoder(dme)
varp = enc.encode(var.sum(axis=-1))
# print(var.sum(axis=-1).reshape((nb_blocks, 4, 4)))

## Decoding
path_keys = '/home/camille/PoolImpHuman/data/main'
keys, vals = get_lookup_arrays(path_keys)

dsd = Design()
dmd = dsd.matrix
dec = SingleBlockDecoder(dmd, keys, vals, 'gp')

vari = np.empty_like(var)
x_shift = dmd.shape[1]  # 16
y_shift = dmd.shape[0]  # 8


def forloop_decoding():
    q = blocks_decoder(nb_blocks, varp.sum(axis=-1), y_shift, keys, vals, 'gp')
    return q


def npapply_decoding():
    vect_decoder = lambda x: single_block_decoder(x, keys, vals, 'gp')
    v2 = varp.sum(axis=-1).reshape((nb_blocks, 8))
    q = np.apply_along_axis(vect_decoder, axis=-1, arr=v2)
    return q


# res = forloop_decoding()
# res = npapply_decoding()
# print(res)

print(timeit.timeit("res = npapply_decoding()",
                    setup="from __main__ import forloop_decoding, npapply_decoding",
                    number=10))

times = [('decoding_method', 'cython_python', 'genotype_format', 'nb_blocks', 'timing_rep', 'time'),
         ('forloop_decoding', 'pure python', 'gp', 12, 10, 0.18752734300051088),
         ('forloop_decoding', 'pure python', 'gp', 120, 10, 1.9894000200001756),
         ('forloop_decoding', 'pure python', 'gp', 1200, 10, 18.713707934000013),
         ('npapply_decoding', 'pure python', 'gp', 12, 10, 0.22564055099974212),
         ('npapply_decoding', 'pure python', 'gp', 120, 10, 2.078016972000114),
         ('npapply_decoding', 'pure python', 'gp', 1200, 10, 22.500562958000046),
         ('forloop_decoding', 'simple cython', 'gp', 12, 10, 0.20502520899935917),
         ('forloop_decoding', 'simple cython', 'gp', 120, 10, 1.7975935060003394),
         ('forloop_decoding', 'simple cython', 'gp', 1200, 10, 19.594972976999998),
         ('npapply_decoding', 'simple cython', 'gp', 12, 10, 0.20081693999964045),
         ('npapply_decoding', 'simple cython', 'gp', 120, 10, 1.9735721830002149),
         ('npapply_decoding', 'simple cython', 'gp', 1200, 10, 20.53694264600017),
         ]

import allel
import os
import time
import warnings
from scipy.stats import bernoulli as bn
import numpy as np
import pandas as pd
from cyvcf2 import VCF, Writer
import collections

warnings.simplefilter('ignore')

print(os.getcwd())


def hello_world(file_data):
    print_lock = threading.Lock()
    with print_lock:
        print("Starting thread : {}".format(threading.current_thread().name))
    print('Hi, here is the object')
    frame = inspect.currentframe()
    with print_lock:
        print("Finished thread : {}".format(threading.current_thread().name))
    return frame.f_locals['file_data']


# idv1_id = allel.read_vcf('./data/ALL.chr20.phase3_GL_subset_1.vcf.gz')['samples']
# head1 = allel.read_vcf('./data/ALL.chr20.phase3_GL_subset_1.vcf.gz', fields=['variants/*'])
data1 = VCF('./data/tests-beagle/ALL.chr20.missing.snps.gt.chunk.vcf.gz')
# meta1 = allel.read_vcf('./data/ALL.chr20.phase3_GL_subset_1.vcf.gz')['headers']
imp = VCF('./data/tests-beagle/IMP.chr20.missing.beagle2.corr.vcf.gz')

# Import existing header
##FORMAT=<ID=GL,Number=3,Type=Float,Description="three log10-scaled likelihoods for RR,RA,AA genotypes">

start = time.time()
stop = time.time()
#print('Processing duration [s]: ', stop-start)

# def f(v,w):
#     w.write_record(v)

# print('SNP #{}'.format(variant.ID))

# for item in itertools.starmap(f, iter(zip(data, itw))):
#     print(item)

# it1 = iter(zip(range(5), 'abcde'))
# it2 = iter(zip(data, itw))
# it3 = list(zip(data, itw))


def minor_allele(arr): # nan or -1?
    if np.sum(arr) == 2:
        return 2
    elif np.isin(arr, 1).any():
        return 1
    else:
        return 0


def compute_maf(vcf_path, verbose=False):
    """

    :param vcf_obj:
    :return:
    """
    print('Computing MAFs from {}'.format(vcf_path).ljust(80, '.'))
    vcf_obj = VCF(vcf_path)
    dic_maf = {}
    for i, var in enumerate(vcf_obj):
        # GT:DS:GP
        if verbose:
            try:
                print(str(i) + '\t' + var.ID + '\t' + var.format('GP'))
            except:
                print(str(i) + '\t', var.ID)
        gt = np.array(var.genotypes)[:, :-1]
        dic_maf[var.ID] = np.sum(np.apply_along_axis(minor_allele, 1, gt))/(len(var.genotypes)*2)
    #     if var.ID == 'rs1000554':
    #         print(str(i) + '\t', gt.shape)
    # print(len(list(dic_maf.items())))
    return dic_maf


def compute_maf_evol(set, df):
    """

    :param set: str, short name of the data set pooled/missing...
    :param df_maf: maf from the original vcf-chunked file, with markers ID as index
    :return:
    """
    print('Set --> ', set.upper())
    o = read_queue(set)
    d = {**o[0], **o[1]}
    postimp_maf = d['postimp']
    preimp_maf = d['preimp']
    # <tests>
    print('postimpmaf --> ', len(list(postimp_maf.items())))
    print('postimpmaf --> ', len(list(preimp_maf.items())))
    print('DF indices coincide?: ', sorted(list(postimp_maf.keys())) == sorted(list(df.index)))
    # </tests>
    dic = dict()
    #for i, k in enumerate(df.index):
        #dic['preimp_' + set], dic['postimp_' + set] = collections.OrderedDict(), collections.OrderedDict()
    dic['preimp_' + set] = preimp_maf # caution, still missing data!
    dic['postimp_' + set] = postimp_maf
    return dic


import threading
from queue import Queue
import inspect
from scripts.poolSNPs.alleles import alleles_tools as alltls

def readlock_vcf(file_data):
    vcf = {}
    k, v = list(file_data.items())[0]
    vcf[k] = compute_maf(v, verbose=False)
    frame = inspect.currentframe()
    return frame.f_locals['vcf']


def process_queue(compress_queue, yield_queue, max_sizeout):
    out = None
    while yield_queue.qsize() < max_sizeout:
        file_data = compress_queue.get()
        out = readlock_vcf(file_data)
        compress_queue.task_done()
        yield_queue.put(out)
        yield_queue.task_done()
    return yield_queue


def read_queue(s): # 2 queues: 1 for pushing input, 1 for getting results
    q_in, q_out = Queue(), Queue()

    out_t = []
    output_names = [{'preimp' : 'IMP.chr20.{}.snps.gt.chunk.vcf.gz'.format(s)},
                    {'postimp': 'IMP.chr20.{}.beagle2.vcf.gz'.format(s)}]

    threads = []
    for i in range(1):
        t = threading.Thread(group=None,
                             target=process_queue,
                             name='thread_' + str(i+1),
                             args=(q_in, q_out, len(output_names),))
        t.setDaemon(True)
        t.start()
        threads.append(t)

    for file_data in output_names:
        q_in.put(file_data)
    q_in.join()

    for file_data in output_names:
        out_t.append(q_out.get())

    for t in threads:
        t.join()

    return out_t

os.chdir('./data/tests-beagle')
set = 'missing'
# o = read_queue(set)
# d = {**o[0], **o[1]}
# print(list(d.keys()))

df_maf = alltls.get_maf('ALL.chr20.snps.gt.chunk.vcf.gz').set_index('id')
dic_maf = compute_maf_evol(set, df_maf)
print(dic_maf)
# for k, v in dic_maf.items():
#     print(k, v)

"""
bcftools query -f '%POS\t%ID\t[ %GT]\n' IMP.chr20.beagle1.vcf.gz | awk '$0="line"NR": "NF' | grep '^252'
Something wrong using reheader after opening the file in the text editor and modifying it by hand.
Re-generate all files with Beagle.
Try sed cmd to edit the header in place: 
$ bcftools view IMP.chr20.missing.beagle2.vcf.gz | sed 's/##FORMAT=<ID=DS,Number=A/##FORMAT=<ID=DS,Number=1/' | sed 's/##FORMAT=<ID=GP,Number=G/##FORMAT=<ID=GP,Number=1/' | bcftools view -Oz -o IMP.chr20.missing.beagle2.vcf.gz
index file then
"""
import os
import threading
from queue import Queue
import inspect
from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs import parameters as prm

chk_sz = prm.CHK_SZ

def readlock_vcf(file_data):
    vcf = {}
    k, v = list(file_data.items())[0]
    vcf[k] = alltls.compute_maf(v, verbose=False)
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
    output_names = [{'preimp' : 'IMP.chr20.{}.snps.gt.chunk{}.vcf.gz'.format(s, str(chk_sz))},
                    {'postimp': 'IMP.chr20.{}.beagle2.chunk{}.corr.vcf.gz'.format(s, str(chk_sz))}]

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

if __name__ == '__main__':
    os.chdir('/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
    set = 'missing'
    o = read_queue(set)
    d = {**o[0], **o[1]}
    print(list(d.keys()))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datasketch import MinHash, MinHashLSH
from collections import Counter
import os

def plot_scurve(b, r):
    xvalues = []
    match = []
    falseP = []
    falseN = []
    for sim in np.arange(0.8, 1.01, 0.01):
        #calculate P(sharing bucket)
        m = 1-((1-sim**r)**b)
        N = (1-sim**r)**b
        P = (sim**r)*b
        xvalues.append(sim)
        match.append(m)
        falseN.append(N)
        if sim<0.98:
            falseP.append(P)
        else:
            falseP.append(None)
    data = pd.DataFrame({"xvalues":xvalues, "P(match)":match, "P(false+)":falseP, "P(false-)":falseN})
    plt.plot( 'xvalues', 'P(match)', data=data)
    plt.plot( 'xvalues', 'P(false+)', data=data)
    plt.plot( 'xvalues', 'P(false-)', data=data)
    plt.legend()
    plt.show()

def optimize_bands(threshold, num_perm, weights=(0.5,0.5), plot=False):
  #use datasketch's optimization function, it's good
  lsh = MinHashLSH(threshold=thresh, num_perm=nh, weights=weight)
  nbands = lsh.b
  nrows = lsh.r
  num_minhash = lsh.b*lsh.r
  if plot:
      plot_scurve(b=nbands, r=nrows)
  return nbands, nrows, num_minhash

 #fn to extract unique kmers from a read/sequence/string
def kmerize(sequence, k):
    kmers = set()
    size = len(sequence)
    for i in range(size-k+1):
        kmers.add(sequence[i:i+k])
    return list(kmers)

#fn to create minhash sketch from list unique kmers
def make_minhash(kmers, num_perm):
    sketch = MinHash(num_perm=num_perm)
    for kmer in kmers:
        sketch.update(kmer.encode('utf8'))
    return sketch

def get_existing_cluster(NNs):
    #if cluster counts are equal counter returns I think the first alphabetically (so lowest number in this case)
    cluster = str(Counter(key[1] for key in NNs).most_common(1)[0][0])
    return cluster

def write_read(cluster, read, header=None, path='clusters/',mode='a'):
    with open(path+cluster+'.clust', mode) as cluster_disk:
        if header:
            cluster_disk.write(header+'\n')
        cluster_disk.write(read+'\n')

if __name__ == "__main__":
    cluster_path='clusters/'
    if not os.path.exists(cluster_path):
        os.makedirs(cluster_path)

    #nbands, nrows, num_minhash = optimize_bands(threshold=0.97, num_perm=200, weights=(0.5,0.5), plot=True)
    #print "optimal number of bands, rows per band, and minhash functions:", nbands, nrows, num_minhash
    nbands = 3
    nrows = 66
    num_minhash = 198

    kmer_size = 31

    #initialize lsh
    lsh = MinHashLSH(num_perm=num_minhash, params=(nbands, nrows))

    test_reads = ['GCTCCGAACTCTCGCGGCCGGCCACGATCTCGGCGATTCGGCGGGCGGTGACGCCCGGATCACCACCGGCCAGTGCCCATCCCCTGGCCCGTACCAGATATGCGGTCGCAGCGGCCTCCGCCGAAACACCCAGGCCGCCTTGGATGTGCACCGCCATGGTGGCCGCCTTCGCGGCCTCTTCGGCCATGAAGACGAACGCC', 'GCTCCGAACTCTCGCGGCCGGCCACGATCTCGGCGATTCGGCGGGCGGTGACGCCCGGATCACCACCGGCCAGTGCCCATCCCCTGGCCCGTACCAGATATGCGGTCGCAGCGGCCTCCGCCGAAACACCCAGGCCGCCTTGGATGTGCACCGCCATGGTGGCCGCCTTCGCGGCCTCTTCGGCCATGAAGACGAACGCC', 'GCTCCGAACGGTCGCGGCCGGTTTGATCTCTTTGATTCGGTATACGGTGACGATCGGATCACCAGGCCAGTCCATATCCTGGCCCTACTATGGCCAGCGGCCTCTTTCGAAACACCCAGCCCGCCTTGGATGTGCACCGCCATGGTGGCCGCCTTCGCGGCCTCTTCGGCCATGAATACGAACGCC']

    maxCluster = 0
    for ix, read in enumerate(test_reads):
        read_key = "read"+str(ix)
        #extract kmers
        ks = kmerize(read, k=kmer_size)
        #create minhash sketch from unique kmers
        sketch = make_minhash(kmers=ks, num_perm=198)
        #query lsh to get nearest neighbors, if any
        NNs = lsh.query(sketch)
        #if no NNs, assign to maxCluster
        if len(NNs)==0:
            clust = str(maxCluster)
            maxCluster += 1 #increment maxCluster so next new cluster will have new assignment
            #insert read in lsh
            lsh.insert((read_key, clust), sketch)
            #write to clust disk as new file
            write_read(cluster=clust, read=read, header='header for '+read_key, path=cluster_path,mode='w')
        if len(NNs)>0:
            #find cluster value as most frequent cluster seen in NNs
            clust = get_existing_cluster(NNs)
            #insert read in lsh
            lsh.insert((read_key, clust), sketch)
            #append to existing clust file
            write_read(cluster=clust, read=read, header='header for '+read_key, path=cluster_path,mode='a')

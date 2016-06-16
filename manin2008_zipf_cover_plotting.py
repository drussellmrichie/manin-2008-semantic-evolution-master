"""
June 8, 2014 -- Russell Richie

This code uses multiprocessing to plot the Zipfian Covering of several Manin 
(2008) simulations and simultaneously.

NOTO BENE -- For some reason, Canopy will not run the model in a pool. Instead,
this script has to be called from the terminal.
"""

import sys, os, itertools, pickle, copy
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
#from matplotlib import rc
from multiprocessing import Pool

current_directory = '/Users/russellrichie/introcompling/FinalProject/ManinExploration/Data_GenMode_IntNum-10000_Delta-relative/'
os.chdir(current_directory)

def getOverlap(a, b):
    """
    Get the overlap between two intervals.
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))    

def unite_intervals(data):
    """
    Compute union of intervals
    """
    data = sorted(data)
    it = iter(data)
    a, b = next(it)
    for c, d in it:
        if b > c:  # Use `if b >= c` if you want (1,2), (2,3) to be
                    # treated as intersection.
            b = max(b, d)
        else:
            yield a, b
            a, b = c, d
    yield a, b    

def overlap_of_all_intervals(interval_edges): #this may be (very?) inefficient...this actually isn't the way manin calculated overlap
    """
    This computes the sum of the overlap between all pairs of intervals, necessary
    for the analysis of Zipfian Covering in Manin (2008).
    """
    totaloverlap = 0
    trimmed_interval_edges = copy.deepcopy(interval_edges)
    for index, first_interval in enumerate(interval_edges):
        #sys.stderr.write('Interval {} of {}...\n'.format(index,len(interval_edges)))
        trimmed_interval_edges.remove(first_interval)
        for second_interval in trimmed_interval_edges:
            totaloverlap += getOverlap(first_interval, second_interval)
    return totaloverlap

def plot_zipfian_covering(filename, rho = 1.1,gap_or_overlap = 'gap',):
    """
    This generates a figure like Fig 6 in Manin (2008), which demonstrates the 
    Zipfian Covering-ness of his generalization and specialization models.
    
    Rho is supposed to be chosen such that the sum of the interval lengths from rank
    k to rho*k is equal to 1.
    
    parameters may be important -- I obtain gap-> 0 as k-> inf, but not replicating 
    power law so far....
    
    Not even getting gap->0 as k-> inf for spec model, and I'm using Manin's exact
    parameters...
    """
    # should I write a function that finds rho such that sum(interval_sizes[k:rho*k]) = 1 ?
    data = pickle.load( open( filename, "rb" ) )
    
    zipped_data = zip(data[0],data[1])
    sorted_data = sorted(zipped_data,key=lambda x:x[1], reverse=True)
    interval_edges = [edges for edges, sizes in sorted_data]
    
    ending_rank = int(len(interval_edges)/rho)
    
    if gap_or_overlap == 'gap':
        for k in range(4,ending_rank): # starting at 4 is, for now, my stupid way of circumventing the fact that 0 is the first index in a list, which messes up k->rho*k
            if k % 50 == 0:
                sys.stderr.write('Interval {} of {}...\n'.format(k,len(interval_edges)))
            curr_union = list(unite_intervals(interval_edges[k:int(k*rho)])) #inspect this....
            curr_gap = 1 - sum(y-x for x, y in curr_union)
            plt.loglog(k,curr_gap, color = 'k', marker='8',markersize=4)
    elif gap_or_overlap == 'overlap':
        for k in range(4,ending_rank): # starting at 4 is, for now, my stupid way of circumventing the fact that 0 is the first index in a list, which messes up k->rho*k
            if k % 50 == 0:
                sys.stderr.write('Interval {} of {}...\n'.format(k,len(interval_edges)))
            curr_overlap = overlap_of_all_intervals(interval_edges[k:int(k*rho)]) #this actually appears not to be the right function for calculating overlap....
            plt.loglog(k,curr_overlap, color = 'k', marker='8',markersize=4)
    title = 'Zipfian_Covering_Manin_GenModel_IntNum-10000_Delta-relative_Rho-1point1_Gap_Or_Overlap-{}'.format(gap_or_overlap)
    plt.title(title)
    #rc('text', usetex=True)
    #plt.xlabel('Starting Rank, \textit{k} (log)') # I think this gives a string index error
    plt.xlabel('Starting Rank, k (log)')
    if gap_or_overlap == 'overlap':
        plt.ylabel('Overlap (log)')
    else:
        plt.ylabel('Gap (log)')
    plt.savefig(title)
    return None

# Set parameters for plotting function

rho = [1.1]
gap_or_overlap = ['gap']
filenames = onlyfiles = [ f for f in listdir(current_directory) if isfile(join(current_directory,f)) and f[-1] == 'p']

graphing_parameter_list = list(itertools.product(filenames,rho,gap_or_overlap))
graphing_parameter_list = [{'filename':a,'rho':b,'gap_or_overlap':c} for a, b, c in graphing_parameter_list]

def graphing_pool_director(kwargs): # kwargs will be a dict of arguments
    return plot_zipfian_covering(**kwargs)

if __name__ == '__main__':
    pool = Pool()
    figure_num = 0
    for plot in pool.imap_unordered(graphing_pool_director,graphing_parameter_list):
        continue
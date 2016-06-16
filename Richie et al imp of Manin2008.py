# -*- coding: utf-8 -*-

"""
This is code reimplementing Manin (2008)'s generalization and specialization
models. Russell Richie originally wrote the reimplementations, and then Stefan
Kaufmann refined them.

First run up to ~line 245 to define all relevant functions. Then run the last 
lines to simulate models and plot results.
"""

import random, sys
#import os
import numpy as np
import matplotlib.pyplot as plt

#os.chdir('/Users/russellrichie/introcompling/FinalProject/')

def getOverlap(a, b):
    """
    Get the overlap between two intervals.
    """
    return max(0.0, min(a[1], b[1]) - max(a[0], b[0]))    

def manin_gen_model(interval_numb = 5000, delta = 'relative'):
    """
    The generalization model from Manin (2008) "Zipf's Law and
    Avoidance of Excessive Synonymy"

    This generates a power-law distribution of interval sizes, and
    supposedly a Zipfian covering as well, but I am having trouble
    replicating this (May 19, 2014).
    
    """    
    if delta == 'relative':
        delta = .01 / interval_numb
    
    unfrozen = [[random.uniform(0,1)]*2 for x in range(interval_numb)]
    unfrozen.sort()
    frozen = []

    # This freezes intervals as soon as it finds them during the loop.
    # That may not be the right approach - those same intervals may
    # have overlap with other non-frozen intervals later in the loop, 
    # which would not be found.
    #
    # while len(unfrozen) > 1:
    #     for i in range(len(unfrozen)):
    #         unfrozen[i][0] = max(unfrozen[i][0] - delta, 0)
    #         unfrozen[i][1] = min(unfrozen[i][1] + delta, 1)
    #     i = 0
    #     while i < len(unfrozen):
    #         j = i+1
    #         while j < len(unfrozen):
    #             if getOverlap(unfrozen[i],unfrozen[j]):
    #                 if random.choice((1,0)):
    #                     frozen.append(unfrozen.pop(i))
    #                     j = i+1
    #                 else:
    #                     frozen.append(unfrozen.pop(j))
    #             else:
    #                 break
    #         i += 1
    # frozen.append(unfrozen.pop(0))
    
    # This freezes intervals only at the end of each cycle, thus 
    # addressing the above problem. The result looks pretty much 
    # the same, though...    
    while len(unfrozen) > 1:
        to_freeze = []
        for i in range(len(unfrozen)):
            unfrozen[i][0] = max(unfrozen[i][0] - delta, 0)
            unfrozen[i][1] = min(unfrozen[i][1] + delta, 1)
        for i in range(len(unfrozen)):
            for j in range(i+1,len(unfrozen)):
                if getOverlap(unfrozen[i],unfrozen[j]):
                    sys.stderr.write('{} intervals frozen.\n'.format(len(frozen)))
                    newly_frozen_interval = random.choice((i,j))
                    if newly_frozen_interval not in to_freeze:
                        to_freeze.append(newly_frozen_interval)
                    #to_freeze.append(random.choice((i,j)))
                else:
                    break
        to_freeze.reverse()
        for i in to_freeze:
            frozen.append(unfrozen.pop(i))
    frozen.append(unfrozen.pop(0))

    # Sort by size
    frozen.sort(key=lambda i : i[1]-i[0], reverse=True)
    return frozen

def manin_spec_model(interval_numb = 5000, gamma = 2.0, cutoff = True):
    """
    
    I think this could be sped up if you only looked at the interval pairs where
    at least one member changed on the previous round? Unless that requires a slew
    of if-thens....
    
    The specialization model from Manin (2008) "Zipf's Law and Avoidance of 
    Excessive Synonymy". Manin used interval_numb = 100000 and gamma = 1.1.

    This generates a power-law distribution of interval sizes, and supposedly a
    Zipfian covering as well.
    
    Manin's description of this model was slightly unclear on the initialization
    of the intervals, so I implemented a cut-off on the interval edges at 1 and 0.
    
    I was also not sure if I should have decreased the lengths on either side of
    the interval, or only on the intersected side...I opted for the former...
    
    My implementation of initialization of intervals, with the cutoff I added, 
    is *not* generating uniform distribution of interval lengths. Omitting cutoff
    makes zipfian covering plot even weirder, it seems, at least for 'gap'. 
    
    Cutting half of overlap from each side of smaller interval appears to lead 
    to never ending while loop, which I think is caused by a Zeno's paradox type 
    situation...?
    """ 
    centers = [random.uniform(0,1) for x in range(interval_numb)]
    lengths = [random.uniform(0,1) for x in range(interval_numb)]
    if cutoff == True:
        intervals = [ [max(0,centers[x]-lengths[x]/2),
                       min(1,centers[x]+lengths[x]/2)] 
                      for x in range(interval_numb)]
    else:
        intervals = [ [centers[x]-lengths[x]/2,centers[x]+lengths[x]/2] 
                      for x in range(interval_numb)]
    intervals.sort()

    overlap = True
    loop_counter = 0
    while overlap == True:
        overlap = False
        loop_counter += 1
        sys.stderr.write('While loop #{}...\n'.format(loop_counter))
        resize_counter = 0

        for i in range(len(intervals)):
            for j in range(i+1,len(intervals)):
                ovl = getOverlap(intervals[i],intervals[j])
                if ovl:
                    lengths = (intervals[i][1]-intervals[i][0],
                               intervals[j][1]-intervals[j][0])
                    if 1.0/gamma < lengths[0]/lengths[1] < gamma:
                        overlap = True
                        if lengths[0] < lengths[1]:
                            intervals[i][0] += ovl/2
                            intervals[i][1] -= ovl/2
                            sys.stderr.write('interval {} resized by {}.\n'.\
                                                 format(i,ovl))
                        else:
                            intervals[j][0] += ovl/2
                            intervals[j][1] -= ovl/2
                            sys.stderr.write('interval {} resized by {}.\n'.\
                                                 format(j,ovl))
                        resize_counter += 1
                # else:
                #     break
        sys.stderr.write('{} intervals resized.\n'.format(resize_counter))
                        
    intervals.sort( key=lambda i:i[1]-i[0], reverse=True)
    return intervals

def overlap(intervals):
    my_intervals = sorted(intervals)
    intersections = []
    for i in range(len(my_intervals)):
        for j in range(i+1, len(my_intervals)):
            max_start = max(my_intervals[i][0],my_intervals[j][0])
            min_end   = min(my_intervals[i][1],my_intervals[j][1])
            if max_start < min_end:
                intersections.append([max_start,min_end])
            else:
                break
    return sum(y-x for x,y in unite_intervals(intersections))

def unite_intervals(intervals):
    """
    Compute union of intervals
    """
    it = iter(sorted(intervals))
    a, b = next(it)
    for c, d in it:
        if b > c:  # Use `if b >= c` if you want (1,2), (2,3) to be
                    # treated as intersection.
            b = max(b, d)
        else:
            yield a, b
            a, b = c, d
    yield a, b    

def plot_interval_sizes(intervals):
    """
    Plot interval size rank vs interval size from Manin (2008) models
    on a log-log plot.
    Assumes that the intervals are sorted by size (descending)
    """
    plt.loglog(range(len(intervals)), [i[1]-i[0] for i in intervals] )
                                    
def plot_zipfian_covering(intervals, rho = 2.0, gap_or_overlap = 'gap'):
    """
    This generates a figure like Fig 6 in Manin (2008), which
    demonstrates the Zipfian Covering-ness of his generalization and
    specialization models.
    
    Rho is supposed to be chosen such that the sum of the interval
    lengths from rank k to rho*k is equal to 1.
    
    """
    # should I write a function that finds rho such that sum(interval_sizes[k:rho*k]) = 1 ?
    
    sums = []
    for k in range(1,int(len(intervals)/rho),10):
        if (k-1) % 500 == 0:
            sys.stderr.write('k = {}, k*rho = {}...\n'.format(k,k*rho))
        sums.append(sum(y-x for x,y in 
                        unite_intervals(intervals[k-1:int(k*rho)-1])))

        if gap_or_overlap == 'gap':
            curr_gap = 1 - sums[-1]
            plt.loglog(k,curr_gap, color = 'k', marker='8',markersize=4)
            
        elif gap_or_overlap == 'overlap':
            curr_overlap = overlap(intervals[k-1:int(k*rho)-1])
            plt.loglog(k,curr_overlap, color = 'k', marker='8',markersize=4)

    sys.stderr.write('Sums: mean = {}; std = {}.\n'.format(
            np.array(sums).mean(), np.array(sums).std()))

    plt.xlabel('Starting Rank, k (log)')
    if gap_or_overlap == 'overlap':
        plt.ylabel('Overlap (log)')
    else:
        plt.ylabel('Gap (log)')
    return plt.show()

"""

After defining functions above, run the blocks below to simulate the
generalization and specialization models, plot their resulting
interval sizes, and plot the resulting (lack of) Zipfian covering.

"""

gen_data = manin_gen_model()
plot_interval_sizes(gen_data)
plot_zipfian_covering(gen_data, gap_or_overlap='gap')

# spec_data = manin_spec_model()
# plot_interval_sizes(spec_data)
# plot_zipfian_covering(spec_data, gap_or_overlap='overlap')

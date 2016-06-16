"""
README

June 7, 2014 -- Russell Richie

This code uses multiprocessing to run several Manin (2008) simulations and 
plotting of said simulations simultaneously

NOTO BENE -- For some reason, Canopy will not run the model in a pool. Instead,
this script has to be called from the terminal.

Another script will use multiprocessing to plot the zipfian coverings of several
data sets simultaneously.
"""
import random, copy, os, sys, pickle, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import itertools
from multiprocessing import Pool

os.chdir('/Users/russellrichie/introcompling/FinalProject/manin-2008-semantic-evolution')

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def manin_gen_model(interval_numb = 1000, delta = 'relative', run_num = None, sim_num = 1):
    print 'Starting simulation {} of {}'.format(run_num,sim_num)
    if delta == 'relative':
        delta = .01 / interval_numb
    unfrozen_interval_edges = [[random.uniform(0,1)]*2 for x in range(interval_numb)]
    finished_interval_edges = []
    while len(unfrozen_interval_edges) > 1:
        for int_index, interval in enumerate(unfrozen_interval_edges):
            unfrozen_interval_edges[int_index][0] = max(interval[0] - delta, 0)
            unfrozen_interval_edges[int_index][1] = min(interval[1] + delta, 1)              
        for interval_edges in unfrozen_interval_edges:
            other_unfrozen_interval_edges = copy.deepcopy(unfrozen_interval_edges)
            other_unfrozen_interval_edges.remove(interval_edges)
            for other_interval_edges in other_unfrozen_interval_edges:
                if getOverlap(interval_edges,other_interval_edges):
                    to_be_frozen_edges = random.choice((interval_edges,other_interval_edges))
                    finished_interval_edges.append(to_be_frozen_edges)
                    unfrozen_interval_edges.remove(to_be_frozen_edges)
                    sys.stderr.write('Interval {} of {} frozen on run_num {}...\n'.format(interval_numb - len(unfrozen_interval_edges),interval_numb,run_num))
                    if to_be_frozen_edges == interval_edges:
                        break
    finished_interval_edges.append(unfrozen_interval_edges.pop(0))
    interval_sizes = [x - y for y,x in finished_interval_edges]
    return (finished_interval_edges, interval_sizes)

def plot_interval_sizes(interval_sizes):
    """
    Plot interval size rank vs interval size from Manin (2008) models on a log-log plot
    """
    plt.loglog(range(len(interval_sizes)),sorted(interval_sizes,reverse=True))            
    
runs_per_setting = range(10)            
interval_numb = [10000]
delta = ['relative']

model_parameter_list = list(itertools.product(interval_numb,delta,runs_per_setting))
model_parameter_list = [{'interval_numb':a,'delta':b,'run_num':c, 'sim_num': len(model_parameter_list)} for a, b, c in model_parameter_list]

def model_pool_director(kwargs): # kwargs will be a dict of arguments
    return manin_gen_model(**kwargs), kwargs

if __name__ == '__main__':  
    pool = Pool() # this will just use all available cores 
    simulation_num = 0 # this uniquely identifies a sim from all the sims run in model_parameter_list, cf. run_num which ID's it for a particular param combo
    for data, kwargs in pool.imap_unordered(model_pool_director, model_parameter_list):            
        simulation_num += 1
        print 'Finishing simulation {} of {}'.format(simulation_num,len(model_parameter_list)) # using counter because couldnt' get enumerate to work here...
        
        interval_numb , delta , run_num = kwargs['interval_numb'] , kwargs['delta'] , kwargs['run_num']
        data_filename = 'Data_GenModel_IntNum-{}_Delta-{}_Run-{}.p'.format(interval_numb,delta,run_num)
        pickle.dump(data,open(data_filename,'wb')) #str object has no attribute 'write'
        
        #plt.figure(1)
        #plot_interval_sizes(data[1])
        #title = 'Interval_Sizes_GenModel_IntNum-{}_Delta-{}_Run-{}'.format(interval_numb,delta,run_num)
        #plt.title(title)
        #plt.xlabel('Interval Size Rank (log)')
        #plt.ylabel('Interval Size (log)')
        #plt.savefig(title)
        #plt.clf()
        
        #plt.figure(2)
        plot_interval_sizes(data[1])
        title = 'Interval_Sizes_GenModel_IntNum-{}_Delta-{}_{}_Runs'.format(interval_numb,delta,len(model_parameter_list))
        plt.title(title)
        plt.xlabel('Interval Size Rank (log)')
        plt.ylabel('Interval Size (log)')
        plt.savefig(title)
        #plt.clf() # deletes the current figure but not the window, since the next runs will open it 
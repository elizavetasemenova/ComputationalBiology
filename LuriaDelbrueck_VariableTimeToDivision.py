"""
        
    The Luria-Delbruck Experiment with variable cell cycle durations.

"""

import numpy as np
from scipy.stats import bernoulli

try:
    import Queue as Q  # Python version < 3.0
except ImportError:
    import queue as Q

# Approximate final population size
N = int(1e+6)
#N = 10

# mutation rate
mu = 1e-6
#mu = 0.5

# number of colonies
n_colonies = 100

# Input: 
#    N - final population size
#    mu - mutation rate
#    n_colonies - number of colonies
# Output: 
#    r_colonies - number of resistant bacteria in the final population for each colony
def cell_cycle_sim(N, mu, n_colonies):

    r_colonies = np.zeros(n_colonies) # number of resistant bacteria in each colony
    
    # Compute number of resistants for each colony
    for colony in range(n_colonies): 

        #=========================================================
        # Create priority queue of tuples (priority, resistance_status).
        # Priority is the time of division of bacteria:
        #     t_division = t_birth + life_time.
        # t_birth is t_division of the parent,
        # life_time is random and distributed as exp(random.lognormal(log(18), log(1.2))).
        # resistance_status is either 0 (wild type), or 1 (mutant/resistant).
        #=========================================================
        
        # Create queue
        q = Q.PriorityQueue()
        
        # Initialize queue
        t = 0 # current clock time
        priority = t + int(round(np.random.lognormal(mean=np.log(18), sigma=np.log(1.2)),2))
        resistance_status = 0
        q.put((priority, resistance_status))
        print q.queue
        
        while q.qsize()<N:
            
            (t, parent_resistance_status) = q.get()
         
            # generate two bacteria out of one
            for _ in range(2):
                priority = t + int(round(np.random.lognormal(mean=np.log(18), sigma=np.log(1.2))))
                if resistance_status == 0:
                    mutate = bernoulli.rvs(mu, size=1)
                    mutate
                    if mutate ==1:
                        child_resistance = 1
                    else:
                        child_resistance = 0
                else:
                    child_resistance = 1
                q.put((priority, child_resistance))
            #print q.queue
            print q.qsize()

        # Count mutants in the progeny
        r = 0
        while not q.empty():
            tup = q.get()
            r += tup[1]
        print r
        r_colonies[colony] = r

    return r_colonies

r_cell_cycle = cell_cycle_sim(N, mu, n_colonies)
                  
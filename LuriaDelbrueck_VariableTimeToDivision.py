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

#=========================================================
# Hypothesis 1 (H1): 
#     ~~~~~~~ Acquired immunity ('adaptation') ~~~~~~~
#     there is no mutation during the growth of bacteria,
#     but there is a small probability, e.g., 10−6, 
#     for the bacteria to become resistant upon
#     contact with the virus.
#=========================================================


# Input: 
#    N - final population size
#    mu - mutation rate
#    n_colonies - number of colonies
# Output: 
#    r_colonies - number of resistant bacteria in the final population for each colony
def adaptation_sim(N, mu, n_colonies):

    # Binomial number of resistant bacteria
    r_colonies = np.random.binomial(N, mu, n_colonies)
    
    return r_colonies


#====================================================
# Hypothesis 2 (H2): 
#     ~~~~~~~ Mutation hypothesis  ~~~~~~~
#     mutations in bacteria take place independenly 
#     of the presence of phage, i.e. before
#     the introduction of phage, as well as after.
#====================================================

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
        
        #print colony

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
        #priority = t + int(round(np.random.lognormal(mean=np.log(18), sigma=np.log(1.2)),2))
        #priority = t + np.exp(np.random.lognormal(mean=np.log(18), sigma=np.log(1.2)))
        priority = t + np.random.lognormal(mean=np.log(18), sigma=np.log(1.2))
        
        resistance_status = 0
        q.put((priority, resistance_status))
        #print q.queue
        
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
            #print q.qsize()

        # Count mutants in the progeny
        r = 0
        while not q.empty():
            tup = q.get()
            r += tup[1]
        #print r
        r_colonies[colony] = r

    return r_colonies

#r_cell_cycle = cell_cycle_sim(N, mu, n_colonies)

#======================================================================================
#  Histogram of the number of resistant bacteria in 100 colonies
#======================================================================================

def compare_hist(N, mu, n_colonies):
    r_h1 = adaptation_sim(N, mu, n_colonies)
    r_h2 = cell_cycle_sim(N, mu, n_colonies)
    
    bins = range(int(max(max(r_h1),max(r_h2))))
    #bins = range(int(max(max(r_h1),100)))
    #plt.style.use('seaborn-deep')
    data = np.vstack([r_h1, r_h2]).T
    plt.hist(data, bins, alpha=0.7, label=['Adaptation Hypothesis', 'Mutation Hypothesis'])
    plt.title("Number of resistant bacteria in " + str(n_colonies) + " colonies")
    plt.xlabel("Number of resistant bacteria")
    plt.ylabel("Frequency")
    plt.legend(loc='upper right')
    plt.show()
    
    bins = range(int(max(max(r_h1),100)))
    plt.hist(data, bins, alpha=0.7, label=['Adaptation Hypothesis', 'Mutation Hypothesis'])
    plt.title("Number of resistant bacteria in " + str(n_colonies) + " colonies")
    plt.xlabel("Number of resistant bacteria")
    plt.ylabel("Frequency")
    plt.legend(loc='upper right')
    plt.show()

print 'Running compare_hist'
compare_hist(N, mu, n_colonies)

#======================================================================================
#  SD to mean ratio
#======================================================================================
mutation_rates = np.array([5*1e-7,1e-6,5*1e-6, 1e-5, 5*1e-5],dtype=float) # mutation rates
mu_n_ratios = np.zeros(shape = (len(mutation_rates),3))

print 'Computing mu_n_ratios'
for i in range(len(mutation_rates)):
    print(i)
    mu = mutation_rates[i]
    r_h1 = adaptation_sim(N, mu, n_colonies)
    #sum(r_h1)
    r_h2 = cell_cycle_sim(N, mu, n_colonies)
    #sum(r_h2)
    mu_n_ratios[i,0] = np.var(r_h1) / np.mean(r_h1)
    mu_n_ratios[i,1] = np.var(r_h2) / np.mean(r_h2)
    mu_n_ratios[i,2] = np.sqrt(n_colonies)/np.log((N)*n_colonies*mu) # formula 12 from the paper
    print('R1='+str(mu_n_ratios[i,0])+'; R2='+str(mu_n_ratios[i,1])+'; Rt='   +str(mu_n_ratios[i,2])) 

''' It took very long to run, so let's save the results

filename = 'C:/Users/semeel/Desktop/Computational Biology 2/solutions/mu_n_ratios.csv'
np.savetxt(filename, mu_n_ratios, fmt='%f', delimiter=',')
from numpy import genfromtxt
mu_n_ratios = genfromtxt(filename, delimiter=',')
print(mu_n_ratios)
'''


plt.plot(mutation_rates, mu_n_ratios[:,0], label='Adaptation Hypothesis')
plt.plot(mutation_rates, mu_n_ratios[:,1], label='Mutation Hypothesis')
plt.plot(mutation_rates, mu_n_ratios[:,2], label='Theoretical')
plt.xlabel('mutation rate', fontsize=18)
plt.ylabel('SD-to-mean ratio', fontsize=16)
plt.title("SD-to-mean ration in " + str(n_colonies) + " colonies")
plt.legend(loc='upper right')
plt.gca().set_xscale('log')
plt.show()
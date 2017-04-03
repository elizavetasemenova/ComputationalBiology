"""
    
Simulate mutations in the population, following Luria Delbrück (1943) paper.   

"""

import numpy as np
import matplotlib.pyplot as plt


# Approximate final population size
N = 10**6

# Compute T (time [number of generations] until population size 
# of a given colony reaches N), i.e. N = 2^T
T = int(round(np.log(N)/np.log(2)))

# mutation rate
mu = 1e-6

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
#    T - final time (equivalently, number of generations)
#    mu - mutation rate
#    n_colonies - number of colonies
# Output: 
#    r_colonies - number of resistant bacteria in the final population for each colony
def adaptation_sim(T, mu, n_colonies):
    # Population size of generation T
    N = 2**T
    
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
#    T - final time (equivalently, number of generations)
#    mu - mutation rate
#    n_colonies - number of colonies
# Output: 
#    r_colonies - number of resistant bacteria in the final population for each colony

def mutation_sim(T, mu, n_colonies):
    r_colonies = np.zeros(n_colonies) # number of resistant bacteria in each colonie
    
    # Compute number of resistants for each colony
    for colony in range(n_colonies):
        r=0
        for t in range(1,T+1):
            r = 2*r + np.random.binomial(2**t - 2 * r, mu)            
        r_colonies[colony] = r
        
    return r_colonies

#======================================================================================
#  Histogram of the number of resistant bacteria in 100 colonies
#======================================================================================

def compare_hist(T, mu, n_colonies):
    r_h1 = adaptation_sim(T, mu, n_colonies)
    r_h2 = mutation_sim(T, mu, n_colonies)
    
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

compare_hist(T, mu, n_colonies)

#======================================================================================
#  SD to mean ratio
#======================================================================================
mutation_rates = np.array([5*1e-7,1e-6,5*1e-6, 1e-5, 5*1e-5],dtype=float) # mutation rates
mu_n_ratios = np.zeros(shape = (len(mutation_rates),3))

for i in range(len(mutation_rates)):
    mu = mutation_rates[i]
    r_h1 = adaptation_sim(T, mu, n_colonies)
    r_h2 = mutation_sim(T, mu, n_colonies)
    mu_n_ratios[i,0] = np.var(r_h1) / np.mean(r_h1)
    mu_n_ratios[i,1] = np.var(r_h2) / np.mean(r_h2)
    mu_n_ratios[i,2] = np.sqrt(n_colonies)/np.log((2**T)*n_colonies*mu) # formula 12 from the paper
    print('R1='+str(mu_n_ratios[i,0])+'; R2='+str(mu_n_ratios[i,1])+'; Rt='   +str(mu_n_ratios[i,2])) 


plt.plot(mutation_rates, mu_n_ratios[:,0], label='Adaptation Hypothesis')
#plt.plot(mutation_rates, mu_n_ratios[:,1], label='Mutation Hypothesis')
plt.plot(mutation_rates, mu_n_ratios[:,2], label='Theoretical')
plt.xlabel('mutation rate', fontsize=18)
plt.ylabel('SD-to-mean ratio', fontsize=16)
plt.title("SD-to-mean ration in " + str(n_colonies) + " colonies")
plt.legend(loc='upper right')
plt.gca().set_xscale('log')
plt.show()

plt.plot(mutation_rates, mu_n_ratios[:,0], label='Adaptation Hypothesis')
plt.plot(mutation_rates, mu_n_ratios[:,1], label='Mutation Hypothesis')
plt.plot(mutation_rates, mu_n_ratios[:,2], label='Theoretical')
plt.xlabel('mutation rate', fontsize=18)
plt.ylabel('SD-to-mean ratio', fontsize=16)
plt.title("SD-to-mean ration in " + str(n_colonies) + " colonies")
plt.legend(loc='upper right')
plt.gca().set_xscale('log')
plt.show()


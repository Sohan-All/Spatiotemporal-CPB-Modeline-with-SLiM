from random import randint
import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np

import pyabc

from pyabc.transition import MultivariateNormalTransition


# Example: use a negative binomial for numSocks (choose r and p as needed).
# loc=11 shifts support so minimum numSocks is 11.
prior = pyabc.Distribution(
    numSocks=pyabc.RV("nbinom", 4.62, 0.133),
    proportionPairs=pyabc.RV("beta", 15, 2)
)

#model to sample 11 socks from a set of x pairs and y singles. returns the number of unpaired socks chosen.
def model(prior):
    #print("-------------------------------model")
    numSocks = prior["numSocks"]
    proportionPairs = prior["proportionPairs"]
    
    numInPairs = int(numSocks * proportionPairs) // 2 * 2
    
    chosen = []
    
    if numSocks < 11:
        return {"singlesChosen": -np.inf}
    
    while len(chosen) < 11:
        #choose an unchosen sock
        while True:
            rand = randint(0, numSocks-1)
            #print("-----------------loop")
            if rand not in chosen:
                #print("------loop")
                chosen.append(rand)
                break
        
    chosen.sort()
    
    singlesChosen = 0
    pairsChosen = 0
    
    for c in chosen:
        #print("-yeezy")
        if c < numInPairs:
            if c % 2 == 0 and (c + 1) in chosen:
                pairsChosen += 1
            elif c % 2 == 1 and (c - 1) in chosen:
                pass
            else:
                singlesChosen += 1
        else:
            singlesChosen += 1
            
    return {"singlesChosen": singlesChosen}
    
    
def distance(x, x0):
    return abs(x["singlesChosen"] - x0["singlesChosen"])

abc = pyabc.ABCSMC(model, prior, distance, population_size=500, transitions=MultivariateNormalTransition())

db_path = os.path.join(tempfile.gettempdir(), "test.db")
observation = {"singlesChosen": 7}
abc.new("sqlite:///" + db_path, observation)

history = abc.run(minimum_epsilon=2, max_nr_populations=5)

df, w = history.get_distribution()
posterior_mean = (df * w[:, None]).sum()
print("Posterior mean:", posterior_mean)
print(df)
print(w)
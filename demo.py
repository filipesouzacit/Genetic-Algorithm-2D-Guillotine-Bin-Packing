#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 20:27:22 2020

@author: Filipe Lucas De Souza
"""
from util import Solution
import matplotlib.pyplot as plt
from random import seed
from GA_C import GA
import numpy as np
    
_popSize            = 100
_maxIterations      = 2300
_widthPlates        = 6000
_heightPlates       = 3210
_crossoverAlgorithm = 'P' # PMX
_mutationAlgorithm  = 'R' # Reciprocal Exchange
_mutationRate       = 0.3
_selectionAlgorithm = 'T' # Tournament Selection 
_localSearchRate    = 1
_eliteRate          = 0.1
instance_name = "data/dataset_A/A6"
instanceSize =  43254870

seed(182510)   
ga = GA(instance_name,_popSize,_mutationRate,_maxIterations,_widthPlates,
   _heightPlates,_crossoverAlgorithm,_mutationAlgorithm,_selectionAlgorithm,
   _localSearchRate,_eliteRate)
ga.search()

ga.getCost(Solution(ga.templateSolution[0:], ga.widthPlates, ga.heightPlates)).printSolution(ga.stacks)

ga.best.printSolution(ga.stacks)

AvgCost = instanceSize/(instanceSize+ np.mean(np.array(ga.output)[:,:-2],axis=1) )
BestCost =  instanceSize/(instanceSize+ np.array(ga.output)[:,-1])

plt.style.use("ggplot")
plt.figure()

plt.plot(np.arange(0,len(AvgCost)),AvgCost,label="AvgCost")
plt.plot(np.arange(0,len(BestCost)),BestCost,label="BestCost")
plt.xlabel("iteration")
plt.ylabel("Occupation Rate")
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))




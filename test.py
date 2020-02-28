#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 20:27:22 2020

@author: filipe
"""
import matplotlib.pyplot as plt
from random import seed
from GA import GA
import cProfile
import pandas as pd
   
_heightPlates = 3210 
_widthPlates  = 6000
_selectionAlgorithm = 'S' # 'R' - Random, and 'S' - Stochastic Universal Sampling
_crossoverAlgorithm = 'P' # 'U' - Uniform, and 'P' - PMX
_mutationAlgorithm  = 'INEW' # 'R' - Reciprocal Exchange, and 'I' - Inversion
_localSearchRate    = 0.2
instance_name = "data/dataset_A/A2"
seed(10)        



#Time 
cProfile.run('GA(instance_name, 1000, 0.4,50,6000,3210,_crossoverAlgorithm,_mutationAlgorithm,_selectionAlgorithm,_localSearchRate).search()')



#
ga = GA(instance_name, 1000, 0.1,50,6000,3210,_crossoverAlgorithm,
                 _mutationAlgorithm,_selectionAlgorithm,_localSearchRate) 
ga.search()
ga.best.printSolution(ga.stacks)
plt.plot( range(len(ga.listBestFit)), ga.listBestFit, marker='', color='olive', linewidth=2 , label="bestFit")
plt.plot( range(len(ga.listBestFit)), ga.listAvgFit, marker='', color='blue', linewidth=2 ,label="AvgFit")
plt.legend()
plt.title('graph of evolution of fitness across the iteration for configuration ')
plt.show() 




#     
instance_names = ["data/dataset_A/A1","data/dataset_A/A2","data/dataset_A/A3","data/dataset_A/A4","data/dataset_A/A5","data/dataset_A/A6","data/dataset_A/A7","data/dataset_A/A8","data/dataset_A/A9","data/dataset_A/A10","data/dataset_A/A11","data/dataset_A/A12","data/dataset_A/A13","data/dataset_A/A14","data/dataset_A/A15","data/dataset_A/A16","data/dataset_A/A17","data/dataset_A/A18","data/dataset_A/A19","data/dataset_A/A20","data/dataset_B/B1","data/dataset_B/B2","data/dataset_B/B3","data/dataset_B/B4","data/dataset_B/B5","data/dataset_B/B6","data/dataset_B/B7","data/dataset_B/B8","data/dataset_B/B9","data/dataset_B/B10","data/dataset_B/B11","data/dataset_B/B12","data/dataset_B/B13","data/dataset_B/B14","data/dataset_B/B15","data/dataset_X/X1","data/dataset_X/X2","data/dataset_X/X3","data/dataset_X/X4","data/dataset_X/X5","data/dataset_X/X6","data/dataset_X/X7","data/dataset_X/X8","data/dataset_X/X9","data/dataset_X/X10","data/dataset_X/X11","data/dataset_X/X12","data/dataset_X/X13","data/dataset_X/X14","data/dataset_X/X15"]
bestCost = []
bestInitialCost = []
for instance_name in instance_names:
    seed(10)        
    ga = GA(instance_name, 1000, 0.1,100,6000,3210,_crossoverAlgorithm,
                 _mutationAlgorithm,_selectionAlgorithm,_localSearchRate) 
    ga.search()
    bestCost.append(ga.best.cost)
    bestInitialCost.append(ga.bestInitialSol)
data =  {'instance': instance_names, 'initialCost': bestInitialCost, 'bestCost':bestCost}
df = pd.DataFrame(data)
df.to_csv('result.csv')









    
    



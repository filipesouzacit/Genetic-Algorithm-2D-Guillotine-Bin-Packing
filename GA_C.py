#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 18:59:47 2020

@author: Filipe Lucas De Souza
"""
import pandas as pd
import numpy as np
from random import choice, randint, uniform, random
import util
from util import Stack, Item, Solution
import time

class GA:
    def __init__(self, _instance, _popSize, _mutationRate, _maxIterations,
                 _widthPlates, _heightPlates,_crossoverAlgorithm,
                 _mutationAlgorithm,_selectionAlgorithm,_localSearchRate, _eliteRate):
        """
        Parameters and general variables
        """
        self.widthPlates        = _widthPlates
        self.heightPlates       = _heightPlates
        self.crossoverAlgorithm = _crossoverAlgorithm
        self.mutationAlgorithm  = _mutationAlgorithm
        self.selectionAlgorithm = _selectionAlgorithm
        self.output         = []
        self.population     = []
        self.elitePopulation= []
        self.newPopulation  = []
        self.eliteRate      = _eliteRate
        self.eliteSize      = _popSize * _eliteRate
        self.newPopulSize   = _popSize - self.eliteSize
        self.matingPool     = []
        self.best           = None
        self.popSize        = _popSize
        self.stacks         = {}
        self.mutationRate   = _mutationRate
        self.localSearchRate= _localSearchRate
        self.maxIterations  = _maxIterations
        self.iteration      = 0
        self.iterationOfBest= 0
        self.instance       = _instance        
        self.listBestFit    = []
        self.listAvgFit     = []
        self.timeIteration  = []
        self.table          = {}
        self.readInstance()
        self.initPopulation()
        self.genSize = len(self.templateSolution)
        
    def readInstance(self):
        batch = self.instance + "_batch.csv"
        batch = pd.read_csv(batch, sep = ";")        
        stack = Stack(batch.STACK[0])
        self.stacks[stack.idStack] = stack
        self.templateSolution = []
        for ix,it in batch.iterrows():
            item = Item(it.ITEM_ID,it.LENGTH_ITEM, it.WIDTH_ITEM)
            if stack.idStack != it.STACK:
                stack = Stack(it.STACK)
                self.stacks[stack.idStack] = stack
            stack.add(item)
            self.templateSolution.append(it.STACK)
    
    def getCost(self, solution):
        [stack.reset() for stack in self.stacks.values()]        
        if random() < self.localSearchRate:
            solution.localSearch(self.stacks)
        else:
            solution.computeCost(self.stacks)
        return solution
            
        
    def initPopulation(self):
        """
        Creating random individuals in the population
        """
        for i in range(0, self.popSize):
            solution = Solution(self.templateSolution[0:], self.widthPlates, self.heightPlates)
            solution = self.getCost(solution)
            self.population.append(solution)
        self.best = self.population[0].copy()
        for sol_i in self.population:
            #elite 
            if len(self.elitePopulation) < self.eliteSize:
                self.elitePopulation.append(sol_i)
            else:
                fit = [i.cost for i in self.elitePopulation]
                fit.append(0)
                maxFit = max(fit)
                if maxFit > sol_i.cost:
                    self.elitePopulation[fit.index(maxFit)] = sol_i
            #best fit    
            if self.best.cost > sol_i.cost:
                self.best = sol_i.copy()
        print ("Best initial sol: ",self.best.cost)
        self.bestInitialSol = self.best.cost
        self.addOutput(0,self.best.cost, self.population)
    
    def updateBest(self, candidate):
        if candidate.cost < self.best.cost:
            self.best = candidate.copy()
            print ("iteration: ",self.iteration, "best: ",self.best.cost)
            self.iterationOfBest = self.iteration
    
    def addOutput(self,time,bestCost, population):
        row = [i.cost for i in population]
        row.append(time)
        row.append(bestCost)
        self.output.append(row)
             
    
    def randomSelection(self):
        """
        Random (uniform) selection of two individuals
        """
        solA = self.matingPool[ randint(0, self.popSize-1) ]
        solB = self.matingPool[ randint(0, self.popSize-1) ]
        return [solA, solB]

    def stochasticUniversalSampling(self):
        """
        stochastic universal sampling Selection Implementation
        """
        solA = self.matingPool[ choice(self.indexs) ]
        solB = self.matingPool[ choice(self.indexs) ]
        return [solA, solB]
    
    def TournamentSelection(self):
        """
        Tournament selection of two individuals
        """
        sol1 = self.matingPool[ randint(0, self.popSize-1) ]
        sol2 = self.matingPool[ randint(0, self.popSize-1) ]
        solA = sol1 if sol1.cost < sol2.cost else sol2

        sol3 = self.matingPool[ randint(0, self.popSize-1) ]
        sol4 = self.matingPool[ randint(0, self.popSize-1) ]
        solB = sol3 if sol3.cost < sol4.cost else sol4
        if solA.cost == solB.cost:
            if solA.cost != sol1.cost:
                solA = sol1
            elif solA.cost != sol2.cost:
                solA = sol2
            elif solA.cost != sol3.cost:
                solA = sol3
            elif solA.cost != sol4.cost:
                solA = sol4
        return [solA, solB]

    def uniformCrossover(self, solA, solB):
        """
        Uniform Crossover Implementation
        """
        child1 = solA.copy()
        child2 = solB.copy()
        tmpIndA = solA.genes[0:]
        tmpIndB = solB.genes[0:]
        tmpIndex= []
        
        for i in range(0, self.genSize):
            if choice([True,False]):
                tmpIndA.remove(child2.genes[i])
                tmpIndB.remove(child1.genes[i])
                
                tmpRotetion = child2.genes[self.genSize+i]
                child2.genes[self.genSize+i] = child1.genes[self.genSize+i]
                child1.genes[self.genSize+i] = tmpRotetion
                
                tmpCut = child2.genes[(self.genSize*2)+i]
                child2.genes[(self.genSize*2)+i] = child1.genes[(self.genSize*2)+i]
                child1.genes[(self.genSize*2)+i] = tmpCut
                                
            else:
                tmpIndex.append(i)
        i=0
        for g in tmpIndex:
            child2.genes[g] = tmpIndA[i]
            child1.genes[g] = tmpIndB[i]
            i +=1
        return (child1, child2)
            
    def updateChild(self,indFix,indComp,index,i ):
        """
        This fuction updates the child with second parent genes based on PMX crossover. 
        """
        if indComp.genes[i] not in indFix.genes[index[0]:index[1]+1]:
            childGene = indComp.genes[i]
        else:
            gene = indComp.genes[indFix.genes.index(indComp.genes[i])]
            while gene in indFix.genes[index[0]:index[1]+1]:
                gene = indComp.genes[indFix.genes.index(gene)]
            childGene = gene
        return childGene


    def pmxCrossover(self, indA, indB):
        """
        PMX Crossover Implementation
        """
        child1 = indA.copy()
        child2 = indB.copy()        
        index = [randint(0, self.genSize-1) for _ in range(2)]        
        index.sort()
        for i in range(0, self.genSize):
            if i >= index[0] and i <= index[1]:
                tmpRotetion = child2.genes[self.genSize+i]
                child2.genes[self.genSize+i] = child1.genes[self.genSize+i]
                child1.genes[self.genSize+i] = tmpRotetion
                
                tmpCut = child2.genes[(self.genSize*2)+i]
                child2.genes[(self.genSize*2)+i] = child1.genes[(self.genSize*2)+i]
                child1.genes[(self.genSize*2)+i] = tmpCut
                
            else:
                 child1.genes[i] = self.updateChild(indA,indB,index,i)
                 child2.genes[i] = self.updateChild(indB,indA,index,i)   
        return (child1, child2)
    
    def reciprocalExchangeMutation(self, ind):
        """
        Reciprocal Exchange Mutation implementation
        """
        if random() < self.mutationRate:
            indexA = randint(0, (self.genSize*3)-1)
            if indexA < self.genSize:
                indexB = randint(0, self.genSize-1)    
                tmp = ind.genes[indexA]
                ind.genes[indexA] = ind.genes[indexB]
                ind.genes[indexB] = tmp
            else:
                ind.genes[indexA] = not(ind.genes[indexA])
            
        ind = self.getCost(ind)
        self.updateBest(ind)
        return ind

    def inversionMutation(self, ind):
        """
        Inversion Mutation implementation
        """
        if random() < self.mutationRate:
            indexA = randint(0, (self.genSize*3)-1)
            if indexA < self.genSize:
                index = [randint(0, self.genSize-1) for _ in range(2)]        
                index.sort()
                ind.genes[index[0]:index[1]+1] = reversed(ind.genes[index[0]:index[1]+1])
            else:
                ind.genes[indexA] = not(ind.genes[indexA])
        ind = self.getCost(ind)
        self.updateBest(ind)
        return ind
    
    def inversionMutationNew(self, ind):
        """
         A variation of Inversion Mutation implementation, where we change the 
         position of a gene with the next gene. I believe that this function works 
         bether when the initial population are create by Nearest neighbor insertion.
        """
        if random() < self.mutationRate:
            indexA = randint(1, (self.genSize*3)-1)
            if indexA < self.genSize:
                ind.genes[indexA-1:indexA+1] = reversed(ind.genes[indexA-1:indexA+1])
            else:
                ind.genes[indexA] = not(ind.genes[indexA])                
        ind = self.getCost(ind)
        self.updateBest(ind)
        return ind
        
    def eliteSurvival(self, ind):
        """
        Ensuring that only the best individuals will be added to the population.
        """
        fit = [i.cost for i in self.elitePopulation]
        fit.append(0)
        maxFit = max(fit)
        if maxFit > ind.cost:
            self.elitePopulation[fit.index(maxFit)] = ind.copy()   
        
        if len(self.newPopulation) < self.popSize:
            self.newPopulation.append(ind)
        else:
            fit = [i.cost for i in self.newPopulation]
            fit.append(0)
            maxFit = max(fit)
            if maxFit > ind.cost:
                self.newPopulation[fit.index(maxFit)] = ind
            
    def updateMatingPool(self):
        """
        Updating the mating pool before creating a new generation
        """
        self.matingPool = [s.copy() for s in self.population]
        
        fit = [i.cost for i in self.matingPool]
        self.listBestFit.append(min(fit)) 
        self.listAvgFit.append(sum(fit)/len(fit))
       

        """
        Updating the indexs for stochastic Universal Sampling before creating a new generation
        """
        if self.selectionAlgorithm == 'S':
            fitnessMinim = [1/i.cost for i in self.matingPool]
            sumFitnessMinim = sum(fitnessMinim)
            fracFitnessMinim = [i/sumFitnessMinim for i in fitnessMinim]
            cumSumFracFitnessMinim = [sum(fracFitnessMinim[:i]) for i in range(1, len(fracFitnessMinim)+1)]
            N = int(len(self.matingPool))
            startPoint = uniform(0, (1/N))
            marks = [startPoint + ((1/N) * i) for i in range(0,N)]
            self.indexs = []
            i = 0
            for point in marks:
                while(cumSumFracFitnessMinim[i]<point):
                    i +=1
                self.indexs.append(i)

    
    def newGeneration(self):
        """
        Creating a new generation
        1. Selection
        2. Crossover
        3. Mutation
        """
        self.newPopulation  = self.elitePopulation[:]
        for i in range(0, int(self.newPopulSize/2)+1):
            """
            Depending of your experiment you need to use the most suitable algorithms for:
            1. Select two candidates
            2. Apply Crossover
            3. Apply Mutation
            """
            if self.selectionAlgorithm == 'S':
                indA, indB = self.stochasticUniversalSampling()
            elif self.selectionAlgorithm == 'T':
                indA, indB = self.TournamentSelection()
            else:
                indA, indB = self.randomSelection()
                
            if self.crossoverAlgorithm == 'P':
                child1,child2 = self.pmxCrossover(indA, indB)
            else:
                child1,child2 = self.uniformCrossover(indA, indB)
            if self.mutationAlgorithm == 'I':
                child1 = self.inversionMutation(child1)
                child2 = self.inversionMutation(child2)
            elif self.mutationAlgorithm == 'INEW':
                child1 = self.inversionMutationNew(child1)
                child2 = self.inversionMutationNew(child2)
            else:
                child1 = self.reciprocalExchangeMutation(child1)
                child2 = self.reciprocalExchangeMutation(child2)
            self.eliteSurvival(child1)
            self.eliteSurvival(child2) 
        self.population = self.newPopulation
      
        
    def filterDuplicate(self):
        sol = np.array([str(i.getSequence(self.stacks,i.localSearchGenes)) for i in self.population])
        indexDupl = [idx for idx, val in enumerate(sol) if val in sol[:idx]]
        for i in indexDupl:
            solution = Solution(self.templateSolution[0:], self.widthPlates, self.heightPlates)
            index = [randint(self.genSize, self.genSize*3) for _ in range(2)]        
            index.sort()
            solution.genes[index[0]:index[1]] = self.best.genes[index[0]:index[1]]
            solution = self.getCost(solution)
            self.population[i] = solution
        self.elitePopulation = []
        for sol_i in self.population:
            #elite 
            if len(self.elitePopulation) < self.eliteSize:
                self.elitePopulation.append(sol_i)
            else:
                fit = [i.cost for i in self.elitePopulation]
                fit.append(0)
                maxFit = max(fit)
                if maxFit > sol_i.cost:
                    self.elitePopulation[fit.index(maxFit)] = sol_i
            #best fit    
            if self.best.cost > sol_i.cost:
                self.best = sol_i.copy()   
            

    def GAStep(self):
        """
        One step in the GA main algorithm
        1. Updating mating pool with current population
        2. Creating a new Generation
        """
        start = time.time()
        self.updateMatingPool()
        self.newGeneration()
        if self.iteration % 50 == 0:
            self.filterDuplicate()
        end = time.time()
        self.timeIteration.append(end-start)
        self.addOutput(end-start,self.best.cost, self.population)

    def search(self):
        """
        General search template.
        Iterates for a given number of steps
        """
        self.iteration = 0
        while self.iteration < self.maxIterations and self.best.cost > 0:
            self.GAStep()
            self.iteration += 1
            if self.iteration% 500 == 0:
                print("Iteration ", self.iteration)

        print ("Total iterations: ",self.iteration)
        print ("Best Solution: ", self.best.cost)

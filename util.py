#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 15:43:35 2020

@author: filipe
"""

from random import choice,sample,random
import matplotlib.pyplot as plt
import cv2
import numpy as np

class Solution:    
    def __init__(self, _templateSolution, _widthPlates, _heightPlates, _genes=None, _cost=0):
        self.widthPlates = _widthPlates
        self.heightPlates = _heightPlates
        self.chromosome = _templateSolution
        self.chromoSize = len(self.chromosome)
        if _genes:
            self.genes = _genes[0:]
            self.cost = _cost
        else:
            self.randomSolution()
    
    def randomSolution(self):
        self.genes = sample(list(range(self.chromoSize)), self.chromoSize)
        self.genes.extend([random() >= .5 for i in range(self.chromoSize*2)]) 

    def copy(self):
        return Solution(self.chromosome,self.widthPlates,self.heightPlates, self.genes, self.cost)            
       
    def computeCost(self,stacks):
        self.cost = 0
        self.idPlate = 0
        x,y = 0,0
        self.branchs = [[self.idPlate,x,y,self.widthPlates,self.heightPlates]]
        self.branch = self.branchs.pop()
        for g in self.genes[0:self.chromoSize]:
            stack = stacks[self.chromosome[g]]
            item = stack.pop()
            if item.width > self.heightPlates:
                self.genes[item.idItem+self.chromoSize] = False
            item = item.get(self.genes[item.idItem+self.chromoSize]) #rotate or not
            while not(self.isFit(item)):
                self.cost +=  self.branch[3]* self.branch[4]
                if self.branchs:
                    self.branch = self.branchs.pop()
                else:
                    self.idPlate += 1
                    self.branch = [self.idPlate,x,y,self.widthPlates,self.heightPlates]                
            if self.genes[item.idItem+(self.chromoSize *2)]:
                self.fistCutVertical(item)
            else:
                self.fistCutHorizontal(item)        
        while len(self.branchs)>0 :
            self.cost +=  self.branch[3]* self.branch[4]
            self.branch = self.branchs.pop()
            
    def isFit(self, _item):
        return (_item.width <= self.branch[3]) and (_item.height <= self.branch[4])
    
    def fistCutVertical(self,_item):
        if _item.width < self.branch[3]:
            branch1,branch2 = self.cutVertical(self.branch, _item.width)
            self.branchs.append(branch1)
            self.branch,self.item = self.cutHorizontal(branch2, _item.height)
    
    def fistCutHorizontal(self,_item):
        if _item.height < self.branch[4]:
            branch1,branch2 = self.cutHorizontal(self.branch, _item.height)
            self.branchs.append(branch1)
            self.branch,self.item = self.cutVertical(branch2, _item.width) 
    
    def cutVertical(self,_branch, _width):
        branch1 = _branch[0:]
        branch1[1] = branch1[1]+_width
        branch1[3] = branch1[3]-_width
        branch2 = _branch[0:]
        branch2[3] = _width
        return (branch1,branch2)
    
    def cutHorizontal(self,_branch, _height):
        branch1 = _branch[0:]
        branch1[2] = branch1[2]+_height
        branch1[4] = branch1[4]-_height
        branch2 = _branch[0:]
        branch2[4] = _height
        return (branch1,branch2)
    
    def localSearch(self,stacks):
        self.cost = 0
        self.idPlate = 0
        x,y = 0,0
        self.branchs = [[self.idPlate,x,y,self.widthPlates,self.heightPlates]]
        self.branch = self.branchs.pop()
        self.tmpSeqGenes = []
        self.tmpIdStacks = []
        for g in self.genes[0:self.chromoSize]:         
            idStack = self.chromosome[g]
            if idStack in self.tmpIdStacks:
                self.tmpIdStacks.remove(idStack)
            else:
                stack = stacks[idStack]
                item = stack.pop()
                if item.width > self.heightPlates:
                    self.genes[item.idItem+self.chromoSize] = False
                item = item.get(self.genes[item.idItem+self.chromoSize]) #rotate or not
                while not(self.isFit(item)):
                    if (self.isFit(item.get(True))):
                        item = item.get(True)
                        self.genes[item.idItem+self.chromoSize] = not(self.genes[item.idItem+self.chromoSize] )
                    else:
                        self.branch = self.findFitItem(stacks,stack.idStack)
                        self.cost +=  self.branch[3]* self.branch[4]
                        if self.branchs:
                            self.branch = self.branchs.pop()
                        else:
                            self.idPlate += 1
                            self.branch = [self.idPlate,x,y,self.widthPlates,self.heightPlates]                
                if self.genes[item.idItem+(self.chromoSize *2)]:
                    self.fistCutVertical(item)
                else:
                    self.fistCutHorizontal(item)
                self.tmpSeqGenes.append(item.idItem)
        self.genes[0:self.chromoSize] = self.tmpSeqGenes[0:]
        while len(self.branchs)>0 :
            self.cost +=  self.branch[3]* self.branch[4]
            self.branch = self.branchs.pop()
    
    def findFitItem(self,stacks,idStack):        
        fitItem = [ (s.get().size,i,False) for i,s in stacks.items() if i != idStack and s.index < s.size and self.isFit(s.get())]
        fitItem.extend([ (s.get().size,i,True) for i,s in stacks.items() if i != idStack and s.index < s.size and self.isFit(s.get().get(True))])
        if len(fitItem)>0:
            size,s,r = max(fitItem)
            item = stacks[s].pop().get(r)
            if self.genes[item.idItem+(self.chromoSize *2)]:
                self.fistCutVertical(item)
            else:
                self.fistCutHorizontal(item)
            self.genes[item.idItem+(self.chromoSize)] = r
            self.tmpSeqGenes.append(item.idItem)
            self.tmpIdStacks.append(s)
            self.branch = self.findFitItem(stacks,idStack)
        return self.branch
        
    def printSolution(self,stacks):
        self.colors = [(128,128,0),(135,206,250),(210,105,30),(65,105,225),(244,164,96),(102,205,170),(119,136,153),(72,61,139),(64,224,208),(255,127,80),(192,192,192),(186,85,211),(128,0,128),(0,255,127),(147,112,219),(255,182,193),(123,104,238),(30,144,255),(112,128,144),(107,142,35),(34,139,34),(184,134,11),(105,105,105),(70,130,180),(124,252,0),(143,188,143),(255,165,0),(255,20,147),(0,128,128),(238,232,170),(127,255,0),(240,128,128),(135,206,235),(0,255,255),(95,158,160),(139,69,19),(0,100,0),(0,191,255),(106,90,205),(75,0,130),(148,0,211),(128,0,0),(153,50,204),(176,224,230),(219,112,147),(0,128,0),(128,128,128),(255,222,173),(72,209,204),(255,105,180),(0,0,255),(205,133,63),(85,107,47),(32,178,170),(218,112,214),(233,150,122),(176,196,222),(205,92,92),(0,0,205),(255,140,0),(138,43,226),(144,238,144),(199,21,133),(255,99,71),(189,183,107),(238,130,238),(0,250,154),(218,165,32),(255,160,122),(46,139,87),(255,192,203),(139,0,0),(255,0,0),(220,20,60),(178,34,34),(216,191,216),(210,180,140),(154,205,50),(0,0,139),(0,255,255),(173,255,47),(127,255,212),(100,149,237),(0,255,0),(47,79,79),(25,25,112),(173,216,230),(60,179,113),(0,139,139),(221,160,221),(255,228,181),(0,206,209),(160,82,45),(0,0,128),(169,169,169),(240,230,140),(139,0,139),(255,215,0),(255,255,0),(250,128,114),(175,238,238),(255,69,0),(188,143,143),(165,42,42),(50,205,50),(255,0,255),(152,251,152),(255, 250, 250),(250,235,215),(245,245,220),(255,235,205),(255,248,220),(220,220,220),(230,230,250)]
        r = 9
        w = 0
        self.font = cv2.FONT_HERSHEY_SIMPLEX
        img = np.zeros((self.heightPlates//r, self.widthPlates//r, 3), np.uint8)
        self.idPlate = 0
        x,y = 0,0
        self.branchs = [[self.idPlate,x,y,self.widthPlates,self.heightPlates]]
        self.branch = self.branchs.pop()
        for g in self.genes[0:self.chromoSize]:
            stack = stacks[self.chromosome[g]]
            item = stack.pop()            
            item = item.get(self.genes[item.idItem+self.chromoSize]) #rotate or not
            while not(self.isFit(item)):
                self.cost +=  self.branch[3]* self.branch[4]
                #self.printWaste(w,r,img, self.branch)
                w+=1
                if self.branchs:
                    self.branch = self.branchs.pop()
                else:
                    self.idPlate += 1
                    self.branch = [self.idPlate,x,y,self.widthPlates,self.heightPlates]
                    plt.show()
                    img = np.zeros((self.heightPlates//r, self.widthPlates//r, 3), np.uint8) 
            if self.genes[item.idItem+(self.chromoSize *2)]:
                self.fistCutVertical(item)
            else:
                self.fistCutHorizontal(item)
            self.printItem(item.idItem,r,img, self.item)
        while len(self.branchs)>0 :
            self.cost +=  self.branch[3]* self.branch[4]
            self.branch = self.branchs.pop()
        x = self.branch[1]
        y = self.branch[2]
        w = self.branch[3]
        h = self.branch[4]
        img = cv2.rectangle(img, (int(x//r), int(y//r)), 
                                (int((x + w)//r), int((y + h)//r)), 
                                (255, 255, 255), -1)
        plt.imshow(img)
        plt.xticks(np.arange(0, self.widthPlates + 1, 1000)//r, np.arange(0, self.widthPlates + 1, 1000))
        plt.yticks(np.arange(0, self.heightPlates + 1, 500)//r, np.arange(0, self.heightPlates + 1, 500)) 
        
        
    
    def printItem(self, idItem,r,img,item):
        x = item[1]
        y = item[2]
        w = item[3]
        h = item[4]
        img = cv2.rectangle(img, (int(x//r), int(y//r)), 
                                (int((x + w)//r), int((y + h)//r)), 
                                self.colors[idItem%100], -1)
        cv2.putText(img, str(idItem), (int((x+w//5)//(r)), int((y + h//1.1)//(r))) ,
                    self.font, 1, (0, 0, 0), 2, cv2.LINE_AA)
        plt.imshow(img)
        plt.xticks(np.arange(0, self.widthPlates + 1, 1000)//r, np.arange(0, self.widthPlates + 1, 1000))
        plt.yticks(np.arange(0, self.heightPlates + 1, 500)//r, np.arange(0, self.heightPlates + 1, 500))
    
    def printWaste(self, idItem,r,img,item):
        x = item[1]
        y = item[2]
        w = item[3]
        h = item[4]
        img = cv2.rectangle(img, (int(x//r), int(y//r)), 
                                (int((x + w)//r), int((y + h)//r)), 
                                self.colors[107+(idItem%7)], -1)
        plt.imshow(img)
        plt.xticks(np.arange(0, self.widthPlates + 1, 1000)//r, np.arange(0, self.widthPlates + 1, 1000))
        plt.yticks(np.arange(0, self.heightPlates + 1, 500)//r, np.arange(0, self.heightPlates + 1, 500))
        
        
    
class Gene:
    def __init__(self, _stack,_rotation,_cut):
        self.stack = _stack
        self.rotation = _rotation
        self.cut = _cut
        
    def get(self):
        return self.rotation, self.cut
    
    def randomSet(self):
        return Gene(self.stack, choice([True, False]),choice([True,False]))
    
    def cutChange(self):
        self.cut = not(self.cut)
    
    def rotationChange(self):
        self.rotation = not(self.rotation)

class Stack:    
    def __init__(self,_idStack):
        self.idStack = _idStack
        self.size = 0
        self.index = 0
        self.items = []
        
    def get(self):
        return self.items[self.index]
            
    def pop(self):
        if self.index < self.size:
            self.index += 1
            return self.items[self.index-1]
        self.reset()
        self.index += 1
        return self.items[self.index-1]
    
    def add(self,_item):
        self.size +=1
        self.items.append(_item)
    
    def reset(self):
        self.index = 0
        
class Item:
    def __init__(self,_idItem,_height,_width):
        self.idItem = _idItem
        self.width = _width
        self.height = _height
        self.size = _height * _width   
    
    def get(self,r):
        if r:
            return Item(self.idItem,self.width,self.height)
        return Item(self.idItem,self.height,self.width)

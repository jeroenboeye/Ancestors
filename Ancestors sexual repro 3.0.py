# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 11:12:58 2012
models a plant species that can adapt to a gradient in deathrate at the cost of reproduction
in a fractal landscape
@author: Jboeye
"""

import numpy as np
from numpy.numarray.random_array import poisson
import fractal2
import Tkinter as tk
import math as math
import random as rnd
CW_SIZE=28
cw_loc=CW_SIZE/2
class Weergave:
    '''This class arranges the visual output. You can draw the landscape in black and white 
    according to suitability,    create the drawing of an individual, delete a drawing and 
    update the canvas.'''

    def __init__(self,max_x,max_y):
        self.zoom=5
        self.max_x=max_x
        self.max_y=max_y
        self.root = tk.Tk()
        self.canvas = tk.Canvas(self.root, width=self.max_x*self.zoom, height=self.max_y*self.zoom)      #create window
        self.canvas.pack()
        self.canvas.config(background='white')   #unsuitable habitat=background
    
    def draw_landscape(self,suitability):
        self.suitability=suitability
        factor=self.zoom/4 
        for x in xrange (self.max_x):
           for z in xrange (self.max_y):
               if suitability[x][z]==1:
                   self.canvas.create_rectangle(self.zoom*x-factor, self.zoom*z-factor, self.zoom*x+factor, self.zoom*z+factor, outline='grey',fill='grey') 
                   
    def update(self):
        self.canvas.update()
        
    def create (self, x, y):
        zx=self.zoom*x
        zy=self.zoom*y
        factor=self.zoom/2        
        self.drawing=self.canvas.create_rectangle(zx-factor, zy-factor, zx+factor, zy+factor, outline='darkgreen',fill='darkgreen')
        return self.drawing
        
    def delete (self,drawing):
        self.drawing=drawing
        self.canvas.delete(self.drawing)
    
    def reset (self):
        self.canvas.delete(tk.ALL)
        
    def draw_line(self,pos1,pos2,gene,extremes):
        minimum=extremes[0]
        maximum=extremes[1]
        color_gene=gene-minimum
        color_gene/=float(maximum-minimum)
        if color_gene<0:
            color_gene=0
        elif color_gene>1:
            color_gene=1        
        r = int(255*color_gene)
        g = 255-r
        b = 0
        rgb = r, g, b     
        Hex = '#%02x%02x%02x' % rgb
        self.canvas.create_line(pos1[0]*self.zoom,pos1[1]*self.zoom,pos2[0]*self.zoom,pos2[1]*self.zoom,fill=str(Hex),width=2)
      


            
class Landscape:
    '''This class is designed to create the landscape (suitability matrix) by use of the fractal.py 
    module. You can check the local suitability and check whether an individual is within the landscape. 
    For visual output a number of functions where added, all of them simply transfer the requested 
    operation to the weergave class.'''
    
    def __init__(self,max_x,max_y,h,p):
        self.max_x=max_x
        self.max_y=max_y
        self.h=h
        self.p=p
        self.suitability=fractal2.create(self.max_x,self.max_y,self.h,self.p,False)
        self.weergave=Weergave(self.max_x,self.max_y)
        #self.weergave.draw_landscape(self.suitability)
        #self.weergave2=Weergave(self.max_x,self.max_y)
        #self.weergave2.draw_landscape(self.suitability)

    def suitability_check(self,x,y):
        return self.suitability[x][y]
        
    def draw_landscape(self):
        self.weergave.draw_landscape(self.suitability)
        
    def boundary_check(self,x,y):
        survival=False
        if (0<x<self.max_x
        and 0<y<self.max_y):      #gradient            
            survival=True
        return survival
        
    def cw_check(self,x):
        survival=False
        if x<=cw_loc+CW_SIZE/2 and x>=cw_loc-CW_SIZE/2:
            survival=True
        return survival        
        
    def update(self):
        self.weergave.update()
        
    def move (self,drawing,x,y):
        self.drawing=drawing
        self.weergave.move(self.drawing,x,y)
       
    def create (self, x, y):
        self.drawing=self.weergave.create(x,y)
        return self.drawing
        
    def delete(self, drawing):
        self.drawing=drawing
        self.weergave.delete(self.drawing)
        
    def reset(self):
        self.weergave.reset()
        
        
class Density:
    '''This class tracks the local densities. You can reset all densities to zero,
        increase the local density, check a local density or check the whole density matrix
        Additionally (and not really related to density) I designed a few functions to 
        keep track of the index of the individual that needs to be processed. I do suspect
        there is a more elegant solution for this.'''
    
    def __init__(self,max_x,max_y):
        self.sum_gendispersal_allele=0
        self.average_gencounter=0
        self.max_x=max_x
        self.max_y=max_y
        self.densloc=[]
        self.male_genes=[]
        self.genindividual=0
        for x in xrange(self.max_x):
            self.densloc.append([0] * self.max_y)           #alles op 0 zetten  
        self.male_genes=np.empty((self.max_x,self.max_y),dtype=list)
            
    def add_male_gene(self,x,y,genes):
        if self.male_genes[x,y]==None:
            self.male_genes[x,y]=[genes]
        else:
            newlist=self.male_genes[x,y]
            newlist.append(genes)
            self.male_genes[x,y]=newlist
        
    def reset_male_genes(self):
        self.male_genes=np.empty((self.max_x,self.max_y),dtype=list)
            
    def get_male_gene(self,x,y):
        if self.male_genes[x,y]!=None:
            if len(self.male_genes[x,y])==1:
                return self.male_genes[x,y][0]
            else:  
                winning_male_gene=rnd.choice(self.male_genes[x,y])
                return winning_male_gene
        else:
            return None
        
    def reset_dispersal_sum(self):
        self.sum_gendispersal_allele=0
        self.average_gencounter=0
        self.sum_specdispersal_allele=0
        self.average_speccounter=0
        
    def add_gendispersal_allele(self,allele):
        self.sum_gendispersal_allele+=allele
        self.average_gencounter+=1    

    def average_gendispersal(self):
        if self.average_gencounter>0:
            average=self.sum_gendispersal_allele/float(self.average_gencounter)
        else: average=0
        return average
        
    def reset_densloc(self):
        del self.densloc[:]        
        for x in xrange(self.max_x):
            self.densloc.append([0] * self.max_y)           #alles op 0 zetten
        
    def densloc_check(self,x,y):
        return self.densloc[x][y]

    def inc_densloc(self,x,y): 
        self.densloc[x][y]+=1

    def kernelsample(self,kernelshape):
        d=int(round(np.random.normal(0,kernelshape)))
        return d
        

class individual:
    def __init__(self, x,
                       y, 
                       kernelshape,
                       male,
                       paternal_gene,
                       maternal_gene):
        self.random_gene = np.random.random()
        self.paternal_gene = paternal_gene
        self.maternal_gene = maternal_gene
        self.male = male
        self.x=x
        self.y=y
        self.kernelshape=kernelshape
        self.drawing=None


class Simulation:
    '''This class coordinates the whole simulation. First a landscape is created, then a number of
    plants are initialized, then for a certain number of generations plants are told to move (disperse)
    and reproduce if they survive'''
    
    def __init__(self,
                 max_x,
                 max_y,
                 maxtime,
                 startpop,
                 equilibrium_density,
                 mutationrate,
                 p,
                 h,
                 sexratio):
        self.sexratio=sexratio
        self.max_x=max_x
        self.max_y=max_y
        self.p=p
        self.h=h
        self.startpop=startpop
        self.maxtime=maxtime
        self.equilibrium_density=equilibrium_density
        self.mutationrate=mutationrate
        self.densloc = Density(self.max_x,self.max_y)
        self.individuals=[]
        Simulation.start(self)
        
    def initialize_individuals(self):
        for c in xrange(self.startpop):
            if np.random.random()>self.sexratio:
                male=True
            else:
                male=False            
            while True:
                x=np.random.randint(0,CW_SIZE)
                y=np.random.randint(0,self.max_y)
                if self.environment.suitability_check(x,y):    #blijven herhalen tot plant in goed habitat gedropt wordt
                    break
            kernelshape=np.random.uniform(1,2)
            paternal_gene, maternal_gene='NA'
            self.individuals.append(individual(x,y,
                                               kernelshape,
                                               male,
                                               paternal_gene, 
                                               maternal_gene)) 
        
    def start(self):            
        global cw_loc
        cw_speed=0.4
        self.history=[] # a list that will store all indviduals that ever live with their location, father 'name' (random number), mother 'name, own 'name
        self.environment=Landscape(self.max_x,self.max_y,self.h,self.p)        
        self.environment.draw_landscape()
        self.densloc.reset_densloc()
        Simulation.initialize_individuals(self) 
        cw_loc=CW_SIZE/2
        self.time=0
        self.survivors=[]
        self.max_kernel = 1
        self.min_kernel = 1
        while (cw_loc<self.max_x-CW_SIZE/2 and len(self.individuals)>0):# or t<270:
            self.history.append([])
            #if cw_loc<self.max_x-CW_SIZE/2:
            cw_loc+=cw_speed
            self.densloc.reset_dispersal_sum()
            del self.survivors[:]
            for ind in self.individuals:
                self.move(ind)
            self.environment.update()
            del self.individuals[:]
            for ind in self.survivors:
                self.reproduction(ind)
            self.densloc.reset_densloc()     
            self.densloc.reset_male_genes()
            print 't=',self.time, "\t #g=",len(self.survivors)  
            #if self.time % 10 == 0:
            #    self.erase_extinct_lineages()   
            self.time+=1                
            #self.environment.weergave.reset()
            #self.environment.draw_landscape()
        self.erase_extinct_lineages()
        n_anchestors = []
        for population in self.history:
            n_anchestors.append(len(population))
        print n_anchestors
        self.draw_history()        
        #self.environment.reset()
        del self.survivors[:]
        del self.individuals[:]
            
        
    def move(self,ind):
        distance=self.densloc.kernelsample(ind.kernelshape)    #sample a distance from a gaussian distribution with mean 0 and a genetically determined variance i.e. dispersiveness
        angle=np.random.uniform(0,2*math.pi)  
        dx=int(round(distance*math.cos(angle)))
        dy=int(round(distance*math.sin(angle)))
        if ind.x+dx<0:
            dx=-dx
        ind.x+=dx
        ind.y=ind.y+dy
        withinlandscape=self.environment.boundary_check(ind.x,ind.y)
        if withinlandscape:
            withincw=self.environment.cw_check(ind.x)
            if withincw:
                suitable=self.environment.suitability_check(ind.x,ind.y)
                if suitable>=1:
                    ind.drawing=self.environment.create(ind.x,ind.y) #only draw individual if it survives
                    self.densloc.inc_densloc(ind.x,ind.y) #increase density if the individual lives
                    self.history[self.time].append((ind.x,
                                                    ind.y,
                                                    ind.kernelshape,
                                                    ind.random_gene,
                                                    ind.paternal_gene, 
                                                    ind.maternal_gene))
                    self.survivors.append(ind)
                    if ind.male:
                        self.densloc.add_male_gene(ind.x,ind.y,[ind.kernelshape,ind.random_gene])
                        
    def reproduction(self,ind):
        self.environment.delete(ind.drawing)   #delete the drawings of all parents (they are about to die)
        if ind.kernelshape > self.max_kernel:
            self.max_kernel = ind.kernelshape 
        elif ind.kernelshape < self.min_kernel:
            self.max_kernel = ind.kernelshape

        if not ind.male:
            male_gene=self.densloc.get_male_gene(ind.x,ind.y) 
            if male_gene!=None:          
                #b=6
                #a=((self.growthrate**(1/float(b))-1)/float(equilibrium_density))
                #mean_n_seeds=(self.growthrate)/(1+(a*(self.density.densloc_check(self.x,self.y)-1))**b)
                self.growthrate = 20                
                a=(self.growthrate-1)/float(self.equilibrium_density)
                mean_n_seeds=self.growthrate/(1+a*(self.densloc.densloc_check(ind.x,ind.y)-1))
                n_seeds=poisson(mean_n_seeds)
                offspring=[]
                new_gene=(ind.kernelshape+male_gene[0])/float(2)
                for r in xrange(n_seeds):   #For all seeds check if they mutate on one of the two genes and transfer the parent traits into the new individual (+ the mutation)
                    mutation=0
                    if np.random.random()<self.mutationrate:
                        mutation=np.random.uniform(-0.5,0.5)
                        if new_gene+mutation<0:
                            mutation=-mutation
                    if np.random.random()>self.sexratio:
                        male=True
                    else:
                        male=False                             
                    offspring.append(individual(ind.x, ind.y,
                                                new_gene+mutation,
                                                male,
                                                male_gene[1],
                                                ind.random_gene))        
                self.individuals.extend(offspring)
                
    def erase_extinct_lineages(self):
        if len(self.individuals)>0:
            old_history = self.history[:]
            del self.history[:]
            n_generation = len(old_history)        
            self.history = [[] for x in xrange(n_generation)]
            self.history[n_generation-1] = old_history[n_generation-1]
            while n_generation > 1:
                n_generation -= 1
                percentage = int(100*(len(old_history)-n_generation)/float(len(old_history)))
                if percentage%10 == 0:                    
                    print percentage,' %'
                for child in self.history[n_generation]:
                    father_found = False
                    mother_found = False    
                    for parent in old_history[n_generation-1]:
                        if child[4] == parent[3]:
                            self.history[n_generation-1].append(parent)
                            old_history[n_generation-1].remove(parent)
                            father_found = True
                        elif child[5] == parent[3]:
                            self.history[n_generation-1].append(parent)
                            old_history[n_generation-1].remove(parent)
                            mother_found = True
                        if father_found and mother_found:
                            break            
                    
    def draw_history(self):
        if len(self.individuals)>0:
            not_yet_drawn=np.zeros((self.max_x,self.max_y),dtype=int)
            #not_yet_drawn=[[True for y in xrange(self.max_y+1)] for x in xrange(self.max_x+1)]   
            n_generation = len(self.history)
            
            while n_generation > 1:
                n_generation -=1
                for child in self.history[n_generation]:
                    not_yet_drawn[child[0],child[1]] += 1
                    if not_yet_drawn[child[0],child[1]]<4:
                        for parent in self.history[n_generation-1]:
                            if (child[4] == parent[3]) or (child[5] == parent[3]):
                                break                                
                        if not((child[0] == parent[0]) and (child[1] == parent[1])):
                            self.environment.weergave.draw_line([child[0],child[1]],
                                                                [parent[0],parent[1]],
                                                                child[2],
                                                                [self.min_kernel,self.max_kernel])
                self.environment.update()

if __name__ == '__main__':                
    simulation=Simulation(max_x=128,
                          max_y=32,
                          maxtime=1000,
                          startpop=10000,
                          equilibrium_density=20,
                          mutationrate=0.001,
                          p=[0.6],
                          h=0.25,
                          sexratio=0.5)
                          
    #import cProfile
    #cProfile.run('simulation.start()')
    print 'Finished !'
tk.mainloop()
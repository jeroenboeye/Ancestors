# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 11:12:58 2012
models a plant species that can adapt to a gradient in deathrate at the cost of reproduction
in a fractal landscape
3.0 Haploid sexual reproduction
3.1 Diploid sexual reproduction
3.2 Selfing of females
3.3 Fixed dispersal, neutral markers for males & females to compare loss in diversity
3.4 plotting only the unique linneages
@author: Jboeye
"""

import numpy as np
from numpy.numarray.random_array import poisson
import fractal2
import Tkinter as tk
import math as math
import random as rnd
CW_SIZE=30
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
                   self.canvas.create_rectangle(self.zoom*x-factor, self.zoom*z-factor+self.max_y*self.zoom*2, self.zoom*x+factor, self.zoom*z+factor+self.max_y*self.zoom*2, outline='grey',fill='grey') 
    def update(self):
        self.canvas.update()
        
    def create (self, x, y, gene):
        zx=self.zoom*x
        zy=self.zoom*y
        factor=self.zoom/2     
        minimum=0
        maximum=1
        color_gene=gene-minimum
        color_gene/=float(maximum-minimum)
        if color_gene<0:
            color_gene=0
        elif color_gene>1:
            color_gene=1        
        r = int(255*color_gene)
        g = r
        b = 255-r
        rgb = r, g, b     
        Hex = '#%02x%02x%02x' % rgb        
        self.drawing=self.canvas.create_rectangle(zx-factor, zy-factor, zx+factor, zy+factor, outline=str(Hex),fill=str(Hex))
        return self.drawing
        
    def delete (self,drawing):
        self.drawing=drawing
        self.canvas.delete(self.drawing)
    
    def reset (self):
        self.canvas.delete(tk.ALL)
        
    def draw_line(self,pos1,pos2,gene,extremes):        
        minimum=extremes[0]
        maximum=extremes[1]
        color_gene=gene[0]-minimum
        if maximum==minimum:
            color_gene/=float(maximum)
        else:
            color_gene/=float(maximum-minimum)
        if color_gene<0:
            color_gene=0
        elif color_gene>1:
            color_gene=1        
        r = int(255*color_gene)
        g = r
        b = 255-r
        rgb = r, g, b     
        Hex = '#%02x%02x%02x' % rgb
        self.canvas.create_line(pos1[0]*self.zoom,
                                pos1[1]*self.zoom,
                                pos2[0]*self.zoom,
                                pos2[1]*self.zoom,
                                fill=str(Hex),width=2)


            
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
       
    def create (self, x, y,gene):
        self.drawing=self.weergave.create(x,y,gene)
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
        self.original_male_linneages = set()
        self.original_female_linneages = set()
        self.individuals=[]
        self.x_dict = dict()
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
            kernelshape1=(np.random.random(),x)
            #kernelshape2=np.random.random()
            kernelshape2=kernelshape1
            paternal_gene, maternal_gene='NA'
            self.individuals.append(individual(x,y,
                                               [kernelshape1,kernelshape2],
                                               male,
                                               paternal_gene, 
                                               maternal_gene))
            if male:
                self.original_male_linneages.add(kernelshape1)                                   
                self.original_male_linneages.add(kernelshape2)
            else:
                self.original_female_linneages.add(kernelshape1)                                   
                self.original_female_linneages.add(kernelshape2)  
            
    def start(self):            
        global cw_loc
        cw_speed=1
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
            if self.time % 10 == 0:
                lin1 = set()    
                for ind in self.individuals:
                    lin1.add(ind.kernelshape[0])
                    lin1.add(ind.kernelshape[1])   
                print 't=',self.time, "\t #g=",len(self.survivors), "unique_lin=",len(lin1)       
            #    self.erase_extinct_lineages()   
            self.time+=1                
            #self.environment.weergave.reset()
            #self.environment.draw_landscape()
            
        lin1 = set()    
        for ind in self.individuals:
            lin1.add(ind.kernelshape[0])
            lin1.add(ind.kernelshape[1])          
        print len(lin1)
        original_x = [0 for x in xrange(CW_SIZE)]
        for anch in lin1:
            original_x[anch[1]]+=1
        print original_x
        self.erase_extinct_lineages(lin1)
        #print len(self.history[0])
        self.erase_unsuccessfull_parents(lin1)
        n_anchestors = []
        for population in self.history:
            n_anchestors.append(len(population))
        print n_anchestors

        
        #print len(self.history[0]) * 2, len(lin)         
        #print 'female % surviving = ', 100*len(lin.intersection(self.original_female_linneages))/float(len(self.original_female_linneages)),' male % surviving = ', 100*len(lin.intersection(self.original_male_linneages))/float(len(self.original_male_linneages))
        self.draw_history(lin1)        
        #self.environment.reset()
        del self.survivors[:]
        del self.individuals[:]
            
        
    def move(self,ind):
        distance=self.densloc.kernelsample(2.0)
        angle=np.random.uniform(0,2*math.pi)  
        dx=int(round(distance*math.cos(angle)))
        dy=int(round(distance*math.sin(angle)))
        if ind.x+dx<0:
            dx=-dx
        ind.x+=dx
        if (ind.y+dy<0) or (ind.y+dy>self.max_y):
            dy = -dy
        ind.y+=dy
        #ind.y = ind.y%self.max_y
        withinlandscape=self.environment.boundary_check(ind.x,ind.y)
        if withinlandscape:
            withincw=self.environment.cw_check(ind.x)
            if withincw:
                suitable=self.environment.suitability_check(ind.x,ind.y)
                if suitable>=1:
                    ind.drawing=self.environment.create(ind.x,ind.y,ind.kernelshape[0][0]) #only draw individual if it survives
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
        if max(ind.kernelshape) > self.max_kernel:
            self.max_kernel = max(ind.kernelshape) 
        elif min(ind.kernelshape) < self.min_kernel:
            self.min_kernel = min(ind.kernelshape)

        if not ind.male:
            male_gene=self.densloc.get_male_gene(ind.x,ind.y) 
            if male_gene!=None:          
                #b=6
                #a=((self.growthrate**(1/float(b))-1)/float(equilibrium_density))
                #mean_n_seeds=(self.growthrate)/(1+(a*(self.density.densloc_check(self.x,self.y)-1))**b)
                self.growthrate = 5               
                a=(self.growthrate-1)/float(self.equilibrium_density)
                mean_n_seeds=self.growthrate/(1+a*(self.densloc.densloc_check(ind.x,ind.y)-1))
                n_seeds=poisson(mean_n_seeds)
                offspring=[]
                for r in xrange(n_seeds):   #For all seeds check if they mutate on one of the two genes and transfer the parent traits into the new individual (+ the mutation)
                    if rnd.random()<0.5:
                        mother_allele = ind.kernelshape[0]
                    else:
                        mother_allele = ind.kernelshape[1]
                    if rnd.random()<0.5:
                        father_allele = male_gene[0][0]
                    else:
                        father_allele = male_gene[0][1]                
                    if np.random.random()>self.sexratio:
                        male=True
                    else:
                        male=False                             
                    offspring.append(individual(ind.x, ind.y,
                                                [mother_allele,
                                                 father_allele],
                                                male,
                                                male_gene[1],
                                                ind.random_gene))        
                self.individuals.extend(offspring)
            else:
                #a=((self.growthrate**(1/float(b))-1)/float(equilibrium_density))
                #mean_n_seeds=(self.growthrate)/(1+(a*(self.density.densloc_check(self.x,self.y)-1))**b)
                self.growthrate = 2.5               
                a=(self.growthrate-1)/float(self.equilibrium_density)
                mean_n_seeds=self.growthrate/(1+a*(self.densloc.densloc_check(ind.x,ind.y)-1))
                n_seeds=poisson(mean_n_seeds)
                offspring=[]
                for r in xrange(n_seeds):   #For all seeds check if they mutate on one of the two genes and transfer the parent traits into the new individual (+ the mutation)
                    mother_allele1 = ind.kernelshape[0]
                    mother_allele2 = ind.kernelshape[1]      
                    if np.random.random()>self.sexratio:
                        male=True
                    else:
                        male=False                             
                    offspring.append(individual(ind.x, ind.y,
                                                [mother_allele1,
                                                 mother_allele2],
                                                male,
                                                ind.random_gene,
                                                ind.random_gene))        
                self.individuals.extend(offspring)                
                
    def erase_extinct_lineages(self,surviving_linneages):
        if len(self.individuals)>0:
            old_history = self.history[:]
            del self.history[:]                          
            n_generation = len(old_history)  
            self.history = [[] for x in xrange(n_generation)]
            for n in xrange(n_generation):
                generation = old_history[n]
                for ind in generation:
                    if (ind[2][0] in surviving_linneages) or (ind[2][1] in surviving_linneages):
                        self.history[n].append(ind)
                    
    def erase_unsuccessfull_parents(self,successfull_linneages):
        if len(self.individuals)>0:
            old_history = self.history[:]
            del self.history[:]
            n_generation = len(old_history)        
            self.history = [[] for x in xrange(n_generation)]
            self.history[n_generation-1] = old_history[n_generation-1]
            while n_generation > 1:
                n_generation -= 1
                percentage = 100*(len(old_history)-n_generation)/float(len(old_history))
                if int(percentage)%10 == 0:                    
                    print percentage,' %'
                for child in self.history[n_generation]:
                    if (child[2][0] not in successfull_linneages) or (child[2][1] not in successfull_linneages):    
                        for parent in old_history[n_generation-1]:
                            if (child[4] == parent[3]) or (child[5] == parent[3]):
                                self.history[n_generation-1].append(parent)
                                old_history[n_generation-1].remove(parent)
                                break      
                    else:
                        father_found = False
                        mother_found = False    
                        for parent in old_history[n_generation-1]:
                            if (child[4] == parent[3]):
                                self.history[n_generation-1].append(parent)
                                old_history[n_generation-1].remove(parent)
                                father_found = True
                            elif child[5] == parent[3]:
                                self.history[n_generation-1].append(parent)
                                old_history[n_generation-1].remove(parent)
                                mother_found = True
                            if father_found and mother_found:
                                break                                
    def draw_history(self,successfull_linneages):
        maximum = -1
        minimum = 2
        for linneage in successfull_linneages:
            if linneage[0]>maximum:
                maximum=linneage[0]
            if linneage[0]<minimum:
                minimum=linneage[0]
        if len(self.individuals)>0:
            not_yet_drawn=np.zeros((self.max_x,self.max_y),dtype=int)
            #not_yet_drawn=[[True for y in xrange(self.max_y+1)] for x in xrange(self.max_x+1)]   
            n_generation = len(self.history)
            while n_generation > 1:
                n_generation -=1
                for child in self.history[n_generation]:
                    not_yet_drawn[child[0],child[1]] += 1
                    if not_yet_drawn[child[0],child[1]]<6:
                        for parent in self.history[n_generation-1]:
                            if (child[4] == parent[3]) or (child[5] == parent[3]):
                                break                                
                        if not((child[0] == parent[0]) and (child[1] == parent[1])):
                            if child[2][0] in successfull_linneages:
                                self.environment.weergave.draw_line([child[0],child[1]],
                                                                    [parent[0],parent[1]],
                                                                    child[2][0],
                                                                    [minimum,maximum])
                            else:
                                self.environment.weergave.draw_line([child[0],child[1]],
                                                                    [parent[0],parent[1]],
                                                                    child[2][1],
                                                                    [minimum,maximum])                                        
#                            elif (child[5] == parent[3]):
#                                self.environment.weergave.draw_line([child[0],child[1]],
#                                                                    [parent[0],parent[1]],
#                                                                    child[2],
#                                                                    [self.min_kernel,self.max_kernel],
#                                                                    False)                                
                self.environment.update()

if __name__ == '__main__':                
    simulation=Simulation(max_x=128,
                          max_y=32,
                          maxtime=1000,
                          startpop=5000,
                          equilibrium_density=15,
                          mutationrate=0.001,
                          p=[1],
                          h=0.25,
                          sexratio=0.5)
                          
    #import cProfile
    #cProfile.run('simulation.start()')
    print 'Finished !'
tk.mainloop()
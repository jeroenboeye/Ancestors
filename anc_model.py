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
3.5 random rgb colors, 4 graphs showing snapshots and traceback of linneages.
4.0 Clean up code, remove fractal landscape, replace selfing by male gene searching in increasingly large square
4.1 Ready for first cluster simulations
4.2 No hard climate window but survival chance exponentially decreasing from optimum
@author: Jboeye
"""

import numpy as np
from numpy.numarray.random_array import poisson
import sys
import math as math
import random as rnd
class Visual:
    '''This class arranges the visual output. You can draw the landscape in black and white 
    according to suitability,    create the drawing of an individual, delete a drawing and 
    update the canvas.'''
    def __init__(self,max_x,max_y):
        import Tkinter as tk
        self.tk = tk
        self.zoom=3
        self.max_x=max_x
        self.max_y=max_y
        self.root = self.tk.Tk()
        self.canvas = self.tk.Canvas(self.root, width=self.max_y*self.zoom*4, height=self.max_x*self.zoom)      #create window
        self.canvas.pack()
        self.canvas.config(background='white')   #unsuitable habitat=background
    
    def update(self):
        self.canvas.update()
        
    def create (self, x, y, gene):
        zx=self.zoom*(self.max_x-x)
        zy=self.zoom*y
        factor=self.zoom/2    
        Hex = '#%02x%02x%02x' % gene        
        self.drawing=self.canvas.create_rectangle(zy-factor+self.max_y*self.zoom*3, zx-factor, zy+factor+self.max_y*self.zoom*3, zx+factor, outline=str(Hex),fill=str(Hex))
        return self.drawing
        
    def delete (self,drawing):
        self.drawing=drawing
        self.canvas.delete(self.drawing)
    
    def reset (self):
        self.canvas.delete(self.tk.ALL)
        
    def save (self,title):
        self.canvas.postscript(file=title)
        
    def draw_line(self,pos1,pos2,gene):        
        Hex = '#%02x%02x%02x' % gene[0]
        self.canvas.create_line(pos1[1]*self.zoom+self.max_y*self.zoom*3,
                                (self.max_x-pos1[0])*self.zoom,
                                pos2[1]*self.zoom+self.max_y*self.zoom*3,
                                (self.max_x-pos2[0])*self.zoom,
                                fill=str(Hex),width=1.5)
                                
    def create_snapshot(self,x,y,gene,count):
        zx=self.zoom*(self.max_x-x)
        zy=self.zoom*y
        factor=self.zoom/2
        Hex = '#%02x%02x%02x' % gene        
        self.canvas.create_rectangle(zy-factor+self.max_y*self.zoom*count,
                                     zx-factor, 
                                     zy+factor+self.max_y*self.zoom*count,
                                     zx+factor, outline=str(Hex),fill=str(Hex))      
        
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
        '''If no male genes are present in the current patch search for them
        in an increasingly large square'''
        genes = []
        if self.male_genes[x,y]!=None:
            genes = self.male_genes[x,y]
        else:
            radius = 1
            while (len(genes)==0) and (radius<30):                
                for dx in xrange(-radius,radius+1):
                    for dy in xrange(-radius,radius+1):
                        if ((0<=(x+dx)<self.max_x)
                        and (0<=(y+dy)<self.max_y)):
                            local_genes = self.male_genes[(x+dx)%self.max_x,(y+dy)%self.max_y]
                            if local_genes!=None:
                                genes.extend(local_genes)
                radius += 1
        if len(genes)==1:
            return genes[0]
        else:
            if len(genes)>1:
                winning_male_gene=rnd.choice(genes)
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
                 max_gen,
                 max_y,
                 cw_size,
                 cw_speed,
                 cw_speed_variance,                 
                 maxtime,
                 startpop,
                 equilibrium_density,
                 growthrate,
                 mutationrate,
                 sexratio,
                 gradient_steepness,
                 repeats,
                 visualisation,
                 asexual):
        self.cw_size = cw_size    
        self.cw_speed = cw_speed
        self.cw_speed_variance = cw_speed_variance
        self.gradient_steepness = gradient_steepness
        self.visualisation=visualisation
        self.growthrate=growthrate
        self.asexual = asexual
        self.max_x = int(max_gen * self.cw_speed + self.cw_size*2) #*2 for margin -> variance
        self.max_y=max_y
        if self.asexual:
            self.sexratio=1
            self.equilibrium_density=equilibrium_density/2.0
            title='asexual_y_%s_cw_%s_dens_%s_r_%s_maxgen_%s.txt'%(str(self.max_y),str(self.cw_size),str(self.equilibrium_density),str(int(self.growthrate)),str(max_gen))
        else:
            self.sexratio=sexratio
            self.equilibrium_density=equilibrium_density
            title='sexual_y_%s_cw_%s_dens_%s_r_%s_maxgen_%s.txt'%(str(self.max_y),str(self.cw_size),str(self.equilibrium_density),str(int(self.growthrate)),str(max_gen))
        self.maxtime=maxtime        
        self.mutationrate=mutationrate
        self.densloc = Density(self.max_x,self.max_y)
        self.original_male_linneages = set()
        self.original_female_linneages = set()
        self.individuals=[]    
        self.output=open(title,'a')        
        self.run(repeats,startpop,max_gen)
        self.output.close()
        
    def initialize_individuals(self,startpop):
        for c in xrange(startpop):
            if np.random.random()>self.sexratio:
                male=True
            else:
                male=False
            x=np.random.randint(0,self.cw_size)
            y=np.random.randint(0,self.max_y)
            r = rnd.randint(0,255)
            g = rnd.randint(0,255)
            b = rnd.randint(0,255)
            kernelshape1=((r,g,b),x)
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
            
    def run(self,repeats,startpop,max_gen):             
        for rep in xrange(repeats):
            self.initialize_individuals(startpop) 
            self.history=[] # a list that will store all indviduals that ever live with their location, father 'name' (random number), mother 'name, own 'name
            if self.visualisation:            
                self.visual=Visual(self.max_x,self.max_y)
            self.densloc.reset_densloc()
            self.cw_loc=self.cw_size/2
            self.time=0
            self.survivors=[]
            self.max_kernel = 1
            self.min_kernel = 1
            init_cw_loc = self.cw_loc
            t1 = 30
            t2 = 70
            t3 = 90
            while (max_gen>self.time and len(self.individuals)>0):            
                self.history.append([])
                #variance on speed of climate change, a minimal factor avoids error
                self.cw_loc=init_cw_loc + np.random.normal(self.cw_speed*self.time,
                                              self.cw_speed_variance+0.000000000001)
                self.densloc.reset_dispersal_sum()
                del self.survivors[:]
                for ind in self.individuals:
                    self.move(ind)
                if self.visualisation:            
                    self.visual.update()
                del self.individuals[:]
                for ind in self.survivors:
                    self.reproduction(ind)
                self.densloc.reset_densloc()     
                self.densloc.reset_male_genes()             
                if (self.visualisation) and (self.time % 10 == 0):
                    lin1 = set()
                    for ind in self.individuals:
                        lin1.add(ind.kernelshape[0])
                        lin1.add(ind.kernelshape[1])
                    print 'rep=', rep,'\t t=',self.time, "\t #g=",len(self.survivors), "\t unique_lin=",len(lin1)   
                if (self.visualisation) and ((self.time == t1) or (self.time == t2) or (self.time == t3)):
                    self.draw_snapshot(t1,t2,t3)
                self.time+=1
            lin1 = set()
            for ind in self.individuals:
                lin1.add(ind.kernelshape[0])
                lin1.add(ind.kernelshape[1])
            #print 'rep=', rep,'\t t=',self.time, "\t #g=",len(self.survivors), "\t unique_lin=",len(lin1)
            original_x = [0 for x in xrange(self.cw_size)]
            for anch in lin1:
                original_x[anch[1]]+=1
            #print original_x
            self.erase_extinct_lineages(lin1)
            self.erase_unsuccessful_parents(lin1)
            n_anchestors = []
            for population in self.history:
                n_anchestors.append(len(population))
            #print n_anchestors        
            #print len(self.history[0]) * 2, len(lin)         
            #print 'female % surviving = ', 100*len(lin.intersection(self.original_female_linneages))/float(len(self.original_female_linneages)),' male % surviving = ', 100*len(lin.intersection(self.original_male_linneages))/float(len(self.original_male_linneages))               
            self.output.write(str(self.cw_speed)+ "\t")
            self.output.write(str(self.cw_speed_variance)+ "\t")
            self.output.write(str(self.gradient_steepness)+ "\t")
            self.output.write(str(len(lin1))+ "\t")
            for location in original_x:
                self.output.write(str(location)+ "\t")
            self.output.write("\n")
            
            
            if self.visualisation:            
                self.draw_history(lin1)        
                self.visual.save('anchestor.eps')
            #self.visual.reset()
            del self.survivors[:]
            del self.individuals[:]
        
    def draw_snapshot(self,t1,t2,t3):
        if self.time == t1:
            t = 0
        elif self.time == t2:
            t = 1
        elif self.time == t3:
            t = 2
        surviving_linneages = set()
        for ind in self.individuals:
            surviving_linneages.add(ind.kernelshape[0])
            surviving_linneages.add(ind.kernelshape[1])
            self.visual.create_snapshot(ind.x,ind.y,ind.kernelshape[0][0],t)
        for anch in self.history[0]:
            if (anch[2][0] in surviving_linneages):
                self.visual.create_snapshot(anch[0],anch[1],anch[2][0][0],t)
            elif (anch[2][1] in surviving_linneages):
                self.visual.create_snapshot(anch[0],anch[1],anch[2][1][0],t)

    def boundary_check(self,x,y):
        survival=False
        if (0<x<self.max_x
        and 0<y<self.max_y):      #gradient            
            survival=True
        return survival
        
    def cw_check(self,x):
        '''procedure for hard climate window (everything or nothing)'''
        survival=False
        if x<=self.cw_loc+self.cw_size/2 and x>=self.cw_loc-self.cw_size/2:
            survival=True
        return survival   
        
    def cw_chec_exponential(self,x):
        '''procedure for survival p exponentially decreasing from optimum'''
        survival=False
        mortality_p=(math.exp((self.cw_loc-x)**2/(float(self.gradient_steepness)))-1)
        if rnd.random()>mortality_p:
            survival=True
        return survival           
        
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
        withinlandscape=self.boundary_check(ind.x,ind.y)
        if withinlandscape:
            withincw=self.cw_chec_exponential(ind.x)
            if withincw:                                    
                if self.visualisation:
                    ind.drawing=self.visual.create(ind.x,ind.y,ind.kernelshape[0][0]) #only draw individual if it survives
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
        if self.visualisation:
            self.visual.delete(ind.drawing)   #delete the drawings of all parents (they are about to die)
        if not ind.male:
            if self.asexual:
                male_gene = ind.kernelshape[0]
            else:
                male_gene=self.densloc.get_male_gene(ind.x,ind.y) 
            if male_gene!=None:          
                #b=6
                #a=((self.growthrate**(1/float(b))-1)/float(equilibrium_density))
                #mean_n_seeds=(self.growthrate)/(1+(a*(self.density.densloc_check(self.x,self.y)-1))**b)            
                a=(self.growthrate-1)/float(self.equilibrium_density)
                mean_n_seeds=self.growthrate/(1+a*(self.densloc.densloc_check(ind.x,ind.y)-1))
                n_seeds=poisson(mean_n_seeds)
                offspring=[]
                for r in xrange(n_seeds):   #For all seeds check if they mutate on one of the two genes and transfer the parent traits into the new individual (+ the mutation)
                    if self.asexual:
                        mother_allele = ind.kernelshape[1]
                        father_allele = ind.kernelshape[0]
                        male=False
                    else:
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
                    
    def erase_unsuccessful_parents(self,successful_linneages):
        if len(self.individuals)>0:
            old_history = self.history[:]
            del self.history[:]
            n_generation = len(old_history)
            self.history = [[] for x in xrange(n_generation)]
            #the last generation are all successful
            self.history[n_generation-1] = old_history[n_generation-1]
            while n_generation > 1:
                n_generation -= 1
                parent_names = set()
                for child in self.history[n_generation]:
                    parent_names.add(child[4]) #add father name
                    parent_names.add(child[5]) #add mother name                
                for parent in old_history[n_generation-1]:
                    if parent[3] in parent_names:
                        self.history[n_generation-1].append(parent)
                                         
    def draw_history(self,successful_linneages):
        if len(self.individuals)>0:
            n_drawn=np.zeros((self.max_x,self.max_y),dtype=int)
            n_generation = len(self.history)
            while n_generation > 1:
                n_generation -=1
                for child in self.history[n_generation]:
                    n_drawn[child[0],child[1]] += 1
                    if n_drawn[child[0],child[1]]<3:
                        for parent in self.history[n_generation-1]:
                            if (child[4] == parent[3]) or (child[5] == parent[3]):
                                break                                
                        if not((child[0] == parent[0]) and (child[1] == parent[1])):
                            if child[2][0] in successful_linneages:
                                self.visual.draw_line([child[0],child[1]],
                                                        [parent[0],parent[1]],
                                                        child[2][0])
                            elif child[2][1] in successful_linneages:
                                self.visual.draw_line([child[0],child[1]],
                                                        [parent[0],parent[1]],
                                                        child[2][1])                                          
                self.visual.update()

if __name__ == '__main__':                
    simulation=Simulation(max_gen=10,
                          max_y=50,
                          cw_size=60,
                          cw_speed=1,
                          cw_speed_variance=0.1,
                          maxtime=1000,
                          startpop=15000,
                          equilibrium_density=15,
                          growthrate=5.0,
                          mutationrate=0.001,
                          sexratio=0.5,
                          gradient_steepness=4000,
                          repeats=3,
                          visualisation=True,
                          asexual=True)
                          
    #import cProfile
    #cProfile.run('simulation.run()')
    #print 'Finished !'
#tk.mainloop()
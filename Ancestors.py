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
        g = r
        b = 255-r        
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
        self.genindividual=0
        for x in xrange(self.max_x):
            self.densloc.append([0] * self.max_y)           #alles op 0 zetten            
            
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

    def get_extremes(self,population):
        '''Calculates the minimal and maximal dispersal values of all anchestors'''
        first=True
        for individual in population:
            for anchestor in individual.anchestor_info:
                if first:
                    minimum=anchestor[2]
                    maximum=anchestor[2]
                    first=False
                else:
                    if anchestor[2]<minimum:
                        minimum=anchestor[2]
                    elif anchestor[2]>maximum:
                        maximum=anchestor[2]
        print 'Minimum d = ',round(minimum,2), 'Maximum d = ',round(maximum,2)
        return [minimum,maximum]
        
    def nr_anchestor_distribution(self,population):
        '''Counts the number off different anchestors per generation
        note that for now identical individuals on the same location 
        are only counted once for all (can be fixed by adding random 
        number to anchestor info tuple)'''
        first=True
        nr_anchestor_distr=[]
        for individual in population:
            if first:
                for anch in individual.anchestor_info:
                    nr_anchestor_distr.append(set([anch]))
                first=False
            else:
                count=0
                for anchestor in individual.anchestor_info:
                    nr_anchestor_distr[count].add(anchestor)
                    count+=1
        return nr_anchestor_distr
        
    def anchestors_per_generation(self,nr_anchestor_distr):
        anch_per_gen=[]
        for n in nr_anchestor_distr:
            anch_per_gen.append(len(n))
        return anch_per_gen
            
    def get_first_generation_xy(self,nr_anchestor_distr):
        first_generation_xy=[]
        for n in nr_anchestor_distr[0]:
            x=n[0]
            y=n[1]
            xy=(x,y)
            first_generation_xy.append(xy)
        return first_generation_xy

class Generalist():
    def __init__(self,x ,y ,kernelshape,environment,density,anchestor_info):
        self.x=x
        self.y=y
        self.kernelshape=kernelshape
        self.environment=environment
        self.density=density
        self.drawing=None
        self.growthrate=5
        self.anchestor_info=anchestor_info
        
    def move(self,survivors):
        dx=self.density.kernelsample(self.kernelshape)    #sample a distance from a gaussian distribution with mean 0 and a genetically determined variance i.e. dispersiveness
        dy=self.density.kernelsample(self.kernelshape)
        self.survivors=survivors
        self.x+=dx
        self.y+=dy
        withinlandscape=self.environment.boundary_check(self.x,self.y)
        if withinlandscape:
            withincw=self.environment.cw_check(self.x)
            if withincw:
                suitable=self.environment.suitability_check(self.x,self.y)
                if suitable>=1:
                    self.drawing=self.environment.create(self.x,self.y) #only draw individual if it survives
                    self.density.inc_densloc(self.x,self.y) #increase density if the individual lives
                    self.anchestor_info.append((self.x,self.y,self.kernelshape))#,random())) #only activate last part if you want to assess nr of ancestors
                    self.survivors.append(self)
       

    def reproduction(self,equilibrium_density,mutationrate,generalists):
        self.generalists=generalists
        self.density.add_gendispersal_allele(self.kernelshape)
        self.environment.delete(self.drawing)   #delete the drawings of all parents (they are about to die)
        #b=6
        #a=((self.growthrate**(1/float(b))-1)/float(equilibrium_density))
        #mean_n_seeds=(self.growthrate)/(1+(a*(self.density.densloc_check(self.x,self.y)-1))**b)
        a=(self.growthrate-1)/float(equilibrium_density)
        mean_n_seeds=self.growthrate/(1+a*(self.density.densloc_check(self.x,self.y)-1))
        n_seeds=poisson(mean_n_seeds)
        offspring=[]
        for r in xrange(n_seeds):   #For all seeds check if they mutate on one of the two genes and transfer the parent traits into the new individual (+ the mutation)
            mutation=0
            if np.random.random()<mutationrate:
                mutation=np.random.uniform(-0.5,0.5)
                if self.kernelshape+mutation<0:
                    mutation=-mutation
            anchestor_info=self.anchestor_info[:]
            offspring.append(Generalist(self.x, self.y,self.kernelshape+mutation,self.environment,self.density,anchestor_info))        
        self.generalists.extend(offspring)
        


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
                 h):                     
        self.max_x=max_x
        self.max_y=max_y
        self.p=p
        self.h=h
        self.startpop=startpop
        self.maxtime=maxtime
        self.equilibrium_density=equilibrium_density
        self.mutationrate=mutationrate
        self.densloc = Density(self.max_x,self.max_y)
        self.generalists=[]
        pratio=str(self.p)
        hurst=str(self.h)
        title='Anchestors pratio %s h %s.out'%(pratio,hurst)
        self.output=open(title,'w')
        Simulation.start(self)
        self.output.close()
        
    def initialize_generalists(self):
        for c in xrange(self.startpop):
            while True:
                x=np.random.randint(0,CW_SIZE)
                y=np.random.randint(0,self.max_y)
                if self.environment.suitability_check(x,y):    #blijven herhalen tot plant in goed habitat gedropt wordt
                    break
            kernelshape=np.random.uniform(1,3)
            anchestor_info=[]
            self.generalists.append(Generalist(x,y,kernelshape,self.environment,self.densloc,anchestor_info)) 
        
    def start(self):            
        global cw_loc
        cw_speed=1
        self.first_generation_xy_list=[]
        replications=1
        self.environment=Landscape(self.max_x,self.max_y,self.h,self.p)
        for r in xrange(replications):            
            self.environment.draw_landscape()
            self.densloc.reset_densloc()
            Simulation.initialize_generalists(self) 
            cw_loc=CW_SIZE/2
            t=0
            self.survivors=[]
            while (cw_loc<self.max_x-CW_SIZE/2 and len(self.generalists)>0):# or t<270:
                #if cw_loc<self.max_x-CW_SIZE/2:
                cw_loc+=cw_speed
                t+=1
                self.densloc.reset_dispersal_sum()  
                del self.survivors[:]
                for ind in self.generalists:
                    ind.move(self.survivors)
                self.environment.update()
                del self.generalists[:]
                for c in self.survivors:
                    c.reproduction(self.equilibrium_density,self.mutationrate,self.generalists)
                self.densloc.reset_densloc()            
                average_gen_disp=round(self.densloc.average_gendispersal(),2)
    
                print 't=',t, "\t #g=",len(self.survivors), "\t d=",average_gen_disp
                #self.environment.weergave.reset()
                #self.environment.draw_landscape()
            self.calculate_anchestor_info(self.survivors)
            self.draw_history(self.survivors)
            #self.environment.reset()
            del self.survivors[:]
            del self.generalists[:]

            
    def calculate_anchestor_info(self,population):
        if len(population)>0:
            nr_anchestor_distribution=self.densloc.nr_anchestor_distribution(population)
            first_generation_xy=self.densloc.get_first_generation_xy(nr_anchestor_distribution)
            self.first_generation_xy_list.extend(first_generation_xy)
            #anch_per_gen=self.densloc.anchestors_per_generation(self.nr_anchestor_distribution)
            
    def draw_history(self,population):
        if len(population)>0:
            not_yet_drawn=[[True for y in xrange(self.max_y+1)] for x in xrange(self.max_x+1)]
            extremes=self.densloc.get_extremes(population)       
            for ind in population:
                first=True
                for ancestor in ind.anchestor_info:
                    if first:
                        previous=ancestor
                        first=False
                    else:
                        if not_yet_drawn[ancestor[0]][ancestor[1]]:
                            self.environment.weergave.draw_line([previous[0],previous[1]],[ancestor[0],ancestor[1]],previous[2],extremes)
                            not_yet_drawn[ancestor[0]][ancestor[1]]=False
                        previous=ancestor

          
            

if __name__ == '__main__':                
    simulation=Simulation(max_x=128,
                          max_y=32,
                          maxtime=1000,
                          startpop=1000,
                          equilibrium_density=10,
                          mutationrate=0.02,
                          p=[1],
                          h=0.3)
                          
    #import cProfile
    #cProfile.run('simulation.start()')
    #print 'Finished !'
tk.mainloop()
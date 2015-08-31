# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 12:44:38 2012

@author: Jboeye
"""

class Generalist():
    def __init__(self,x ,y ,kernelshape,anchestor_info):
        self.x=x
        self.y=y
        self.kernelshape=kernelshape
        self.anchestor_info=anchestor_info
        
    def move(self):
        self.x+=self.kernelshape
        self.y+=self.kernelshape
        self.anchestor_info.append([self.x,self.y,self.kernelshape])

   
individual=Generalist(0,0,1,[])         
for t in xrange(20):
    individual.move()
print individual.anchestor_info
print len(individual.anchestor_info)

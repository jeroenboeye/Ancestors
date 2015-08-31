# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 14:23:08 2011

@author: Jboeye
"""
from numpy import *
from numpy.random import *
import math
#import Tkinter as tk
#root = tk.Tk()
seedvalue=10
variance = 1

def create(max_x=1024,max_z=128,h=0.5,p=[0.5],gradient=False): #max_x and max_z must be a power of two (e.g. 2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)
    """
    Enter max_x, max_z, hurst exponent, proportion of habitat (List!) and whether
    there is a gradient in autocorrelation. The p value must be given in a list 
    format since this allows to select several qualities of habitat. The lowest 
    value of p stands for the proportion of highest quality patches and vice versa. 
    max_x and max_z must both be a power of 2. Max_x can be larger than max_z in
    case you want a long landscape. Although, max_x must be a power of 2 times
    longer than max_z (2,4,8..). The Hurst exponent and proportions of habitat
    must be values between 0 and 1.
    """
#    canvas = tk.Canvas(root, width=max_x, height=max_z)      #create window
#    canvas.pack()
#    canvas.config(background='white')   #good landscape=background
    
    def diamondstep (iteration, x, z, standev):
          global matrix
          noemer=2**iteration
          rechts=x+((max_z)/noemer)
          links =x-((max_z)/noemer)
          boven =z+((max_z)/noemer)
          onder =z-((max_z)/noemer)
          matrix[x,z]=normal(((matrix[rechts,boven]+matrix[links,boven]+matrix[rechts,onder]+matrix[links,onder])/4),standev)
    
    
    def squarestep (iteration, x, z, standev,extra_x):
          global matrix
          noemer=2**(iteration-1)
          rechts=x+(((max_z)/2)/noemer)
          links =x-(((max_z)/2)/noemer)
          boven =z+(((max_z)/2)/noemer)
          onder =z-(((max_z)/2)/noemer)
          if (extra_x < x < extra_x+max_z) and (0<z<max_z):  # voor als ze ergens in het midden liggen
              matrix[x,z]=normal(((matrix[rechts,z]+matrix[links,z]+matrix[x,boven]+matrix[x,onder])/4),standev)
          elif (z==0):
              matrix[x,z]=normal(((matrix[rechts,z]+matrix[links,z]+matrix[x,boven])/3),standev)
          elif (z==max_z):
              matrix[x,z]=normal(((matrix[rechts,z]+matrix[links,z]+matrix[x,onder])/3),standev)
          elif (x==0):
              matrix[x,z]=normal(((matrix[rechts,z]+matrix[x,boven]+matrix[x,onder])/3),standev)   # deze 4 voor als ze aan de kant liggen
          elif (x==extra_x+max_z):
              matrix[x,z]=normal(((matrix[links,z]+matrix[x,boven]+matrix[x,onder])/3),standev)
    
    #create matrix and set all to zero
    global matrix, variance
    matrix=zeros((max_x+1,max_z+1),dtype=float)
    #set corners to seedvalue
    for x in [0,max_x]:
        for z in [0,max_z]:
            matrix[x,z]=seedvalue
            
    #Divide long landscape in smaller squares (actually you only set "seed" values here, the actual division happens lower)
    n_iterations=int(math.log(max_x/max_z,2))
    for i in xrange(n_iterations):
        iteration=i+1
        noemer=2**(iteration)
        v=variance/float(4)
        sd=sqrt(v)
        for z in  [0,max_z]:
            x=max_x/noemer
            while x<max_x:
                matrix[x,z]=normal(((matrix[x-max_x/noemer,z]+matrix[x+max_x/noemer,z])/2),sd)
                x+=2*max_x/noemer
                
    #Diamond-square algorithm for each individual square
    maxiterations=int(math.log(max_z,2))
    squareiteration=max_x/max_z
    extra_x=0
    for s in xrange (squareiteration):
        v=variance
        for i in xrange(maxiterations):
            iteration=i+1                       #Because i starts at 0    
            noemer = 2**(iteration)             #2,4,8,16....
            v=v/(2**(2*h))                      #variantie berekenen adhv hurst exponent, v wordt kleiner per iteratie
            sd=sqrt(v)                          #standaard deviatie (nodig voor formule normal())
            r = (max_z)/noemer           #beginnen op eerste rij
            while r<max_z:                      #beginnen op eerste kolom en ook op de volgende rij op de eerste kolom beginnen
                 c = extra_x+(max_z)/noemer      #Terug naar eerste kolom (op een andere rij)
                 while c<extra_x+max_z:                 #herhaal alles hierbinnen over alle rijen
                    if iteration>=maxiterations-2 and gradient==True and c>0: 
                        sd=sqrt(v)*(c/max_x) #voor de laatste stappen wordt de stochastisiteit adhv gradient berekend (let op maxiterations-2 = laatste 3 stappen)
                    diamondstep(iteration,c,r,sd)
                    #vanaf hier squarestep, 4 maal, naast de locatie van de diamondstep
                    boven =r+(max_z)/noemer
                    squarestep(iteration,c,boven,sd,extra_x)
                    onder =r-(max_z)/noemer
                    squarestep(iteration,c,onder,sd,extra_x)
                    rechts=c+(max_z)/noemer
                    squarestep(iteration,rechts,r,sd,extra_x)
                    links =c-(max_z)/noemer
                    squarestep(iteration,links,r,sd,extra_x)
                    # einde squarestep
                    c += max_z/(noemer/2)  # volgende kolom (rekening houdend met iteratie)
                 r += max_z/(noemer/2)      # volgende rij (rekening houdend met iteratie)
        extra_x+=max_z
   
    proportionmatrix=p   #check difference between [] and ()
    proportionmatrix.sort()
    proportionmatrix.reverse()
    sortedmatrix=sort(matrix,axis=None)
    old_treshold=-10000
    integermatrix=[[len(proportionmatrix) for z in xrange (max_z+1)] for x in xrange (max_x+1)]
    suitability_value=0

    for c in xrange (len(proportionmatrix)):
        index=round((1-proportionmatrix[c])*(len(sortedmatrix)))-1
        if index<0: index=0
        treshold=sortedmatrix[index]
        for x in xrange (max_x+1): 
            for z in xrange (max_z+1):
                if old_treshold<matrix[x][z]<treshold:
                    integermatrix[x][z]=suitability_value
        old_treshold=treshold
        suitability_value+=1
    #outputmatrix[x,0:max_z+1]=integercolumn

#    draw the landscape
#    for x in range (max_x+1):
#        for z in range (max_z+1):
#            if booleanmatrix[x,z]==True:
#                canvas.create_rectangle(x-0.5, z-0.5, x+0.5, z+0.5, outline='black',fill='black')
#    canvas.update()
#    root.mainloop()

    #print outputmatrix
    return(integermatrix)
#import cProfile
#cProfile.run('create()')
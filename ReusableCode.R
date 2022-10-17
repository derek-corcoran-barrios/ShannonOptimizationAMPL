## ---- LoadPackages --------

library(janitor)
library(terra)
library(magrittr)
library(geodata)
library(sf)
library(tidyverse)
library(mregions)
library(tidyterra)

## ---- TestCode --------

mtcars %>% 
  ggplot(aes(x = hp, y = mpg)) +
  geom_point()

## ---- AMPLModel --------

set V;   # vertex set of the spatial graph
set SP; #Set of Species
param L; # Landuse types
param budget >= 0;
param u {SP, V,1..L} ; # node presence absence pred for species SP, in landuse in time L
param c {V,1..L} ; # node carbon content, in landuse in time L
param m {V}; # node cost (all 1) this is for calculating the budget
var y {v in V,l in 1..L} >= 0; # Desition on which landuse to use for node V
var z {v in V} >= 0; # Desition to preserve or not node

maximize ShanonDiv: sum {s in SP, v in V, l in 0..L} m[v]*y[s,v,l]*c[v,l];

subj to MaxBudget {v in V, l in 1..L}:
  budget <= sum {(v) in V} m[v]*z[v];

## ---- AMPTemplate --------

set V;   # vertex set of the spacial graph
set E within {V,V};   # edge set of spacial graph. It should contain all loops and the pair (u,v) and (v,u) for an edge uv in the graph.
set SP; #Set of Species

param T; # time horizon

param nchains{SP} >= 0;
param c {V}; #Costo
param u {SP, V,0..T} ; # node biomass capacity for species SP, in node V in time T

var y {s in SP,v in V,t in 0..T} >= 0; # flow on sites
var r {SP,E,1..T} >= 0; # inner flow on transitions

###########OPTION 1###############

minimize Quad_Cost: sum {s in SP, v in V, t in 0..T} c[v]*y[s,v,t]*c[v]*y[s,v,t];

subj to Initial_Flow{s in SP}:
  sum{v in V} y[s,v,0] = nchains[s];
  
  subj to Final_Flow{s in SP}:
    sum{v in V} y[s,v,7] = nchains[s];
    
    ###########OPTION 2###############
    
    #maximize Flow: sum {v in V} y[v,0];
    
    ################################
    
    subj to Flow_ConservationIN {s in SP, w in V, t in 1..T}:
      y[s,w,t] = sum {(v,w) in E} r[s,v,w,t];
      
      subj to Flow_ConservationOUT {s in SP,w in V, t in 0..T-1}:
        y[s,w,t] = sum {(w,v) in E} r[s,w,v,t+1];
        
        subj to Flow_Capacity {s in SP, w in V, t in 0..T}:
          y[s,w,t] <= u[s,w,t];
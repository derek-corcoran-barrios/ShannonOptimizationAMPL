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

## ---- AMPLModel1IdealWorld --------

set Cells;   # vertex set of the spatial graph
set Species; #Set of Species
param Landuses; # Landuse types
param LanduseSuitability {Cells,Landuses}; #If cell is suitable or not in each cell
param SpeciesSuitability {Species, Cells,Landuses}; # node presence absence pred for species SP, in landuse in time L
param Cost {Cells}; # node cost (all 1) this is for calculating the budget
var LanduseDecision {c in Cells,l in Landuses} binary; #else >= 0  decision on which landuse to use for cell Cell

minimize InvShanonDiv: sum{s in Species, c in Cells, l in Landuses} LanduseDecision[c,l]*LanduseSuitability[v,l]*SpeciesSuitability[s,v,l]*log(LanduseDecision[c,l]*LanduseSuitability[v,l]*SpeciesSuitability[s,v,l]);

subj to PropotionalUse{c in Cells}:
  sum{l in Landuses} LanduseDecision[c,l] = 1;
## ---- AMPLModel2TransitionConstrained --------  
  
set Cells;   # vertex set of the spatial graph
set Species; #Set of Species
param Landuses; # Landuse types
param LanduseSuitability {Cells,Landuses}; #If cell is suitable or not in each cell
param TransitionCost {Cells,Landuses}; #Cost of changing the cell landuse, set to 0 
                                      #if it stays the same, 1 if it is a change or infitinte if we want to set landuse
param MaxTransitionCost >= 0;
param SpeciesSuitability {Species, Cells,Landuses}; # node presence absence pred for species SP, in landuse in time L
param Cost {Cells}; # node cost (all 1) this is for calculating the budget
var LanduseDecision {c in Cells,l in Landuses} binary; #else >= 0  decision on which landuse to use for cell Cell
  
minimize InvShanonDiv: sum{s in Species, c in Cells, l in Landuses} LanduseDecision[c,l]*LanduseSuitability[v,l]*SpeciesSuitability[s,v,l]*log(LanduseDecision[c,l]*LanduseSuitability[v,l]*SpeciesSuitability[s,v,l]);
  
subj to PropotionalUse{c in Cells}:
  sum{l in Landuses} LanduseDecision[c,l] = 1;
  
subj to MaxTransition:
  sum{l in Landuses, c in Cells} LanduseDecision[c,l]*TransitionCost[c,l] <= MaxTransitionCost;


## ---- AMPLModel --------

set Cells;   # vertex set of the spatial graph
set Species; #Set of Species
param Landuses; # Landuse types
param budget >= 0;
param LanduseSuitability {Cells,Landuses}; #If cell is suitable or not in each cell
param SpeciesSuitability {Species, Cells,Landuses}; # node presence absence pred for species SP, in landuse in time L
param CarbonContent {Cells,Landuses} ; # node carbon content, in landuse in time L
param Cost {Cells}; # node cost (all 1) this is for calculating the budget
var LanduseDecision {c in Cells,l in Landuses} binary; #else >= 0  decision on which landuse to use for cell Cell
var ConervationDecision {c in Cells} >= 0; # Desition to preserve or not node

minimize InvShanonDiv: sum{s in Species, c in Cells, l in Landuses} LanduseDecision[c,l]*LanduseSuitability[v,l]*SpeciesSuitability[s,v,l]*log(LanduseDecision[c,l]*LanduseSuitability[v,l]*SpeciesSuitability[s,v,l]);

subj to MaxBudget {c in Cells, l in Landuses}:
  budget <= sum {c in Cells} Cost[c]*ConervationDecision[c];

subj to PropotionalUse{c in Cells}:
  sum{l in Landuses} LanduseDecision[c,l] = 1;

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
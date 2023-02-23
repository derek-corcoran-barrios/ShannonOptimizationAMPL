Optimization for Shannon diversity in AMPL
================
Derek Corcoran
23/02, 2023

- <a href="#1-models" id="toc-1-models">1 Models</a>
  - <a href="#11-model-1-ideal-world" id="toc-11-model-1-ideal-world">1.1
    Model 1 ideal world</a>
  - <a href="#12-model-2-add-cost-of-transition"
    id="toc-12-model-2-add-cost-of-transition">1.2 Model 2 add “Cost” of
    transition</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

The goal of this repository is to explain and prepare models for the
optimization of geographical areas in order to maximize diversity
(Shannon or other indexes) using AMPL.

# 1 Models

## 1.1 Model 1 ideal world

In this model we try to decipher the ideal landuse configuration of
Denmark taking only into account landuse suitability and diversity, we
just want to maximize Shannon’s diversity and descide on the best
composition of landuse for denmark given that.

``` r
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
```

The only constrain we have here is the landuse composition of Denmark to
maximize the diversity as defined by an index, in this example, it is
Shannon’s diversity index, but this can be easily be changed.

We can make some smart and easy “Constraints” such as not making
possible to change the landuse from city to another use.

The mathematical equations for this would be

$$\begin{equation*}
\begin{aligned}
& \underset{x}{\text{minimize}}
& & -H=\sum_{s=1}^{Species}\sum_{c=1}^{Cells}\sum_{l=1}^{Landuses} LanduseDecision_{v,l}\times LanduseSuitability{c,l}\times SpeciesSuitability{s,v,l} log_2\,(LanduseDecision_{v,l}\times LanduseSuitability{c,l}\times SpeciesSuitability{s,v,l}) \\
& \text{subject to Propotional Use}
& &
  \sum_{1=l}^{Landuses} LanduseDecision_l = 1
\end{aligned}
\end{equation*}$$

## 1.2 Model 2 add “Cost” of transition

One of the main issues with model 1 is that this assumes that we could
change all the cells in Denmark, this could be easily changed by adding
a parameter called transition cost, which could be set to 0 if the
landuse stays the same and 1 if it changes for sake of simplicity, but
more information could be add to it by using actual transition costs or
even mantainance costs. That model looks as follows

``` r
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
```

We also add a maximum transition cost, which is the maximum number of
cells that can change their landuse, this could be set to the number of
cells that are 10% of Denmark for example

## Using the inventory data 

1. What’s the difference between NL%IED_INIT_MODE = 5 or = 6
    > - it's not important
    > - only use 0/6

2. How will the inventory plot data links with the patch in size/area 
    > - Only the relivant fraction of stem density will be considered in inventory data. In such case, it can be easily transfered into an average value in the  simulition: `site[fraction of every patch] -> patch[fraction for every cohort] - cohort[every individual tree will be set as one individual cohort]`
    > - The results in the simulation of ED2 will be averaged into squre meter level like: $GPP = \frac{3kgC}{1year \times ha}$. 
    > - It can be scaled up to plot/regional level by multipling area size.

3. The accurancy of latitude and longtitude.
    >  - It doesn't matter.
    >  - Check the climate driver informtion, and make sure it is close to the site.

## TLS_ED2.2

1. TLS_ED model: what’s the difference between closed canopy and finite crown.
    > - Closed crowns are evenly spread throughout the patch area and cohorts are stacked on top of each other: 所有的chort的冠层都在养地上平铺开来，最上层的chort会接收到所有的直接辐射；而下层cohort冠层也是铺满整个样地但是只能接收间接辐射。这样传送到林下植物的光通常较少。
    > - Finite Cohorts have a finite radius and are stacked on top of each other (Dietze et al., 2008): 冠层的面积依据 `DBH` 来确定，一般不会占满整个样地，这样就还是会有直射光抵达下层的植物。这样传送到林下植物的光通常较多。

2. How the TLS-drived parameters works in plant allometry?
    > - allometry parameters will be set in the `config.xml` file

4. How to let SLA and Vcmax change across light level? (4787)
    > - the setting `TRAIT_PLASTICITY_SCHEME`, let the traits can be updated monthly
    
6. Have you compared the difference between one PFT and two PFTs?
    > - one PFT is more like to be the average of two PFTs;
    > - but have not compare the results with the real data pairly.

7. Is there some codes to convert the PFT parameters into .xml file.
    > - write.config.xml

## PEcAn - to see whether the change improves the robustness of model

1. How to calculate the posterior distributions of parameters
    > - there is a peice of codes that extract the main function of pecan `generate_SA.R`.
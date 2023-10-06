# Inputs of model

## Initial state
> Every simulaitonrequires an initial state for forest structure and compositon. 
> 
> To initialize a plant community from inventory data, one must have either the diameter at breast height of every individual or the stem density of different diameter size classes, along with plant functional type identification and location; in addition, necromass from the litter layer, woody debris, and soil organic carbon are needed. 
> Alternatively, initial conditions can be obtained from airborne lidar measurements (Antonarakis et al., 2011, 2014) or a prescribed near-bare-ground condition may be used for long-term spin up simulations. Previous simulations can be used as initial conditions as well.

## Edaphic conditions
> THe user must also provide soil characteistics such as:
> - total soil depth
> - total number of soil layers
> - the thickness of each layer
> - soil texture: can be read from standard data sets (e.g., Global Soil Data Task, 2000; Hengl et al., 2017) or provided directly by the user.
> - color, the bottom soil boundary condition (bedrock, reduced drainage, free drainage, or permanent water table): must be provided directly by the user as of ED-2.2.
>
> In addition, simulations with multiple sites per polygon also need to provide the fractional areas of each site and the mean soil texture class, soil depth, slope, aspect, elevation, and topographic moisture index of each site.




## Atmospheric conditions
> 
> Meteorological conditions needed to drive ED2.2 including:
> - temperature
> - specific humidity
> - CO2 molar fraction: it should be provide comparable temporal and spatial resolution to other meteorlogical drivers. Alternatively the meteorological forcing (including CO2)  may be provided directly by BRAMS
> - pressure of the air above the conopy
> - precipitation
> - incoming solar (shortwave) irradiance (radiation flux)
> - incomming thermal (longwave) irradiance (a few meters above the canopy)
> 
> beter measurements 0.5-6h

1. single location drivers
2. grided meteorological drivers


## Plant functional types
> The user must specify which PFTs are allowed to occur in any given simulation.
> ED-2.2 has a list of default PFTs, with parameters described in [github wiki](https://github.com/EDmodel/ED2/wiki/Plant-functional-types).
> The user can modify the [parameters](https://github.com/EDmodel/ED2/wiki/PFT-parameters) of existing PFTs or define new PFTs through an extensible markup language ([XML](https://github.com/EDmodel/ED2/wiki/Model-parameters-and-xml-parameter-files)) file, which is read during the model initilization ([ED2IN](https://github.com/EDmodel/ED2/wiki/ED2IN-namelist)).

## Example of input file for BIC
1. VEG_DATABASE: land/water mask
    > vegetation database, used only to determine the land/water mask. Fill with the path and the prefix. 
    ```Shell
    NL%VEG_DATABASE = '/user/data/gent/gvo000/gvo00074/ED_common_data/veg_oge/OGE2_'
    ```

2. SOIL_DATABASE: soil type data
    > soil database, used to determine the soil type.  Fill with the path and the prefix. 
    ```Shell
    NL%SOIL_DATABASE = '/user/data/gent/gvo000/gvo00074/ED_common_data/soils/FAO/FAO_'
    ```
3. LU_DATABASE: land-use change disturbance rates
    > land-use change disturbance rates database, used only when IANTH_DISTURB is set to 1.  Fill with the path and the prefix. 
    ```Shell
    NL%LU_DATABASE = ''
    ```
    
4. PLANTATION_FILE: forest plantation fraction
    > Character string for the path to the forest plantation fraction file.  This is used only when IANTH_DISTURB is set to 1 and the user wants to simulate forest plantations. Otherwise, leave it empty (PLANTATION_FILE='')   
    ```Shell
    NL%PLANTATION_FILE = ''
    ```

5. THSUMS_DATABASE: phenology about chilling-degree and growing-degree days
    > input directory with dataset to initialise chilling-degree and growing-degree days, which is used to drive the cold-deciduous phenology (you must always provide this, even when your PFTs are not cold deciduous). 
    ```Shell
    NL%THSUMS_DATABASE = '/user/data/gent/gvo000/gvo00074/ED_common_data/ed_inputs/'
    ```

6. ED_MET_DRIVER_DB: meteorological driver instructions
    > File containing information for meteorological driver instructions (the "header" file).
    ```Shell
    NL%ED_MET_DRIVER_DB = '/data/gent/vo/000/gvo00074/felicien/R/climate.site/site.lat9.25N.lon79.75W/ED2/ED_MET_DRIVER_HEADER'
    ```

7. OBSTIME_DB: 
    > File containing times of desired IOOUTPUT Reference file: /ED/run/obstime_template.time
    ```Shell
    NL%SOILSTATE_DB = ''
    ```

8. SOILSTATE_DB: soil moisture and temperature infromation if ISOILSTATEINIT=1
    > If ISOILSTATEINIT=1, this variable specifies the full path of the file that contains soil moisture and temperature information. 
    ```Shell
    
    ```

9. SOILDEPTH_DB: soil moisture and temperature information if ISOILDEPTHFLG=1
    > If ISOILDEPTHFLG=1, this variable specifies the full path of the file that contains soil moisture and temperature information. 
    ```Shell
    NL%SOILDEPTH_DB = ''
    ```


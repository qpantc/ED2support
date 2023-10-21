# Meteorology driver

## Using Fluxtower as drivers data

1. required variables
    | Variables| Description and unit|
    | ------------- | ------------------- |
    | co2 (optional)|  surface CO2 mixing ratio (μmolCO2 molAir-1)|
    | dlwrf|  downward long wave radiation(Wm-2)|
    | nbdsf|  near infrared beam downward solar radiation (W m-2)|
    | nddsf|  near IR diffuse downward solar radiation (W m-2)|
    | vbdsf|  visible beam downward solar radiation (W m-2)|
    | vddsf|  visible diffuse downward solar radiation (W m-2)|
    | prate|  precipitation rate (kgH2O m-2 s-1)|
    | pres |  atmospheric pressure (Pa)|
    | hgt  |  geopotential height (m)|
    | ugrd |  zonal wind (m s-1)|
    | vgrd |  meridional wind (m s-1)|
    | sh   |  specific humidity (kgH2O kgAir-1)|
    | tmp  |  air temperature (K)|


2. Conversion
    the script `support/tower/tower.driver.R` can be used to convert the fluxnet data into HDF5 and generate the `ED_MET_DRIVER_DB`
    ```shell
    1
    /Users/quan/projects/ED2.2/EDsupport/tower/sites/FLX_GF-Guy2004_2020/FLX_GF-Guy2004_2020_
    2 2 1 1 -53 5
    12
    'hgt'   'tmp'  'pres'    'sh'  'ugrd'  'vgrd' 'prate' 'dlwrf' 'nbdsf' 'nddsf' 'vbdsf' 'vddsf'
    80 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800 1800
    4 1 1 1 1 1 0 1 1 1 1 1
    ```

## For generating regional drivers


# Phenology

There are three options to prescribe the phenology in ED:

0 – Calculate based on meteorology
1 – User prescribed
2 – MODIS derived
For the first of these cases (0), phenology is not a direct boundary condition, but it based on the phenological model of Botta et al. (2000) (cold-deciduous), Longo (2014) (drought-deciduous) or Kim et al. (2012) (light-driven). The second case (1) represents a hard-wired phenological timing based on observations. Contact David Medvigy or Alexander Antonarakis for additional information on how to prepare the files. The third case represents phenology driven by MODIS satellite data. Phenology data is assumed to be in the “phenology” folder of the directory given by IOPTINPT (i.e. Optimizer inputs) and assumes that files are names <year></year>.dat where year is hardwired to cover 2001-2004. MODIS phenology is driven by the function “modis_phenology”.

# Land use and plantation rotation

## Land use transition matrix 

## Forest plantation rotation 

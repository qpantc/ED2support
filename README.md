# Tutorial for ED2 at Ugent

## Compile ED2

1. download the model
    ```git
    git clone git@github.com:EDmodel/ED2.git
    ```

2. change the makefile in `ED2/ED/build/make`
    ```shell
    cd ED2/ED/build/make
    cp include.mk.intel include.mk.intel_hpc
    ```

3. Edit include.mk.intel_hpc: `input and output lib`
    >  mainly for input and output lib
    ```Makefile
    HDF5_INCS=-I/apps/gent/CO7/skylake-ib/software/HDF5/1.12.1-iimpi-2021b/include
    HDF5_LIBS=-lm -lz -L/apps/gent/CO7/skylake-ib/software/HDF5/1.12.1-iimpi-2021b/bin -lhdf5 -lhdf5_fortran -lhdf5_hl
    ```
4. load the module to complie
    ```shell
    ml purge; ml intel-compilers/2021.4.0 HDF5/1.12.1-iimpi-2021b UDUNITS/2.2.28-GCCcore-11.2.0; ulimit -s unlimited
    ```
5. cd to the build dir and compile
    > #-k E fast and efficiently run/ -k A will give some clue when it broken.
    ```
    cd .. 
    ./install.sh -k E -p intel_hpc
    ```

6. see what you get
    > 

## Preparing the driving data

1. download data prepare tool
    ```git
    git clone git@github.com:qpantc/ED2support.git
    ```
2. load R module, cd to script dir and run R
    > `ml av R` see module suitable for R
    ```shell
    ml purge; 
    ml R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0 
    
    cd ED2support/R
    R
    ```

3. run the data download and convert script
    ```R
    source("download.and.convert.input.ED2.R")
    ```
4. Edit ED2IN file
   - location of the simulation
   ```shell
   NL%N_POI = 1
   NL%POI_LAT = 9.25
   NL%POI_LON = -79.75
   NL%POI_RES = 1
   ```
   - start of the simulation
   ```shell
   NL%IMONTHA = 1
   NL%IDATEA = 1
   NL%IYEARA = 1901
   NL%ITIMEA = 0
   ```
   - end of the simulation
   ```shell
   NL%IMONTHZ = 1
   NL%IDATEZ = 1
   NL%IYEARZ = 1903
   NL%ITIMEZ = 0
   ```

    - driver location
    ```shell
    NL%VEG_DATABASE = '/user/data/gent/gvo000/gvo00074/ED_common_data/veg_oge/OGE2_'
    NL%SOIL_DATABASE = '/user/data/gent/gvo000/gvo00074/ED_common_data/soils/FAO/FAO_'
    NL%LU_DATABASE = ''
    NL%PLANTATION_FILE = ''
    NL%THSUMS_DATABASE = '/user/data/gent/gvo000/gvo00074/ED_common_data/ed_inputs/'
    NL%ED_MET_DRIVER_DB = '/user/gent/442/vsc44253/site.lat9.25N.lon79.75W/ED2/ED_MET_DRIVER_HEADER'
    NL%SOILSTATE_DB = ''
    NL%SOILDEPTH_DB = ''
    ```
    - output folder
    ```shell
    NL%FFILOUT = '/scratch/gent/vo/000/gvo00074/vsc44253/ED_results/BCI/analy/analysis'
    NL%SFILOUT = '/scratch/gent/vo/000/gvo00074/vsc44253/ED_results/BCI/hist/history'
   ```
   - empty the IEDCNFGF config
   ```shell
   NL%IEDCNFGF   = '/mypath/config.xml'
   NL%EVENT_FILE = '/mypath/event.xml'
   ```

   - SFILIN
   ```shell
   NL%SFILIN = ''
   ```




5. run the model
   > Job.sh file Or inside an interactive job
   > After load the needed modules
   ```
   ../build/ed_2.2-opt-master-8d4c3aff -f ED2IN
   ```

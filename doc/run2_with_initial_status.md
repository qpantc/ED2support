# NBG: Run the model from near-bare ground
## spinup: run the model initially from 1500 to 1900
1. preparing the crucnep drivers 
    > Follow the README.md

2. set start and end time
    ```R
    NL%IMONTHA = 1
    NL%IDATEA = 1
    NL%IYEARA = 1500
    NL%ITIMEA = 0


    NL%IMONTHZ = 1
    NL%IDATEZ = 1
    NL%IYEARZ = 1901
    NL%ITIMEZ = 0
    ```
3. set the recycled drivers
    ```R
    NL%METCYC1 = 1981
    NL%METCYCF = 1999
    ```

## Run the model historically from 1900:
1. set NL%IED_INIT_MODE
    `NL%IED_INIT_MODE=6`
    >- 0: 表示植被将用几棵幼苗(近裸地)初始化
    >- 6: 用于单个站点时,可以使用历史记录重新启动,常用于使用森林清单数据进行运行的情况

2. set historical file path in ED2IN
    > to the dir of historical file
    ```R
    NL%SFILIN = '/scratch/gent/vo/000/gvo00074/vsc44253/ED2.2/EDsupport/run/GF_Guy/outputs/histo/history'
    ```
3. set the continue date
    >  the continue date must be the start date that you have your drivers
    ```R
    NL%IYEARH = 2004
    NL%IMONTHH = 1
    NL%IDATEH = 1
    NL%ITIMEH = 0
    ```


# Census: Run the model with inventory data
模型的初始化就是给每个模拟单元中的sites,patch和cohort都设置相应的属性参数.
将森林调查数据作为输入,就是将每个地块对应上一个patch, patch上的每棵树对应一个cohort.

1. set NL%IED_INIT_MODE
    `NL%IED_INIT_MODE=6`
    >- 0: 表示植被将用几棵幼苗(近裸地)初始化
    >- 6: 用于单个站点时,可以使用历史记录重新启动,常用于使用森林清单数据进行运行的情况

2. preparing pss and css file
    - Patch file: *.pss
        > 每行该site里包含的一个patch,记录了patch的基本信息
        > area表明该patch占整体site的比例

        ```R
        "time" "patch" "trk" "age" "area" "water" "fsc" "stsc" "stsl" "ssc" "lai" "msn" "fsn" "nep" "gpp" "rh"
        2015 1 2 0 1 0 0.15 5 5 4.546 0 0 0 0 0 0
        ```
    - Cohort file: *.css
        > 每行中,每棵树都被当作一个cohort
        > n: Local (within-patch) stem density.

        ```R
        "time" "patch" "cohort" "dbh" "hite" "pft" "n" "bdead" "balive" "lai"
        2015 1 1 88 33.6156234911455 4 0.0025 0 0 0
        2015 1 2 45 24.0580642031917 3 0.0025 0 0 0
        2015 1 3 22 16.0251000958566 3 0.0025 0 0 0

        ......

        2015 1 119 17 13.7182437752212 3 0.0025 0 0 0
        2015 1 120 20 15.1374071779411 2 0.0025 0 0 0
        2015 1 121 20 15.1374071779411 2 0.0025 0 0 0
        ```
3. set file path in ED2IN
    only the shared part except .css/.pss
    ```
    NL%SFILIN  = './csspss/lianas_s83_default.lat5.000lon-52.000'
    ```


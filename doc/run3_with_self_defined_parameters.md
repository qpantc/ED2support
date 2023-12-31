
## Plant functional types



> The user must specify which PFTs are allowed to occur in any given simulation.
> ED-2.2 has a list of default PFTs, with parameters described in [github wiki](https://github.com/EDmodel/ED2/wiki/Plant-functional-types).
> The user can modify the [parameters](https://github.com/EDmodel/ED2/wiki/PFT-parameters) of existing PFTs or define new PFTs through an extensible markup language ([XML](https://github.com/EDmodel/ED2/wiki/Model-parameters-and-xml-parameter-files)) file, which is read during the model initilization ([ED2IN](https://github.com/EDmodel/ED2/wiki/ED2IN-namelist)).


- config file
    config file for PFTs
    ```xml
    <?xml version="1.0"?>
        <!DOCTYPE config SYSTEM "ed.dtd">
            <config>
                <pft>
                    <num>3</num>
                    <is_tropical>1</is_tropical>
                    <Vm0>12</Vm0>
                    <wood_Kmax>0.06</wood_Kmax>
                </pft>
                <pft>
                    <num>4</num>
                    <is_tropical>1</is_tropical>
                </pft>
            </config>
    ```

## set the file path in ED2In

```R
NL%IEDCNFGF = '/data/gent/vo/000/gvo00074/pecan/output/other_runs/Wytham/growth_storage_resp/run/near_bare_ground/config.xml'
```

## PEcAn in ED2.2
Using PEcAn to help generate single or multiple ed2 running files.
PEcAn aims to synthesize plant trait data to estimate model parameters. It based on prior knowledge that is refined by species-level data using Bayesian meta-analysis. PEcAn allows to use a rigorous meta-analysis to inform the parameters of a mechanistic ecosystem model.

PEcAn的基本思路就是根据已有的数据,拟合一个最佳的先验的参数分布概率的模型,然后基于这个先验的概率分布,并在此基础上:
- 进行Meta分析: 一般是基于一种混合效应模型估计在当前条件下最优的参数的均值和标准差,并假设真实的参数为满足该条件的正态分布.
- 直接基于现有的先验分布取直进行模拟.
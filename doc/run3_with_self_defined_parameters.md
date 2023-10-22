
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
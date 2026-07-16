# Weekly Meeting on July 16, 2026

## Test Run of `CaSR` forcing field

**Something wrong?**

    Hit the time limit: slurmstepd: error: JOB 17724174 ON c698 CANCELLED AT 2026-07-16T01:22:18 DUE TO TIME LIMIT 

`-> info : Register new Context : nemo`

`-> info : CClientBuffer: allocated 2 x 13717553 bytes for server 0 with a maximum of 124 buffered events`

### What changed?

**nibi-example.yaml**

    - $HOME/MEOPAR/SS-run-sets/v202111/namelist.atmos_rivers.hrdps
    
    - /home/junqiqiu/MEOPAR/SS-run-sets/SalishSea/jqiu/v202111_CaSR_test_run/namelist.atmos_rivers.hrdps
    
    - $HOME/MEOPAR/SS-run-sets/v202111/namelist.light
    
    - /home/junqiqiu/MEOPAR/SS-run-sets/SalishSea/jqiu/v202111_CaSR_test_run/namelist.light

**namelist.atmos_rivers.hrdps**

    file name `hrdps` $\to$ `Casr`
    weight file `grid/weights-continental2.5-hrdps_202108_23feb23onward.nc` $\to$ `grid/weights-CaSR_202108.nc`

**Command used:**

    `command: pixi run -m $HOME/MEOPAR/SalishSeaCmd salishsea run /home/junqiqiu/MEOPAR/SS-run-sets/SalishSea/jqiu/v202111_CaSR_test_run/nibi-example.yaml  $SCRATCH/MEOPAR/junqi_test_results/junqi_test_results_CaSR`


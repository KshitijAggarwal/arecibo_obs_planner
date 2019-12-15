# Arecibo Observation Planner
Script to generate observation tables for Arecibo Observations

# Requirements
* Astropy
* Astroplan
* Pandas
* Numpy
* Tabulate

# Usage
Create a source list file (`source_list`) of the following format: 
    
    Source, ra, dec
    B0301+19, 03h04m33.115s, +19d32m51.4s
    
The observation tables for each observation for a proposal (say P3231) can be generated using: 

    python arecibo_obs_planner.py -l www.naic.edu/~gomathi/cgi-bin/pdbsearchct.cgi?pid=P1234 -s source_list
    
A table (eg: `2019-12-11_22:30:00.000.tab`) will be generated for each observation:

    Observations from 2019-12-11 22:30:00.000 to 2019-12-11 23:59:00.000
    Length of observations: 1.48hr
    +---------------+---------+----------+----------------+-------------------------+-------------------------+
    | target name   |      RA |      DEC |   uptime (min) | tstart_obs (AST)        | tend_obs (AST)          |
    |---------------+---------+----------+----------------+-------------------------+-------------------------|
    
Here `tstart_obs` and `tend_obs` will indicate the time window within which the source is observable (within the observation duration). Output table will only contain sources with uptime greater than `min_uptime`. These tables can be used to make the `.cmd` files for arecibo observations.

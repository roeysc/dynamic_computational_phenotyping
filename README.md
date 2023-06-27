# dynamic_computational_phenotyping
Code and data for: Schurr, Reznik, Hillman, Bhui, Gershman (submitted). Dynamic computational phenotyping of human cognition.

The code depositiry description:
- Stan: Stan code containing the computational models of each of the seven tasks. Each task has two model versions, both of them are hierarchical models - the 'independent' version (which assumes no temporal structure of the phenotype) and the 'dynamic' version, which embodies specific assumptions about the temporal evolution of the phenotype. This folder contains a short tutorial on how to run these scripts via python interface. Some of our Stan models were adopted from the hBayesDM Package (https://github.com/CCS-Lab/hBayesDM) created by Woo-Young, Haines, and Zhang 2017.
- Python and R (required libraries are indicated in the code): Mostly code for extracting posterior estimates from the posterior samples files (the 'chains'), submitting the Stan models, and generating required data structures. 
- Matlab (required libraries are indicated in the code): Analysis code at the level of behavior and parameter estimates.
- Tasks: All the scripts that were used for running the tasks. The tasks were coded by Hanna Hillman with the assistance of Jasmine Zhou.

Data: Data accompanying this masncuript will be made publicly availble in the future.

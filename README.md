# dynamic_computational_phenotyping
Code and data for: Schurr, Reznik et al. (submitted). Dynamic computational phenotyping of human cognition.

Code consists of several groups:
- Stan: A Stan file of the computational model of each task. Each model has two versions, both of them are hierarchical models - the 'independent' version (which assumes no temporal structure of the phenotype) and the 'dynamic' version, which embodies specific assumptions about the temporal evolution of the phenotype.
- Python and R (required libraries are indicated in the code): Mostly code for extracting posterior estimates from the posterior samples files (the 'chains') 
- Matlab (required libraries are indicated in the code): Analysis code at the level of behavior and parameter estimates.
- Tasks: All the scripts that were used for running the tasks. The tasks were coded by Hanna Hillman with the assistance of Jasmine Zhou

Data: Data will be added in the near future.

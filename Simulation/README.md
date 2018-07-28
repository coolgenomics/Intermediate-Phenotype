# Simulation

This directory contains the scripts for simulation. This simulates a cases and controls for a custom number of subgroups based on a logistic model that operates on intermediate phenotypes.

To create a simulation, first run 

```
python3 calc_simulation_param.py [weight_path] [genos_path] [result_path] [...gene_props]
```

This creates a randomized set of initial parameters for a trial run. Once the parameters are calculated, run

```
python3 simulation_logistic.py [weight_path] [param_path] [result_path] [heritability] [control_size] [...case_sizes]
```

This uses the previously generated parameters and simulates cases and controls. The length of case_sizes and gene_props must be the same.

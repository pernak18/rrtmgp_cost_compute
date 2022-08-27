# Preamble 

Originally designed as part of RRTMGP [_g_-point reduction](https://github.com/pernak18/g-point-reduction), the code in this repository calculates a user-defined cost function for a given formulation of RRTMGP, then compares the total cost and its components to other formulations. Different formulations exist because of, e.g., new RRTMGP library versions (code changes) or modified _k_-distributions or profile specifications (model input). Users define the cost function with:

1. component name (broadband flux or heating rate, band fluxes or heating rates, forcing)
2. levels for each component
3. weights for each component

# Comparing Costs

As of this writing, only the longwave (LW) has been tested. More work needs to be done for the shortwave (SW). Consequently, we will only go into detail with the LW.

Initially, the code was developed in the [`abs_val`](https://github.com/pernak18/g-point-reduction/tree/abs_val) branch of the _g_-point reduction repository. Eventually, the two code bases (reduction and cost computation) were decoupled.

## Inputs

### YAML

Most of the configuration for this code is given in YAML files, and templates for the [LW](https://github.com/pernak18/rrtmgp_cost_compute/blob/main/costCalc_wrapper_LW.yaml) and [SW](https://github.com/pernak18/rrtmgp_cost_compute/blob/main/costCalc_wrapper_SW.yaml) are in version control. The only argument into the script is the path to the YAML file, and this is optional (by default, it is the [LW template](https://github.com/pernak18/rrtmgp_cost_compute/blob/main/costCalc_wrapper_LW.yaml)).

Some documentation exists for each field in the YAML files, but we can elaborate:

| YAML Field | Description |
| :--- | :--- |
| `ref_path` | single string; path to LBLRTM flux netCDF file that follows the RRTMGP convention, an example of which is provided in the output of the `rrtmgp_garand_atmos` RRTMGP driver in [_g_-point reduction](https://github.com/pernak18/g-point-reduction/tree/master/garand_atmos) |
| `test_path` | single string; path to RRTMGP flux file (same convention as `ref_path`) with full _k_-distribution results (used for normalizing other model formulations) |
| `others` | list of paths (strings) to flux files for other RRTMGP formulations; can be any number of files; again, RRTMGP flux file convention |
| ``do_lw` | LW or SW switch (boolean) |
| `boundaries` | list of strings; YAML variable that is just reused (many times) in the `components` field -- we commonly only focus on the cost at the surface, tropopause, and top-of-atmosphere boundaries, and in this field we just explicitly say that and use the boundaries for many components; this field likely will only need to be changed if users only want to focus on, e.g., the surface for many components |
| `level_indices` | array indices used in the code to represent the surface, tropopause, and top-of-atmosphere boundaries |
| `components` | can be any number of components; each field name is the component number/index; each value is a list of 3 elements -- [component name](#naming), levels at which to compute the cost for the component, and the weight of the component in the total cost calculation |

## Component Naming Convention <a name="naming"></a>

Component names start with a variable name that can be found in the RRTMGP flux files -- `flux_net`, `flux_dn`, `flux_up`, or `heating_rate` for broadband, band the same 4 fields with `band_` prepended for by-band (SW options include `flux_dir_dn` and ``flux_dif_dn` as well). By themselves, these names are considered present-day Garand specifications. They all can have a suffix that represents a "record" or "experiment" name, of which there are 19 (following the [_g_-point combination naming convention](https://github.com/pernak18/g-point-reduction/wiki/LW-Forcing-Number-Convention#g-point-reduction-convention-) for all of the Garand "experiments" or "records" that are used by AER). The mapping of (unit-offset) record index (e.g., in the `flux_net` netCDF array) to record name is:

1. 'garand' (present day)
2. 'preind'
3. '2xch4' 
4. '2xco2'
5. 'pi_pd-no2' 
6. '4xch4'
7. '4xco2'
8. 'xs-all'
9. 'ccl4'
10. 'cfc11'
11. 'cfc12'
12. 'cfc22'
13. 'hfc143a'
14. 'hfc125'
15. 'hfc23' 
16. 'hfc32'
17. 'hfc134a'
18. 'cf4'
19. 'no2xs'

Example component names could be `flux_net_ccl4` or `band_heating_rate_4xch4`. Again, `flux_*` and `heating_rate` with no suffix is assumed to be present day ("garand").

Additionally, forcing components can be provided. Using an example of forcing due to double the CO<sub>2</sub> with respect to preindustrial specifications, we would want scenario 4 to be subtracted from scenario 2. This compuation happens in the [`flux_cost_compute.py`](https://github.com/pernak18/rrtmgp_cost_compute/blob/main/flux_cost_compute.py) module automatically (baselines are defined in the `costCalc` function), but the user must provide the "forcing" string in the component name. So in our example, one cost component name could be `flux_net_forcing_2xco2`, which can be read as "forcing in net flux due to double CO<sub>2</sub>."

## Flux and Heating Rate Calculation

To be fleshed out with standard output and possibly commentary...

```
./costCalc_wrapper.py
```

For the shortwave:

```
./costCalc_wrapper.py -i costCalc_wrapper_SW.yaml
```
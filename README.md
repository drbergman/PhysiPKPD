# PhysiPKPD
## Getting Started...
### ...by downloading PhysiPKPD add adding to a working PhysiCell directory
1. Download the repository and unzip the file.
2. Move the folder `PhysiPKPD/addons/PhysiPKPD` into `PhysiCell/addons/`
3. Move the folder `PhysiPKPD/sample_projects_phsyipkpd` into `PhysiCell`
4. Open `PhysiCell/sample_projects/Makefile-default` (the one that `make reset` will will put in the main PhysiCell directory)
5. Add the text from `Makefile-PhysiPKPD_Addendum` to `PhysiCell/sample_projects/Makefile-default` (anywhere should work, perhaps best around line 195 at the end of the other sample projects)
6. Replace `PhysiCell/BioFVM/BioFVM_microenvironment.cpp` with `PhysiPKPD/BioFVM_microenvironment_robust_dcs.cpp` (and similarly with the `.h` file). Make sure the new files are both named `BioFVM_microenvironment` with the proper extension.
7. Delete `PhysiCell/BioFVM_microenvironment.o` so the next `make` call will compile the new `BioFVM_microenvironment` code.

### ...by cloning the repository
1. Fork this repository to your own GitHub account.
2. Clone the resulting forked repository onto your machine.
3. Copy all the PhysiCell files in your PhysiCell directory **except addons**
4. Copy the subfolders in `PhysiCell/addons` into your cloned directory's `addons` folder
5. Open `PhysiCell/sample_projects/Makefile-default` (the one that `make reset` will will put in the main PhysiCell directory)
6. Add the text from `Makefile-PhysiPKPD_Addendum` to `PhysiCell/sample_projects/Makefile-default` (anywhere should work, perhaps best around line 195 at the end of the other sample projects)
7. Replace `PhysiCell/BioFVM/BioFVM_microenvironment.cpp` with `PhysiPKPD/BioFVM_microenvironment_robust_dcs.cpp` (and similarly with the `.h` file). Make sure the new files are both named `BioFVM_microenvironment` with the proper extension.
8. Delete `PhysiCell/BioFVM_microenvironment.o` so the next `make` call will compile the new `BioFVM_microenvironment` code.

Congratulations! You're ready to try out PhysiPKPD!

## Running the samples
There are 5 sample projects currently distributed with PhysiPKPD.
There is one for each supported Mechanism of Action (MOA) and one combination treatment.
To run one of these samples, do the following:

1. `make reset` to make sure you have the newly edited Makefile in your top directory
2. Make your preferred project:
    * `make moa_proliferation` 
    * `make moa_apoptosis` 
    * `make moa_necrosis` 
    * `make moa_motility` 
    * `make combo`
3. Compile your project: `make`
4. Run your project: `./project ./config/mymodel.xml`
5. Look at the snapshots in `output/` and the living cell counts in `output/cell_counts.csv`

## Reconfiguring, editing, and re-running
**Note:** *While the below functionality is present, it is discouraged because it is likely to inadvertently affect the work of others.
Instead, it is recommended to instead save any changes to these files in a non-tracked directory and manually copy them into their proper places after `make`-ing the sample project.*

Instead of editing the configuration file copied into `PhysiCell/config/mymodel.xml`, you can choose to edit the original in `PhysiCell/sample_projects_physipkpd/[project_name]/config/mymodel.xml` to save the changes for future runs.
The command `make rc` will reconfigure from the original `mymodel.xml` to facilitate editing in the latter fashion.

Similarly, you can edit the custom modules in `PhysiCell/sample_projects_physipkpd/[project_name]/custom_modules/` to save changes for future runs.
After making these changes, you can run `make redo` and this will automatically move those changes to their proper places and recompile the project.

## Setting parameters
PhysiPKPD parameters are largely concentrated in two areas in `mymodel.xml`: PK parameters are at the bottom in `user_parameters` and PD parameters are in `cell_definitions` in the `custom_data` for each cell type.
PhysiPKPD can add PK and/or PD dynamics to any substrates in a PhysiCell simulation.
PK dynamics must be set for each PK substrate and PD dynamics determined for each cell type affected by a particular substrate.

In what follows, `S` stands for the name of a substrate.

### PK and PD substrates
You can specify which substrates you want to include PK dynamics.
You can do the same for PD dynamics.
These lists need not have *any* relationships.
These are the names PhysiPKPD will use to look through all the other `user_parameters` and `custom_data` to implement the PKPD dynamics.
Add the following to `user_parameters`:

| Parameter | Description |
|---|---|
| PKPD_pk_substrate_names | Comma-separated list of substrates also following PK dynamics |
| PKPD_pd_substrate_names | Comma-separated list of substrates also producing PD dynamics in some cells |

PhysiPKPD will issue a warning if it finds a name in either of these lists that is not a substrate.
Watch the initial output closely!
Example:
```
<PKPD_pk_substrate_names type="string">myDrug,myDrug_no_PD</PKPD_pk_substrate_names>
<PKPD_pd_substrate_names type="string">myDrug,myDrug_no_PK</PKPD_pd_substrate_names>
```
In the above example, `myDrug` and `myDrug_no_PD` will follow PK dynamics.
In addition, `myDrug` and `myDrug_no_PK` will lead to PD effects for any cell type with `S_moa_is_X` set to `1.0` (see [below](#pd_pars)).

**Note:** Any spaces in these lists will cause PhysiPKPD to look for substrate names with those spaces.

### PK parameters <a name="pk_pars"></a>
PK dynamics in PhysiPKPD currently follows a 2-compartment model[^oldpk]:

[^oldpk]: For those using the simplified 2-compartment model that PhysiPKPD used to use, see the [`S_flux_across_capillaries`](#old_flux_par) entry in the table.

$$
\begin{aligned}
C' & = \frac{k_{21}}{R}P - k_{12}C - \lambda C \\
P' & = k_{12}RC - k_{21}P
\end{aligned}
$$

A 1-compartment model can be modeled by setting $k_{12}=0$[^divby0].

[^divby0]: The implemented analytical solution, in this one case, requires $k_{21}\neq\lambda$ to avoid a divide-by-zero error. Since the periphery is effectively excluded in this case, the value of $k_{21}$ is irrelevant, so we just add one and move on.

We are also working on including a way for users to implement even more complex PK dynamics.

For each drug, you can set the following parameters in `user_parameters`.
For example, a substrate called `myDrug` with 10 doses administered would have `myDrug_max_number_doses` set to `10`.

| Parameter | Type | Description |
| ---  | --- | --- |
| `S_number_loading_doses` | `int` | Number of loading doses to give before switching to regular doses |
| `S_max_number_doses` | `int` | Total number of doses to give including loading doses |
| `S_dose_interval` | `double` | Time between successive doses, loading or regular (in minutes) |
| `S_set_first_dose_time` | `bool` | Boolean determining if the first dose time is fixed or if a confluence condition will be used to determine the first dose time |
| `S_first_dose_time` | `double` | Time of first dose if given at fixed time (in minutes) |
| `S_confluence_condition` | `double` | Proportion of microenvironment filled with cells at which to give first dose; confluence calculated by sum of cross-sectional area of all cells divided by area of microenvironment
| `d1_color_ec50` | `double` | If `damage_coloring` is used for plotting, this sets the damage from drug 1 that causes half the maximum redshift in the cell cytoplasm |
| `d1_color_hp` | `double` | If `damage_coloring` is used for plotting, this is the Hill coefficient used to calculate the amount of redshift in the cytoplasm |
| `d2_color_ec50` | `double` | If `damage_coloring` is used for plotting, this sets the damage from drug 2 that causes half the maximum blueshift in the cell nucleus |
| `d2_color_hp` | `double` | If `damage_coloring` is used for plotting, this is the Hill coefficient used to calculate the amount of blueshift in the nucleus |
| `S_central_increase_on_loading_dose` | `double` | Increase in concentration in central compartment after a loading dose |
| `S_central_increase_on_dose` | `double` | Increase in concentration in central compartment after a regular dose |
| `S_central_elimination_rate` $(\lambda)$ | `double` | Linear elimination rate in central compartment (in mintues<sup>-1</sup>) |
| `S_flux_across_capillaries`<a name="old_flux_par"></a> | `double` | **While this is still allowed, consider using the following two parameters to quantify intercompartmental clearance rates.**[^1] Rate of change in concentration in central compartment due to distribution and redistribution (in minutes<sup>-1</sup>) |
| `S_central_to_periphery_clearance_rate` $(k_{12})$ | `double` | Rate of change in concentration in central compartment due to distribution (in minutes<sup>-1</sup>) |
| `S_periphery_to_central_clearance_rate` $(k_{21})$ | `double` | Rate of change in concentration in periphery compartment due to redistribution (in minutes<sup>-1</sup>) |
| `S_biot_number` | `double` | Ratio of drug concentration on boundary of microenvironment (Dirichlet condition) and concentration in systemic circulation |
|`S_central_to_periphery_volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `double` | Ratio of central compartment to periphery compartment |
|`central_to_periphery_volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `double` | Ratio of central compartment to periphery compartment *for any substrates without a specific volume ratio as above* |

[^1]: To use these new parameters, you will want to set $k_{12}$ as your original flux rate and $k_{21}$ as `S_flux_across_capillaries * S_central_to_periphery_volume_ratio`.

You can also set the following parameters in `microenvironment_setup` for each drug:
| Parameter | Description |
| ---| --- |
| `diffusion_coefficient` | Diffusion rate in the microenvironment |
| `decay_rate` | Rate of decay in the microenvironment |

### PD parameters <a name="pd_pars"></a>
For each cell type, all of the PD parameters are in `custom_data`.
In the table below, `X` can stand for any one of `prolif`, `apop`, `necrosis`, or `motility`.
For example, a substrate called `myDrug` causing a proliferation effect would have `myDrug_moa_is_prolif` set to `1.0`.
|Parameter|Description|
|---|---|
| `S_moa_is_X` | Used as boolean to determine which effects to apply to this cell type based on the damage from drug 1; values > 0.5 will apply the effect |
| `S_X_saturation_rate` | Rate of `X` as damage from drug 1 approaches infinity |
| `S_X_EC50` | Damage from drug 1 at which the rate of `X` is halfway between the base and saturation rates (in damage) |
| `S_X_hill_power` | Hill coefficient for calculating the effect of drug 1 on the rate of `X` |
| `S_damage` | Not a parameter; data that tracks the current damage to the cell |
| `S_repair_rate` | **While this is still allowed, consider using the following two parameters to quantify repair rates instead.** Zero-order elimination rate of damage from drug 1 (in damage per minute) |
| `S_repair_rate_constant` | Zero-order elimination rate of damage from drug 1 (in damage per minute) |
| `S_repair_rate_linear` | First-order elimination rate of damage from drug 1 (in minutes<sup>-1</sup>) |
| `S_metabolism_rate` | Rate of elimination of drug 1 from inside a cell (in minutes<sup>-1</sup>) |

### Miscellaneous parameters
The following are user parameters that provide some control over how the dynamics are solved.
Currently, the only PK and PD models implemented in PhysiPKPD are readily solved using analytic techniques.
By informal observation, the analytic methods are not slower than numerical methods.
If anything the analtyic methods are faster.
Therefore, all simulations use analytic solutions.

To maximize the efficiency of these analytic solutions, many terms are pre-computed ahead of time.
**However**, if the parameters governing your PD dynamics vary--e.g., due to heterogeneity within the affected cell type--then you will not want to do pre-computations.
These pre-computations are an all-or-nothing for every (substrate, cell type) pairing.
You can include the following as Booleans in `user_parameters` to control this.
In the following table, `C` stands for a cell type name.
For example, if a substrate `myDrug` affects cell type `tumor` and you want to pre-compute the PD quantities, set `myDrug_precompute_pd_for_tumor` to `True`.

|Parameter|Type|Description|
|---|---|---|
| `PKPD_precompute_all_pd_quantities` | 'bool' | Boolean to determine whether to pre-compute values used to analytically solve *all* PD dynamics; **turn this off if PD parameters can vary within a cell type OR if your `mechanics_dt` is not a multiple of your `diffusion_dt`**; defaults to `False` if not set |
| `S_precompute_pd_for_C` | 'bool' | **Only affects simulations if `PKPD_precompute_all_pd_quantities==False`** Boolean to determine whether to pre-compute values used to analytically solve PD dynamics for a single (substrate, cell type) pairing; **turn this off if PD parameters can vary within this cell type OR if your `mechanics_dt` is not a multiple of your `diffusion_dt`**; defaults to `False` if not set |

## Making your own project using PhysiPKPD
If you wish to make your own project that uses PhysiPKPD (and not just one of the pre-built sample projects), this is how you can proceed.
1. Make the PKPD template project: `make template_pkpd`
2. Edit the configuration file to set the Dirichlet conditions, [PK Parameters](#pk_pars), and [PD Parameters](#pd_pars) for the two PKPD drugs and the default cell type `cell`.
3. Add additional substrates as normal (using the Model Builder for this is untested)
4. Add additional cell types as normal (using the Model Builder for this is untested).
Copy the `custom_data` block for `cell` into any newly created cell types, setting these as desired.
5. By default, each cell type is assigned the same `update_phenotype` function, which is `cell_phenotype` found in the `custom.cpp` file.
Add new phenotype functions as desired for each cell type.
6. For each phenotype function, make sure to uncomment the line resetting the mechanism of action to its base value.
7. If the mechanism of action is motility, then uncomment the line setting the `update_migration_bias` or add that line for each cell type that undergoes a motility effect.

**Note:** You must manually put any any chemotactic signals in the `update_migration_bias` function if you use a motility effect.
See `PhysiCell/core/PhysiCell_standard_models.cpp` for the `chemotaxis_function`, `advanced_chemotaxis_function`, and `advanced_chemotaxis_function_normalized`.
Hopefully, this will not be necessary in the future.

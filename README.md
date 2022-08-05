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

## Varying parameters
PhysiPKPD parameters are largely concentrated in two areas in `mymodel.xml`: PK parameters are at the bottom in `user_parameters` and PD parameters are in `cell_definitions` in the `custom_data` for each cell type.
PhysiPKPD comes hardcoded with two drugs and neither can be excluded. Of course, you can set them so that there are no doses or that doses result in no increase to the drug concentration.
PK dynamics must be set for each drug and PD dynamics determined for each cell type for each drug.

### PK parameters <a name="pk_pars"></a>
For each drug, you can set the following parameters in `user_parameters`:

| Parameter | Description |
| ---  | --- |
| `PKPD_D1_number_loading_doses` | Number of loading doses to give before switching to regular doses |
| `PKPD_D1_max_number_doses` | Total number of doses to give including loading doses |
| `PKPD_D1_dose_interval` | Time between successive doses, loading or regular (in minutes) |
| `PKPD_D1_set_first_dose_time` | Boolean determining if the first dose time is fixed or if a confluence condition will be used to determine the first dose time |
| `PKPD_D1_first_dose_time` | Time of first dose if given at fixed time (in minutes) |
| `PKPD_D1_confluence_condition` | Proportion of microenvironment filled with cells at which to give first dose; confluence calculated by sum of cross-sectional area of all cells divided by area of microenvironment
| `d1_color_ec50` | If `damage_coloring` is used for plotting, this sets the damage from drug 1 that causes half the maximum redshift in the cell cytoplasm |
| `d1_color_hp` | If `damage_coloring` is used for plotting, this is the Hill coefficient used to calculate the amount of redshift in the cytoplasm |
| `d2_color_ec50` | If `damage_coloring` is used for plotting, this sets the damage from drug 2 that causes half the maximum blueshift in the cell nucleus |
| `d2_color_hp` | If `damage_coloring` is used for plotting, this is the Hill coefficient used to calculate the amount of blueshift in the nucleus |
| `PKPD_D1_central_increase_on_loading_dose` | Increase in concentration in central compartment after a loading dose |
| `PKPD_D1_central_increase_on_dose` | Increase in concentration in central compartment after a regular dose |
| `PKPD_D1_central_elimination_rate` | Linear elimination rate in central compartment (in mintues<sup>-1</sup>) |
| `PKPD_D1_flux_across_capillaries` | **While this is still allowed, consider using the following two parameters to quantify intercompartmental clearance rates.**[^1] Rate of change in concentration in central compartment due to distribution and redistribution (in minutes<sup>-1</sup>) |
| `PKPD_D1_central_to_periphery_clearance_rate` $(k_{12})$ | Rate of change in concentration in central compartment due to distribution (in minutes<sup>-1</sup>) |
| `PKPD_D1_periphery_to_central_clearance_rate` $(k_{21})$ | Rate of change in concentration in periphery compartment due to redistribution (in minutes<sup>-1</sup>) |
| `PKPD_D1_biot_number` | Ratio of drug concentration on boundary of microenvironment (Dirichlet condition) and concentration in systemic circulation |
|`central_to_periphery_volume_ratio` $(V_1/V_2 \text{ or } V_C/V_P)$ | Ratio of central compartment to periphery compartment |

[^1]: To use these new parameters, you will want to set $k_{12}$ as your original flux rate and $k_{21}$ as `PKPD_D1_flux_across_capillaries * central_to_periphery_volume_ratio`.

You can also set the following parameters in `microenvironment_setup` for each drug:
| Parameter | Description |
| ---| --- |
| `diffusion_coefficient` | Diffusion rate in the microenvironment |
| `decay_rate` | Rate of decay in the microenvironment |

### PD parameters <a name="pd_pars"></a>
For each cell type, all of the PD parameters are in `custom_data` for each cell type.
In the table below, `X` can stand for any one of `prolif`, `apop`, `necrosis`, or `motility`.
|Parameter|Description|
|---|---|
| `PKPD_D1_moa_is_X` | Used as boolean to determine which effects to apply to this cell type based on the damage from drug 1; values > 0.5 will apply the effect |
| `PKPD_D1_X_saturation_rate` | Rate of `X` as damage from drug 1 approaches infinity |
| `PKPD_D1_X_EC50` | Damage from drug 1 at which the rate of `X` is halfway between the base and saturation rates (in damage) |
| `PKPD_D1_X_hill_power` | Hill coefficient for calculating the effect of drug 1 on the rate of `X` |
| `PKPD_D1_damage` | Not a parameter; data that tracks the current damage to the cell |
| `PKPD_D1_repair_rate` | **While this is still allowed, consider using the following two parameters to quantify repair rates instead.** Zero-order elimination rate of damage from drug 1 (in damage per minute) |
| `PKPD_D1_repair_rate_constant` | Zero-order elimination rate of damage from drug 1 (in damage per minute) |
| `PKPD_D1_repair_rate_linear` | First-order elimination rate of damage from drug 1 (in minutes<sup>-1</sup>) |
| `PKPD_D1_metabolism_rate` | Rate of elimination of drug 1 from inside a cell (in minutes<sup>-1</sup>) |

### Miscellaneous parameters
The following are user parameters that provide some control over how the dynamics are solved.
Currently, the only PK and PD models implemented in PhysiPKPD are readily solved using analytic techniques.
By informal observation, the analytic solutions are not slower than the numerical solutions.
If anything the analtyic solutions are faster.
But you do you!
|Parameter|Description|
|---|---|
| `PKPD_use_analytic_pk_solutions` | Boolean to determine whether to use analytic (matrix exponential) or numerical (explicit Euler) to solve PK dynamics; defaults to `True` |
| `PKPD_use_analytic_pd_solutions` | Boolean to determine whether to use analytic (matrix exponential) or numerical (explicit Euler) to solve PD dynamics; defaults to `False` |
| `PKPD_precompute_pd_quantities` | Boolean to determine whether to pre-compute values used to analytically solve PD dynamics; **turn this off if PD parameters can vary within a cell type OR if your `mechanics_dt` is not a multiple of your `diffusion_dt`; defaults to `False` |

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

**Note:** You must manually put any any chemotactic signals here.
See `PhysiCell/core/PhysiCell_standard_models.cpp` for the `chemotaxis_function`, `advanced_chemotaxis_function`, and `advanced_chemotaxis_function_normalized`.
Hopefully, this will not be necessary in the future.

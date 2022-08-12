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
PhysiPKPD parameters are found in two areas in `mymodel.xml`: PK parameters are at the bottom in `user_parameters` and PD parameters are in `cell_definitions` in the `custom_data` for each cell type.
PhysiPKPD can add PK and/or PD dynamics to any substrates in a PhysiCell simulation.
PK dynamics must be set for each PK substrate and PD dynamics determined for each cell type affected by a particular substrate.

In what follows, `S` stands for the name of a substrate.

### PK and PD substrates
You can specify which substrates you want to include PK dynamics.
You can do the same for PD dynamics.
These lists need not have *any* relationships.
The parameters in the following table are those PhysiPKPD will use to look through all the other `user_parameters` and `custom_data` to implement the PKPD dynamics.
Add them to `user_parameters`.
Omitting `PKPD_pk_substrate_names` will result in no PK dynamics and similarly for the PD dynamics.

| Parameter | Description |
|:--|---|
| `PKPD_pk_substrate_names` | Comma-separated list of substrates also following PK dynamics |
| `PKPD_pd_substrate_names` | Comma-separated list of substrates also producing PD dynamics in some cells |
<p align="center">
    <b>Table:</b> Identifying PK and PD Substrates
</p>

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

[^divby0]: The implemented analytical solution, in the case $k_{12}=0$, requires $k_{21}\neq\lambda$ to avoid a divide-by-zero error.
Since the periphery is effectively excluded in this case, the value of $k_{21}$ is irrelevant, so PhysiPKPD updates `k_{21}+=1` and carries on.

We are also working on including a way for users to implement even more complex PK dynamics.

For each substrate, you can set the following parameters in `user_parameters`.
For example, a substrate called `myDrug` with 10 doses administered would have the following in `user_parameters`: 

```
<myDrug_max_number_doses type="int">10</myDrug_max_number_doses>
```

| Parameter | Type | Description | If Missing |
| :-- | --- | :-: | --- |
| `S_number_loading_doses` | `int` | Number of loading doses to give before switching to regular doses | Set to `0` |
| `S_max_number_doses` | `int` | Total number of doses to give including loading doses | Set to `0` |
| `S_dose_interval` | `double` | Time between successive doses, loading or regular (in minutes) | If `S_max_number_doses`>0, throws an error |
| `S_set_first_dose_time` | `bool` | Boolean determining if the first dose time is fixed or if a confluence condition will be used to determine the first dose time | Set to `False` |
| `S_first_dose_time` | `double` | Time of first dose if given at fixed time (in minutes) | Set to `current_time` |
| `S_confluence_condition` | `double` | Proportion of microenvironment filled with cells at which to give first dose; confluence calculated by sum of cross-sectional area of all cells divided by area of microenvironment | Never check confluence condition |
| `S_central_increase_on_loading_dose` | `double` | Increase in concentration in central compartment after a loading dose | If `S_number_loading_doses`>0, throws an error |
| `S_central_increase_on_dose` | `double` | Increase in concentration in central compartment after a regular dose | If `S_max_number_doses`>`S_number_loading_doses`, throws an error |
| `S_central_elimination_rate` $(\lambda)$ | `double` | Linear elimination rate in central compartment (in mintues<sup>-1</sup>) | Set to `0` |
| `S_central_to_periphery_volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `double` | Ratio of central compartment to periphery compartment | Set to `central_to_periphery_volume_ratio` |
| `central_to_periphery_volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `double` | Ratio of central compartment to periphery compartment *for any substrates without a specific volume ratio as above* | Set to `1` |
| `S_central_to_periphery_clearance_rate` $(k_{12})$ | `double` | Rate of change in concentration in central compartment due to distribution (in minutes<sup>-1</sup>) | Set to `S_flux_across_capillaries`, if present. Otherwise, set to `0` |
| `S_periphery_to_central_clearance_rate` $(k_{21})$ | `double` | Rate of change in concentration in periphery compartment due to redistribution (in minutes<sup>-1</sup>) | Set to `S_flux_across_capillaries` $\times$ `S_central_to_periphery_volume_ratio`, if present. Otherwise, set to `0` |
| `S_flux_across_capillaries`<a name="old_flux_par"></a> | `double` | **While this is still allowed, consider using the above two parameters to quantify intercompartmental clearance rates.**[^1] Rate of change in concentration in central compartment due to distribution and redistribution (in minutes<sup>-1</sup>) | See above |
| `S_biot_number` | `double` | Ratio of substrate concentration on boundary of microenvironment (Dirichlet condition) and concentration in systemic circulation | Set to `1` |
<p align="center">
    <b>Table:</b> PK Parameters
</p>

[^1]: To use these new parameters, you will want to set $k_{12}$ as your original flux rate and $k_{21}$ as `S_flux_across_capillaries * S_central_to_periphery_volume_ratio`.

You can also set the following parameters in `microenvironment_setup` for each substrate:
| Parameter | Description |
| :--| --- |
| `diffusion_coefficient` | Diffusion rate in the microenvironment |
| `decay_rate` | Rate of decay in the microenvironment |
<p align="center">
    <b>Table:</b> PK Parameters in PhysiCell
</p>

### PD parameters <a name="pd_pars"></a>
For each cell type, all of the PD parameters are in `custom_data`.


#### Damage Accumulation Parameters <a name="dam_pars"></a>
Each cell affected by a substrate `S` accumulates damage, $D$, based on the *amount*, $A$ of the substrate internalized.
This accumulation of damage obeys the following differential equation:

$$
\begin{aligned}
A' & = -mA \\
D' & = A - r_1D - r_0
\end{aligned}
$$

Note that since damage is an abstract quantity, we do not include a rate parameter as a coefficient for $A$ in the equation for $D'$.
Each of these parameters **must** be set[^old_repair].
If they are not set, PhysiPKPD will throw an error.
By default, PhysiPKPD uses the `mechanics_dt` set in the configuration file.
You can change this by adding a user parameter of the format `S_dt_C`.
For example, if you wish to use a timestep of `0.01` for updating the damage of `myDrug` on `tumor`, then add the following to `user_parameters`:

```
<myDrug_dt_tumor type="double">0.01</myDrug_dt_tumor>
```

[^old_repair]: Old versions of PhysiPKPD only had a constant repair rate. For backwards compatibility, `S_repair_rate` can be set instead.
If the repair rates in the [table](#dam_pars_table) are not set, then PhysiPKPD will use this value to set a constant repair rate.

|Parameter|Description|
|:--|---|
| `S_damage` $(D)$ | Not a parameter; data that tracks the current damage to the cell |
| `S_repair_rate_constant` $(r_0)$ | Zero-order elimination rate of damage from `S` (in damage per minute) |
| `S_repair_rate_linear` $(r_1)$ | First-order elimination rate of damage from `S` (in minutes<sup>-1</sup>) |
| `S_metabolism_rate` $(m)$ | Rate of elimination of `S` from inside a cell (in minutes<sup>-1</sup>) |

<p align="center">
    <b>Table:</b> Damage Accumulation Parameters <a name="dam_pars_table"></a>
</p>

#### Cell Effect Parameters
The accumulated damage then affects the cell through one or more mechanisms of action (MOA).
These MOAs are specified by setting `S_moa_is_X` to `1.0` in `custom_data`.
Currently, PhysiPKPD supports the following MOAs: `prolif`, `apop`, `necrosis`, or `motility`.
Replace `X` with the corresponding one of these you wish for your MOA.
For example, a substrate called `myDrug` causing a proliferation effect would have the following in `custom_data` for the affected cell:

```
<myDrug_moa_is_prolif>1.0</myDrug_moa_is_prolif>
```

MOAs that are not present in `custom_data` are assumed off.
*Beware of typos!*

For each MOA, the damage is input into a Hill-type function so that the effect eventually saturates.
These parameters also go in `custom_data`, but if any are omitted for an active MOA, PhysiCell will relentlessly spam the standard output.
The effect is applied to whatever the *current* relevant rate for the cell is.
This way, you can easily add these PD effects on top of other effects already present in your model.
Just take care that before calling `pd_function` within `update_phenotype` that the affected rate is reset to the base value lest PhysiPKPD compound these effects each time.
See the sample projects targeting `proliferation`, `apop`, and `necrosis`.
If you target `motility`, do the same but within `update_migration_bias`.

|Parameter|Description|
|:--|:--|
| `S_moa_is_X` | Used as boolean to determine which effects to apply to this cell type based on the damage from `S`; values > 0.5 will apply the effect |
| `S_X_saturation_rate` | Rate of `X` as damage from `S` approaches infinity |
| `S_X_EC50` | Damage from `S` at which the rate of `X` is halfway between the base and saturation rates (in damage) |
| `S_X_hill_power` | Hill coefficient for calculating the effect of `S` on the rate of `X` |
<p align="center">
    <b>Table:</b> Cell Effect Parameters
</p>

In the event that multiple substrates have the same MOA for a given cell type, PhysiPKPD combines their effects multiplicatively.
That is, PhysiPKPD internally uses the `S_X_saturation_rate` and the associated base rate for the given `Cell_Definition` to compute a saturation factor.
The value of the Hill-type function then determines the magnitude of the applied factor for that substrate.
All of these factors are multiplied and finally multiplied to the *current* rate for the given cell.
In the case of necrosis, however, it is our experience that most models assume a base necrosis rate of `0`, which would result in divide-by-zero errors in this formulation.
Instead, we choose to add up the effects rather than computing and then multiplying the factors.
This can be illustrated with an example:

##### Example: MOA targeted by multiple substrates
Cell `C` of type `CD` is targeted by `S1` and `S2`, both with an `apop` MOA on `CD`.
Currently, `C` has damages `S1_damage` and `S2_damage` from `S1` and `S2`, respectively.
Using the distinct EC50s and Hill powers for `S1` and `S2` on `C`, these Hill-type functions output $h_1 = 0.25$ and $h_2=0.5$ for `S1` and `S2`, respectively.
These values can be interpreted as `S1` causing a quarter of its maximal effect on `C` and `S2` causing half its maximal effect on `C`.

Now, suppose that the saturation effects on `C` are `1E-3` $\text{min}^{-1}$ and `1E-2` $\text{min}^{-1}$ for `S1` and `S2`, respectively.
The base apoptosis rate for `CD` is `1E-5` $\text{min}^{-1}$.
Then, `S1` and `S2` have saturation factors of $f_{\text{sat},1} = 10^{-3} / 10^{-5} = 100$ and $f_{\text{sat},2} = 10^{-2} / 10^{-5} = 1000$, respectively.
These can be interpreted as the substrates causing up to a 100-fold and 1000-fold, respectively, increase in the apoptosis rate of `C`.

Thus, the factor change from `S1` is $f_{\text{sat},1} \cdot h_1 = 25$; from `S2` it is $f_{\text{sat},2} \cdot h_2 = 500$.
The total factor change applied to the apoptosis rate for `C` is then $25\cdot 500 = 12,500$.
Finally, this is multplied to the *current* apoptosis rate for `C`.
If other aspects of the model cause `C` to deviate from the base apoptosis rate of `CD`, then PhysiPKPD honors these.
If, for example the current apoptosis rate of `C` is not `1E-5` but instead `1E-6` perhaps because it is in a quiescent niche, then PhysiPKPD will update its apoptosis rate to $12,500 \cdot 10^{-6} = 0.0125\ \text{min}^{-1}$.

### Miscellaneous parameters
The following are user parameters that provide some control over how the dynamics are solved.
Currently, the only PK and PD models implemented in PhysiPKPD are readily solved using analytic techniques.
By informal observation, the analytic methods are not slower than numerical methods.
If anything the analtyic methods are faster.
Therefore, all simulations use analytic solutions.

To maximize the efficiency of these analytic solutions, many terms are pre-computed ahead of time.
**However**, if the parameters governing your PD dynamics vary--e.g., due to heterogeneity within the affected cell type--then you will not want to do pre-computations.
These pre-computations are an all-or-nothing for every (substrate, cell type) pairing.
In other words, you cannot specify that some of the pre-computations can be done but not others.
You can include the following as Booleans in `user_parameters` to control this.
In the following table, `C` stands for a cell type name.
For example, if a substrate `myDrug` affects cell type `tumor` and you want to pre-compute the PD quantities, include the following in `user_parameters`:

```
<myDrug_precompute_pd_for_tumor type="bool">True</myDrug_precompute_pd_for_tumor>
```

|Parameter|Type|Description| If Missing |
|---|---|---|:--|
| `PKPD_precompute_all_pd_quantities` | 'bool' | Boolean to determine whether to pre-compute values used to analytically solve *all* PD dynamics; **turn this off if PD parameters can vary within a cell type OR if your `mechanics_dt` is not a multiple of your `diffusion_dt`** | Set to `False` |
| `S_precompute_pd_for_C` | 'bool' | **Only affects simulations if `PKPD_precompute_all_pd_quantities==False`** Boolean to determine whether to pre-compute values used to analytically solve PD dynamics for a single (substrate, cell type) pairing; **turn this off if PD parameters can vary within this cell type OR if your `mechanics_dt` is not a multiple of your `diffusion_dt`** | Set to `False` |
<p align="center">
    <b>Table:</b> Miscellaneous Parameters
</p>

## Making your own project using PhysiPKPD
If you wish to make your own project that uses PhysiPKPD (and not just one of the pre-built sample projects), this is how you can proceed.
1. Make the PKPD template project: `make template_pkpd`
2. Edit the configuration file to set the Dirichlet conditions, [PK Parameters](#pk_pars), and [PD Parameters](#pd_pars) for the two PKPD substrates and the default cell type `cell`.
3. Add additional substrates as normal (using the Model Builder for this is untested)
4. Add additional cell types as normal (using the Model Builder for this is untested).
5. By default, each cell type is assigned the same `update_phenotype` function, which is `cell_phenotype` found in the `custom.cpp` file.
Add new phenotype functions as desired for each cell type.
6. For each phenotype function, make sure to uncomment the line resetting the mechanism of action to its base value.
7. If the mechanism of action is motility, then uncomment the line setting the `update_migration_bias` or add that line for each cell type that undergoes a motility effect[^mot].

**Note:** The `custom.cpp` file that is loaded with the `template_pkpd` project has two substrates hardcoded to use for the motility MOA.
If you wish to add additional substrates that have a motility MOA (or change the names of the default substrates), you will need to change the `motility_rule` to reflect this for these substrates to affect cell migration speed.

[^mot]: You must manually put any any chemotactic signals in the `update_migration_bias` function if you use a motility effect.
See `PhysiCell/core/PhysiCell_standard_models.cpp` for the `chemotaxis_function`, `advanced_chemotaxis_function`, and `advanced_chemotaxis_function_normalized`.
Hopefully, this will not be necessary in the future.

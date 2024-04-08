# PhysiPKPD
## Getting started...
### ...with a new project
1. Download (or copy contents of) `get_physipkpd.py` to a local file.
2. Run 
```
python get_physipkpd.py --dir path/to/project/dir --studio
```
   * The `--studio` flag downloads a PhysiPKPD-enabled version of PhysiCell-Studio to help build the model
3. Follow output for next steps

### ...adding PhysiPKPD to an existing project
1. Download (or copy contents of) `add_physipkpd.py` to a local file.
2. Run
```
python add_physipkpd.py path/to/project/dir --studio
```
   * The `--studio` flag downloads a PhysiPKPD-enabled version of PhysiCell-Studio to help build the model
3. Follow output for next steps

### ...by downloading PhysiPKPD and adding to a working PhysiCell directory <a name="DL"></a>
1. Download the repository and unzip the file.
2. Move the folder `PhysiPKPD/addons/PhysiPKPD` into `PhysiCell/addons/`
3. Move the folder `PhysiPKPD/sample_projects_phsyipkpd` into `PhysiCell`
4. Open `PhysiCell/sample_projects/Makefile-default` (the one that `make reset` will put it in the main PhysiCell directory)
5. Add the text from `PhysiCell/addons/PhysiPKPD/Makefile-PhysiPKPD_Addendum` to `PhysiCell/sample_projects/Makefile-default` (anywhere should work, perhaps best around line 195 at the end of the other sample projects)

### ...by cloning the repository
1. Fork this repository to your own GitHub account.
2. Clone the resulting forked repository onto your machine.
3. Copy all the PhysiCell files in your PhysiCell directory **except addons**
4. Copy the subfolders in `PhysiCell/addons` into your cloned directory's `addons` folder
5. Continue from #4 [above](#DL).

Congratulations! You're ready to try out PhysiPKPD!

## Running the sample projects
There are 5 sample projects currently distributed with PhysiPKPD.
They highlight some of the possible behaviors achievable with PhysiPKPD.
To run one of these samples, do the following:

1. `make reset` to make sure you have the newly edited Makefile in your top directory
2. Make your preferred sample project:
    * `make pkpd-proliferation-sample` 
    * `make pkpd-apoptosis-sample` 
    * `make pkpd-necrosis-sample` 
    * `make pkpd-motility-sample` 
    * `make pkpd-combo-sample`
    * `make pkpd-confluence-start-sample` (currently deprecated)
3. Compile your model: `make`
4. Run your model
   * On Mac: `./pkpd_sample ./config/pkpd_model.xml`
   * On Windows: `pkpd_sample.exe .\config\pkpd_model.xml`
5. Look at the snapshots in `output/` and the living cell counts in `output/cell_counts.csv`

## Making your own project
There are two template projects available.
Load one with `make pkpd-template` or `make pkpd-template-sbml`.
The former sets up a model with two PKPD substrates, one with a dosing schedule defined by a CSV file.
The latter sets up a model with two PKPD substrates, one defined by an SBML model.

## Adding PhysiPKPD to an existing project
To add PhysiPKPD to an existing project, take the following steps:
1. In `main.cpp`...
   * add `setup_pharmacodynamics();` as a new line after `setup_tissue();`
   * add `PK_model( PhysiCell_globals.current_time );` immediately *before* `microenvironment.simulate_diffusion_decay( diffusion_dt );`
   * add `PD_model( PhysiCell_globals.current_time );` immediately *after* `microenvironment.simulate_diffusion_decay( diffusion_dt );`
2. In `custom_modules/custom.h`...
   * add `#include "../addons/PhysiPKPD/src/PhysiPKPD.h"` to the top
3. In the configuration file, e.g., `config/PhysiCell_settings.xml`...
   * add a `PK` element to every PK substrate (`microenvironment_setup//variable`) (see [PK Templates](#pk_templates))
   * add a `PD` element to every `cell_definition` that will be affected by a PD substrate; add a `substrate` block within this `PD` element for each PD substrate affecting the given cell definition (see [Templates](#pd_templates))
   * add `S_damage` to *every* `cell_definition//custom_data` for every PD substrate `S`
4. In `Makefile`...
   * add a line `PhysiPKPD_OBJECTS := PhysiPKPD_PK.o PhysiPKPD_PD.o`
   * add `$(PhysiPKPD_OBJECTS)` to the list of `PhysiCell_OBJECTS`
   * define `PhysiPKPD_PK.o` and `PhysiPKPD_PD.o` with the following lines
```
PhysiPKPD_PK.o: ./addons/PhysiPKPD/src/PhysiPKPD_PK.cpp
	$(COMPILE_COMMAND) -c ./addons/PhysiPKPD/src/PhysiPKPD_PK.cpp

PhysiPKPD_PD.o: ./addons/PhysiPKPD/src/PhysiPKPD_PD.cpp
	$(COMPILE_COMMAND) -c ./addons/PhysiPKPD/src/PhysiPKPD_PD.cpp
```
5. Ensure Dirichlet conditions are enabled for any PK substrates.
6. Ensure cell types that are affected by a PD substrate have a nonzero uptake rate for that substrate.
7. Ensure rules are enabled and set in the rules CSV file for any PD substrates.
<!-- ## Reconfiguring, editing, and re-running
**Note:** *While the below functionality is present, it is discouraged because it is likely to inadvertently affect the work of others.
Instead, it is recommended to instead save any changes to these files in a non-tracked directory and manually copy them into their proper places after `make`-ing the sample project.*

Instead of editing the configuration file copied into `PhysiCell/config/pkpd_model.xml`, you can choose to edit the original in `PhysiCell/sample_projects_physipkpd/[project_name]/config/pkpd_model.xml` to save the changes for future runs.
The command `make rc` will reconfigure from the original `pkpd_model.xml` to facilitate editing in the latter fashion.

Similarly, you can edit the custom modules in `PhysiCell/sample_projects_physipkpd/[project_name]/custom_modules/` to save changes for future runs.
After making these changes, you can run `make redo` and this will automatically move those changes to their proper places and recompile the project. -->

## Pharmacokinetics
PK models are added as `microenvironment_setup//variable//PK` element, i.e., at the same level as `Dirichlet_options`, of a substrate, `S`.
An `enabled="true"` attribute is required to turn on the PK model.

### PK models
There are three available `model`'s: `1C`, `2C`, and `SBML`[^sbml].
An example template for each is provided [below](#pk_templates).
For each `model`, the `circulation_concentration` of the PK model is what interacts with the microenvironment.
Specifically, any Dirichlet node for `S` is set as `circulation_concentration x biot_number`.
The `biot_number` is thus just the factor dropoff in substrate concentration from the blood to the perivascular niche.

[^sbml]: To use an SBML-defined PK model, you must have libRoadRunner installed.
If you can run the `sample_projects_intracellular/ode` projects, you are ready to run these PhysiPKPD models.

| Model | Description | Specification |
| :-- | :-- | :-: |
| 1-compartment | Circulation compartment with linear elimination | `1C` |
| 2-compartment | `1C` plus a periphery compartment with linear intercompartmental clearance rates | `2C` |
| SBML-defined | Any SBML-defined model. Place the file in the `./config/` folder. PhysiPKPD will look for the species named `circulation_concentration` to use for updating Dirichlet nodes | `SBML` |
<p align="center">
    <b>Table:</b> PK model specifications
</p>

| Parameter | Description | If Missing |
| :-- | :-- | :-- |
| `model` | Defines the PK model used for the substrate | Error
| `biot_number` | Ratio of substrate concentration on boundary of microenvironment (Dirichlet condition) and concentration in systemic circulation | Set to `1.0` |
<p align="center">
    <b>Table:</b> XML elements for all PK models
</p>

#### `1C` model <a name="1C"></a>
The `1C` model is a 1-comparment model defined by 
$$
\begin{aligned}
C' & = - \lambda C
\end{aligned}
$$
where $\lambda$ is the elimination rate of the `circulation_concentration`, $C$.

#### `2C` model
The `2C` model is a 2-comparment model defined by 
$$
\begin{aligned}
C' & = \frac{k_{21}}{R}P - k_{12}C - \lambda C \\
P' & = k_{12}RC - k_{21}P
\end{aligned}
$$
where $\lambda$ and $C$ are as [above](#1C), $P$ is the concentration in a peripheral compartment, though not the PhysiCell microenvironment.
The intercompartmental clearance rates, $k_{12}$ and $k_{21}$, define the exchange rate between these compartments while $R$ is the ratio of the compartmental volumes: $V_C/V_P$.

| Parameter | Models | Description | If Missing |
| :-- | :-: | :-- | :-- |
| `elimination_rate` $(\lambda)$ | `1C`,`2C` | Linear elimination rate in central compartment (in mintues<sup>-1</sup>) | Error |
| `volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `2C` | Ratio of central compartment to periphery compartment | Error for `2C` |
| `k12` $(k_{12})$ | `2C` | Rate of change in concentration in central compartment due to distribution (in minutes<sup>-1</sup>) | Error for `2C` |
| `k21` $(k_{21})$ | `double` | Rate of change in concentration in periphery compartment due to redistribution (in minutes<sup>-1</sup>) | Error for `2C` |
<!-- | `central_to_periphery_volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `double` | Ratio of central compartment to periphery compartment *for any substrates without a specific volume ratio as above* | Set to `1` | -->
<!-- | `S_flux_across_capillaries`<a name="old_flux_par"></a> | `double` | **Consider using the above parameters to quantify intercompartmental clearance rates.**[^1] Rate of change in concentration in central compartment due to distribution and redistribution (in minutes<sup>-1</sup>) | See above | -->
<p align="center">
    <b>Table:</b> PK parameters for `1C` and `2C` models
</p>

#### `SBML` model
Any SBML-defined ODE with one state variable named `circulation_concentration` can be used for a PK model as well.
Place the file in the `./config/` folder.
If `sbml_filename` is not set, PhysiPKPD will look for `./config/PK_default.xml` to apply.

| Parameter | Description | If Missing |
| :-- | :-- | :-- |
| `sbml_filename` | Filename of SBML file, e.g. `PK_default.xml` | Set to `PK_default.xml` |
<p align="center">
    <b>Table:</b> PK parameters for `SBML` models
</p>

### Dosing schedules

Every PK model requires a dosing schedule to achieve nonzero `circulation_concentration`.
For `SBML` models, these must be set as events in the SBML file.
For `1C` and `2C` models, these can be set either with parameters of a CSV file.
This is set with the `schedule` element and its attribute `format`, set to either `parameters` or `csv`.
See the [`1C` template](#1c_template) for an example using `parameters`, and see the [`2C` template](#2c_template) for an example using `csv`.

#### Dosing parameters

Two types of doses can be set this way, loading doses and regular doses.
The number of loading doses, `loading_doses`, determines how many loading doses are given before regular doses.
The total number of doses, loading *and* regular, is set by `total_doses`.
Loading doses are given first, then regular doses.
If `loading_doses` $>$ `total_doses`, then `total_doses` of loading doses are given and no regular doses.
The time between any consecutive doses, loading or regular, is `dose_interval`.
The `units` attribute of `dose_interval` can be set to `min`, `hours`, or `days`.
If none is supplied, `min` is used.
The time of the first dose is given by `first_dose_time` which can also use the `units` attribute exactly like `dose_interval`.
The dose given on a loading dose is set by `loading_dose` and on a regular dose by `regular_dose`.

| Parameter | Description | If Missing |
| :-- | :-- | :-- |
| `total_doses`| Total number of doses to give including loading doses | Error |
| `loading_doses` | Number of loading doses to give before switching to regular doses | Set to `0` |
| `regular_dose` | Increase in concentration in central compartment after a regular dose | If `total_doses`>`loading_doses`, throws an error |
| `loading_dose` | Increase in concentration in central compartment after a loading dose | If `loading_doses`>0, throws an error |
| `first_dose_time` | Time of first dose if given at fixed time (in `units="units"`) | Set to `current_time` |
| `dose_interval` | Time between successive doses, loading or regular (in `units="units"`) | If `total_doses`>1, throws an error |
<p align="center">
    <b>Table:</b> Dosing parameters
</p>

#### Dosing by CSV

If the `format` attribute in `schedule` is set to `"csv"`, PhysiPKPD will look for the file `./config/S_dose_schedule.csv` to define the dosing schedule.
This file consists of two columns: time of dose (min), dose amount.
No header row is needed.
<!-- | `S_read_dose_schedule_from_csv`[^dose_csv] | `bool` | If true, ignore the following parameters and read in a dosing schedule from `./config/S_dose_schedule.csv` | Set to `False` | -->
<!-- | `S_set_first_dose_time` | Boolean determining if the first dose time is fixed or if a confluence condition will be used to determine the first dose time | Set to `False` | -->
<!-- | `S_confluence_condition` | Proportion of microenvironment filled with cells at which to give first dose; confluence calculated by sum of cross-sectional area of all cells divided by area of microenvironment | Never check confluence condition | -->

<!-- [^dose_csv]: The CSV file for a dosing schedule must have the dosing times (in minutes) in the first column and the dose amounts in the second column.
No header row should be included. -->

<!-- ## Setting parameters
PhysiPKPD parameters are found in two areas in `pkpd_model.xml`: PK parameters are at the bottom in `user_parameters` and PD parameters are in `cell_definitions` in the `custom_data` for each cell type.
PhysiPKPD can add PK and/or PD dynamics to any substrates in a PhysiCell simulation.
PK dynamics must be set for each PK substrate.
PD dynamics must be set for each pairing `(S,C)` where substrate `S` acts on cell type `C`.
PhysiPKPD will attempt to use defaults if none are supplied and issue warnings/errors until you have fully specified everything.

In what follows, `S` stands for the name of a substrate, `C` stands for the name of a cell type, and `X` stands for a MOA.

### PK and PD substrates
You can specify which substrates you want to include PK dynamics.
You can do the same for PD dynamics.
These lists can be set independently of one another so that every PhysiCell substrate can have PK, PD, both, or neither.
The parameters in the following table are those PhysiPKPD will use to define which substrates follow PK and PD dynamics.
Add them to `user_parameters`.
Omitting `PKPD_pk_substrate_names` will result in no PK dynamics and similarly for the PD dynamics.
PhysiPKPD will use the comma-separated values in these parameters to search for the relevant PK/PD parameter names.

| Parameter | Description |
|:--|---|
| `PKPD_pk_substrate_names` | Comma-separated list of substrates also following PK dynamics |
| `PKPD_pd_substrate_names` | Comma-separated list of substrates also producing PD dynamics in some cells |
<p align="center">
    <b>Table:</b> Identifying PK and PD substrates
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

**Note:** Avoid spaces! -->

<!-- ### PK parameters <a name="pk_pars"></a>
There are three PK models to choose from.
To specify which one `S` uses, set `S_pk_model` according to Specification below.

If `S_pk_model` is not set, PhysiPKPD will default to a 2-compartment model[^oldpk]: -->

<!-- [^oldpk]: For those using the simplified 2-compartment model that PhysiPKPD used to use, see the [`S_flux_across_capillaries`](#old_flux_par) entry in the table. -->

<!-- $$
\begin{aligned}
C' & = \frac{k_{21}}{R}P - k_{12}C - \lambda C \\
P' & = k_{12}RC - k_{21}P
\end{aligned}
$$

For all PK models, a Biot number should be specified.
This describes the ratio of substrate concentration just outside capillary walls to that inside the blood vessel.
It can be related to a reflection coefficient. -->
<!-- 
#### Non-SBML models

For each substrate following `1C` or `2C`, you can set the following parameters in `user_parameters`.
For example, a substrate called `myDrug` with 10 doses administered would have the following in `user_parameters`: 

```
<myDrug_max_number_doses type="int">10</myDrug_max_number_doses>
``` -->

<!-- 
| Parameter | Type | Description | If Missing |
| :-- | :-: | :-- | :-- |
| `S_central_elimination_rate` $(\lambda)$ | `double` | Linear elimination rate in central compartment (in mintues<sup>-1</sup>) | Set to `0` |
| `S_central_to_periphery_volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `double` | Ratio of central compartment to periphery compartment | Set to `central_to_periphery_volume_ratio` |
| `S_central_to_periphery_clearance_rate` $(k_{12})$ | `double` | Rate of change in concentration in central compartment due to distribution (in minutes<sup>-1</sup>) | Set to `S_flux_across_capillaries`, if present. Otherwise, set to `0` |
| `S_periphery_to_central_clearance_rate` $(k_{21})$ | `double` | Rate of change in concentration in periphery compartment due to redistribution (in minutes<sup>-1</sup>) | Set to `S_flux_across_capillaries` $\times$ `S_central_to_periphery_volume_ratio`, if present. Otherwise, set to `0` |
| `central_to_periphery_volume_ratio` $(R = V_1/V_2 = V_C/V_P)$ | `double` | Ratio of central compartment to periphery compartment *for any substrates without a specific volume ratio as above* | Set to `1` |
<!-- | `S_flux_across_capillaries`<a name="old_flux_par"></a> | `double` | **Consider using the above parameters to quantify intercompartmental clearance rates.**[^1] Rate of change in concentration in central compartment due to distribution and redistribution (in minutes<sup>-1</sup>) | See above |
<p align="center">
    <b>Table:</b> PK parameters
</p> 

[^1]: To use these two parameters, you will want to set $k_{12}$ as your original flux rate and $k_{21}$ as `S_flux_across_capillaries * S_central_to_periphery_volume_ratio`. 
-->

<!-- #### SBML-defined models

You must supply a user parameter with the name of the SBML file:

Place this file in the `./config/` folder.
You must specify dosing events yourself in SBML files. -->

### PhysiCell PK parameters

You can also set the following parameters in `microenvironment_setup` for each substrate:
| Parameter | Description |
| :--| --- |
| `diffusion_coefficient` | Diffusion rate in the microenvironment |
| `decay_rate` | Rate of decay in the microenvironment |
<p align="center">
    <b>Table:</b> PK parameters in PhysiCell
</p>

## Pharmacodynamics
PD models are added as `cell_definition//PD` element, i.e.,  at the same level as `custom_data`.
Within this element, multiple `cell_definition//PD//substrate` elements define which substrates affect the given cell type.
Each `substrate` includes the attribute `name="S"` to set the substrate.

### PD models
Two similar models are supplied for PD, both using the principle of Area Under the Curve or AUC.
One uses the AUC of the internalized *concentration*, `AUC`, and the other uses the AUC of the internalized *amount*, `AUC_amount`.
PhysiPKPD will turn on `track_internalized_substrates` for you if it detects a PD model.
When using `AUC`, PhysiPKPD divides the internalized amount that PhysiCell tracks by the cell's current volume to determine the initial condition for $A$.
Both update according to the following differential equation:

$$
\begin{aligned}
A' & = -mA \\
D' & = A - r_1D - r_0
\end{aligned}
$$

where $A$ is the internalized quantity, amount or concentration, and $D$ is the accumulated damage or AUC[^auc_no_par].
The internalized quantity, $A$, is metabolized at a rate $m$, damage is repaired with linear rate $r_1$ and constant rate $r_0$.

[^auc_no_par]: Because $D$ is measuring the AUC, the coefficient for $A$ in $D'$ is fixed to $1$. That is, $D$ has units amount x min or concentration x min, depending on the model used.

| Parameter | Description | If Missing |
| :-- | :-- | :-- |
| `model` | Defines the PD model used for the substrate | Set to `AUC` |
| `metabolism_rate` $(m)$ | Rate of elimination of `S` from inside a cell (in minutes<sup>-1</sup>) | Error |
| `linear_repair_rate` $(r_1)$ | First-order elimination rate of damage from `S` (in minutes<sup>-1</sup>) | Error |
| `constant_repair_rate` $(r_0)$ | Zero-order elimination rate of damage from `S` (in damage per minute) | Error |
<a name="tab:pd_pars">
<p align="center">
    <b>Table:</b> XML elements for all PD models
</p>
</a>

Let us know if you would like to see a different PD model included.
Unlike PK dynamics, integration with SBML solvers is not implemented.

### The damage variable, `S_damage`
If a substrate, `S`, affects any cell type via a PD model, **every** cell type **must** include `custom_data//S_damage` in its `custom_data`.
This is where the damage, $D$, above is stored and accessed for affecting cell behavior.

<!-- 
### PD parameters <a name="pd_pars"></a>
The PD dynamics in PhysiPKPD work by the target cell type accumulating damage over time due to internalized substrate.
Thus, for each (substrate, cell type)-pairing, a PD model must be specified.
There are two models for this currently implemented in PhysiPKPD: damage accumulates based on either concentration or amount of the internalized substrate.
By default, PhysiPKPD will assume concentration is being used.
To specify this model, create a string parameter called `S_on_C_pd_model` with value selected from the following table:

| Model | Description | Specification |
| :-- | :-- | :-: |
| Concentration-based | Damage accumulates based on internal concentration | `AUC`[^auc] |
| Amount-based | Damage accumulates based on internal amount | `AUC_amount` | -->

<!-- <p align="center">
    <b>Table:</b> PD model specifications
</p> -->

<!-- [^auc]: Area under the curve (AUC) is a common metric in pharmaceutical sciences measuring the extent of exposure to a substance. -->


<!-- #### Damage Accumulation Parameters <a name="dam_pars"></a>
For both of the above models, the cell type, `C`, affected by substrate `S` accumulates damage, $D$, based on the internalized substrate, $A$.
The equation is unchanged if $A$ is measured in amount or concentration.
This accumulation of damage obeys the following differential equation:

$$
\begin{aligned}
A' & = -mA \\
D' & = A - r_1D - r_0
\end{aligned}
$$

Note that since damage is an abstract quantity, we do not include a rate parameter as a coefficient for $A$ in the equation for $D'$.
Each of these parameters **should** be set.<!--[^old_repair]. -->
<!-- If they are not set, PhysiPKPD will issue a warning and set them to 0.
These parameters go in the `custom_data` of `C`. -->

<!-- [^old_repair]: For backwards compatibility, `S_repair_rate` can be set instead of `S_repair_rate_constant`. -->

<!-- |Parameter|Description|
|:--|---|
| `S_damage` $(D)$ | Not a parameter; data that tracks the current damage to the cell |
| `S_repair_rate_constant` $(r_0)$ | Zero-order elimination rate of damage from `S` (in damage per minute) |
| `S_repair_rate_linear` $(r_1)$ | First-order elimination rate of damage from `S` (in minutes<sup>-1</sup>) |
| `S_metabolism_rate` $(m)$ | Rate of elimination of `S` from inside a cell (in minutes<sup>-1</sup>) |

<p align="center">
    <b>Table:</b> Damage accumulation parameters <a name="dam_pars_table"></a>
</p> -->

### Mechanisms of action
The effect of the PD model on the cell type is no longer defined by PhysiPKPD.
Instead, the user must supply rules that define how each cell type is affected by each PD substrate, if at all.
The intention of PhysiPKPD is that users will use the signal `custom:S_damage` to alter cell phenotype via the rules.
We recommend using PhysiCell-Studio to help generate a rules file.

### Optional PD parameters

#### Setting PD timesteps
By default, PhysiPKPD uses the `mechanics_dt` set in the configuration file to determine how often to update each PD model.
You can change this by adding a `PD//substrate//dt` element with the desired time in minutes.
For precomputation purposes, see [below](#precomputation), this time step *must* be a multiple of `diffusion_dt`.

#### Precomputation <a name="precomputation"></a>
PhysiPKPD currently requires all PD models to use precomputation, which speeds up calculations.
This relies on two assumptions: homogeneity of response to `S` within a given cell type and that the same amount of simulation time passes between updates.
The latter is assured by having the time step be a multiple of `diffusion_dt`.
The former is currently imposed by how PhysiPKPD uses the [PD parameters](#tab:pd_pars) uniformly for all cells within a cell type.
If you want these parameters to vary by cells *within* a cell type, that is not currently supported.

| Parameter | Description | If Missing |
| :-- | :-- | :-- |
| `dt` | Sets the time interval between PD updates (in minutes) | Set to `mechanics_dt` |
| `precompute` $(m)$ | Boolean for precomputations | **Must be** `true` |
<a name="tab:pd__optional_pars">
<p align="center">
    <b>Table:</b> Optional XML elements for all PD models
</p>
</a>

<!-- For example, if you wish to use a timestep of `0.01` for updating the damage of `myDrug` on `tumor`, then add the following to `user_parameters`:

```
<myDrug_dt_tumor type="double">0.01</myDrug_dt_tumor>
``` -->

<!-- 
The accumulated damage then affects the cell through one or more mechanisms of action (MOA) using the rules grammar.
Each PD substrate will have an associated damage variable called `custom:S_damage` that is designed to be used as the signal for the rule.
As an example, if a PD substrate called `myDrug` affects cell type `tumor` by increasing its apoptosis, add the following row to the rules CSV:
```
tumor,custom:myDrug_damage,increases,apoptosis,1,1e5,2,0
```
The `1,1e5,2,0` are the saturation apoptosis rate, damage EC50, hill coefficient, and Boolean for applying the rule to dead cells, respectively. -->

<!-- #### Cell Effect Parameters
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
Just take care that before calling `pd_phenotype_function` within `update_phenotype` that the affected rate is reset to the base value lest PhysiPKPD compound these effects each time.
See the sample projects targeting `proliferation`, `apop`, and `necrosis`.
If you target `motility`, do the same but within `custom_cell_rule` before calling `pd_custom_function`.

|Parameter|Description|
|:--|:--|
| `S_moa_is_X` | Used as boolean to determine which effects to apply to this cell type based on the damage from `S`; values > 0.5 will apply the effect |
| `S_X_saturation_rate` | Rate of `X` as damage from `S` approaches infinity |
| `S_X_EC50` | Damage from `S` at which the rate of `X` is halfway between the base and saturation rates (in damage) |
| `S_X_hill_power` | Hill coefficient for calculating the effect of `S` on the rate of `X` |

<p align="center">
    <b>Table:</b> Cell effect parameters
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
If, for example the current apoptosis rate of `C` is not `1E-5` but instead `1E-6` perhaps because it is in a quiescent niche, then PhysiPKPD will update its apoptosis rate to $12,500 \cdot 10^{-6} = 0.0125\ \text{min}^{-1}$. -->
<!-- 
### Precomputation
Currently, all (non-SBML) PK and PD models implemented in PhysiPKPD are readily solved using analytic techniques.
By informal observation, the analytic methods are not slower than numerical methods.
If anything, the analtyic methods are faster.
Therefore, all simulations use analytic solutions.

To maximize the efficiency of these analytic solutions, many terms are pre-computed ahead of time.
**However**, if the parameters governing your PD dynamics vary&mdash;e.g., due to heterogeneity within the affected cell type&mdash;then you will not want to do pre-computations.
You can include the following as Booleans in `user_parameters` to control this.
For example, if a substrate `myDrug` affects cell type `tumor` and you want to pre-compute the PD quantities, include the following in `user_parameters`:

```
<myDrug_precompute_pd_for_tumor type="bool">True</myDrug_precompute_pd_for_tumor>
```

|Parameter|Type|Description| If Missing |
|---|:-:|---|:--|
| `PKPD_precompute_all_pd_quantities` | `bool` | Boolean to determine the default behavior for pre-computations for *all* PD dynamics | Set to `True` |
| `S_precompute_pd_for_C` | `bool` | Boolean to override default behavior for a single (substrate, cell type) pairing; **turn this off if PD parameters can vary within this cell type OR if your PD time step is not a multiple of your `diffusion_dt`** | No effect |
<p align="center">
    <b>Table:</b> Miscellaneous parameters
</p> -->

## Templates <a name="template"></a>

### PK templates <a name="pk_templates"></a>
Place any of the following within the `variable` element:
```
<microenvironment_setup>
   <variable name="PKPD_D1" units="dimensionless" ID="0">
      ...
      </Dirichlet_options>

      <--! Place PK model here -->

   </variable>
   ...
</microenvironment_setup>
```

#### `1C` template <a name="1c_template"></a>


```
<PK enabled="true">
    <model>1C</model>
    <schedule format="parameters">
        <total_doses>4</total_doses>
        <loading_doses>0</loading_doses>
        <first_dose_time units="min">180</first_dose_time>
        <dose_interval units="min">360</dose_interval>
        <regular_dose>500</regular_dose>
        <loading_dose>1000</loading_dose>
    </schedule>
    <elimination_rate units="1/min">0.0027</elimination_rate>
    <biot_number>1</biot_number>
</PK>
```

#### `2C` template <a name="2c_template"></a>

```
<PK enabled="true">
    <model>2C</model>
    <schedule format="csv" />
    <elimination_rate units="1/min">0.0027</elimination_rate>
    <k12 units="1/min">0.0048</k12>
    <k21 units="1/min">.0048</k21>
    <volume_ratio>1</volume_ratio>
    <biot_number>1</biot_number>
</PK>
```

#### `SBML` template

```
<PK enabled="true">
    <model>SBML</model>
    <sbml_filename>PK_default.xml</sbml_filename>
    <biot_number>1</biot_number>
</PK>
```

### PD template <a name="pd_templates"></a>
Place the following within the `cell_definition` element:
```
<cell_definitions>
   <cell_definition name="PKPD_cell" ID="0">
      ...
      </phenotype>

      <--! Place PD model(s) here -->

      </cell_definition>
   ...
</cell_definitions>
```

```
<PD>
    <substrate name="PKPD_D1">
        <model>AUC</model>
        <metabolism_rate>0.05</metabolism_rate>
        <constant_repair_rate>5e3</constant_repair_rate>
        <linear_repair_rate>2e-2</linear_repair_rate>
        <precompute>true</precompute>
        <dt>0.1</dt>
    </substrate>
    <substrate name="PKPD_D2">
        <model>AUC</model>
        <metabolism_rate>0.05</metabolism_rate>
        <constant_repair_rate>5e3</constant_repair_rate>
        <linear_repair_rate>2e-2</linear_repair_rate>
        <precompute>true</precompute>
        <dt>0.1</dt>
    </substrate>
</PD>
```
<!-- 
## `user_parameters` template <a name="up_template"></a>
Add the following block to `user_parameters` one time to include PK and PD substrates.
Fill in the names of PK and PD substrates in the corresponding parameter lists.

```
<PKPD_pk_substrate_names type="string" description="names of substrates following pk dynamics">S1,S2</PKPD_pk_substrate_names>
<PKPD_pd_substrate_names type="string" description="names of substrates following pk dynamics">S1,S3</PKPD_pd_substrate_names>
<PKPD_precompute_all_pd_quantities type="bool" units="dimensionless" description="If using analytic solutions (see PKPD_use_analytic_pd_solutions), whether to precompute or not. precompute should not be done if pd parameters vary by agent/time">True</PKPD_precompute_all_pd_quantities>
<csv_data_interval type="double" units="min" description="time between writing data to CSV">30</csv_data_interval>
```
For each PK substrate `S`, add the following parameters.
Replace `S` with the name of the substrate. and `C` should be replaced with the name of the cell definition.
```
<S_pk_model type="string" description="What PK model to use. Some options: 2C">1C</S_pk_model>
<S_read_dose_schedule_from_csv type="bool" units="none" description="if true, reads in a dosing schedule from a csv">False</S_read_dose_schedule_from_csv>
<S_biot_number type="double" units="dimensionless" description="ratio of drug concentration outside capillary to inside capillary (or just all systemic circulation)">1.0</S_biot_number>

<S_max_number_doses type="int" units="none" description="max doses, including loading doses">4</S_max_number_doses>
<S_number_loading_doses type="int" units="none" description="number of loading doses">0</S_number_loading_doses>
<S_dose_interval type="double" units="min" description="time between dosings">360.0</S_dose_interval>
<S_central_increase_on_dose type="double" units="mmHg" description="concentration added to systemic circulation on dose">5e2</S_central_increase_on_dose>
<S_central_increase_on_loading_dose type="double" units="mmHg" description="concentration added to systemic circulation on dose">400.0</S_central_increase_on_loading_dose>

<S_set_first_dose_time type="bool" units="none" description="user sets first dose time to start therapy for drug 1; false=user gives confluence condition">True</S_set_first_dose_time>
<S_first_dose_time type="double" units="min" description="time of first dose for drug 1">0.0</S_first_dose_time>
<S_confluence_condition type="double" units="none" description="confluence as proportion at which to start therapy for drug 1">0.6</S_confluence_condition>

<S_central_elimination_rate type="double" units="1/min" description="decay rate in systemic circulation">0.0027</S_central_elimination_rate>
<S_central_to_periphery_clearance_rate type="double" units="1/min" description="rate of central concentration change due to distribution">0.0048</S_central_to_periphery_clearance_rate>
<S_periphery_to_central_clearance_rate type="double" units="1/min" description="rate of periphery concentration change due to redistribution">0.0048</S_periphery_to_central_clearance_rate>
<S_central_to_periphery_volume_ratio type="double" units="dimensionless" description="ratio of volume in systemic circulation to the volume of the periphery">1.0</S_central_to_periphery_volume_ratio>
```
If `S` is a PD substrate affecting cell type `C`, add the following to `user_parameters`, replacing `S` with the name of the substrate and `C` with the name of the cell type.
```
<S_on_PKPD_C_pd_model type="string">AUC</S_on_PKPD_C_pd_model>
<S_dt_PKPD_C type="double">0.1</S_dt_PKPD_C>
<S_precompute_pd_for_C type="bool">True</S_precompute_pd_for_C>
```

```
<S_damage>0</S_damage>
<S_repair_rate_constant units="damage/min">0</S_repair_rate_constant>
<S_repair_rate_linear units="damage/min">2e-2</S_repair_rate_linear>
<S_metabolism_rate units="1/min">0.05</S_metabolism_rate>
``` -->

<!-- ## Making your own project using PhysiPKPD
If you wish to make your own project that uses PhysiPKPD (and not just one of the pre-built sample projects), this is how you can proceed.
1. Make one of the PKPD template projects: `make pkpd-template` or `make pkpd-template-sbml`
2. Edit the configuration file to set the Dirichlet conditions, [PK Parameters](#pk_pars), and [PD Parameters](#pd_pars) for the two PKPD substrates and the default cell type `cell`.
3. Add additional substrates as normal (using the Model Builder for this is untested)
4. Add additional cell types as normal (using the Model Builder for this is untested). -->
<!-- 5. By default, each cell type is assigned the same `update_phenotype` function, which is `cell_phenotype` found in the `custom.cpp` file. -->
<!-- Add new phenotype functions as desired for each cell type. -->
<!-- 5. Notice that the `update_phenotype` function resets the three MOA-targeted rates before calling `pd_phenotype_function`. -->
<!-- 7. By default, each cell type is assigned the same `custom_cell_rule` function, which is `custom_function` found in the `custom.cpp` file. -->
<!-- 8. Notice that the `custom_cell_rule` function resets the migration speed before calling `pd_custom_function`. -->

<!-- **Note:** The order in which PhysiCell updates the various aspects of a cell matters!
For example, intracellular updates happen first.
So, if those intracellular updates affect proliferation, apoptosis, or necrosis rates, then resetting those rates in your `update_phenotype` function will erase those changes.
Similarly, `custom_cell_rule` is called after `update_phenotype` and prior to `update_velocity`, and `update_velocity` in turn calls `update_migration_bias`.
Thus, changes to migration speed in `update_phenotype` or from an intracellular model would be erased if `custom_cell_rule` resets migration speed.
Resetting the migration speed in `update_migration_bias` erases the PD effects calculated in `custom_cell_rule`. -->

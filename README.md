# PhysiPKPD
## Getting started...
### ...with a new project
1. Download (or copy contents of) `get_physipkpd.py` and `helpers_for_get_and_add.py` to local files in the same directory.
2. Run 
```
python get_physipkpd.py --dir path/to/project/dir --studio
```
   * The `--studio` flag downloads a PhysiPKPD-enabled version of PhysiCell-Studio to help build the model
3. Follow output for next steps

### ...adding PhysiPKPD to an existing project
1. Download (or copy contents of) `add_physipkpd.py` and `helpers_for_get_and_add.py` to local files in the same directory.
2. Run
```
python add_physipkpd.py path/to/project/dir --studio
```
   * The `--studio` flag downloads a PhysiPKPD-enabled version of PhysiCell-Studio to help build the model
3. Follow output for next steps

**Note**: This script will perform the steps [below](#AddingManually), except updates to the configuration file. Specifically, it will...
*  add lines to your project's `main.cpp`, `custom.h`, `custom.cpp`, and `Makefile`
*  add PhysiPKPD to your project's `addons` if it does not exist already
*  create a new copy of the PhysiPKPD-enabled studio at the same directory level as your project

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

## Adding PhysiPKPD to an existing project <a name="AddingManually"></a>
To add PhysiPKPD to an existing project, take the following steps:
1. In `main.cpp`...
   * add `setup_pharmacodynamics();` as a new line after `setup_tissue();`
   * add `PK_model( PhysiCell_globals.current_time );` immediately *before* `microenvironment.simulate_diffusion_decay( diffusion_dt );`
   * add `PD_model( PhysiCell_globals.current_time );` immediately *after* `microenvironment.simulate_diffusion_decay( diffusion_dt );`
2. In `custom_modules/custom.h`...
   * add `#include "../addons/PhysiPKPD/src/PhysiPKPD.h"` to the top
3. In the configuration file, e.g., `config/PhysiCell_settings.xml`...
   * add a `PK` element to every PK substrate (`microenvironment_setup//variable`) (see [PK Templates](#pk_templates))
   * add a `PD` element to every `cell_definition` that will be affected by a PD substrate; add a `substrate` block within this `PD` element for each PD substrate affecting the given cell definition (see [PD Templates](#pd_templates))
   * add `S_damage` to *every* `cell_definition//custom_data` for every PD substrate `S` (even those unaffected by `S` because PhysiCell requires all cells to have the same custom data)
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

## Pharmacokinetics
PK models are added as `microenvironment_setup//variable//PK` element, i.e., at the same level as `Dirichlet_options`, of a substrate, `S`.
An `enabled="true"` attribute is required to turn on the PK model.

### PK models
There are four available `model`'s: `Constant`, `1C`, `2C`, and `SBML`[^sbml].
An example template for each is provided [below](#pk_templates).
For each `model`, the `circulation_concentration` of the PK model is what interacts with the microenvironment.
Specifically, any Dirichlet node for `S` is set as `circulation_concentration x biot_number`.
The `biot_number` is thus just the factor dropoff in substrate concentration from the blood to the perivascular niche.

[^sbml]: To use an SBML-defined PK model, you must have libRoadRunner installed.
If you can run the `sample_projects_intracellular/ode` projects, you are ready to run these PhysiPKPD models.

| Model | Description | Specification |
| :-- | :-- | :-: |
| Constant | Piecewise constant circulation compartment | `Constant` |
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

#### `Constant` model <a name="Constant"></a>
The `Constant` model fixes the circulation compartment concentration at user-specified times and values.
The times and values must be supplied in a CSV (see [below](#csv_dosing)).

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

For `Constant` models, use a CSV in the same format as the `1C` and `2C` models, i.e., times in column one and values in column two.
In this case, the values are what the circulation compartment concentration is set to rather than the amount added to that compartment.
It will remain at that value indefinitely, so to turn "off" the substrate, add a (time, value) pair of (t<sub>off</sub>, 0).

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

#### Dosing by CSV <a name="csv_dosing"></a>

If the `format` attribute in `schedule` is set to `"csv"`, PhysiPKPD will look for the file `./config/S_dose_schedule.csv` to define the dosing schedule.
Alternatively, you can specify the path similar to the `cell_positions` and `cell_rules` files by supplying a `folder` and `filename` (see the [2C template](#2c_template)).
This file consists of two columns: time of dose (min), dose amount.
No header row is needed.

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
**See templates [here](#pd_templates).**

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

### Mechanisms of action
The effect of the PD model on the cell type is no longer defined by PhysiPKPD.
Instead, you must supply rules that define how each cell type is affected by each PD substrate, if at all.
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
| `precompute` | Boolean for precomputations | **Must be** `true` |
<a name="tab:pd__optional_pars">
<p align="center">
    <b>Table:</b> Optional XML elements for all PD models
</p>
</a>

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

#### `Constant` template <a name="constant_template"></a>

```
<PK enabled="true">
    <model>Constant</model>
    <schedule format="csv" />
    <biot_number>1</biot_number>
</PK>
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
    <schedule format="csv">
        <folder>./config</folder>
        <filename>S_dose_schedule.csv</filename>
    </schedule>
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
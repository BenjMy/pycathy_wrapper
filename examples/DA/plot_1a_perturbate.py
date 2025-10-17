"""
Creating and Managing Parameter Perturbation Scenarios for Data Assimilation
============================================================================

This notebook illustrates how to create, read, and prepare parameter perturbation scenarios
(using placeholders or suggested values) for Data Assimilation (DA) in pyCATHY.

The notebook explains how to:

1. **Initialize a DA simulation**
2. **Create scenario files** with placeholders or suggested values
3. **Read and inspect scenario files**
4. **Apply perturbations to the simulation**

*Estimated time to run the notebook: ~2 min*
"""
#%%
# Initialize the DA simulation
# ----------------------------
# Create a DA object that will manage your simulations, parameter ensembles, and perturbations.
from pyCATHY.DA.cathy_DA import DA

simu = DA(
    dirName='.',
    prj_name='test_parameters_perturbation'
)

#%%
# Create a scenario file using placeholders
# -----------------------------------------
# This creates a generic scenario where parameters are indicated with placeholders
# like <value>, <sigma>, etc. It is useful for template purposes.
from pyCATHY import cathy_utils

cathy_utils.create_scenario_file(
    scenario_name="scenario_placeholders",
    param_names=["ic", "Ks"],
    use_common_values=False,  # Use placeholders instead of suggested values
    filetype="json"
)

#%%
# Create a scenario file using suggested values
# ---------------------------------------------
# Here we create a scenario with realistic suggested values for parameters.
cathy_utils.create_scenario_file_single_control(
    "scenario",
    parameters=["ic", "Ks"],
    use_suggested=True,
    control_type=None  # 'layers', 'zone', 'root_map', or None
)

# Read the scenario file and convert it to a DataFrame for inspection
scenario = cathy_utils.read_scenario_file('scenario.json')
scenario_df = cathy_utils.scenario_dict_to_df_list(scenario, 'scenario')
print(scenario_df)

#%%
# Scenario with perturbations by zones
# -----------------------------------
# Generate a scenario where perturbations are controlled per zone.
cathy_utils.create_scenario_file_single_control(
    "scenario_zone",
    parameters=["ic", "Ks"],
    use_suggested=True,
    control_type='zone',  # Control perturbations per zone
    nzones=2
)

scenario_multipleZones = cathy_utils.read_scenario_file('scenario_zone.json')
scenario_df = cathy_utils.scenario_dict_to_df_list(scenario_multipleZones, 'scenario_zone')
print(scenario_df)

#%%
# Scenario with perturbations by root map
# --------------------------------------
# Create a scenario where perturbations are controlled per root map element.
cathy_utils.create_scenario_file_single_control(
    "scenario_root_map",
    parameters=["ic", "porosity"],
    use_suggested=False,  # Use placeholders
    control_type='root_map',
    nveg=2
)

scenario_multipleLayers = cathy_utils.read_scenario_file('scenario_root_map.json')
scenario_df = cathy_utils.scenario_dict_to_df_list(scenario_multipleLayers, 'scenario_root_map')
print(scenario_df)

#%%
# Apply perturbations to the DA simulation
# ----------------------------------------
# This step generates ensembles of perturbed parameters according to the scenario.
from pyCATHY.DA import perturbate

simu.MAXVEG = 1  # Example configuration
list_pert = perturbate.perturbate(
    simu,
    scenario_multipleZones['scenario_zone'],
    256  # Ensemble size
)

var_per_dict_stacked = {}
for dp in list_pert:
    var_per_dict_stacked = perturbate.perturbate_parm(
        var_per_dict_stacked,
        parm=dp,
        type_parm=dp['type_parm'],  # Can also be VAN GENUCHTEN parameters
        mean=dp['mean'],
        sd=dp['sd'],
        sampling_type=dp['sampling_type'],
        ensemble_size=dp['ensemble_size'],  # Size of the ensemble
        per_type=dp['per_type'],
        savefig=dp['savefig']
    )

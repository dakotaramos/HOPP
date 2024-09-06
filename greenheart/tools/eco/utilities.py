import os
import os.path
import yaml
import copy

import numpy as np
import numpy_financial as npf
import pandas as pd
import openpyxl

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker

import ORBIT as orbit

from hopp.simulation.technologies.resource.wind_resource import WindResource
from hopp.simulation.technologies.resource.cambium_data import CambiumData

from hopp.simulation import HoppInterface

from hopp.utilities import load_yaml

from hopp.simulation.technologies.dispatch import plot_tools

from .finance import adjust_orbit_costs

from hopp import ROOT_DIR

"""
This function returns the ceiling of a/b (rounded to the nearest greater integer). 
The function was copied from https://stackoverflow.com/a/17511341/5128616
"""
def ceildiv(a, b):
    return -(a // -b)

# Function to load inputs
def get_inputs(
    filename_hopp_config,
    filename_greenheart_config,
    filename_orbit_config,
    filename_turbine_config,
    filename_floris_config=None,
    verbose=False,
    show_plots=False,
    save_plots=False,
):

    ############### load turbine inputs from yaml

    # load turbine inputs
    turbine_config = load_yaml(filename_turbine_config)

    # load hopp inputs
    hopp_config = load_yaml(filename_hopp_config)

    # load eco inputs
    greenheart_config = load_yaml(filename_greenheart_config)

    ################ load plant inputs from yaml
    if filename_orbit_config != None:
        orbit_config = orbit.load_config(filename_orbit_config)

        # print plant inputs if desired
        if verbose:
            print("\nPlant configuration:")
            for key in orbit_config.keys():
                print(key, ": ", orbit_config[key])

        # check that orbit and hopp inputs are compatible
        if (
            orbit_config["plant"]["capacity"]
            != hopp_config["technologies"]["wind"]["num_turbines"]
            * hopp_config["technologies"]["wind"]["turbine_rating_kw"]
            * 1e-3
        ):
            raise (
                ValueError("Provided ORBIT and HOPP wind plant capacities do not match")
            )

    # update floris_config file with correct input from other files
    # load floris inputs
    if (
        hopp_config["technologies"]["wind"]["model_name"] == "floris"
    ):  # TODO replace elements of the file
        if filename_floris_config is None:
            raise (ValueError("floris input file must be specified."))
        else:
            floris_config = load_yaml(filename_floris_config)
            floris_config.update({"farm": {"turbine_type": turbine_config}})
    else:
        floris_config = None

    # print turbine inputs if desired
    if verbose:
        print("\nTurbine configuration:")
        for key in turbine_config.keys():
            print(key, ": ", turbine_config[key])

    ############## provide custom layout for ORBIT and FLORIS if desired
    if filename_orbit_config != None:
        if orbit_config["plant"]["layout"] == "custom":
            # generate ORBIT config from floris layout
            for i, x in enumerate(floris_config["farm"]["layout_x"]):
                floris_config["farm"]["layout_x"][i] = x + 400

            layout_config, layout_data_location = convert_layout_from_floris_for_orbit(
                floris_config["farm"]["layout_x"],
                floris_config["farm"]["layout_y"],
                save_config=True,
            )

            # update orbit_config with custom layout
            # orbit_config = orbit.core.library.extract_library_data(orbit_config, additional_keys=layout_config)
            orbit_config["array_system_design"]["location_data"] = layout_data_location

    # if hybrid plant, adjust hybrid plant capacity to include all technologies
    total_hybrid_plant_capacity_mw = 0.0
    for tech in hopp_config["technologies"].keys():
        if tech == "grid":
            continue
        elif tech == "wind":
            total_hybrid_plant_capacity_mw += (
                hopp_config["technologies"][tech]["num_turbines"]
                * hopp_config["technologies"][tech]["turbine_rating_kw"]
                * 1e-3
            )
        elif tech == "pv":
            total_hybrid_plant_capacity_mw += (
                hopp_config["technologies"][tech]["system_capacity_kw"] * 1e-3
            )
        elif tech == "wave":
            total_hybrid_plant_capacity_mw += (
                hopp_config["technologies"][tech]["num_devices"]
                * hopp_config["technologies"][tech]["device_rating_kw"]
                * 1e-3
            )

    # initialize dict for hybrid plant
    if filename_orbit_config != None:
        if total_hybrid_plant_capacity_mw != orbit_config["plant"]["capacity"]:
            orbit_hybrid_electrical_export_config = copy.deepcopy(orbit_config)
            orbit_hybrid_electrical_export_config["plant"][
                "capacity"
            ] = total_hybrid_plant_capacity_mw
            orbit_hybrid_electrical_export_config["plant"].pop(
                "num_turbines"
            )  # allow orbit to set num_turbines later based on the new hybrid capacity and turbine rating
        else:
            orbit_hybrid_electrical_export_config = {}

    if verbose:
        print(
            f"Total hybrid plant rating calculated: {total_hybrid_plant_capacity_mw} MW"
        )

    if filename_orbit_config is None:
        orbit_config = None
        orbit_hybrid_electrical_export_config = {}

    ############## return all inputs

    return (
        hopp_config,
        greenheart_config,
        orbit_config,
        turbine_config,
        floris_config,
        orbit_hybrid_electrical_export_config,
    )


def convert_layout_from_floris_for_orbit(turbine_x, turbine_y, save_config=False):

    turbine_x_km = (np.array(turbine_x) * 1e-3).tolist()
    turbine_y_km = (np.array(turbine_y) * 1e-3).tolist()

    # initialize dict with data for turbines
    turbine_dict = {
        "id": list(range(0, len(turbine_x))),
        "substation_id": ["OSS"] * len(turbine_x),
        "name": list(range(0, len(turbine_x))),
        "longitude": turbine_x_km,
        "latitude": turbine_y_km,
        "string": [0] * len(turbine_x),  # can be left empty
        "order": [0] * len(turbine_x),  # can be left empty
        "cable_length": [0] * len(turbine_x),
        "bury_speed": [0] * len(turbine_x),
    }
    string_counter = -1
    order_counter = 0
    for i in range(0, len(turbine_x)):
        if turbine_x[i] - 400 == 0:
            string_counter += 1
            order_counter = 0

        turbine_dict["order"][i] = order_counter
        turbine_dict["string"][i] = string_counter

        order_counter += 1

    # initialize dict with substation information
    substation_dict = {
        "id": "OSS",
        "substation_id": "OSS",
        "name": "OSS",
        "longitude": np.min(turbine_x_km) - 200 * 1e-3,
        "latitude": np.average(turbine_y_km),
        "string": "",  # can be left empty
        "order": "",  # can be left empty
        "cable_length": "",
        "bury_speed": "",
    }

    # combine turbine and substation dicts
    for key in turbine_dict.keys():
        # turbine_dict[key].append(substation_dict[key])
        turbine_dict[key].insert(0, substation_dict[key])

    # add location data
    file_name = "osw_cable_layout"
    save_location = "./input/project/plant/"
    # turbine_dict["array_system_design"]["location_data"] = data_location
    if save_config:
        if not os.path.exists(save_location):
            os.makedirs(save_location)
        # create pandas data frame
        df = pd.DataFrame.from_dict(turbine_dict)

        # df.drop("index")
        df.set_index("id")

        # save to csv
        df.to_csv(save_location + file_name + ".csv", index=False)

    return turbine_dict, file_name


def visualize_plant(
    hopp_config,
    greenheart_config,
    turbine_config,
    wind_cost_outputs,
    hopp_results,
    platform_results,
    desal_results,
    h2_storage_results,
    electrolyzer_physics_results,
    design_scenario,
    colors,
    plant_design_number,
    show_plots=False,
    save_plots=False,
    output_dir="./output/",
):
    # save plant sizing to dict
    component_areas = {}

    plt.rcParams.update({"font.size": 7})

    if hopp_config["technologies"]["wind"]["model_name"] != "floris":
        raise (
            NotImplementedError(
                f"`visualize_plant()` only works with the 'floris' wind model, `model_name` \
                                  {hopp_config['technologies']['wind']['model_name']} has been specified"
            )
        )

    # set colors
    turbine_rotor_color = colors[0]
    turbine_tower_color = colors[1]
    pipe_color = colors[2]
    cable_color = colors[8]
    electrolyzer_color = colors[4]
    desal_color = colors[9]
    h2_storage_color = colors[6]
    substation_color = colors[7]
    equipment_platform_color = colors[1]
    compressor_color = colors[0]
    if hopp_config["site"]["solar"]:
        solar_color = colors[2]
    if hopp_config["site"]["wave"]:
        wave_color = colors[8]
    battery_color = colors[8]

    # set hatches
    solar_hatch = "//"
    wave_hatch = "\\\\"
    battery_hatch = "+"
    electrolyzer_hatch = "///"
    desalinator_hatch = "xxxx"

    # Views
    # offshore plant, onshore plant, offshore platform, offshore turbine

    # get plant location

    # get shore location

    # get cable/pipe locations
    if design_scenario["wind_location"] == "offshore":
        cable_array_points = (
            wind_cost_outputs.orbit_project.phases["ArraySystemDesign"].coordinates
            * 1e3
        )  # ORBIT gives coordinates in km, convert to m
        pipe_array_points = (
            wind_cost_outputs.orbit_project.phases["ArraySystemDesign"].coordinates
            * 1e3
        )  # ORBIT gives coordinates in km, convert to m

        # get turbine tower base diameter
        tower_base_diameter = wind_cost_outputs.orbit_project.config["turbine"][
            "tower"
        ]["section_diameters"][
            0
        ]  # in m
        tower_base_radius = tower_base_diameter / 2.0

        # get turbine locations
        turbine_x = (
            wind_cost_outputs.orbit_project.phases[
                "ArraySystemDesign"
            ].turbines_x.flatten()
            * 1e3
        )  # ORBIT gives coordinates in km, convert to m
        turbine_x = turbine_x[~np.isnan(turbine_x)]
        turbine_y = (
            wind_cost_outputs.orbit_project.phases[
                "ArraySystemDesign"
            ].turbines_y.flatten()
            * 1e3
        )  # ORBIT gives coordinates in km, convert to m
        turbine_y = turbine_y[~np.isnan(turbine_y)]

        # get offshore substation location and dimensions
        substation_x = (
            wind_cost_outputs.orbit_project.phases["ArraySystemDesign"].oss_x * 1e3
        )  # ORBIT gives coordinates in km, convert to m (treated as center)
        substation_y = (
            wind_cost_outputs.orbit_project.phases["ArraySystemDesign"].oss_y * 1e3
        )  # ORBIT gives coordinates in km, convert to m (treated as center)
        substation_side_length = 20  # [m] just based on a large substation (https://www.windpowerengineering.com/making-modern-offshore-substation/) since the dimensions are not available in ORBIT

        # get equipment platform location and dimensions
        equipment_platform_area = platform_results["toparea_m2"]
        equipment_platform_side_length = np.sqrt(equipment_platform_area)
        equipment_platform_x = (
            substation_x - substation_side_length - equipment_platform_side_length / 2
        )  # [m] (treated as center)
        equipment_platform_y = substation_y  # [m] (treated as center)

        # get platform equipment dimensions
        if design_scenario["electrolyzer_location"] == "turbine":
            desal_equipment_area = desal_results[
                "per_turb_equipment_footprint_m2"
            ]  # equipment_footprint_m2
        elif design_scenario["electrolyzer_location"] == "platform":
            desal_equipment_area = desal_results["equipment_footprint_m2"]
        else:
            desal_equipment_area = 0

        desal_equipment_side = np.sqrt(desal_equipment_area)

        # get pipe points
        pipe_x = np.array([substation_x - 1000, substation_x])
        pipe_y = np.array([substation_y, substation_y])

        # get cable points
        cable_x = pipe_x
        cable_y = pipe_y

    else:
        turbine_x = np.array(
            hopp_config["technologies"]["wind"]["floris_config"]["farm"]["layout_x"]
        )
        turbine_y = np.array(
            hopp_config["technologies"]["wind"]["floris_config"]["farm"]["layout_y"]
        )
        cable_array_points = []
    
    # wind farm area
    turbine_length_x = np.max(turbine_x)-np.min(turbine_x)
    turbine_length_y = np.max(turbine_y)-np.min(turbine_y)
    turbine_area = turbine_length_x * turbine_length_y

    # compressor side # not sized
    compressor_area = 25
    compressor_side = np.sqrt(compressor_area)
    ## create figure
    fig, ax = plt.subplots(2, 2, figsize=(10, 6))

    # get turbine rotor diameter
    rotor_diameter = turbine_config["rotor_diameter"]  # in m
    rotor_radius = rotor_diameter / 2.0

    # set onshore substation dimensions
    onshore_substation_x_side_length = 127.25  # [m] based on 1 acre area https://www.power-technology.com/features/making-space-for-power-how-much-land-must-renewables-use/
    onshore_substation_y_side_length = 31.8  # [m] based on 1 acre area https://www.power-technology.com/features/making-space-for-power-how-much-land-must-renewables-use/
    onshore_substation_area = onshore_substation_x_side_length * onshore_substation_y_side_length

    if greenheart_config["h2_storage"]["type"] == "pressure_vessel":
        h2_storage_area = h2_storage_results["tank_footprint_m2"]
        h2_storage_side = np.sqrt(h2_storage_area)
    else:
        h2_storage_side = 0
        h2_storage_area = 0

    electrolyzer_area = electrolyzer_physics_results["equipment_footprint_m2"]
    if design_scenario["electrolyzer_location"] == "turbine":
        electrolyzer_area /= hopp_config["technologies"]["wind"]["num_turbines"]

    electrolyzer_side = np.sqrt(electrolyzer_area)

    # set onshore origin
    onshorex = 50
    onshorey = 50

    wind_buffer = np.min(turbine_x) - (
        onshorey + 2 * rotor_diameter + electrolyzer_side
    )
    if "pv" in hopp_config["technologies"].keys():
        wind_buffer -= np.sqrt(hopp_results["hybrid_plant"].pv.footprint_area)
    if "battery" in hopp_config["technologies"].keys():
        wind_buffer -= np.sqrt(hopp_results["hybrid_plant"].battery.footprint_area)
    if wind_buffer < 50:
        onshorey += wind_buffer - 50

    if design_scenario["wind_location"] == "offshore":
        origin_x = substation_x
        origin_y = substation_y
    else:
        origin_x = 0.0
        origin_y = 0.0

    ## create figure
    if design_scenario["wind_location"] == "offshore":
        fig, ax = plt.subplots(2, 2, figsize=(10, 6))
        ax_index_plant = (0, 0)
        ax_index_detail = (1, 0)
        ax_index_wind_plant = (0, 1)
        ax_index_turbine_detail = (1, 1)
    else:
        fig, ax = plt.subplots(1, 2, figsize=(10, 6))
        ax_index_plant = 0
        ax_index_wind_plant = 0
        ax_index_detail = 1
        ax_index_turbine_detail = False

    # plot the stuff

    # onshore plant | offshore plant
    # platform/substation | turbine

    ## add turbines
    def add_turbines(ax, turbine_x, turbine_y, radius, color):
        i = 0
        for x, y in zip(turbine_x, turbine_y):
            if i == 0:
                rlabel = "Wind Turbine Rotor"
                tlabel = "Wind Turbine Tower"
                i += 1
            else:
                rlabel = None
                tlabel = None
            turbine_patch = patches.Circle(
                (x, y), radius=radius, color=color, fill=False, label=rlabel, zorder=10,
            )
            ax.add_patch(turbine_patch)

    add_turbines(
        ax[ax_index_wind_plant], turbine_x, turbine_y, rotor_radius, turbine_rotor_color
    )
    component_areas["turbine_area_m2"] = turbine_area
    # turbine_patch01_tower = patches.Circle((x, y), radius=tower_base_radius, color=turbine_tower_color, fill=False, label=tlabel, zorder=10)
    # ax[0, 1].add_patch(turbine_patch01_tower)
    if design_scenario["wind_location"] == "onshore":
        add_turbines(
            ax[ax_index_detail], turbine_x, turbine_y, rotor_radius, turbine_rotor_color
        )

    if ax_index_turbine_detail:
        # turbine_patch11_rotor = patches.Circle((turbine_x[0], turbine_y[0]), radius=rotor_radius, color=turbine_rotor_color, fill=False, label=None, zorder=10)
        tlabel = "Wind Turbine Tower"
        turbine_patch11_tower = patches.Circle(
            (turbine_x[0], turbine_y[0]),
            radius=tower_base_radius,
            color=turbine_tower_color,
            fill=False,
            label=tlabel,
            zorder=10,
        )
        # ax[1, 1].add_patch(turbine_patch11_rotor)
        ax[ax_index_turbine_detail].add_patch(turbine_patch11_tower)

    # add pipe array
    if design_scenario["transportation"] == "hvdc+pipeline" or (
        design_scenario["h2_storage_location"] != "turbine"
        and design_scenario["electrolyzer_location"] == "turbine"
    ):
        i = 0
        for point_string in pipe_array_points:
            if i == 0:
                label = "Array Pipes"
                i += 1
            else:
                label = None
            ax[0, 1].plot(
                point_string[:, 0],
                point_string[:, 1] - substation_side_length / 2,
                ":",
                color=pipe_color,
                zorder=0,
                linewidth=1,
                label=label,
            )
            ax[1, 0].plot(
                point_string[:, 0],
                point_string[:, 1] - substation_side_length / 2,
                ":",
                color=pipe_color,
                zorder=0,
                linewidth=1,
                label=label,
            )
            ax[1, 1].plot(
                point_string[:, 0],
                point_string[:, 1] - substation_side_length / 2,
                ":",
                color=pipe_color,
                zorder=0,
                linewidth=1,
                label=label,
            )

    ## add cables
    if (len(cable_array_points) > 1) and (
        design_scenario["h2_storage_location"] != "turbine"
        or design_scenario["transportation"] == "hvdc+pipeline"
    ):
        i = 0
        for point_string in cable_array_points:
            if i == 0:
                label = "Array Cables"
                i += 1
            else:
                label = None
            ax[0, 1].plot(
                point_string[:, 0],
                point_string[:, 1] + substation_side_length / 2,
                "-",
                color=cable_color,
                zorder=0,
                linewidth=1,
                label=label,
            )
            ax[1, 0].plot(
                point_string[:, 0],
                point_string[:, 1] + substation_side_length / 2,
                "-",
                color=cable_color,
                zorder=0,
                linewidth=1,
                label=label,
            )
            ax[1, 1].plot(
                point_string[:, 0],
                point_string[:, 1] + substation_side_length / 2,
                "-",
                color=cable_color,
                zorder=0,
                linewidth=1,
                label=label,
            )

    ## add offshore substation
    if design_scenario["wind_location"] == "offshore" and (
        design_scenario["h2_storage_location"] != "turbine"
        or design_scenario["transportation"] == "hvdc+pipeline"
    ):
        substation_patch01 = patches.Rectangle(
            (
                substation_x - substation_side_length,
                substation_y - substation_side_length / 2,
            ),
            substation_side_length,
            substation_side_length,
            fill=True,
            color=substation_color,
            label="Substation*",
            zorder=11,
        )
        substation_patch10 = patches.Rectangle(
            (
                substation_x - substation_side_length,
                substation_y - substation_side_length / 2,
            ),
            substation_side_length,
            substation_side_length,
            fill=True,
            color=substation_color,
            label="Substation*",
            zorder=11,
        )
        ax[0, 1].add_patch(substation_patch01)
        ax[1, 0].add_patch(substation_patch10)

        component_areas['offshore_substation_area_m2'] = substation_side_length ** 2

    ## add equipment platform
    if design_scenario["wind_location"] == "offshore" and (
        design_scenario["h2_storage_location"] == "platform"
        or design_scenario["electrolyzer_location"] == "platform"
    ):  # or design_scenario["transportation"] == "pipeline":
        equipment_platform_patch01 = patches.Rectangle(
            (
                equipment_platform_x - equipment_platform_side_length / 2,
                equipment_platform_y - equipment_platform_side_length / 2,
            ),
            equipment_platform_side_length,
            equipment_platform_side_length,
            color=equipment_platform_color,
            fill=True,
            label="Equipment Platform",
            zorder=1,
        )
        equipment_platform_patch10 = patches.Rectangle(
            (
                equipment_platform_x - equipment_platform_side_length / 2,
                equipment_platform_y - equipment_platform_side_length / 2,
            ),
            equipment_platform_side_length,
            equipment_platform_side_length,
            color=equipment_platform_color,
            fill=True,
            label="Equipment Platform",
            zorder=1,
        )
        ax[0, 1].add_patch(equipment_platform_patch01)
        ax[1, 0].add_patch(equipment_platform_patch10)

        component_areas['equipment_platform_area_m2'] = equipment_platform_area

    ## add hvdc cable
    if (
        design_scenario["transportation"] == "hvdc"
        or design_scenario["transportation"] == "hvdc+pipeline"
    ):
        ax[0, 0].plot(
            [onshorex + onshore_substation_x_side_length, 1000],
            [48, 48],
            "--",
            color=cable_color,
            label="HVDC Cable",
        )
        ax[0, 1].plot(
            [-5000, substation_x],
            [substation_y - 100, substation_y - 100],
            "--",
            color=cable_color,
            label="HVDC Cable",
            zorder=0,
        )
        ax[1, 0].plot(
            [-5000, substation_x],
            [substation_y - 2, substation_y - 2],
            "--",
            color=cable_color,
            label="HVDC Cable",
            zorder=0,
        )

    ## add onshore substation
    if (
        design_scenario["transportation"] == "hvdc"
        or design_scenario["transportation"] == "hvdc+pipeline"
    ):
        onshore_substation_patch00 = patches.Rectangle(
            (
                onshorex + 0.2 * onshore_substation_y_side_length,
                onshorey - onshore_substation_y_side_length * 1.2,
            ),
            onshore_substation_x_side_length,
            onshore_substation_y_side_length,
            fill=True,
            color=substation_color,
            label="Substation*",
            zorder=11,
        )
        ax[0, 0].add_patch(onshore_substation_patch00)

        component_areas['onshore_substation_area_m2'] = onshore_substation_area

    ## add transport pipeline
    if design_scenario["transportation"] == "colocated":
        # add hydrogen pipeline to end use
        linetype = "-."
        label = "Pipeline to Storage/End-Use"
        linewidth = 1.0

        ax[ax_index_plant].plot(
            [onshorex, -10000],
            [onshorey, onshorey],
            linetype,
            color=pipe_color,
            label=label,
            linewidth=linewidth,
            zorder=0,
        )

        ax[ax_index_detail].plot(
            [onshorex, -10000],
            [onshorey, onshorey],
            linetype,
            color=pipe_color,
            label=label,
            linewidth=linewidth,
            zorder=0,
        )
    if (
        design_scenario["transportation"] == "pipeline"
        or design_scenario["transportation"] == "hvdc+pipeline"
        or (
            design_scenario["transportation"] == "hvdc"
            and design_scenario["h2_storage_location"] == "platform"
        )
    ):
        linetype = "-."
        label = "Transport Pipeline"
        linewidth = 1.0

        ax[ax_index_plant].plot(
            [onshorex, 1000],
            [onshorey + 2, onshorey + 2],
            linetype,
            color=pipe_color,
            label=label,
            linewidth=linewidth,
            zorder=0,
        )

        if design_scenario["wind_location"] == "offshore":
            ax[ax_index_wind_plant].plot(
                [-5000, substation_x],
                [substation_y + 100, substation_y + 100],
                linetype,
                linewidth=linewidth,
                color=pipe_color,
                label=label,
                zorder=0,
            )
            ax[ax_index_detail].plot(
                [-5000, substation_x],
                [substation_y + 2, substation_y + 2],
                linetype,
                linewidth=linewidth,
                color=pipe_color,
                label=label,
                zorder=0,
            )

            if (
                design_scenario["transportation"] == "hvdc"
                or design_scenario["transportation"] == "hvdc+pipeline"
            ) and design_scenario["h2_storage_location"] == "platform":
                h2cx = onshorex - compressor_side
                h2cy = onshorey - compressor_side + 2
                h2cax = ax[ax_index_plant]
            else:
                h2cx = substation_x - substation_side_length
                h2cy = substation_y
                h2cax = ax[ax_index_detail]

        if design_scenario["wind_location"] == "onshore":
            compressor_patch01 = patches.Rectangle(
                (origin_x, origin_y),
                compressor_side,
                compressor_side,
                color=compressor_color,
                fill=None,
                label="Transport Compressor*",
                hatch="+++",
                zorder=20,
            )
            ax[ax_index_plant].add_patch(compressor_patch01)

        compressor_patch10 = patches.Rectangle(
            (h2cx, h2cy),
            compressor_side,
            compressor_side,
            color=compressor_color,
            fill=None,
            label="Transport Compressor*",
            hatch="+++",
            zorder=20,
        )
        h2cax.add_patch(compressor_patch10)

        component_areas['compressor_area_m2'] = compressor_area

    ## add plant components
    if design_scenario["electrolyzer_location"] == "onshore":
        electrolyzer_x = onshorex
        electrolyzer_y = onshorey
        electrolyzer_patch = patches.Rectangle(
            (electrolyzer_x, electrolyzer_y),
            electrolyzer_side,
            electrolyzer_side,
            color=electrolyzer_color,
            fill=None,
            label="Electrolyzer",
            zorder=20,
            hatch=electrolyzer_hatch,
        )
        ax[ax_index_plant].add_patch(electrolyzer_patch)
        component_areas['electrolyzer_area_m2'] = electrolyzer_area

        if design_scenario["wind_location"] == "onshore":
            electrolyzer_patch = patches.Rectangle(
                (onshorex - h2_storage_side, onshorey + 4),
                electrolyzer_side,
                electrolyzer_side,
                color=electrolyzer_color,
                fill=None,
                label="Electrolyzer",
                zorder=20,
                hatch=electrolyzer_hatch,
            )
            ax[ax_index_detail].add_patch(electrolyzer_patch)

    elif design_scenario["electrolyzer_location"] == "platform":
        dx = equipment_platform_x - equipment_platform_side_length / 2
        dy = equipment_platform_y - equipment_platform_side_length / 2
        e_side_y = equipment_platform_side_length
        e_side_x = electrolyzer_area / e_side_y
        d_side_y = equipment_platform_side_length
        d_side_x = desal_equipment_area / d_side_y
        electrolyzer_x = dx + d_side_x
        electrolyzer_y = dy

        electrolyzer_patch = patches.Rectangle(
            (electrolyzer_x, electrolyzer_y),
            e_side_x,
            e_side_y,
            color=electrolyzer_color,
            fill=None,
            zorder=20,
            label="Electrolyzer",
            hatch=electrolyzer_hatch,
        )
        ax[ax_index_detail].add_patch(electrolyzer_patch)
        desal_patch = patches.Rectangle(
            (dx, dy),
            d_side_x,
            d_side_y,
            color=desal_color,
            zorder=21,
            fill=None,
            label="Desalinator",
            hatch=desalinator_hatch,
        )
        ax[ax_index_detail].add_patch(desal_patch)
        component_areas['desalination_area_m2'] = desal_equipment_area

    elif design_scenario["electrolyzer_location"] == "turbine":
        electrolyzer_patch11 = patches.Rectangle(
            (turbine_x[0], turbine_y[0] + tower_base_radius),
            electrolyzer_side,
            electrolyzer_side,
            color=electrolyzer_color,
            fill=None,
            zorder=20,
            label="Electrolyzer",
            hatch=electrolyzer_hatch,
        )
        ax[ax_index_turbine_detail].add_patch(electrolyzer_patch11)
        desal_patch11 = patches.Rectangle(
            (turbine_x[0] - desal_equipment_side, turbine_y[0] + tower_base_radius),
            desal_equipment_side,
            desal_equipment_side,
            color=desal_color,
            zorder=21,
            fill=None,
            label="Desalinator",
            hatch=desalinator_hatch,
        )
        ax[ax_index_turbine_detail].add_patch(desal_patch11)
        component_areas['desalination_area_m2'] = desal_equipment_area
        i = 0
        for x, y in zip(turbine_x, turbine_y):
            if i == 0:
                elable = "Electrolyzer"
                dlabel = "Desalinator"
            else:
                elable = None
                dlabel = None
            electrolyzer_patch01 = patches.Rectangle(
                (x, y + tower_base_radius),
                electrolyzer_side,
                electrolyzer_side,
                color=electrolyzer_color,
                fill=None,
                zorder=20,
                label=elable,
                hatch=electrolyzer_hatch,
            )
            desal_patch01 = patches.Rectangle(
                (x - desal_equipment_side, y + tower_base_radius),
                desal_equipment_side,
                desal_equipment_side,
                color=desal_color,
                zorder=21,
                fill=None,
                label=dlabel,
                hatch=desalinator_hatch,
            )
            ax[ax_index_wind_plant].add_patch(electrolyzer_patch01)
            ax[ax_index_wind_plant].add_patch(desal_patch01)
            i += 1

    h2_storage_hatch = "\\\\\\"
    if design_scenario["h2_storage_location"] == "onshore" and (
        greenheart_config["h2_storage"]["type"] != "none"
    ):
        h2_storage_patch = patches.Rectangle(
            (onshorex - h2_storage_side, onshorey - h2_storage_side - 2),
            h2_storage_side,
            h2_storage_side,
            color=h2_storage_color,
            fill=None,
            label="H$_2$ Storage",
            hatch=h2_storage_hatch,
        )
        ax[ax_index_plant].add_patch(h2_storage_patch)
        component_areas["h2_storage_area_m2"] = h2_storage_area

        if design_scenario["wind_location"] == "onshore":
            h2_storage_patch = patches.Rectangle(
                (onshorex - h2_storage_side, onshorey - h2_storage_side - 2),
                h2_storage_side,
                h2_storage_side,
                color=h2_storage_color,
                fill=None,
                label="H$_2$ Storage",
                hatch=h2_storage_hatch,
            )
            ax[ax_index_detail].add_patch(h2_storage_patch)
            component_areas["h2_storage_area_m2"] = h2_storage_area
    elif design_scenario["h2_storage_location"] == "platform" and (
        greenheart_config["h2_storage"]["type"] != "none"
    ):
        s_side_y = equipment_platform_side_length
        s_side_x = h2_storage_area / s_side_y
        sx = equipment_platform_x - equipment_platform_side_length / 2
        sy = equipment_platform_y - equipment_platform_side_length / 2
        if design_scenario["electrolyzer_location"] == "platform":
            sx += equipment_platform_side_length - s_side_x

        h2_storage_patch = patches.Rectangle(
            (sx, sy),
            s_side_x,
            s_side_y,
            color=h2_storage_color,
            fill=None,
            label="H$_2$ Storage",
            hatch=h2_storage_hatch,
        )
        ax[ax_index_detail].add_patch(h2_storage_patch)
        component_areas["h2_storage_area_m2"] = h2_storage_area

    elif design_scenario["h2_storage_location"] == "turbine":

        if greenheart_config["h2_storage"]["type"] == "turbine":
            h2_storage_patch = patches.Circle(
                (turbine_x[0], turbine_y[0]),
                radius=tower_base_diameter / 2,
                color=h2_storage_color,
                fill=None,
                label="H$_2$ Storage",
                hatch=h2_storage_hatch,
            )
            ax[ax_index_turbine_detail].add_patch(h2_storage_patch)
            component_areas["h2_storage_area_m2"] = h2_storage_area
            i = 0
            for x, y in zip(turbine_x, turbine_y):
                if i == 0:
                    slable = "H$_2$ Storage"
                else:
                    slable = None
                h2_storage_patch = patches.Circle(
                    (x, y),
                    radius=tower_base_diameter / 2,
                    color=h2_storage_color,
                    fill=None,
                    label=None,
                    hatch=h2_storage_hatch,
                )
                ax[ax_index_wind_plant].add_patch(h2_storage_patch)
        elif greenheart_config["h2_storage"]["type"] == "pressure_vessel":
            h2_storage_side = np.sqrt(
                h2_storage_area / greenheart_config["plant"]["num_turbines"]
            )
            h2_storage_patch = patches.Rectangle(
                (
                    turbine_x[0] - h2_storage_side - desal_equipment_side,
                    turbine_y[0] + tower_base_radius,
                ),
                width=h2_storage_side,
                height=h2_storage_side,
                color=h2_storage_color,
                fill=None,
                label="H$_2$ Storage",
                hatch=h2_storage_hatch,
            )
            ax[ax_index_turbine_detail].add_patch(h2_storage_patch)
            component_areas["h2_storage_area_m2"] = h2_storage_area
            i = 0
            for x, y in zip(turbine_x, turbine_y):
                if i == 0:
                    slable = "H$_2$ Storage"
                else:
                    slable = None
                h2_storage_patch = patches.Rectangle(
                    (
                        turbine_x[i] - h2_storage_side - desal_equipment_side,
                        turbine_y[i] + tower_base_radius,
                    ),
                    width=h2_storage_side,
                    height=h2_storage_side,
                    color=h2_storage_color,
                    fill=None,
                    label=slable,
                    hatch=h2_storage_hatch,
                )
                ax[ax_index_wind_plant].add_patch(h2_storage_patch)
                i += 1

    ## add battery
    if "battery" in hopp_config["technologies"].keys():
        component_areas['battery_area_m2'] = hopp_results["hybrid_plant"].battery.footprint_area
        if design_scenario["battery_location"] == "onshore":
            battery_side_y = np.sqrt(
                hopp_results["hybrid_plant"].battery.footprint_area
            )
            battery_side_x = battery_side_y

            batteryx = electrolyzer_x

            batteryy = electrolyzer_y + electrolyzer_side + 10

            battery_patch = patches.Rectangle(
                (batteryx, batteryy),
                battery_side_x,
                battery_side_y,
                color=battery_color,
                fill=None,
                label="Battery Array",
                hatch=battery_hatch,
            )
            ax[ax_index_plant].add_patch(battery_patch)

            if design_scenario["wind_location"] == "onshore":

                battery_patch = patches.Rectangle(
                    (batteryx, batteryy),
                    battery_side_x,
                    battery_side_y,
                    color=battery_color,
                    fill=None,
                    label="Battery Array",
                    hatch=battery_hatch,
                )
                ax[ax_index_detail].add_patch(battery_patch)

        elif design_scenario["battery_location"] == "platform":
            battery_side_y = equipment_platform_side_length
            battery_side_x = (
                hopp_results["hybrid_plant"].battery.footprint_area / battery_side_y
            )

            batteryx = equipment_platform_x - equipment_platform_side_length / 2
            batteryy = equipment_platform_y - equipment_platform_side_length / 2

            battery_patch = patches.Rectangle(
                (batteryx, batteryy),
                battery_side_x,
                battery_side_y,
                color=battery_color,
                fill=None,
                label="Battery Array",
                hatch=battery_hatch,
            )
            ax[ax_index_detail].add_patch(battery_patch)   

    else:
        battery_side_y = 0.0
        battery_side_x = 0.0   
    
    ## add solar
    if hopp_config["site"]["solar"]:
        component_areas['pv_area_m2'] = hopp_results["hybrid_plant"].pv.footprint_area
        if design_scenario["pv_location"] == "offshore":
            solar_side_y = equipment_platform_side_length
            solar_side_x = hopp_results["hybrid_plant"].pv.footprint_area / solar_side_y

            solarx = equipment_platform_x - equipment_platform_side_length / 2
            solary = equipment_platform_y - equipment_platform_side_length / 2

            solar_patch = patches.Rectangle(
                (solarx, solary),
                solar_side_x,
                solar_side_y,
                color=solar_color,
                fill=None,
                label="Solar Array",
                hatch=solar_hatch,
            )
            ax[ax_index_detail].add_patch(solar_patch)
        else:
            solar_side_y = np.sqrt(hopp_results["hybrid_plant"].pv.footprint_area)
            solar_side_x = hopp_results["hybrid_plant"].pv.footprint_area / solar_side_y

            solarx = electrolyzer_x

            solary = electrolyzer_y + electrolyzer_side + 10

            if "battery" in hopp_config["technologies"].keys():
                solary += battery_side_y + 10

            solar_patch = patches.Rectangle(
                (solarx, solary),
                solar_side_x,
                solar_side_y,
                color=solar_color,
                fill=None,
                label="Solar Array",
                hatch=solar_hatch,
            )

            ax[ax_index_plant].add_patch(solar_patch)

            solar_patch = patches.Rectangle(
                (solarx, solary),
                solar_side_x,
                solar_side_y,
                color=solar_color,
                fill=None,
                label="Solar Array",
                hatch=solar_hatch,
            )

            ax[ax_index_detail].add_patch(solar_patch)
    else:
        solar_side_x = 0.0
        solar_side_y = 0.0

    ## add wave
    if hopp_config["site"]["wave"]:
        # get wave generation area geometry
        num_devices = hopp_config["technologies"]["wave"]["num_devices"]
        distance_to_shore = (
            hopp_config["technologies"]["wave"]["cost_inputs"]["distance_to_shore"]
            * 1e3
        )
        number_rows = hopp_config["technologies"]["wave"]["cost_inputs"]["number_rows"]
        device_spacing = hopp_config["technologies"]["wave"]["cost_inputs"][
            "device_spacing"
        ]
        row_spacing = hopp_config["technologies"]["wave"]["cost_inputs"]["row_spacing"]

        # calculate wave generation area dimenstions
        wave_side_y = device_spacing * np.ceil(num_devices / number_rows)
        wave_side_x = row_spacing * (number_rows)
        wave_area = wave_side_x * wave_side_y
        component_areas['wave_area_m2'] = wave_area

        # generate wave generation patch
        wavex = substation_x - wave_side_x
        wavey = substation_y + distance_to_shore
        wave_patch = patches.Rectangle(
            (wavex, wavey),
            wave_side_x,
            wave_side_y,
            color=wave_color,
            fill=None,
            label="Wave Array",
            hatch=wave_hatch,
            zorder=1,
        )
        ax[ax_index_wind_plant].add_patch(wave_patch)

        # add electrical transmission for wave
        wave_export_cable_coords_x = [substation_x, substation_x]
        wave_export_cable_coords_y = [substation_y, substation_y + distance_to_shore]

        ax[ax_index_wind_plant].plot(
            wave_export_cable_coords_x,
            wave_export_cable_coords_y,
            cable_color,
            zorder=0,
        )
        ax[ax_index_detail].plot(
            wave_export_cable_coords_x,
            wave_export_cable_coords_y,
            cable_color,
            zorder=0,
        )

    if design_scenario["wind_location"] == "offshore":
        allpoints = cable_array_points.flatten()
    else:
        allpoints = turbine_x

    allpoints = allpoints[~np.isnan(allpoints)]

    if design_scenario["wind_location"] == "offshore":
        roundto = -2
        ax[ax_index_plant].set(
            xlim=[
                round(np.min(onshorex - 100), ndigits=roundto),
                round(
                    np.max(
                        onshorex
                        + onshore_substation_x_side_length
                        + electrolyzer_side
                        + 200
                    ),
                    ndigits=roundto,
                ),
            ],
            ylim=[
                round(np.min(onshorey - 100), ndigits=roundto),
                round(
                    np.max(
                        onshorey
                        + battery_side_y
                        + electrolyzer_side
                        + solar_side_y
                        + 100
                    ),
                    ndigits=roundto,
                ),
            ],
        )
        ax[ax_index_plant].set(aspect="equal")
    else:
        roundto = -3
        ax[ax_index_plant].set(
            xlim=[
                round(np.min(allpoints - 6000), ndigits=roundto),
                round(np.max(allpoints + 6000), ndigits=roundto),
            ],
            ylim=[
                round(np.min(onshorey - 1000), ndigits=roundto),
                round(np.max(turbine_y + 4000), ndigits=roundto),
            ],
        )
        ax[ax_index_plant].autoscale()
        ax[ax_index_plant].set(aspect="equal")
        ax[ax_index_plant].xaxis.set_major_locator(ticker.MultipleLocator(2000))
        ax[ax_index_plant].yaxis.set_major_locator(ticker.MultipleLocator(1000))

    roundto = -3
    ax[ax_index_wind_plant].set(
        xlim=[
            round(np.min(allpoints - 6000), ndigits=roundto),
            round(np.max(allpoints + 6000), ndigits=roundto),
        ],
        ylim=[
            round((np.min([np.min(turbine_y), onshorey]) - 1000), ndigits=roundto),
            round(np.max(turbine_y + 4000), ndigits=roundto),
        ],
    )
    # ax[ax_index_wind_plant].autoscale()
    ax[ax_index_wind_plant].set(aspect="equal")
    ax[ax_index_wind_plant].xaxis.set_major_locator(ticker.MultipleLocator(5000))
    ax[ax_index_wind_plant].yaxis.set_major_locator(ticker.MultipleLocator(1000))

    if design_scenario["wind_location"] == "offshore":
        roundto = -2
        ax[ax_index_detail].set(
            xlim=[
                round(origin_x - 400, ndigits=roundto),
                round(origin_x + 100, ndigits=roundto),
            ],
            ylim=[
                round(origin_y - 200, ndigits=roundto),
                round(origin_y + 200, ndigits=roundto),
            ],
        )
        ax[ax_index_detail].set(aspect="equal")
    else:
        roundto = -2

        if "pv" in hopp_config["technologies"].keys():
            xmax = round(
                np.max([onshorex + 510, solarx + solar_side_x + 100]), ndigits=roundto
            )
            ymax = round(solary + solar_side_y + 100, ndigits=roundto)
        else:
            xmax = round(np.max([onshorex + 510, 100]), ndigits=roundto)
            ymax = round(100, ndigits=roundto)
        ax[ax_index_detail].set(
            xlim=[round(onshorex - 10, ndigits=roundto), xmax,],
            ylim=[round(onshorey - 200, ndigits=roundto), ymax,],
        )
        ax[ax_index_detail].set(aspect="equal")

    if design_scenario["wind_location"] == "offshore":
        tower_buffer0 = 10
        tower_buffer1 = 10
        roundto = -1
        ax[ax_index_turbine_detail].set(
            xlim=[
                round(
                    turbine_x[0] - tower_base_radius - tower_buffer0 - 50,
                    ndigits=roundto,
                ),
                round(
                    turbine_x[0] + tower_base_radius + 3 * tower_buffer1,
                    ndigits=roundto,
                ),
            ],
            ylim=[
                round(
                    turbine_y[0] - tower_base_radius - 2 * tower_buffer0,
                    ndigits=roundto,
                ),
                round(
                    turbine_y[0] + tower_base_radius + 4 * tower_buffer1,
                    ndigits=roundto,
                ),
            ],
        )
        ax[ax_index_turbine_detail].set(aspect="equal")
        ax[ax_index_turbine_detail].xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax[ax_index_turbine_detail].yaxis.set_major_locator(ticker.MultipleLocator(10))
        # ax[0,1].legend(frameon=False)
        # ax[0,1].axis('off')

    if design_scenario["wind_location"] == "offshore":
        labels = [
            "(a) Onshore plant",
            "(b) Offshore plant",
            "(c) Equipment platform and substation",
            "(d) NW-most wind turbine",
        ]
    else:
        labels = ["(a) Full plant", "(b) Non-wind plant detail"]
    for axi, label in zip(ax.flatten(), labels):
        axi.legend(frameon=False, ncol=2)  # , ncol=2, loc="best")
        axi.set(xlabel="Easting (m)", ylabel="Northing (m)")
        axi.set_title(label, loc="left")
        # axi.spines[['right', 'top']].set_visible(False)

    ## save the plot
    plt.tight_layout()
    savepaths = [
            output_dir + "figures/layout/",
            output_dir + "data/",
        ]
    if save_plots:
        for savepath in savepaths:
            if not os.path.exists(savepath):
                os.makedirs(savepath)
        plt.savefig(
            savepaths[0] + "plant_layout_%i.png" % (plant_design_number), transparent=True
        )
        
        df = pd.DataFrame([component_areas])
        df.to_csv(savepaths[1] + "component_areas_layout_%i.csv" % (plant_design_number), index=False)

    if show_plots:
        plt.show()
    return 0


def save_energy_flows(
    hybrid_plant: HoppInterface.system, 
    electrolyzer_physics_results, 
    solver_results, 
    hours, 
    h2_storage_results,
    ax=None, 
    simulation_length=8760, 
    output_dir="./output/",
):

    

    if ax == None:
        fig, ax = plt.subplots(1)

    output = {}
    if hybrid_plant.pv:
        solar_plant_power = np.array(
            hybrid_plant.pv.generation_profile[0:simulation_length]
        )
        output.update({"pv generation [kW]": solar_plant_power})
    if hybrid_plant.wind:
        wind_plant_power = np.array(
            hybrid_plant.wind.generation_profile[0:simulation_length]
        )
        output.update({"wind generation [kW]": wind_plant_power})
    if hybrid_plant.wave:
        wave_plant_power = np.array(
            hybrid_plant.wave.generation_profile[0:simulation_length]
        )
        output.update({"wave generation [kW]": wave_plant_power})
    if hybrid_plant.battery:
        battery_power_out_mw = hybrid_plant.battery.outputs.P 
        output.update({"battery discharge [kW]": [(int(p>0))*p*1E3 for p in battery_power_out_mw]}) # convert from MW to kW and extract only discharging
        output.update({"battery charge [kW]": [-(int(p<0))*p*1E3 for p in battery_power_out_mw]}) # convert from MW to kW and extract only charging
        output.update({"battery state of charge [%]": hybrid_plant.battery.outputs.dispatch_SOC})

    output.update({"total renewable energy production hourly [kW]": [solver_results[0]]*simulation_length})
    output.update({"grid energy usage hourly [kW]": [solver_results[1]]*simulation_length})
    output.update({"desal energy hourly [kW]": [solver_results[2]]*simulation_length})
    output.update({"electrolyzer energy hourly [kW]": electrolyzer_physics_results["power_to_electrolyzer_kw"]})
    output.update({"transport compressor energy hourly [kW]": [solver_results[3]]*simulation_length})
    output.update({"storage energy hourly [kW]": [solver_results[4]]*simulation_length})
    output.update({"h2 production hourly [kg]": electrolyzer_physics_results["H2_Results"]["Hydrogen Hourly Production [kg/hr]"]})
    if "hydrogen_storage_soc" in h2_storage_results:
        output.update({"hydrogen storage SOC [kg]": h2_storage_results["hydrogen_storage_soc"]})
    
    df = pd.DataFrame.from_dict(output)

    filepath = os.path.abspath(output_dir + "data/production/")

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    df.to_csv(os.path.join(filepath, "energy_flows.csv"))

    return output

###NOTE / TODO:
    # Add all config for LCA to greenheart_config yaml as input
    # update greenheart_simulation.py 
        # call of post_process_simulation() with flags for run_lca()
    # update post_process_simulation() to call run_lca()
#QUESTIONS:
    # if cambium_year does not align with 2025:2055:5 intervals should we extrapolate back? ie cambium_year = 2023, extrapolate for 2023 and 2024?
    # Are capex_EI values over stated (capex EI calculations for every year, do those value include startup / construction capex?)
#NOTE: Suedocode:
    # x 1. define conversions
    # x 2. pull greet values (start- pull all, optimal- based on hopp/system config)
    # x 3. define variables to each years LCA calculation data
    # x 4. logic to convert atb_year to cambium_year (+5 yrs)
    # x 5. define lists to hold data for all LCA calculations / cambium years
    # x 6. read in hopp data as df (energy to electrolyzer kwh, energy from grid kwh, energy from renewables kwh, total energy kwh)
    # x 7. loop through cambium files, read in data, concat with hopp data, perform calculations based on grid case, append data to lists from 5
    # x 8. after looping through all cambium files, create dataframe with lists from 5/7 (emissions_intensities_df)
    # 9. calculate endoflife_year = cambium_year + system_life
    # 10. define lists for interpolated data
    # 11. loop through each year between cambium_year and endoflife_year
        # if year <= max year interpolate values for that year
        # else append last value of emissions_intensities_df (from 8) to interpolated lists
    # 12. sum interpolated lists * annual h2 prod sum / h2prod_life_sum to calculate each _LCA value
    # 13. put all cumulative metrics into dictionary and then dataframe 
    # 14. save as csv
    # grid_case put in greenheart_config as input

def run_lca(
    hopp_results,
    electrolyzer_physics_results,
    greenheart_config,
    hopp_config,
    ):

    # Load HOPP Data:  
    project_lifetime = greenheart_config['project_parameters']['project_lifetime']                                  # system lifetime (years)
    wind_annual_energy_kwh = hopp_results['annual_energies']['wind']                                                # annual energy from wind (kWh)
    solar_pv_annual_energy_kwh = hopp_results['annual_energies']['pv']                                              # annual energy from solar (kWh)
    battery_system_capacity_kwh = hopp_results['hybrid_plant'].battery.system_capacity_kwh                          # battery rated capacity (kWh)
    energy_to_electrolyzer = np.array(electrolyzer_physics_results['power_to_electrolyzer_kw'])                     # total power to the electrolyzer (kW*1hr = kWh)
    energy_to_electrolyzer_from_renewables = np.array(hopp_results['combined_hybrid_power_production_hopp'])        # power to the electrolyzer from renewable sources (kW*1hr = kWh)
    energy_to_electrolyzer_from_grid = energy_to_electrolyzer - energy_to_electrolyzer_from_renewables              # power to the electrolyzer from grid (kW*1hr = kWh)
    # Sets energy from grid to zero if negative (some cases where energy to electrolyzer = 0 but energy from renewables > 0, other cases where floating point precision results in negative value from grid)
    energy_to_electrolyzer_from_grid = np.maximum(energy_to_electrolyzer_from_grid,0)
    h2_hourly_prod_kg = np.array(electrolyzer_physics_results['H2_Results']['Hydrogen Hourly Production [kg/hr]'])  # Hourly H2 production (kg/hr * 1hr = kg)

    # Create dataframe for electrolyzer power profiles
    electrolyzer_profiles_data_dict = {'Energy to electrolyzer (kWh)': energy_to_electrolyzer,
                                       'Energy from grid (kWh)': energy_to_electrolyzer_from_grid,
                                       'Energy from renewables (kWh)': energy_to_electrolyzer_from_renewables,
                                       'Hydrogen Hourly production (kg)': h2_hourly_prod_kg}
    electrolyzer_profiles_df = pd.DataFrame(data=electrolyzer_profiles_data_dict)
    electrolyzer_profiles_df = electrolyzer_profiles_df.reset_index().rename(columns={'index':'Interval'})
    electrolyzer_profiles_df['Interval'] = electrolyzer_profiles_df['Interval']+1
    electrolyzer_profiles_df = electrolyzer_profiles_df.set_index('Interval')

    # Instantiate object to hold EI values per year
    smr_Scope3_EI = 'NA'
    smr_Scope2_EI = 'NA'
    smr_Scope1_EI = 'NA'
    smr_total_EI  = 'NA'
    smr_ccs_Scope3_EI = 'NA'
    smr_ccs_Scope2_EI = 'NA'
    smr_ccs_Scope1_EI = 'NA'
    smr_ccs_total_EI  = 'NA'
    NH3_smr_Scope3_EI = 'NA'
    NH3_smr_Scope2_EI = 'NA'
    NH3_smr_Scope1_EI = 'NA'
    NH3_smr_total_EI  = 'NA'
    NH3_smr_ccs_Scope3_EI = 'NA'
    NH3_smr_ccs_Scope2_EI = 'NA'
    NH3_smr_ccs_Scope1_EI = 'NA'
    NH3_smr_ccs_total_EI  = 'NA'
    steel_smr_Scope3_EI = 'NA'
    steel_smr_Scope2_EI = 'NA'
    steel_smr_Scope1_EI = 'NA'
    steel_smr_total_EI  = 'NA'
    steel_smr_ccs_Scope3_EI = 'NA'
    steel_smr_ccs_Scope2_EI = 'NA'
    steel_smr_ccs_Scope1_EI = 'NA'
    steel_smr_ccs_total_EI  = 'NA'
    electrolysis_Scope3_EI = 'NA'
    electrolysis_Scope2_EI = 'NA'
    electrolysis_Scope1_EI = 'NA'
    electrolysis_total_EI  = 'NA'
    NH3_electrolysis_Scope3_EI = 'NA'
    NH3_electrolysis_Scope2_EI = 'NA'
    NH3_electrolysis_Scope1_EI = 'NA'
    NH3_electrolysis_total_EI  = 'NA'
    steel_electrolysis_Scope3_EI = 'NA'
    steel_electrolysis_Scope2_EI = 'NA'
    steel_electrolysis_Scope1_EI = 'NA'
    steel_electrolysis_total_EI  = 'NA'

    # Instantiate lists to hold data for all LCA calculations for all cambium years
    electrolysis_emission_intensity = []
    electrolysis_Scope3_emission_intensity = []
    electrolysis_Scope2_emission_intensity = []
    smr_Scope3_emission_intensity = []
    smr_Scope2_emission_intensity = []
    smr_emission_intensity = []
    smr_ccs_Scope3_emission_intensity = []
    smr_ccs_Scope2_emission_intensity = []
    smr_ccs_emission_intensity = []
    NH3_electrolysis_Scope3_emission_intensity = []
    NH3_electrolysis_Scope2_emission_intensity = []
    NH3_electrolysis_emission_intensity = []
    steel_electrolysis_Scope3_emission_intensity = []
    steel_electrolysis_Scope2_emission_intensity = []
    steel_electrolysis_emission_intensity = []
    NH3_smr_Scope3_emission_intensity = []
    NH3_smr_Scope2_emission_intensity = []
    NH3_smr_emission_intensity = []
    steel_smr_Scope3_emission_intensity = []
    steel_smr_Scope2_emission_intensity = []
    steel_smr_emission_intensity = []
    NH3_smr_ccs_Scope3_emission_intensity = []
    NH3_smr_ccs_Scope2_emission_intensity = []
    NH3_smr_ccs_emission_intensity = []
    steel_smr_ccs_Scope3_emission_intensity = []
    steel_smr_ccs_Scope2_emission_intensity = []
    steel_smr_ccs_emission_intensity = []
    
    ## GREET Data
    # Define conversions
    g_to_kg  = 0.001            # 1 g = 0.001 kg
    MT_to_kg = 1000             # 1 metric tonne = 1000 kg
    kWh_to_MWh = 0.001          # 1 kWh = 0.001 MWh
    MWh_to_kWh = 1000           # 1 MWh = 1000 kWh
    gal_H2O_to_MT = 0.00378541  # 1 US gallon of H2O = 0.00378541 metric tonnes (1 gal = 3.78541 liters, 1 liter H2O = 1 kg, 1000 kg = 1 metric tonne)

    # Instantiate GreetData class object, parse greet if not already parsed, return class object and load data dictionary
    greet_data = GreetData(greet_year=2023)
    greet_data_dict = greet_data.data

    #------------------------------------------------------------------------------
    # Renewable infrastructure embedded emission intensities
    #------------------------------------------------------------------------------
    #TODO: add electrolyzer type (PEM, alkaline, SOEC) to lca_config yaml
    if greenheart_config['lca_config']['electrolyzer_type'] == 'PEM':
        ely_stack_capex_EI = greet_data_dict['pem_ely_stack_capex_EI']                                                  # PEM electrolyzer CAPEX emissions (kg CO2e/kg H2)
        ely_stack_and_BoP_capex_EI = greet_data_dict['pem_ely_stack_and_BoP_capex_EI']                                  # PEM electrolyzer stack CAPEX + Balance of Plant emissions (kg CO2e/kg H2)
    elif greenheart_config['lca_config']['electrolyzer_type'] == 'Alkaline':
        ely_stack_capex_EI = greet_data_dict['alk_ely_stack_capex_EI']                                                  # Alkaline electrolyzer CAPEX emissions (kg CO2e/kg H2)
        ely_stack_and_BoP_capex_EI = greet_data_dict['alk_ely_stack_and_BoP_capex_EI']                                  # Alkaline electrolyzer stack CAPEX + Balance of Plant emissions (kg CO2e/kg H2)
    elif greenheart_config['lca_config']['electrolyzer_type'] == 'SOEC':
        ely_stack_capex_EI = greet_data_dict['soec_ely_stack_capex_EI']                                                 # SOEC electrolyzer CAPEX emissions (kg CO2e/kg H2)
        ely_stack_and_BoP_capex_EI = greet_data_dict['soec_ely_stack_and_BoP_capex_EI']                                 # SOEC electrolyzer stack CAPEX + Balance of Plant emissions (kg CO2e/kg H2)
    wind_capex_EI = greet_data_dict['wind_capex_EI']                                                                    # Wind CAPEX emissions (g CO2e/kWh)
    solar_pv_capex_EI = greet_data_dict['solar_pv_capex_EI']                                                            # Solar PV CAPEX emissions (g CO2e/kWh)
    #TODO: add battery / PV install type (residential vs commercial) to lca_config yaml
    if greenheart_config['lca_config']['battery_install_type'] == 'Residential':
        battery_EI = greet_data_dict['battery_LFP_residential_EI'] * (project_lifetime/battery_system_capacity_kwh)     # Battery embodied emissions for residential solar PV applications (g CO2e/kWh)
    elif greenheart_config['lca_config']['battery_install_type'] == 'Commercial':
        battery_EI = greet_data_dict['battery_LFP_commercial_EI'] * (project_lifetime/battery_system_capacity_kwh)      # Battery embodied emissions for commercial solar PV applications (g CO2e/kWh)
    #TODO: add nuclear type (BWR vs PWR) to lca_config yaml
    if greenheart_config['lca_config']['nuclear_type'] == 'BWR':
        nuclear_capex_EI = greet_data_dict['nuclear_BWR_capex_EI']                                                      # Nuclear Boiling Water Reactor (BWR) CAPEX emissions (g CO2e/kWh)
    elif greenheart_config['lca_config']['nuclear_type'] == 'PWR':
        nuclear_capex_EI = greet_data_dict['nuclear_PWR_capex_EI']                                                      # Nuclear Pressurized Water Reactor (PWR) CAPEX emissions (g CO2e/kWh)
    coal_capex_EI = greet_data_dict['coal_capex_EI']                                                                    # Coal CAPEX emissions (g CO2e/kWh)
    gas_capex_EI = greet_data_dict['gas_capex_EI']                                                                      # Natural Gas Combined Cycle (NGCC) CAPEX emissions (g CO2e/kWh)
    hydro_capex_EI = greet_data_dict['hydro_capex_EI']                                                                  # Hydro CAPEX emissions (g CO2e/kWh)
    bio_capex_EI = greet_data_dict['bio_capex_EI']                                                                      # Biomass CAPEX emissions (g CO2e/kWh)
    #TODO: add geothermal type (EGS, binary, flash) to lca_config yaml
    if greenheart_config['lca_config']['geothermal_type'] == 'EGS':
        geothermal_capex_EI = greet_data_dict['geothermal_EGS_capex_EI']                                                # Geothermal EGS CAPEX emissions (g CO2e/kWh)
    elif greenheart_config['lca_config']['geothermal_type'] == 'Binary':
        geothermal_capex_EI = greet_data_dict['geothermal_binary_capex_EI']                                             # Geothermal Binary CAPEX emissions (g CO2e/kWh)
    elif greenheart_config['lca_config']['geothermal_type'] == 'Flash':
        geothermal_capex_EI = greet_data_dict['geothermal_flash_capex_EI']                                              # Geothermal Flash CAPEX emissions (g CO2e/kWh)

    #------------------------------------------------------------------------------
    # Steam methane reforming (SMR) and Autothermal Reforming (ATR) - Incumbent H2 production processes
    #------------------------------------------------------------------------------
    #TODO: add ATR vs SMR to lca_config and add logic to set properly
    smr_HEX_eff = greet_data_dict['smr_HEX_eff']                                                                        # SMR Heat exchange efficiency (%)
    smr_NG_supply = greet_data_dict['smr_NG_supply']                                                                    # Natural gas extraction and supply to SMR plant assuming 2% CH4 leakage rate (g CO2e/MJ)
    smr_NG_combust =  greet_data_dict['smr_NG_combust']                                                                 # SMR Natural gas combustion emission factor (g CO2e/MJ)
    ccs_PO_consume = greet_data_dict['ccs_PO_consume']                                                                  # Power consumption for CCS (kWh/kg CO2)
    smr_ccs_steam_prod = greet_data_dict['smr_ccs_steam_prod']                                                          # SMR Steam exported w/ CCS (MJ/kg H2)
    smr_ccs_perc_capture = greet_data_dict['smr_ccs_perc_capture']                                                      # CCS rate for SMR (%)
    atr_ccs_steam_prod = greet_data_dict['atr_ccs_steam_prod']                                                          # ATR Steam exported w/CCS (MJ/kg H2)
    atr_ccs_perc_capture = greet_data_dict['atr_ccs_perc_capture']                                                      # CCS rate for Autothermal Reforming (%)
    smr_steam_prod = greet_data_dict['smr_steam_prod']                                                                  # SMR Steam exported w/out CCS (MJ/kg H2)
    atr_steam_prod = greet_data_dict['atr_steam_prod']                                                                  # ATR Steam exported w/out CCS (MJ/kg H2)
    if greenheart_config['lca_config']['H2_prod_NG_type'] == 'NG':
        smr_ccs_NG_consume = greet_data_dict['smr_ccs_NG_consume']                                                      # SMR w/ CCS Well to Gate (WTG) Natural Gas (NG) consumption (MJ-LHV/kg H2)
        smr_ccs_PO_consume = greet_data_dict['smr_ccs_NG_PO_consume']                                                   # SMR via NG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        atr_ccs_NG_consume = greet_data_dict['atr_ccs_NG_consume']                                                      # ATR w/ CCS WTG NG consumption (MJ-LHV/kg H2)
        atr_ccs_PO_consume = greet_data_dict['atr_ccs_NG_PO_consume']                                                   # ATR via NG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        smr_NG_consume = greet_data_dict['smr_NG_consume']                                                              # SMR w/out CCS WTG NG consumption (MJ-HHV/kg H2)
        smr_PO_consume = greet_data_dict['smr_NG_PO_consume']                                                           # SMR via NG w/out CCS WTG Total Energy consumption (kWh/kg H2)
        atr_NG_consume = greet_data_dict['atr_NG_consume']                                                              # ATR w/out CCS WTG NG consumption (MJ-HHV/kg H2)
        atr_PO_consume = greet_data_dict['atr_NG_PO_consume']                                                           # ATR via NG w/out CCS WTG Total Energy consumption (kWh/kg H2)
    elif greenheart_config['lca_config']['H2_prod_NG_type'] == 'RNG':
        smr_ccs_NG_consume = greet_data_dict['smr_ccs_RNG_consume']                                                     # SMR w/ CCS WTG Renewable Natural Gas (RNG) consumption (MJ-LHV/kg H2)
        smr_ccs_PO_consume = greet_data_dict['smr_ccs_RNG_PO_consume']                                                  # SMR via RNG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        atr_ccs_NG_consume = greet_data_dict['atr_ccs_RNG_consume']                                                     # ATR w/ CCS WTG RNG consumption (MJ-LHV/kg H2)
        atr_ccs_PO_consume = greet_data_dict['atr_ccs_RNG_PO_consume']                                                  # ATR via RNG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        smr_NG_consume = greet_data_dict['smr_RNG_consume']                                                             # SMR w/out CCS WTG RNG consumption (MJ-HHV/kg H2)
        smr_PO_consume = greet_data_dict['smr_RNG_PO_consume']                                                          # SMR via RNG w/out CCS WTG Total Energy consumption (kWh/kg H2)
        atr_NG_consume = greet_data_dict['atr_RNG_consume']                                                             # ATR w/out CCS WTG RNG consumption (MJ-HHV/kg H2)
        atr_PO_consume = greet_data_dict['atr_RNG_PO_consume']                                                          # ATR via RNG w/out CCS WTG Total Energy consumption (kWh/kg H2)

    #------------------------------------------------------------------------------
    # Hydrogen production via water electrolysis
    #------------------------------------------------------------------------------
    grid_trans_losses = greet_data_dict['grid_trans_losses']                                                            # Grid losses of 5% are assumed (-)
    fuel_to_grid_curr = greet_data_dict['fuel_to_grid_curr']                                                            # Fuel mix emission intensity for current power grid (g CO2e/kWh)
    fuel_to_grid_futu = greet_data_dict['fuel_to_grid_futu']                                                            # Fuel mix emission intensity for future power grid (g CO2e/kWh)
    #NOTE: GREET 2023 does not provide power consumption for SOEC or Alkaline electrolyzers
    #TODO: Add SOEC and Alkaline power consumption values to greet_data.py (hardcoded, from GREET, or from other model)
    ely_PO_consume = greet_data_dict['pem_ely_PO_consume']                                                              # PEM Electrolysis power consumption per kg h2 (kWh/kg H2)

    #------------------------------------------------------------------------------
    # Ammonia (NH3)
    #------------------------------------------------------------------------------
    #TODO: add green vs blue vs conventional NH3 production type in lca_config
    NH3_boiler_EI = greet_data_dict['NH3_boiler_EI']                                                                    # Boiler combustion of methane for Ammonia (kg CO2e/kg NH3)
    # Values for green NH3
    if greenheart_config['lca_config']['NH3_type'] == 'Green':
        NH3_H2_consume = greet_data_dict['green_NH3_H2_consume']                                                        # Green Ammonia production Hydrogen consumption (kg H2/kg NH3)
        NH3_PO_consume = greet_data_dict['green_NH3_PO_consume']                                                        # Green Ammonia production Total Energy consumption (kWh/kg NH3)
    
    # Values for blue NH3
    elif greenheart_config['lca_config']['NH3_type'] == 'Blue':
        NH3_H2_consume = greet_data_dict['blue_NH3_H2_consume']                                                         # Blue Ammonia production Hydrogen consumption (kg H2/kg NH3)
        NH3_PO_consume = greet_data_dict['blue_NH3_PO_consume']                                                         # Blue Ammonia production Total Energy consumption (kWh/kg NH3)

    # Values for conventional NH3
    elif greenheart_config['lca_config']['NH3_type'] == 'Conventional':
        NH3_H2_consume = greet_data_dict['conventional_NH3_H2_consume']                                                 # Conventional Ammonia production Hydrogen consumption (kg H2/kg NH3)
        NH3_PO_consume = greet_data_dict['conventional_NH3_PO_consume']                                                 # Conventional Ammonia production Total Energy consumption (kWh/kg NH3)

    #------------------------------------------------------------------------------
    # Steel
    #------------------------------------------------------------------------------
    # Values agnostic of DRI-EAF config
    steel_H2O_EI = greet_data_dict['steel_H2O_EI']                                                                      # Water consumption emissions for use in DRI-EAF Steel production (kg CO2e/gal H20) 
    steel_NG_supply_EI = greet_data_dict['steel_NG_supply_EI']                                                          # Upstream Natural Gas emissions for DRI-EAF Steel production (g CO2e/MJ)
    steel_iron_ore_EI = greet_data_dict['steel_iron_ore_EI']                                                            # Iron ore production emissions for use in DRI-EAF Steel production (kg CO2e/kg iron ore)
    steel_iron_ore_consume = greet_data_dict['steel_iron_ore_consume']                                                  # Iron ore consumption for EAF and LRF Steel production from DRI (metric tonne iron ore/metric tonne steel production)
    steel_lime_EI = greet_data_dict['steel_lime_EI']                                                                    # Lime production emissions for use in DRI-EAF Steel production (kg CO2e/kg lime)
    steel_lime_consume = greet_data_dict['steel_lime_consume']                                                          # Lime consumption for EAF and LRF Steel production from DRI (metric tonne lime/metric tonne steel production)

    # TODO: add DRI-EAF configuration to lca_config and add logic to set properly if desired
    steel_CH4_prod = greet_data_dict['steel_CH4_prod']                                                                  # CH4 emissions for DRI-EAF Steel production w/ 83% H2 and 0% scrap (kg CO2e/metric tonne annual steel lab production)
    steel_CO2_prod = greet_data_dict['steel_CO2_prod']                                                                  # CO2 emissions for DRI-EAF Steel production w/ 83% H2 and 0% scrap (kg CO2e/metric tonne annual steel lab production)
    steel_H2O_consume = greet_data_dict['steel_H2O_consume']                                                            # H2O consumption for DRI-EAF Steel production w/ 83% H2 and 0% scrap (metric tonne H2O/metric tonne steel production)
    steel_H2_consume = greet_data_dict['steel_H2_consume']                                                              # Hydrogen consumption for DRI-EAF Steel production w/ 83% H2 regardless of scrap (metric tonnes H2/metric tonne steel production)
    steel_NG_consume = greet_data_dict['steel_NG_consume']                                                              # Natural gas consumption for DRI-EAF Steel production (GJ/ton steel)
    steel_PO_consume = greet_data_dict['steel_PO_consume']                                                              # Total Energy consumption for DRI-EAF Steel production w/ 83% H2 and 0% scrap (MWh/metric tonne steel production)
    
    ## Cambium
    # Define cambium_year
    cambium_year = (greenheart_config['project_parameters']['atb_year'] + 5)            # NOTE: current hopp logic for LCOH = atb_year + 2yr + install_period(3yrs) = 5 years
    # Pull / download cambium data files
    cambium_data = CambiumData(lat = hopp_config["site"]["data"]["lat"],
                               lon = hopp_config["site"]["data"]["lon"],
                               year = cambium_year,
                               project_uuid = greenheart_config["cambium"]["project_uuid"],
                               scenario = greenheart_config["cambium"]["scenario"],
                               location_type = greenheart_config["cambium"]["location_type"],
                               time_type = greenheart_config["cambium"]["time_type"],
                               )

    # Read in Cambium data and combine with hopp electrolyzer data
    for resource_file in cambium_data.resource_files:
        cambium_data_df = pd.read_csv(resource_file,
                                      index_col= None,
                                      header = 0, 
                                      usecols = ['lrmer_co2_c','lrmer_ch4_c','lrmer_n2o_c','lrmer_co2_p','lrmer_ch4_p','lrmer_n2o_p','lrmer_co2e_c','lrmer_co2e_p','lrmer_co2e',\
                                                 'generation','battery_MWh','biomass_MWh','beccs_MWh','canada_MWh','coal_MWh','coal-ccs_MWh','csp_MWh','distpv_MWh',\
                                                 'gas-cc_MWh','gas-cc-ccs_MWh','gas-ct_MWh','geothermal_MWh','hydro_MWh','nuclear_MWh','o-g-s_MWh','phs_MWh,upv_MWh','wind-ons_MWh','wind-ofs_MWh']
                                    )
        cambium_data_df = cambium_data_df.reset_index().rename(columns = {'index':'Interval',
                                                                          'lrmer_co2_c':'LRMER CO2 combustion (kg-CO2/MWh)','lrmer_ch4_c':'LRMER CH4 combustion (g-CH4/MWh)',
                                                                          'lrmer_n2o_c':'LRMER N2O combustion (g-N2O/MWh)','lrmer_co2_p':'LRMER CO2 production (kg-CO2/MWh)',
                                                                          'lrmer_ch4_p':'LRMER CH4 production (g-CH4/MWh)','lrmer_n2o_p':'LRMER N2O production (g-N2O/MWh)',
                                                                          'lrmer_co2e_c':'LRMER CO2 equiv. combustion (kg-CO2e/MWh)','lrmer_co2e_p':'LRMER CO2 equiv. production (kg-CO2e/MWh)',
                                                                          'lrmer_co2e':'LRMER CO2 equiv. total (kg-CO2e/MWh)'})
        cambium_data_df['Interval'] = cambium_data_df['Interval']+1
        cambium_data_df = cambium_data_df.set_index('Interval')

        combined_data_df = pd.concat([electrolyzer_profiles_df, cambium_data_df], axis=1)

        # Calculate hourly grid emissions factors (kg CO2e)
        combined_data_df['Total grid emissions (kg-CO2e)'] = (combined_data_df['Energy from grid (kWh)'] / 1000) * combined_data_df['LRMER CO2 equiv. total (kg-CO2e/MWh)']
        combined_data_df['Scope 2 (combustion) grid emissions (kg-CO2e)'] = (combined_data_df['Energy from grid (kWh)'] / 1000) * combined_data_df['LRMER CO2 equiv. combustion (kg-CO2e/MWh)']
        combined_data_df['Scope 3 (production) grid emissions (kg-CO2e)'] = (combined_data_df['Energy from grid (kWh)'] / 1000) * combined_data_df['LRMER CO2 equiv. production (kg-CO2e/MWh)']

        # Calculate annual and lifetime Hydrogen production (kg H2)
        h2prod_annual_sum = combined_data['Hydrogen Hourly production (kg)'].sum()
        h2prod_life_sum = combined_data['Hydrogen Hourly production (kg)'].sum() * project_lifetime

        # Sum total grid emissions
        total_grid_emissions_annual_sum = combined_data_df['Total grid emissions (kg-CO2e)'].sum()                      # (kg CO2e)
        scope2_grid_emissions_annual_sum = combined_data_df['Scope 2 (combustion) grid emissions (kg-CO2e)'].sum()      # (kg CO2e)
        scope3_grid_emissions_annual_sum = combined_data_df['Scope 3 (production) grid emissions (kg-CO2e)'].sum()      # (kg CO2e)
        ren_annual_sum_MWh= combined_data_df['Energy from renewables (kWh)'].sum() / 1000                               # (MWh)
        grid_annual_sum_MWh = combined_data_df['Energy from grid (kWh)'].sum() / 1000                                   # (MWh)
        grid_emission_intensity_annual_average = combined_data_df['LRMER CO2 equiv. total (kg-CO2e/MWh)'].mean()        # (kg CO2e/MWh)

        # Calculate annual percentages of solar, wind and fossil in grid mix (%)
        generation_annual_total_MWh = cambium_data_df['generation'].sum()
        generation_annual_nuclear_fraction = cambium_data_df['nuclear_MWh'].sum() / generation_annual_total_MWh
        generation_annual_coal_oil_fraction = (cambium_data_df['coal_MWh'].sum() + cambium_data_df['coal-ccs_MWh'].sum() + cambium_data_df['o-g-s_MWh'].sum()) / generation_annual_total_MWh
        generation_annual_gas_fraction = (cambium_data_df['gas-cc_MWh'].sum() + cambium_data_df['gas-cc-ccs_MWh'].sum() + cambium_data_df['gas-ct_MWh'].sum()) / generation_annual_total_MWh
        generation_annual_bio_fraction = (cambium_data_df['biomass_MWh'].sum() + cambium_data_df['beccs_MWh'].sum()) / generation_annual_total_MWh
        generation_annual_geothermal_fraction = cambium_data_df['geothermal_MWh'].sum() / generation_annual_total_MWh
        generation_annual_hydro_fraction = (cambium_data_df['hydro_MWh'].sum() + cambium_data_df['phs_MWh'].sum()) / generation_annual_total_MWh
        generation_annual_wind_fraction = (cambium_data_df['wind-ons_MWh'].sum() + cambium_data_df['wind-ofs_MWh'].sum()) / generation_annual_total_MWh
        generation_annual_solar_fraction = (cambium_data_df['upv_MWh'].sum() + cambium_data_df['distpv_MWh'].sum() + cambium_data_df['csp_MWh'].sum()) / generation_annual_total_MWh
        generation_annual_battery_fraction = (cambium_data_df['battery_MWh'].sum()) / generation_annual_total_MWh

        grid_generation_fraction = {'Nuclear':generation_annual_nuclear_fraction,
                                    'Coal & Oil':generation_annual_coal_oil_fraction,
                                    'Gas':generation_annual_gas_fraction,
                                    'Bio':generation_annual_bio_fraction,
                                    'Geothermal':generation_annual_geothermal_fraction,
                                    'Hydro':generation_annual_hydro_fraction,
                                    'Wind':generation_annual_wind_fraction,
                                    'Solar':generation_annual_solar_fraction,
                                    'Battery':generation_annual_battery_fraction}

        grid_imbedded_EI = (generation_annual_nuclear_fraction * nuclear_capex_EI) + (generation_annual_coal_oil_fraction * coal_capex_EI) + (generation_annual_gas_fraction * gas_capex_EI) + (generation_annual_bio_fraction * bio_capex_EI)\
                         + (generation_annual_geothermal_fraction * geothermal_capex_EI) + (generation_annual_hydro_fraction * hydro_capex_EI) + (generation_annual_wind_fraction * wind_capex_EI) + (generation_annual_solar_fraction * solar_pv_capex_EI)\
                         + (generation_annual_battery_fraction * battery_EI) #(g CO2e/kwh)

        if 'hybrid-grid' in grid_case:
            # Calculate grid-connected electrolysis emissions (kg CO2e/kg H2), future cases should reflect targeted electrolyzer electricity usage
            electrolysis_Scope3_EI = ely_stack_capex_EI + ((scope3_grid_emissions_annual_sum + (wind_capex_EI * g_to_kg * wind_annual_energy_kwh) + (solar_pv_capex_EI * g_to_kg * solar_pv_annual_energy_kwh) + (grid_imbedded_EI * g_to_kg * grid_annual_sum_MWh * MWh_to_kWh)) / h2prod_annual_sum)
            electrolysis_Scope2_EI = scope2_grid_emissions_annual_sum / h2prod_annual_sum 
            electrolysis_Scope1_EI = 0
            electrolysis_total_EI  = electrolysis_Scope1_EI + electrolysis_Scope2_EI + electrolysis_Scope3_EI 
            electrolysis_total_EI_policy_grid = electrolysis_total_EI
            #TODO: Masha, shouldn't electrolysis_total_EI_policy_offgrid still include ely_stack_capex_EI?
            #NOTE: electrolysis_total_EI_policy_offgrid not used in subsequent calculations
            electrolysis_total_EI_policy_offgrid = 0 
            # Calculate ammonia emissions via hybrid grid electrolysis (kg CO2e/kg NH3)
            NH3_electrolysis_Scope3_EI = (NH3_H2_consume * electrolysis_total_EI) + (NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean())
            NH3_electrolysis_Scope2_EI = NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean()
            NH3_electrolysis_Scope1_EI = NH3_boiler_EI
            NH3_electrolysis_total_EI  = NH3_electrolysis_Scope1_EI + NH3_electrolysis_Scope2_EI + NH3_electrolysis_Scope3_EI
            # Calculate steel emissions via hybrid grid electrolysis (kg CO2e/metric tonne steel)
            steel_electrolysis_Scope3_EI = (steel_H2_consume * MT_to_kg * electrolysis_total_EI) + (steel_lime_EI * steel_lime_consume * MT_to_kg) + (steel_iron_ore_EI * steel_iron_ore_consume  * MT_to_kg) + (steel_NG_supply_EI * steel_NG_consume) + (cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean() * steel_PO_consume) + ((steel_H2O_EI / gal_H2O_to_MT) * steel_H2O_consume)
            steel_electrolysis_Scope2_EI = steel_PO_consume * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean()  
            steel_electrolysis_Scope1_EI = steel_CH4_prod + steel_CO2_prod
            steel_electrolysis_total_EI  = steel_electrolysis_Scope1_EI + steel_electrolysis_Scope2_EI + steel_electrolysis_Scope3_EI

        if 'grid-only' in grid_case:
            # Calculate SMR emissions. SMR and SMR + CCS are always grid-connected (kg CO2e/kg H2)
            smr_Scope3_EI = (smr_NG_supply * g_to_kg * (smr_NG_consume - smr_steam_prod/smr_HEX_eff)) + (smr_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean())
            smr_Scope2_EI = smr_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean()
            smr_Scope1_EI = smr_NG_combust * g_to_kg * (smr_NG_consume - smr_steam_prod/smr_HEX_eff)
            smr_total_EI  = smr_Scope1_EI + smr_Scope2_EI + smr_Scope3_EI
            electrolysis_total_EI_policy_grid = electrolysis_total_EI
            #TODO: Masha, shouldn't electrolysis_total_EI_policy_offgrid still include ely_stack_capex_EI?
            #NOTE: electrolysis_total_EI_policy_offgrid not used in subsequent calculations
            electrolysis_total_EI_policy_offgrid = 0 
            
            # Calculate ammonia emissions via SMR process (kg CO2e/kg NH3)
            NH3_smr_Scope3_EI = (NH3_H2_consume * smr_total_EI) + (NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean())
            NH3_smr_Scope2_EI = NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean()
            NH3_smr_Scope1_EI = NH3_boiler_EI
            NH3_smr_total_EI = NH3_smr_Scope1_EI + NH3_smr_Scope2_EI + NH3_smr_Scope3_EI   
            
            # Calculate steel emissions via SMR process (kg CO2e/metric tonne steel)
            steel_smr_Scope3_EI = (smr_total_EI * steel_H2_consume * MT_to_kg) + (steel_lime_EI * steel_lime_consume * MT_to_kg) + (steel_iron_ore_EI  * steel_iron_ore_consume  * MT_to_kg) + (steel_NG_supply_EI * steel_NG_consume) + (cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean() * steel_PO_consume) + ((steel_H2O_EI / gal_H2O_to_MT) * steel_H2O_consume)
            steel_smr_Scope2_EI = cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean() * steel_PO_consume 
            steel_smr_Scope1_EI = steel_CH4_prod + steel_CO2_prod
            steel_smr_total_EI  = steel_smr_Scope1_EI + steel_smr_Scope2_EI + steel_smr_Scope3_EI
            
            # Calculate SMR + CCS emissions (kg CO2e/kg H2)
            # TODO: Masha, original formula use smr_steam_prod, GREET shows smr_ccs_steam_prod = 0, validate correct value 
            # TODO: Masha, ccs_PO_consume units (kWh/kg CO2), if value !=0 below equations don't make sense (kWh/kg H2 + kWh/kg CO2) * (kg CO2e/MWh), units dont align
            smr_ccs_Scope3_EI = (smr_NG_supply * g_to_kg * (smr_ccs_NG_consume - smr_steam_prod/smr_HEX_eff)) + (((smr_ccs_PO_consume +  ccs_PO_consume) * kWh_to_MWh) * cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean())
            smr_ccs_Scope2_EI = (((smr_ccs_PO_consume +  ccs_PO_consume) *  kWh_to_MWh) * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean())
            smr_ccs_Scope1_EI = (1-smr_ccs_perc_capture) * smr_NG_combust * g_to_kg * (smr_NG_consume_CCS - smr_steam_prod/smr_HEX_eff)
            smr_ccs_total_EI  = smr_ccs_Scope1_EI + smr_ccs_Scope2_EI + smr_ccs_Scope3_EI    
            
            # Calculate ammonia emissions via SMR with CCS process (kg CO2e/kg NH3)
            NH3_smr_ccs_Scope3_EI = (NH3_H2_consume * smr_ccs_total_EI) + (NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean())
            NH3_smr_ccs_Scope2_EI = NH3_smr_Scope2_EI
            NH3_smr_ccs_Scope1_EI = NH3_smr_Scope1_EI
            NH3_smr_ccs_total_EI = NH3_smr_ccs_Scope1_EI + NH3_smr_ccs_Scope2_EI + NH3_smr_ccs_Scope3_EI   
            
            # Calculate steel emissions via SMR with CCS process (kg CO2e/metric tonne steel)
            steel_smr_ccs_Scope3_EI = (smr_ccs_total_EI * steel_H2_consume * MT_to_kg) + (steel_lime_EI * steel_lime_consume * MT_to_kg) + (steel_iron_ore_EI  * steel_iron_ore_consume  * MT_to_kg) + (steel_NG_supply_EI * steel_NG_consume)  + (cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean() * steel_PO_consume) + ((steel_H2O_EI / gal_H2O_to_MT) * steel_H2O_consume)  
            steel_smr_ccs_Scope2_EI = steel_smr_Scope2_EI 
            steel_smr_ccs_Scope1_EI = steel_smr_Scope1_EI 
            steel_smr_ccs_total_EI  = steel_smr_ccs_Scope1_EI + steel_smr_ccs_Scope2_EI + steel_smr_ccs_Scope3_EI    

            # Calculate grid-connected electrolysis emissions (kg CO2e/kg H2)
            electrolysis_Scope3_EI = ely_stack_capex_EI + ((scope3_grid_emissions_annual_sum + (grid_imbedded_EI * g_to_kg * grid_annual_sum_MWh * MWh_to_kWh))/h2prod_annual_sum)
            electrolysis_Scope2_EI = scope2_grid_emissions_annual_sum / h2prod_annual_sum 
            electrolysis_Scope1_EI = 0
            electrolysis_total_EI = electrolysis_Scope1_EI + electrolysis_Scope2_EI + electrolysis_Scope3_EI

            # Calculate ammonia emissions via grid only electrolysis (kg CO2e/kg NH3)
            NH3_electrolysis_Scope3_EI = (NH3_H2_consume * electrolysis_total_EI) + (NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean())
            NH3_electrolysis_Scope2_EI = NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean()
            NH3_electrolysis_Scope1_EI = NH3_boiler_EI
            NH3_electrolysis_total_EI  = NH3_electrolysis_Scope1_EI + NH3_electrolysis_Scope2_EI + NH3_electrolysis_Scope3_EI

            # Calculate steel emissions via grid only electrolysis (kg CO2e/metric tonne steel)
            steel_electrolysis_Scope3_EI = (steel_H2_consume * MT_to_kg * electrolysis_total_EI) + (steel_lime_EI * steel_lime_consume * MT_to_kg) + (steel_iron_ore_EI  * steel_iron_ore_consume * MT_to_kg) + (steel_NG_supply_EI * steel_NG_consume) + (cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean() * steel_PO_consume) + ((steel_H2O_EI / gal_H2O_to_MT) * steel_H2O_consume)  
            steel_electrolysis_Scope2_EI = steel_PO_consume * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean()  
            steel_electrolysis_Scope1_EI = steel_CH4_prod + steel_CO2_prod
            steel_electrolysis_total_EI  = steel_electrolysis_Scope1_EI + steel_electrolysis_Scope2_EI + steel_electrolysis_Scope3_EI
        if 'off-grid' in grid_case:
            # Calculate renewable only electrolysis emissions (kg CO2e/kg H2)       
            electrolysis_Scope3_EI = ely_stack_capex_EI + (((wind_capex_EI * g_to_kg * wind_annual_energy_kwh) + (solar_pv_capex_EI * g_to_kg * solar_pv_annual_energy_kwh)) /h2prod_annual_sum)
            electrolysis_Scope2_EI = 0
            electrolysis_Scope1_EI = 0
            electrolysis_total_EI = electrolysis_Scope1_EI + electrolysis_Scope2_EI + electrolysis_Scope3_EI
            electrolysis_total_EI_policy_offgrid = electrolysis_total_EI
            #TODO: Masha, shouldn't electrolysis_total_EI_policy_grid still include ely_stack_capex_EI?
            #NOTE: electrolysis_total_EI_policy_grid not used in subsequent calculations
            electrolysis_total_EI_policy_grid = 0

            # Calculate ammonia emissions via renewable electrolysis (kg CO2e/kg NH3)
            NH3_electrolysis_Scope3_EI = (NH3_H2_consume * electrolysis_total_EI) + (NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean())
            NH3_electrolysis_Scope2_EI = NH3_PO_consume * kWh_to_MWh * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean()
            NH3_electrolysis_Scope1_EI = NH3_boiler_EI
            NH3_electrolysis_total_EI = NH3_electrolysis_Scope1_EI + NH3_electrolysis_Scope2_EI + NH3_electrolysis_Scope3_EI

            # Calculate steel emissions via renewable electrolysis (kg CO2e/metric tonne steel)
            steel_electrolysis_Scope3_EI = (steel_H2_consume * MT_to_kg * electrolysis_total_EI) + (steel_lime_EI * steel_lime_consume * MT_to_kg) + (steel_iron_ore_EI * steel_iron_ore_consume * MT_to_kg) + (steel_NG_supply_EI * steel_NG_consume) + (cambium_data['LRMER CO2 equiv. combustion (kg-CO2e/MWh)'].mean() * steel_PO_consume) + ((steel_H2O_EI / gal_H2O_to_MT) * steel_H2O_consume)
            steel_electrolysis_Scope2_EI = steel_PO_consume * cambium_data['LRMER CO2 equiv. production (kg-CO2e/MWh)'].mean() 
            steel_electrolysis_Scope1_EI = steel_CH4_prod + steel_CO2_prod
            steel_electrolysis_total_EI  = steel_electrolysis_Scope1_EI + steel_electrolysis_Scope2_EI + steel_electrolysis_Scope3_EI
        
        # Append emission intensity values for each year to lists
        electrolysis_Scope3_emission_intensity.append(electrolysis_Scope3_EI)
        electrolysis_Scope2_emission_intensity.append(electrolysis_Scope2_EI)
        electrolysis_emission_intensity.append(electrolysis_total_EI)
        smr_Scope3_emission_intensity.append(smr_Scope3_EI)
        smr_Scope2_emission_intensity.append(smr_Scope2_EI)
        smr_emission_intensity.append(smr_total_EI)
        smr_ccs_Scope3_emission_intensity.append(smr_Scope3_EI)
        smr_ccs_Scope2_emission_intensity.append(smr_Scope2_EI)
        smr_ccs_emission_intensity.append(smr_ccs_total_EI)
        NH3_electrolysis_Scope3_emission_intensity.append(NH3_electrolysis_Scope3_EI)
        NH3_electrolysis_Scope2_emission_intensity.append(NH3_electrolysis_Scope2_EI)
        NH3_electrolysis_emission_intensity.append(NH3_electrolysis_total_EI)
        steel_electrolysis_Scope3_emission_intensity.append(steel_electrolysis_Scope3_EI)
        steel_electrolysis_Scope2_emission_intensity.append(steel_electrolysis_Scope2_EI)
        steel_electrolysis_emission_intensity.append(steel_electrolysis_total_EI)
        NH3_smr_Scope3_emission_intensity.append(NH3_smr_Scope3_EI)
        NH3_smr_Scope2_emission_intensity.append(NH3_smr_Scope2_EI)
        NH3_smr_emission_intensity.append(NH3_smr_total_EI)
        steel_smr_Scope3_emission_intensity.append(steel_smr_Scope3_EI)
        steel_smr_Scope2_emission_intensity.append(steel_smr_Scope2_EI)
        steel_smr_emission_intensity.append(steel_smr_total_EI)
        NH3_smr_ccs_Scope3_emission_intensity.append(NH3_smr_ccs_Scope3_EI)
        NH3_smr_ccs_Scope2_emission_intensity.append(NH3_smr_ccs_Scope2_EI)
        NH3_smr_ccs_emission_intensity.append(NH3_smr_ccs_total_EI)
        steel_smr_ccs_Scope3_emission_intensity.append(steel_smr_ccs_Scope3_EI)
        steel_smr_ccs_Scope2_emission_intensity.append(steel_smr_ccs_Scope2_EI)
        steel_smr_ccs_emission_intensity.append(steel_smr_ccs_total_EI)
    
    # Instantiate dataframe from dictionary of emission intensity lists
    emission_intensities_df = pd.DataFrame({'Year':cambium_data.cambium_years,
                                            'electrolysis Scope3 EI (kg CO2e/kg H2)':electrolysis_Scope3_emission_intensity, 
                                            'electrolysis Scope2 EI (kg CO2e/kg H2)':electrolysis_Scope2_emission_intensity, 
                                            'electrolysis EI (kg CO2e/kg H2)':electrolysis_emission_intensity, 
                                            'smr Scope3 EI (kg CO2e/kg H2)': smr_Scope3_emission_intensity, 
                                            'smr Scope2 EI (kg CO2e/kg H2)': smr_Scope2_emission_intensity, 
                                            'smr EI (kg CO2e/kg H2)': smr_emission_intensity, 
                                            'smr ccs Scope3 EI (kg CO2e/kg H2)': smr_ccs_Scope3_emission_intensity, 
                                            'smr ccs Scope2 EI (kg CO2e/kg H2)': smr_ccs_Scope2_emission_intensity, 
                                            'smr ccs EI (kg CO2e/kg H2)': smr_ccs_emission_intensity,      
                                            'NH3 electrolysis Scope3 EI (kg CO2e/kg H2)': NH3_electrolysis_Scope3_emission_intensity, 
                                            'NH3 electrolysis Scope2 EI (kg CO2e/kg H2)': NH3_electrolysis_Scope2_emission_intensity, 
                                            'NH3 electrolysis EI (kg CO2e/kg H2)': NH3_electrolysis_emission_intensity, 
                                            'steel electrolysis Scope3 EI (kg CO2e/kg H2)': steel_electrolysis_Scope3_emission_intensity, 
                                            'steel electrolysis Scope2 EI (kg CO2e/kg H2)': steel_electrolysis_Scope2_emission_intensity, 
                                            'steel electrolysis EI (kg CO2e/kg H2)': steel_electrolysis_emission_intensity,
                                            'NH3 smr Scope3 EI (kg CO2e/kg H2)': NH3_smr_Scope3_emission_intensity, 
                                            'NH3 smr Scope2 EI (kg CO2e/kg H2)': NH3_smr_Scope2_emission_intensity, 
                                            'NH3 smr EI (kg CO2e/kg H2)': NH3_smr_emission_intensity, 
                                            'steel smr Scope3 EI (kg CO2e/kg H2)': steel_smr_Scope3_emission_intensity, 
                                            'steel smr Scope2 EI (kg CO2e/kg H2)': steel_smr_Scope2_emission_intensity, 
                                            'steel smr EI (kg CO2e/kg H2)': steel_smr_emission_intensity,
                                            'NH3 smr ccs Scope3 EI (kg CO2e/kg H2)': NH3_smr_ccs_Scope3_emission_intensity, 
                                            'NH3 smr ccs Scope2 EI (kg CO2e/kg H2)': NH3_smr_ccs_Scope2_emission_intensity, 
                                            'NH3 smr ccs EI (kg CO2e/kg H2)': NH3_smr_ccs_emission_intensity, 
                                            'steel smr ccs Scope3 EI (kg CO2e/kg H2)': steel_smr_ccs_Scope3_emission_intensity, 
                                            'steel smr ccs Scope2 EI (kg CO2e/kg H2)': steel_smr_ccs_Scope2_emission_intensity, 
                                            'steel smr ccs EI (kg CO2e/kg H2)': steel_smr_ccs_emission_intensity,
                                            })
    ## Interpolation of emission intensities for years not captured by cambium (cambium 2023 offers 2025-2050 in 5 year increments)
    endoflife_year = cambium_year + project_lifetime

    # Instantiate lists to hold interpolated data
    electrolysis_Scope3_EI_interpolated = []
    electrolysis_Scope2_EI_interpolated = []
    electrolysis_EI_interpolated = []
    smr_Scope3_EI_interpolated = []
    smr_Scope2_EI_interpolated = []
    smr_EI_interpolated = []
    smr_ccs_Scope3_EI_interpolated = []
    smr_ccs_Scope2_EI_interpolated = []
    smr_ccs_EI_interpolated = []
    NH3_electrolysis_Scope3_EI_interpolated = []
    NH3_electrolysis_Scope2_EI_interpolated = []
    NH3_electrolysis_EI_interpolated = []
    steel_electrolysis_Scope3_EI_interpolated = []
    steel_electrolysis_Scope2_EI_interpolated = []
    steel_electrolysis_EI_interpolated = []
    NH3_smr_Scope3_EI_interpolated = []
    NH3_smr_Scope2_EI_interpolated = []
    NH3_smr_EI_interpolated = []
    steel_smr_Scope3_EI_interpolated = []
    steel_smr_Scope2_EI_interpolated = []
    steel_smr_EI_interpolated = []
    NH3_smr_ccs_Scope3_EI_interpolated = []
    NH3_smr_ccs_Scope2_EI_interpolated = []
    NH3_smr_ccs_EI_interpolated = []
    steel_smr_ccs_Scope3_EI_interpolated = []
    steel_smr_ccs_Scope2_EI_interpolated = []
    steel_smr_ccs_EI_interpolated = []

    # Loop through years between cambium_year and endoflife_year, interpolate values
    for year in range(cambium_year,endoflife_year):
        # thoughts on logic for extrapolation
        # if year < min(cambium_data.cambium_years):
            # logic to extrapolate, likely linear?
        if year <= max(emission_intensities_df['Year']):
            electrolysis_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['electrolysis Scope3 EI (kg CO2e/kg H2)']))
            electrolysis_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['electrolysis Scope2 EI (kg CO2e/kg H2)']))
            electrolysis_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['electrolysis EI (kg CO2e/kg H2)']))
            smr_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['smr Scope3 EI (kg CO2e/kg H2)']))
            smr_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['smr Scope2 EI (kg CO2e/kg H2)']))
            smr_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['smr EI (kg CO2e/kg H2)']))
            smr_ccs_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['smr ccs Scope3 EI (kg CO2e/kg H2)']))
            smr_ccs_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['smr ccs Scope2 EI (kg CO2e/kg H2)']))
            smr_ccs_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['smr ccs EI (kg CO2e/kg H2)']))
            NH3_electrolysis_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 electrolysis Scope3 EI (kg CO2e/kg H2)']))
            NH3_electrolysis_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 electrolysis Scope2 EI (kg CO2e/kg H2)']))
            NH3_electrolysis_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 electrolysis EI (kg CO2e/kg H2)']))
            steel_electrolysis_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel electrolysis Scope3 EI (kg CO2e/kg H2)']))
            steel_electrolysis_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel electrolysis Scope2 EI (kg CO2e/kg H2)']))
            steel_electrolysis_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel electrolysis EI (kg CO2e/kg H2)']))
            NH3_smr_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 smr Scope3 EI (kg CO2e/kg H2)']))
            NH3_smr_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 smr Scope2 EI (kg CO2e/kg H2)']))
            NH3_smr_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 smr EI (kg CO2e/kg H2)']))
            steel_smr_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel smr Scope3 EI (kg CO2e/kg H2)']))
            steel_smr_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel smr Scope2 EI (kg CO2e/kg H2)']))
            steel_smr_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel smr EI (kg CO2e/kg H2)']))  
            NH3_smr_ccs_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 smr ccs Scope3 EI (kg CO2e/kg H2)']))
            NH3_smr_ccs_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 smr ccs Scope2 EI (kg CO2e/kg H2)']))
            NH3_smr_ccs_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['NH3 smr ccs EI (kg CO2e/kg H2)']))
            steel_smr_ccs_Scope3_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel smr ccs Scope3 EI (kg CO2e/kg H2)']))
            steel_smr_ccs_Scope2_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel smr ccs Scope2 EI (kg CO2e/kg H2)']))
            steel_smr_ccs_EI_interpolated.append(np.interp(year,emission_intensities_df['Year'],emission_intensities_df['steel smr ccs EI (kg CO2e/kg H2)']))  
        else:
            electrolysis_Scope3_EI_interpolated.append(emission_intensities_df['electrolysis Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            electrolysis_Scope2_EI_interpolated.append(emission_intensities_df['electrolysis Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            electrolysis_EI_interpolated.append(emission_intensities_df['electrolysis EI (kg CO2e/kg H2)'].values[-1:][0])
            smr_Scope3_EI_interpolated.append(emission_intensities_df['smr Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            smr_Scope2_EI_interpolated.append(emission_intensities_df['smr Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            smr_EI_interpolated.append(emission_intensities_df['smr EI (kg CO2e/kg H2)'].values[-1:][0])
            smr_ccs_Scope3_EI_interpolated.append(emission_intensities_df['smr ccs Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            smr_ccs_Scope2_EI_interpolated.append(emission_intensities_df['smr ccs Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            smr_ccs_EI_interpolated.append(emission_intensities_df['smr ccs EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_electrolysis_Scope3_EI_interpolated.append(emission_intensities_df['NH3 electrolysis Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_electrolysis_Scope2_EI_interpolated.append(emission_intensities_df['NH3 electrolysis Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_electrolysis_EI_interpolated.append(emission_intensities_df['NH3 electrolysis EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_electrolysis_Scope3_EI_interpolated.append(emission_intensities_df['steel electrolysis Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_electrolysis_Scope2_EI_interpolated.append(emission_intensities_df['steel electrolysis Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_electrolysis_EI_interpolated.append(emission_intensities_df['steel electrolysis EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_smr_Scope3_EI_interpolated.append(emission_intensities_df['NH3 smr Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_smr_Scope2_EI_interpolated.append(emission_intensities_df['NH3 smr Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_smr_EI_interpolated.append(emission_intensities_df['NH3 smr EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_smr_Scope3_EI_interpolated.append(emission_intensities_df['steel smr Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_smr_Scope2_EI_interpolated.append(emission_intensities_df['steel smr Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_smr_EI_interpolated.append(emission_intensities_df['steel smr EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_smr_ccs_Scope3_EI_interpolated.append(emission_intensities_df['NH3 smr ccs Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_smr_ccs_Scope2_EI_interpolated.append(emission_intensities_df['NH3 smr ccs Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            NH3_smr_ccs_EI_interpolated.append(emission_intensities_df['NH3 smr ccs EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_smr_ccs_Scope3_EI_interpolated.append(emission_intensities_df['steel smr ccs Scope3 EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_smr_ccs_Scope2_EI_interpolated.append(emission_intensities_df['steel smr ccs Scope2 EI (kg CO2e/kg H2)'].values[-1:][0])
            steel_smr_ccs_EI_interpolated.append(emission_intensities_df['steel smr ccs EI (kg CO2e/kg H2)'].values[-1:][0])

    # Calculate lifetime LCA values
    electrolysis_Scope3_LCA = sum(np.asarray(electrolysis_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    electrolysis_Scope2_LCA = sum(np.asarray(electrolysis_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    electrolysis_total_LCA = sum(np.asarray(electrolysis_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    smr_Scope3_LCA = sum(np.asarray(smr_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    smr_Scope2_LCA = sum(np.asarray(smr_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    smr_total_LCA = sum(np.asarray(smr_EI_interpolated) * h2prod_annual_sum)/h2prod_life_sum
    smr_ccs_Scope3_LCA = sum(np.asarray(smr_ccs_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    smr_ccs_Scope2_LCA = sum(np.asarray(smr_ccs_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    smr_ccs_total_LCA = sum(np.asarray(smr_ccs_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_electrolysis_Scope3_LCA = sum(np.asarray(NH3_electrolysis_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_electrolysis_Scope2_LCA = sum(np.asarray(NH3_electrolysis_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_electrolysis_total_LCA = sum(np.asarray(NH3_electrolysis_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_electrolysis_Scope3_LCA = sum(np.asarray(steel_electrolysis_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_electrolysis_Scope2_LCA = sum(np.asarray(steel_electrolysis_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_electrolysis_total_LCA = sum(np.asarray(steel_electrolysis_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_smr_Scope3_LCA = sum(np.asarray(NH3_smr_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_smr_Scope2_LCA = sum(np.asarray(NH3_smr_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_smr_total_LCA = sum(np.asarray(NH3_smr_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_smr_Scope3_LCA = sum(np.asarray(steel_smr_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_smr_Scope2_LCA = sum(np.asarray(steel_smr_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_smr_total_LCA = sum(np.asarray(steel_smr_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_smr_ccs_Scope3_LCA = sum(np.asarray(NH3_smr_ccs_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_smr_ccs_Scope2_LCA = sum(np.asarray(NH3_smr_ccs_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    NH3_smr_ccs_total_LCA = sum(np.asarray(NH3_smr_ccs_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_smr_ccs_Scope3_LCA = sum(np.asarray(steel_smr_ccs_Scope3_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_smr_ccs_Scope2_LCA = sum(np.asarray(steel_smr_ccs_Scope2_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum
    steel_smr_ccs_total_LCA = sum(np.asarray(steel_smr_ccs_EI_interpolated) * h2prod_annual_sum) /h2prod_life_sum

    # Put all cumulative metrics into a dictionary, then dataframe, then save results to csv
    lca_dict = {'Total Life Cycle H2 Production (tonnes-H2/MW)': [h2prod_annual_sum],
                'Grid Total Scope 2 (Combustion) GHG Emissions (tonnes-CO2e/MW)': [scope2_grid_emissions_annual_sum],
                'Grid Total Scope 3 (Production) GHG Emissions (tonnes-CO2e/MW)': [scope3_grid_emissions_annual_sum],
                'Grid Total Life Cycle Emissions (tonnes-CO2e/MW)' : [total_grid_emissions_annual_sum],
                'Annaul Average Grid Emission Intensity (kg-CO2/MWh)': [grid_emission_intensity_annual_average],
                'SMR Scope 3 GHG Emissions (kg-CO2e/kg-H2)': [smr_Scope3_LCA],
                'SMR Scope 2 GHG Emissions (kg-CO2e/kg-H2)': [smr_Scope2_LCA],
                'SMR Scope 1 GHG Emissions (kg-CO2e/kg-H2)': [smr_Scope1_EI],
                #'SMR Total GHG Emissions (kg-CO2e/kg-H2)': [smr_total_EI],  
                'SMR Total GHG Emissions (kg-CO2e/kg-H2)': [smr_total_LCA], 
                'Ammonia SMR Scope 3 GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_Scope3_LCA],
                'Ammonia SMR Scope 2 GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_Scope2_LCA], 
                'Ammonia SMR Scope 1 GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_Scope1_EI],
                #'Ammonia SMR Total GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_total_EI], 
                'Ammonia SMR Total GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_total_LCA], 
                'Steel SMR Scope 3 GHG Emissions (kg-CO2e/MT steel)': [steel_smr_Scope3_LCA],
                'Steel SMR Scope 2 GHG Emissions (kg-CO2e/MT steel)': [steel_smr_Scope2_LCA],
                'Steel SMR Scope 1 GHG Emissions (kg-CO2e/MT steel)': [steel_smr_Scope1_EI],
                #'Steel SMR Total GHG Emissions (kg-CO2e/MT steel)': [steel_smr_total_EI],    
                'Steel SMR Total GHG Emissions (kg-CO2e/MT steel)': [steel_smr_total_LCA],  
                'SMR with CCS Scope 3 GHG Emissions (kg-CO2e/kg-H2)': [smr_ccs_Scope3_LCA],
                'SMR with CCS Scope 2 GHG Emissions (kg-CO2e/kg-H2)': [smr_ccs_Scope2_LCA],
                'SMR with CCS Scope 1 GHG Emissions (kg-CO2e/kg-H2)': [smr_ccs_Scope1_EI],
                #'SMR with CCS Total GHG Emissions (kg-CO2e/kg-H2)': [smr_ccs_total_EI],     
                'SMR with CCS Total GHG Emissions (kg-CO2e/kg-H2)': [smr_ccs_total_LCA],   
                'Ammonia SMR with CCS Scope 3 GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_ccs_Scope3_LCA],
                'Ammonia SMR with CCS Scope 2 GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_ccs_Scope2_LCA], 
                'Ammonia SMR with CCS Scope 1 GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_ccs_Scope1_EI],
                #'Ammonia SMR with CCS Total GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_ccs_total_EI],    
                'Ammonia SMR with CCS Total GHG Emissions (kg-CO2e/kg-NH3)': [NH3_smr_ccs_total_LCA], 
                'Steel SMR with CCS Scope 3 GHG Emissions (kg-CO2e/MT steel)': [steel_smr_ccs_Scope3_LCA],
                'Steel SMR with CCS Scope 2 GHG Emissions (kg-CO2e/MT steel)': [steel_smr_ccs_Scope2_LCA],
                'Steel SMR with CCS Scope 1 GHG Emissions (kg-CO2e/MT steel)': [steel_smr_ccs_Scope1_EI],
                # 'Steel SMR with CCS Total GHG Emissions (kg-CO2e/MT steel)': [steel_smr_ccs_total_EI],  
                'Steel SMR with CCS Total GHG Emissions (kg-CO2e/MT steel)': [steel_smr_ccs_total_LCA],                  
                'Electrolysis Scope 3 GHG Emissions (kg-CO2e/kg-H2)':[electrolysis_Scope3_LCA],
                'Electrolysis Scope 2 GHG Emissions (kg-CO2e/kg-H2)':[electrolysis_Scope2_LCA],
                'Electrolysis Scope 1 GHG Emissions (kg-CO2e/kg-H2)':[electrolysis_Scope1_EI],   
                #'Electrolysis Total GHG Emissions (kg-CO2e/kg-H2)':[electrolysis_total_EI],            
                'Electrolysis Total GHG Emissions (kg-CO2e/kg-H2)':[electrolysis_total_LCA],
                'Ammonia Electrolysis Scope 3 GHG Emissions (kg-CO2e/kg-NH3)':[NH3_electrolysis_Scope3_LCA],
                'Ammonia Electrolysis Scope 2 GHG Emissions (kg-CO2e/kg-NH3)':[NH3_electrolysis_Scope2_LCA],
                'Ammonia Electrolysis Scope 1 GHG Emissions (kg-CO2e/kg-NH3)':[NH3_electrolysis_Scope1_EI],   
                #'Ammonia Electrolysis Total GHG Emissions (kg-CO2e/kg-NH3)':[NH3_electrolysis_total_EI],     
                'Ammonia Electrolysis Total GHG Emissions (kg-CO2e/kg-NH3)':[NH3_electrolysis_total_LCA],                              
                'Steel Electrolysis Scope 3 GHG Emissions (kg-CO2e/MT steel)':[steel_electrolysis_Scope3_LCA],
                'Steel Electrolysis Scope 2 GHG Emissions (kg-CO2e/MT steel)':[steel_electrolysis_Scope2_LCA],
                'Steel Electrolysis Scope 1 GHG Emissions (kg-CO2e/MT steel)':[steel_electrolysis_Scope1_EI],   
                #'Steel Electrolysis Total GHG Emissions (kg-CO2e/MT steel)':[steel_electrolysis_total_EI]
                'Steel Electrolysis Total GHG Emissions (kg-CO2e/MT steel)':[steel_electrolysis_total_LCA]
                }
    emissions_and_h2_df = pd.DataFrame(data=lca_dict)

# set up function to post-process HOPP results
def post_process_simulation(
    lcoe,
    lcoh,
    pf_lcoh,
    pf_lcoe,
    hopp_results,
    electrolyzer_physics_results,
    hopp_config,
    greenheart_config,
    orbit_config,
    turbine_config,
    h2_storage_results,
    capex_breakdown,
    opex_breakdown,
    wind_cost_results,
    platform_results,
    desal_results,
    design_scenario,
    plant_design_number,
    incentive_option,
    solver_results=[],
    show_plots=False,
    save_plots=False,
    verbose=False,
    output_dir="./output/",
):  # , lcoe, lcoh, lcoh_with_grid, lcoh_grid_only):
    # colors (official NREL color palette https://brand.nrel.gov/content/index/guid/color_palette?parent=61)
    colors = [
        "#0079C2",
        "#00A4E4",
        "#F7A11A",
        "#FFC423",
        "#5D9732",
        "#8CC63F",
        "#5E6A71",
        "#D1D5D8",
        "#933C06",
        "#D9531E",
    ]
    # load saved results

    # post process results
    if verbose:
        print("LCOE: ", round(lcoe * 1e3, 2), "$/MWh")
        print("LCOH: ", round(lcoh, 2), "$/kg")
        print(
            "hybrid electricity plant capacity factor: ",
            round(
                np.sum(hopp_results["combined_hybrid_power_production_hopp"])
                / (hopp_results["hybrid_plant"].system_capacity_kw.hybrid * 365 * 24),
                2,
            ),
        )
        print(
            "electrolyzer capacity factor: ",
            round(
                np.sum(electrolyzer_physics_results["power_to_electrolyzer_kw"])
                * 1e-3
                / (greenheart_config["electrolyzer"]["rating"] * 365 * 24),
                2,
            ),
        )
        print(
            "Electorlyzer CAPEX installed $/kW: ",
            round(
                capex_breakdown["electrolyzer"]
                / (greenheart_config["electrolyzer"]["rating"] * 1e3),
                2,
            ),
        )

    if show_plots or save_plots:
        visualize_plant(
            hopp_config,
            greenheart_config,
            turbine_config,
            wind_cost_results,
            hopp_results,
            platform_results,
            desal_results,
            h2_storage_results,
            electrolyzer_physics_results,
            design_scenario,
            colors,
            plant_design_number,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )
    savepaths = [
        output_dir + "data/",
        output_dir + "data/lcoe/",
        output_dir + "data/lcoh/",
        output_dir + "data/lca/"
    ]
    for sp in savepaths:
        if not os.path.exists(sp):
            os.makedirs(sp)

    pf_lcoh.get_cost_breakdown().to_csv(
        savepaths[2]
        + "cost_breakdown_lcoh_design%i_incentive%i_%sstorage.csv"
        % (
            plant_design_number,
            incentive_option,
            greenheart_config["h2_storage"]["type"],
        )
    )
    pf_lcoe.get_cost_breakdown().to_csv(
        savepaths[1]
        + "cost_breakdown_lcoe_design%i_incentive%i_%sstorage.csv"
        % (
            plant_design_number,
            incentive_option,
            greenheart_config["h2_storage"]["type"],
        )
    )

    # create dataframe for saving all the stuff
    greenheart_config["design_scenario"] = design_scenario
    greenheart_config["plant_design_number"] = plant_design_number
    greenheart_config["incentive_options"] = incentive_option

    # save power usage data
    if len(solver_results) > 0:
        hours = len(hopp_results["combined_hybrid_power_production_hopp"])
        annual_energy_breakdown = {
            "electricity_generation_kwh": sum(
                hopp_results["combined_hybrid_power_production_hopp"]
            ),
            "electrolyzer_kwh": sum(
                electrolyzer_physics_results["power_to_electrolyzer_kw"]
            ),
            "renewable_kwh": solver_results[0] * hours,
            "grid_power_kwh": solver_results[1] * hours,
            "desal_kwh": solver_results[2] * hours,
            "h2_transport_compressor_power_kwh": solver_results[3] * hours,
            "h2_storage_power_kwh": solver_results[4] * hours,
        }


    ######################### save detailed ORBIT cost information
    if wind_cost_results.orbit_project:
        _, orbit_capex_breakdown, wind_capex_multiplier = adjust_orbit_costs(
            orbit_project=wind_cost_results.orbit_project,
            greenheart_config=greenheart_config,
        )

        # orbit_capex_breakdown["Onshore Substation"] = orbit_project.phases["ElectricalDesign"].onshore_cost
        # discount ORBIT cost information
        for key in orbit_capex_breakdown:
            orbit_capex_breakdown[key] = -npf.fv(
                greenheart_config["finance_parameters"]["costing_general_inflation"],
                greenheart_config["project_parameters"]["cost_year"]
                - greenheart_config["finance_parameters"]["discount_years"]["wind"],
                0.0,
                orbit_capex_breakdown[key],
            )

        # save ORBIT cost information
        ob_df = pd.DataFrame(orbit_capex_breakdown, index=[0]).transpose()
        savedir = output_dir + "data/orbit_costs/"
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        ob_df.to_csv(
            savedir
            + "orbit_cost_breakdown_lcoh_design%i_incentive%i_%sstorage.csv"
            % (
                plant_design_number,
                incentive_option,
                greenheart_config["h2_storage"]["type"],
            )
        )
        ###############################

        ###################### Save export system breakdown from ORBIT ###################

        _, orbit_capex_breakdown, wind_capex_multiplier = adjust_orbit_costs(
            orbit_project=wind_cost_results.orbit_project,
            greenheart_config=greenheart_config,
        )

        onshore_substation_costs = (
            wind_cost_results.orbit_project.phases["ElectricalDesign"].onshore_cost
            * wind_capex_multiplier
        )

        orbit_capex_breakdown["Export System Installation"] -= onshore_substation_costs

        orbit_capex_breakdown[
            "Onshore Substation and Installation"
        ] = onshore_substation_costs

        # discount ORBIT cost information
        for key in orbit_capex_breakdown:
            orbit_capex_breakdown[key] = -npf.fv(
                greenheart_config["finance_parameters"]["costing_general_inflation"],
                greenheart_config["project_parameters"]["cost_year"]
                - greenheart_config["finance_parameters"]["discount_years"]["wind"],
                0.0,
                orbit_capex_breakdown[key],
            )

        # save ORBIT cost information
        ob_df = pd.DataFrame(orbit_capex_breakdown, index=[0]).transpose()
        savedir = output_dir + "data/orbit_costs/"
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        ob_df.to_csv(
            savedir
            + "orbit_cost_breakdown_with_onshore_substation_lcoh_design%i_incentive%i_%sstorage.csv"
            % (
                plant_design_number,
                incentive_option,
                greenheart_config["h2_storage"]["type"],
            )
        )

    ##################################################################################
    if (
        hasattr(hopp_results["hybrid_plant"], "dispatch_builder")
        and hopp_results["hybrid_plant"].battery
    ):
        savedir = output_dir + "figures/production/"
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        plot_tools.plot_generation_profile(
            hopp_results["hybrid_plant"],
            start_day=0,
            n_days=10,
            plot_filename=os.path.abspath(savedir + "generation_profile.pdf"),
            font_size=14,
            power_scale=1 / 1000,
            solar_color="r",
            wind_color="b",
            wave_color="g",
            discharge_color="b",
            charge_color="r",
            gen_color="g",
            price_color="r",
            show_price=False,
        )
    else:
        print(
            "generation profile not plotted because HoppInterface does not have a "
            "'dispatch_builder'"
        )

    # save production information
    hourly_energy_breakdown = save_energy_flows(
        hopp_results["hybrid_plant"],
        electrolyzer_physics_results,
        solver_results,
        hours,
        h2_storage_results,
        output_dir=output_dir
    )

    # save hydrogen information
    key = "Hydrogen Hourly Production [kg/hr]"
    np.savetxt(
        output_dir + "h2_usage",
        electrolyzer_physics_results["H2_Results"][key],
        header="# " + key
    )

    print("*******************TESTING******************")
    print("IN POST_PROCESS_SIMULATION()")
    # print("power_to_electrolyzer_kw")
    # print(np.array(electrolyzer_physics_results['power_to_electrolyzer_kw']))
    # print(np.array(electrolyzer_physics_results['power_to_electrolyzer_kw']).shape)
    # print("energy_to_electrolyzer_from_renewables")
    # print(np.array(hopp_results['combined_hybrid_power_production_hopp']))
    # print(np.array(hopp_results['combined_hybrid_power_production_hopp']).shape)
    # print(np.array(electrolyzer_physics_results['H2_Results']['Hydrogen Hourly Production [kg/hr]']))
    # print(np.array(electrolyzer_physics_results['H2_Results']['Hydrogen Hourly Production [kg/hr]']).shape)
    # energy_to_electrolyzer = np.array(electrolyzer_physics_results['power_to_electrolyzer_kw'])
    # energy_to_electrolyzer_from_renewables = np.array(hopp_results['combined_hybrid_power_production_hopp'])
    # energy_to_electrolyzer_from_grid = energy_to_electrolyzer - energy_to_electrolyzer_from_renewables
    # energy_to_electrolyzer_from_grid = np.maximum(energy_to_electrolyzer_from_grid,0)
    # h2_hourly_prod_kg = np.array(electrolyzer_physics_results['H2_Results']['Hydrogen Hourly Production [kg/hr]'])

    # electrolyzer_profiles_data_dict = {'Energy to electrolyzer (kWh)': energy_to_electrolyzer,
    #                                    'Energy from grid (kWh)': energy_to_electrolyzer_from_grid,
    #                                    'Energy from renewables (kWh)': energy_to_electrolyzer_from_renewables,
    #                                    'Hydrogen Hourly production (kg)': h2_hourly_prod_kg}

    # electrolyzer_profiles_df = pd.DataFrame(data=electrolyzer_profiles_data_dict)
    # electrolyzer_profiles_df = electrolyzer_profiles_df.reset_index().rename(columns={'index':'Interval'})
    # electrolyzer_profiles_df['Interval'] = electrolyzer_profiles_df['Interval']+1
    # electrolyzer_profiles_df = electrolyzer_profiles_df.set_index('Interval')
    # print(electrolyzer_profiles_df)
    
    return annual_energy_breakdown, hourly_energy_breakdown

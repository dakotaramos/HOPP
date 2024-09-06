import os
import csv
from pathlib import Path
from typing import Union
import requests
import json
import yaml
import openpyxl

from hopp.utilities.log import hybrid_logger as logger
from hopp.simulation.technologies.resource.resource import Resource
from hopp import ROOT_DIR

class GREETData:
    """
    Class to manage GREET data

    Args:
        year: (int) version / vintage of GREET to pull data from
        path_resource: directory where to save data files
        filepath: file path of resource file to load
        preprocess_greet: Make an API call even if there is an existing data file. Default == False
        kwargs: additiona keyword arguments

    """

    def __init__(
        self,
        greet_year: int = 2023,
        path_resource: Union[str, Path] = ROOT_DIR / "simulation" / "resource_files",
        filepath: Union[str, Path] ="",
        preprocess_greet: bool = False,
        **kwargs
    ):

        self.year = greet_year

        self.__dict__.update(kwargs)

        # Check if path_resource is a directory, if yes define as self.path_resource attribute
        if os.path.isdir(path_resource):
            self.path_resource = path_resource
        
        # update path with GREET directory
        self.path_resource = os.path.join(self.path_resource, 'greet', str(self.year))

        # Force override internal definitions if kwargs passed in
        self.__dict__.update(kwargs)

        # Define the filepath and filename for the resource file
        if filepath == "":
            filepath = os.path.join(self.path_resource, "greet" + "_" + str(self.year) + "_" + "processed.yaml")
        self.filename = filepath

        # Check if the download directory exists (HOPP/hopp/simulation/resource_files/greet/<year>), if not make the directory
        self.check_download_dir()

        # If a processed resource file does not already exist or preprocess_greet flag == True, process the data
        if not os.path.isfile(self.filename) or preprocess_greet:
            self.preprocess_greet()

        # Check if greet_X_processed.yaml exists, if not error, if yes load yaml to dictionary and save to self.data
        self.format_data()

        # logger.info("CambiumData: {}".format(self.filename))

    def check_download_dir(self):
        if not os.path.isdir(os.path.dirname(self.filename)):
            os.makedirs(os.path.dirname(self.filename))

    def preprocess_greet(self):
        ## Define Conversions for GREET data points
        # Unit conversions
        btu_to_kWh = 0.00029307107                          # 1 btu = 0.00029307107 kWh
        btu_to_MJ = 0.00105505585262                        # 1 btu = 0.00105505585262 MJ 
        mmbtu_to_kWh = 293.07107                            # 1 MMbtu = 293.07107 kWh 
        mmbtu_to_MWh = 0.29307107                           # 1 MMbtu = 0.29307107 MWh
        mmbtu_to_MJ = 1055.05585262                         # 1 MMbtu = 1055.05585262 MJ
        mmbtu_to_GJ = 1.0550585262                          # 1 MMbtu = 1.05505585262 GJ
        kWh_to_MJ = 3.6                                     # 1 kWh = 3.6 MJ
        MJ_to_kWh = (1/3.6)                                 # 1 MJ = (1/3.6) kWh ~= 0.2777777777777778 kWh
        ton_to_kg = 907.18474                               # 1 ton = 907.18474 kg
        ton_to_MT = 0.90718474                              # 1 ton = 0.90718474 metric tonne
        g_to_kg  = 0.001                                    # 1 g = 0.001 kg
        kg_to_MT = 0.001                                    # 1 kg = 0.001 metric tonne
        MT_to_kg = 1000                                     # 1 metric tonne = 1000 kg
        kWh_to_MWh = 0.001                                  # 1 kWh = 0.001 MWh
        MWh_to_kWh = 1000                                   # 1 MWh = 1000 kWh

        # Chemical properties
        mmbtuhhv_per_kg_h2 = 0.134                          # 1 kg H2 = 0.134 MMbtu-hhv h2
        kg_h2_per_mmbtuhhv = 7.462                          # 1 MMbtu-hhv H2 = 7.462 kg H2
        mmbtulhv_per_kg_h2 = 0.114                          # 1 kg H2 = 0.114 MMbtu-lhv H2
        kg_h2_per_mmbtulhv = 8.772                          # 1 MMbtu-lhv H2 = 8.772 kg H2
        kg_CH4_to_kg_CO2e = 29.8                            # 1 kg CH4 = 29.8 kg CO2e 
        MJHHV_per_kg_h2 = 141.88                            # Higher Heating Value of hydrogen = 141.88 MJ-HHV/kg H2
        MJLHV_per_kg_h2 = 119.96                            # Lower Heating Value of hydrogen = 119.96 MJ-LHV/kg H2
        kWh_per_kg_h2_LHV = (MJLHV_per_kg_h2 * MJ_to_kWh)   # kWh per kg of hydrogen using LHV, ~= 33.3222222 kWh/kg H2
        kWh_per_kg_h2_HHV = (MJHHV_per_kg_h2 * MJ_to_kWh)   # kWh per kg of hydrogen using HHV, ~= 39.4111111 kWh/kg H2
        gal_H2O_to_MT = 0.00378541                          # 1 US gallon of H2O = 0.00378541 metric tonnes (1 gal = 3.78541 liters, 1 liter H2O = 1 kg, 1000 kg = 1 metric tonne)
        MMBTU_NG_to_kg_CO2e = 53                            # 1 MMBTU Natural Gas = ~53kg CO2e when combusted

        # Chemical conversion formulas for greenhouse gases (GHGs) to CO2e emissions intensities (EI) with GWP and Carbon ratios
        # CO2 (VOC, CO, CO2) = CO2 + (VOC*0.85/0.27) + (CO*0.43/0.27)
        # GHGs = CO2 (VOC, CO, CO2) + (CH4*29.8) + (N2O)*273 + (VOC*0) + (CO*0) + (NOx*0) + (BC*0) + (OC*0)
        # GHGs = CO2 + (VOC*0.85/0.27) + (CO*0.43/0.27) + (CH4*29.8) + (N20*273)
        # Carbon Ratios
        VOC_to_CO2e = (0.85/0.272727)
        CO_to_CO2e = (0.4285710/0.272727)
        CH4_to_CO2e = (.75/.272727)
        # Global Warming Potential (relative to CO2)
        CH4_gwp_to_CO2e = 29.8          # 1g CH4 == 29.8 g CO2e (not combusted)
        N2O_gwp_to_CO2e = 273           # 1g N2O == 273 g CO2e (not combusted)

        ## Define hardcoded values for efficiencies, emissions intensities (EI), combustion, consumption, and production processes
        #TODO: In future, update to pull hardcoded values from GREET or other models programmatically if possible
        # Following values determined through communications with GREET / ANL team
        NH3_boiler_EI = 0.5             # Boiler combustion of methane for Ammonia (kg CO2e/kg NH3)
        smr_NG_combust = 56.2           # Natural gas combustion emission factor (g CO2e/MJ)
        smr_HEX_eff = 0.9               # Heat exchange efficiency (%)
        smr_NG_supply = 9               # Natural gas extraction and supply to SMR plant assuming 2% CH4 leakage rate (g CO2e/MJ)
        ccs_PO_consume = 0              # Power consumption for CCS (kWh/kg CO2)
        atr_steam_prod = 0              # No steam exported during ATR (MJ/kg H2)
        atr_ccs_steam_prod = 0          # No steam exported during ATR (MJ/kg H2)

        # Following values determined through communications with LBNL
        # TODO: confirm / validate we can pull these values from Greet with Masha (compare original values with GREET values)
        # steel_H2_consume = 0.06596      # NOTE: possibly pull from Greet2 > Steel > AK66, Metric tonnes of H2 per tonne of steel (ton H2/ton steel), from comms with LBNL
        # steel_NG_consume = 0.71657      # NOTE: possibly pull from Greet2 > Steel > AE63+AG63+AK63+AM63 OR AE78+AG78+AK78+AM78, GJ-LHV per tonne of steel (GJ-LHV/ ton steel), from comms with LBNL
        # steel_lime_consume = 0.01812    # NOTE: possibly pull from Greet2 > Steel > AM68, metric tonne of lime per tonne of steel (ton lime / ton steel), from comms with LBNL
        # steel_iron_ore_consume = 1.629  # NOTE: possilby pull from Greet2 > Steel > AM69, metric tonnes of iron ore per metric tonne of steel, from comms with LBNL
        # steel_PO_consume = 0.5502       # NOTE: possibly pull from Greet2 > Steel > Y108, MWh per metric tonne of steel, from comms with LBNL
        # steel_H2O_consume = 0.8037      # NOTE: possibly pull from Greet2 > Steel > Y113, metric tonnes of H2O per tonne of steel, from comms with LBNL
        
        # Following values determined from prior NREL knowledge of electrolysis and assumptions about future grid mix
        # TODO: confirm if fuel_to_grid_curr and _futu can be determined through GREET / Cambium
        grid_trans_losses = 0.05    # Grid losses of 5% are assumed (-)
        fuel_to_grid_curr = 48      # Fuel mix emission intensity for current power grid (g CO2e/kWh)
        fuel_to_grid_futu = 14      # Fuel mix emission intensity for future power grid (g CO2e/kWh) #TODO: define future power grid
        # ely_PO_consume = 55         # NOTE: possibly pull from Greet1 > Hydrogen > J239:L239, Electrolysis power consumption per kg h2 (kWh/kg H2)

        ## Pull GREET Values
        # NOTE / TODO: following logic / cells to pull data from was created for GREET 2023, in future versions may need to add if statements for year / update cells to pull from
        # Define GREET filepaths
        greet1_blue_NH3 = os.path.join(self.path_resource, "blue_NH3_prod", "GREET1_2023_Rev1.xlsm")
        greet2_blue_NH3 = os.path.join(self.path_resource, "blue_NH3_prod", "GREET2_2023_Rev1.xlsm")
        greet1_green_NH3 = os.path.join(self.path_resource, "green_NH3_prod", "GREET1_2023_Rev1.xlsm")
        greet2_green_NH3 = os.path.join(self.path_resource, "green_NH3_prod", "GREET2_2023_Rev1.xlsm")
        greet1_conventional_NH3 = os.path.join(self.path_resource, "conventional_NH3_prod", "GREET1_2023_Rev1.xlsm")
        greet2_conventional_NH3 = os.path.join(self.path_resource, "conventional_NH3_prod", "GREET2_2023_Rev1.xlsm")
        greet1_ccs_central_h2 = os.path.join(self.path_resource, "ccs_central_h2_prod", "GREET1_2023_Rev1.xlsm")
        greet2_ccs_central_h2 = os.path.join(self.path_resource, "ccs_central_h2_prod", "GREET2_2023_Rev1.xlsm")
        greet1_no_ccs_central_h2 = os.path.join(self.path_resource,"no_ccs_central_h2_prod", "GREET1_2023_Rev1.xlsm")
        greet2_no_ccs_central_h2 = os.path.join(self.path_resource,"no_ccs_central_h2_prod", "GREET2_2023_Rev1.xlsm")

        #------------------------------------------------------------------------------
        # Renewable infrastructure embedded emission intensities, Hydrogen production via water electrolysis, and Steel
        #------------------------------------------------------------------------------
        # NOTE: Capex EI, Electrolysis, and Steel GREET values agnostic of ccs/no_ccs and NH3 production methods, ie: can use any version of greet to pull
        # NOTE: For Steel, alternative DRI-EAF configurations (w/ and w/out scrap, H2 vs NG) found in greet2 > Steel > W107:Z136
                # Iron ore vs scrap % controlled by B24:C24 values
                # greet2 > Steel > B17:B18 controls NG vs RNG, changing these values drastically changes steel_NG_supply_EI
                    # May require hosting of different steel config GREET versions if desired
                # Values below are for DRI-EAF 83% H2, 100% DRI 0% Scrap
        greet1 = openpyxl.load_workbook(greet1_ccs_central_h2, data_only=True)
        # Renewable Infrastructure
        wind_capex_EI = (greet1['ElecInfra']['G112'].value / mmbtu_to_kWh)                          # NOTE: original value = 10, greet value = 9.77, Wind CAPEX emissions (g CO2e/kWh)
        solar_pv_capex_EI = (greet1['ElecInfra']['H112'].value / mmbtu_to_kWh)                      # NOTE: original value = 37, greet value = 35.87, Solar PV CAPEX emissions (g CO2e/kWh)
        nuclear_PWR_capex_EI = (greet1['ElecInfra']['D112'].value / mmbtu_to_kWh)                   # NOTE: original value = 0.3, greet value = 0.22, Nuclear Pressurized Water Reactor (PWR) CAPEX emissions (g CO2e/kWh)
        nuclear_BWR_capex_EI = (greet1['ElecInfra']['E112'].value / mmbtu_to_kWh)                   # NOTE: original value = 0.3, greet value = 0.21, Nuclear Boiling Water Reactor (BWR) CAPEX emissions (g CO2e/kWh)
        coal_capex_EI = (greet1['ElecInfra']['B112'].value / mmbtu_to_kWh)                          # NOTE: original value = 0.8, greet value = 0.78, Coal CAPEX emissions (g CO2e/kWh)
        gas_capex_EI = (greet1['ElecInfra']['C112'].value / mmbtu_to_kWh)                           # NOTE: original value = 0.42, greet value = 0.42, Natural Gas Combined Cycle (NGCC) CAPEX emissions (g CO2e/kWh)
        hydro_capex_EI = (greet1['ElecInfra']['F112'].value / mmbtu_to_kWh)                         # NOTE: original value = 7.22, greet value = 7.12, Hydro CAPEX emissions (g CO2e/kWh)
        bio_capex_EI = (greet1['ElecInfra']['L112'].value / mmbtu_to_kWh)                           # NOTE: original value = 0.81, greet value = 0.80, Biomass CAPEX emissions (g CO2e/kWh)
        geothermal_EGS_capex_EI = (greet1['ElecInfra']['I112'].value / mmbtu_to_kWh)                # NOTE: original value = 20.71, greet value = 20.52, Geothermal EGS CAPEX emissions (g CO2e/kWh)
        geothermal_flash_capex_EI = (greet1['ElecInfra']['J112'].value / mmbtu_to_kWh)              # NOTE: original value = 20.71, greet value = 4.61, Geothermal Flash CAPEX emissions (g CO2e/kWh)
        geothermal_binary_capex_EI = (greet1['ElecInfra']['K112'].value / mmbtu_to_kWh)             # NOTE: original value = 20.71, greet value = 20.44, Geothermal Binary CAPEX emissions (g CO2e/kWh)
        # Electrolysis
        pem_ely_PO_consume = (greet1['Hydrogen']['J239'].value * btu_to_kWh * mmbtulhv_per_kg_h2)   # NOTE: original value = 55, greet value = 58.47, PEM water electrolysis energy consumption (kWh/kg H2)
        # Steel
        steel_lime_EI = ((greet1['Chemicals']['BA247'].value +                                      # NOTE: original value = 1.28, greet value = 1.28, Lime production emissions for use in DRI-EAF Steel production (kg CO2e/kg lime)
                            (greet1['Chemicals']['BA237'].value * VOC_to_CO2e) +             
                            (greet1['Chemicals']['BA238'].value * CO_to_CO2e) +
                            (greet1['Chemicals']['BA245'].value * CH4_gwp_to_CO2e) +
                            (greet1['Chemicals']['BA246'].value * N2O_gwp_to_CO2e)
                            ) * g_to_kg * (1/ton_to_kg))
        greet1.close()

        greet2 = openpyxl.load_workbook(greet2_ccs_central_h2, data_only=True)
        # Renewable Infrastructure
        battery_LFP_residential_EI = (greet2['Solar_PV']['DI289'].value)                            # NOTE: original value = 20 (no indication of residential or commercial), greet value = 172.92 (30 yrs, 100 MWh battery), Battery embodied emissions for residential solar PV applications (g CO2e/kWh), assumed LFP batteries (can update battery chemistry in cell K155)
                                                                                                        # value must be multiplied by factor = (project_lifetime / battery_system_capacity_kwh) before use in LCA calculations, opted not to divide here so value remains agnostic of project lifetime and system battery size        
        battery_LFP_commercial_EI = (greet2['Solar_PV']['DJ289'].value)                             # NOTE: original value = 20 (no indication of residential or commercial), greet value = 106362.41 (30 yrs, 100 MWh battery), Battery embodied emissions for commercial solar PV applications (g CO2e/kWh), assumed LFP batteries (can update battery chemistry in cell W155)
                                                                                                        # value must be multiplied by factor = (project_lifetime / battery_system_capacity_kwh) before use in LCA calculations, opted not to divide here so value remains agnostic of project lifetime and system battery size
        # Electrolysis TODO: confirm with Masha stack or stack+BOP
        pem_ely_stack_capex_EI = (greet2['Electrolyzers']['I257'].value * g_to_kg)                  # NOTE: original value = 0.019, greet value = 0.01357 PEM electrolyzer stack CAPEX emissions (kg CO2e/kg H2)
        pem_ely_stack_and_BoP_capex_EI = (greet2['Electrolyzers']['L257'].value * g_to_kg)          # NOTE: original value = 0.019, greet value = 0.03758 PEM electrolyzer stack CAPEX + Balance of Plant emissions (kg CO2e/kg H2)
        alk_ely_stack_capex_EI = (greet2['Electrolyzers']['O257'].value * g_to_kg)                  # Alkaline electrolyzer stack CAPEX emissions (kg CO2e/kg H2)
        alk_ely_stack_and_BoP_capex_EI = (greet2['Electrolyzers']['R257'].value * g_to_kg)          # Alkaline electrolyzer stack CAPEX + Balance of Plant emissions (kg CO2e/kg H2)
        soec_ely_stack_capex_EI = (greet2['Electrolyzers']['C257'].value * g_to_kg)                 # SOEC electrolyzer stack CAPEX emissions (kg CO2e/kg H2)
        soec_ely_stack_and_BoP_capex_EI = (greet2['Electrolyzers']['F257'].value * g_to_kg)         # SOEC electrolyzer stack CAPEX + Balance of Plant emissions (kg CO2e/kg H2)
        # Steel
        steel_CH4_prod = (greet2['Steel']['Y123'].value * CH4_gwp_to_CO2e * g_to_kg / ton_to_MT)                    # NOTE: original value = 39.29, greet value = 67.21, CH4 emissions for DRI-EAF Steel production w/ 83% H2 and 0% scrap (kg CO2e/metric tonne annual steel lab production)
        steel_CO2_prod = (greet2['Steel']['Y125'].value * g_to_kg / ton_to_MT)                                      # NOTE: original value = 174.66, greet value = 1043.87, CO2 emissions for DRI-EAF Steel production w/ 83% H2 and 0% scrap (kg CO2e/metric tonne annual steel lab production)
        steel_NG_supply_EI = ((greet2['Steel']['B260'].value +                                                      # NOTE: original value = 13, greet value = 12.67, these are upstream emissions of NG, not exact emissions of NG used in steel process, Upstream Natural Gas emissions for DRI-EAF Steel production (g CO2e/MJ)
                                (greet2['Steel']['B250'].value * VOC_to_CO2e) +                  
                                (greet2['Steel']['B251'].value * CO_to_CO2e) + 
                                (greet2['Steel']['B258'].value * CH4_gwp_to_CO2e) + 
                                (greet2['Steel']['B259'].value * N2O_gwp_to_CO2e)
                                ) / mmbtu_to_MJ
                                )
        steel_iron_ore_EI = ((greet2['Steel']['B92'].value +                                                        # NOTE: original value = 0.048, greet value = 0.045 Iron ore production emissions for use in DRI-EAF Steel production (kg CO2e/kg iron ore)
                                (greet2['Steel']['B82'].value * VOC_to_CO2e) +                     
                                (greet2['Steel']['B83'].value * CO_to_CO2e) +
                                (greet2['Steel']['B90'].value * CH4_gwp_to_CO2e) + 
                                (greet2['Steel']['B91'].value * N2O_gwp_to_CO2e)
                                ) * (g_to_kg / ton_to_kg)
                            )
        # TODO: confirm / validate we can pull these values from Greet with Masha (compare original values with GREET values)
        steel_H2O_EI = (MMBTU_NG_to_kg_CO2e / greet2['Steel']['B249'].value)                                        # TODO: check conversion, original value = 0.00013, greet value = 15.33, Water consumption emissions for use in DRI-EAF Steel production (kg CO2e/gal H20)      
        steel_H2O_consume = (greet2['Steel']['Y113'].value * (gal_H2O_to_MT/ton_to_MT))                             # NOTE: original value = 0.8037, greet value = 6.296, H2O consumption for DRI-EAF Steel production w/ 83% H2 and 0% scrap (metric tonne H2O/metric tonne steel production)
        steel_H2_consume = (greet2['Steel']['AK66'].value * (mmbtu_to_MJ/MJLHV_per_kg_h2) * (kg_to_MT/ton_to_MT))   # NOTE: original value = 0.06596, greet value = 0.08184, Hydrogen consumption for DRI-EAF Steel production w/ 83% H2 regardless of scrap (metric tonnes H2/metric tonne steel production)
        steel_NG_consume = ((greet2['Steel']['AE63'].value +                                                        # NOTE: original value = 0.71657, greet value = 4.5729, Natural gas consumption for DRI-EAF Steel production (GJ/ton steel)
                                greet2['Steel']['AG63'].value +                                    
                                greet2['Steel']['AK63'].value + 
                                greet2['Steel']['AM63'].value
                                ) * (mmbtu_to_GJ / ton_to_MT)
                            )
        steel_lime_consume = (greet2['Steel']['AM68'].value)                                                        # NOTE: original value = 0.01812, greet value = 0.01269, Lime consumption for DRI-EAF Steel production (metric tonne lime/metric tonne steel production)
        steel_iron_ore_consume = (greet2['Steel']['AM69'].value)                                                    # NOTE: original value = 1.629, greet value = 1.82333, Iron ore consumption for DRI-EAF Steel production (metric tonne iron ore/metric tonne steel production)
        steel_PO_consume = (greet2['Steel']['Y108'].value * (mmbtu_to_MWh/ton_to_MT))                               # NOTE: original value = 0.5502, greet value = 12.259, Total Energy consumption for DRI-EAF Steel production w/ 83% H2 and 0% scrap (MWh/metric tonne steel production)

        greet2.close()

        #------------------------------------------------------------------------------
        # Steam methane reforming (SMR) and Autothermal Reforming (ATR) - Incumbent H2 production processes
        #------------------------------------------------------------------------------
        #TODO: Validate new greet values compared to original values used with Masha
        #TODO; Confirm proper use of HHV vs LHV with Masha
        # Values with Carbon Capture Sequestration (CCS)
        greet1 = openpyxl.load_workbook(greet1_ccs_central_h2, data_only=True)
        smr_ccs_NG_consume = (greet1['Hydrogen']['C392'].value * btu_to_MJ * mmbtulhv_per_kg_h2)            # NOTE: original value = 177.1, greet value = 81.93 SMR, w/ CCS Well to Gate (WTG) Natural Gas (NG) consumption (MJ-LHV/kg H2)
        smr_ccs_RNG_consume = (greet1['Hydrogen']['L392'].value * btu_to_MJ * mmbtulhv_per_kg_h2)           # NOTE: greet value = 20.451, SMR w/ CCS WTG Renewable Natural Gas (RNG) consumption (MJ-LHV/kg H2)
        smr_ccs_NG_PO_consume = (greet1['Hydrogen']['C389'].value * btu_to_kWh * mmbtulhv_per_kg_h2)        # NOTE: original value = 1.5, greet value = 24.5, SMR via NG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        smr_ccs_RNG_PO_consume = (greet1['Hydrogen']['L389'].value * btu_to_kWh * mmbtulhv_per_kg_h2)       # NOTE: greet value = -39.92, SMR via RNG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        smr_ccs_steam_prod = (greet1['Inputs']['I1105'].value * btu_to_MJ * mmbtulhv_per_kg_h2)             # NOTE: value = 0, no steam exported w/ CCS? SMR Steam exported w/ CCS (MJ/kg H2)
        smr_ccs_perc_capture = (greet1['Hydrogen']['B11'].value)                                            # CCS rate for SMR (%)
        atr_ccs_NG_consume = (greet1['Hydrogen']['N392'].value * btu_to_MJ * mmbtulhv_per_kg_h2)            # NOTE: greet value = 200.1 ATR w/ CCS WTG NG consumption (MJ-LHV/kg H2)
        atr_ccs_RNG_consume = (greet1['Hydrogen']['O392'].value * btu_to_MJ * mmbtulhv_per_kg_h2)           # NOTE: greet value = 26.69 ATR w/ CCS WTG RNG consumption (MJ-LHV/kg H2)
        atr_ccs_NG_PO_consume = (greet1['Hydrogen']['N389'].value * btu_to_kWh * mmbtulhv_per_kg_h2)        # NOTE: greet value = 59.55 ATR via NG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        atr_ccs_RNG_PO_consume = (greet1['Hydrogen']['O389'].value * btu_to_kWh * mmbtulhv_per_kg_h2)       # NOTE: greet value = -33.78 ATR via RNG w/ CCS WTG Total Energy consumption (kWh/kg H2)
        atr_ccs_perc_capture = (greet1['Hydrogen']['B15'].value)                                            # CCS rate for Autothermal Reforming (%)
        greet1.close()

        # Values without CCS
        greet1 = openpyxl.load_workbook(greet1_no_ccs_central_h2, data_only=True)
        smr_NG_consume = (greet1['Hydrogen']['C392'].value * btu_to_MJ * mmbtuhhv_per_kg_h2)                # NOTE: original value = 167, greet value = 35.15, SMR w/out CCS WTG NG consumption (MJ-HHV/kg H2)
        smr_RNG_consume = (greet1['Hydrogen']['L392'].value * btu_to_MJ * mmbtuhhv_per_kg_h2)               # NOTE: greet value = -24.37, SMR w/out CCS WTG RNG consumption (MJ-HHV/kg H2)
        smr_NG_PO_consume = (greet1['Hydrogen']['C389'].value * btu_to_kWh * mmbtuhhv_per_kg_h2)            # NOTE: original value = 0.13, greet value = 10.02, SMR via NG w/out CCS WTG Total Energy consumption (kWh/kg H2)
        smr_RNG_PO_consume = (greet1['Hydrogen']['L389'].value * btu_to_kWh * mmbtuhhv_per_kg_h2)           # NOTE: greet value = -58.83, SMR via RNG w/out CCS WTG Total Energy consumption (kWh/kg H2)
        smr_steam_prod = (greet1['Inputs']['I1105'].value * btu_to_MJ * mmbtulhv_per_kg_h2)                 # NOTE: original value = 25.6, greet value = 25.66 SMR Steam exported w/out CCS (MJ/kg H2)
        atr_NG_consume = (greet1['Hydrogen']['N392'].value * btu_to_MJ * mmbtuhhv_per_kg_h2)                # NOTE: greet value = 235.21, same value as w/ CCS enabled, no pathway for ATR w/out CCS, ATR w/out CCS WTG NG consumption (MJ-HHV/kg H2)
        atr_RNG_consume = (greet1['Hydrogen']['O392'].value * btu_to_MJ * mmbtuhhv_per_kg_h2)               # NOTE: greet value = 31.37, same value as w/ CCS enabled, no pathway for ATR w/out CCS, ATR w/out CCS WTG RNG consumption (MJ-HHV/kg H2)
        atr_NG_PO_consume = (greet1['Hydrogen']['N389'].value * btu_to_kWh * mmbtuhhv_per_kg_h2)            # NOTE: greet value = 70.0, same value as w/ CCS enabled, no pathway for ATR w/out CCS, ATR via NG w/out CCS WTG Total Energy consumption (kWh/kg H2)
        atr_RNG_PO_consume = (greet1['Hydrogen']['O389'].value * btu_to_kWh * mmbtuhhv_per_kg_h2)           # NOTE: greet value = -39.72, same value as w/ CCS enabled, no pathway for ATR w/out CCS # ATR via RNG w/out CCS WTG Total Energy consumption (kWh/kg H2)
        greet1.close()

        #------------------------------------------------------------------------------
        # Ammonia (NH3)
        #------------------------------------------------------------------------------
        #NOTE: Hydrogen consumption for conventional / blue NH3 production not explicitly given in GREET. 
            #  However, value is constant for NH3 from coal/poplar gasification, PEM electrolysis, and GREEN Ammonia ~= 0.197 (kg h2/kg NH3)
            #  Thus, using this value for all NH3 production pathways
        # Values for Green NH3
        greet1 = openpyxl.load_workbook(greet1_green_NH3, data_only=True)
        green_NH3_H2_consume = (greet1['Ag_Inputs']['AM54'].value)                                  # NOTE: original value = 0.2, greet value = 0.197, Green Ammonia production Hydrogen consumption (kg H2/kg NH3)
        green_NH3_PO_consume = (greet1['Hydrogen']['B359'].value * btu_to_kWh)                      # NOTE: original value = 0.1207, greet value = 11.87, Green Ammonia production Total Energy consumption (kWh/kg NH3)
        greet1.close()

        # Values for Blue NH3
        greet1 = openpyxl.load_workbook(greet1_blue_NH3, data_only=True)
        blue_NH3_H2_consume = (greet1['Ag_Inputs']['AM54'].value)                                   # Blue Ammonia production Hydrogen consumption (kg H2/kg NH3)
        blue_NH3_PO_consume = (greet1['Hydrogen']['B359'].value * btu_to_kWh)                       # NOTE: original value = 0.1207, greet value = 10.51, Blue Ammonia production Total Energy consumption (kWh/kg NH3)
        greet1.close()

        # Values for Conventional NH3
        greet1 = openpyxl.load_workbook(greet1_conventional_NH3, data_only=True)
        conventional_NH3_H2_consume = (greet1['Ag_Inputs']['AM54'].value)                           # Conventional Ammonia production Hydrogen consumption (kg H2/kg NH3)
        conventional_NH3_PO_consume = (greet1['Hydrogen']['B359'].value * btu_to_kWh)               # NOTE: original value = 0.1207, greet value = 10.23, Conventional Ammonia production Total Energy consumption (kWh/kg NH3)
        greet1.close()

        #------------------------------------------------------------------------------
        # Steel
        #------------------------------------------------------------------------------
        # NOTE: Alternative DRI-EAF configurations (w/ and w/out scrap, H2 vs NG) found in greet2 > Steel > W107:Z136
                # Iron or vs scrap % controlled by B24:C24 values
        # NOTE: greet2 > Steel > B17:B18 controls NG vs RNG, changing these values drastically changes steel_NG_supply_EI
            # May require hosting of different steel config GREET versions if desired ^
        # Values for DRI-EAF 83% H2, 100% DRI 0% Scrap
        # greet1 = openpyxl.load_workbook(greet1_ccs_central_h2, data_only=True)
        # steel_lime_EI = ((greet1['Chemicals']['BA247'].value +                                                      # NOTE: original value = 1.28, greet value = 1.28, Lime production emissions for use in DRI-EAF Steel production (kg CO2e/kg lime)
        #                     (greet1['Chemicals']['BA237'].value * VOC_to_CO2e) +             
        #                     (greet1['Chemicals']['BA238'].value * CO_to_CO2e) +
        #                     (greet1['Chemicals']['BA245'].value * CH4_gwp_to_CO2e) +
        #                     (greet1['Chemicals']['BA246'].value * N2O_gwp_to_CO2e)
        #                     ) * g_to_kg * (1/ton_to_kg))
        # greet1.close()

        # greet2 = openpyxl.load_workbook(greet2_ccs_central_h2, data_only=True)
        # steel_CH4_prod = (greet2['Steel']['Y123'].value * CH4_gwp_to_CO2e * g_to_kg / ton_to_MT)                    # NOTE: original value = 39.29, greet value = 67.21, CH4 emissions for DRI-EAF Steel production w/ 83% H2 and 0% scrap (kg CO2e/metric tonne annual steel lab production)
        # steel_CO2_prod = (greet2['Steel']['Y125'].value * g_to_kg / ton_to_MT)                                      # NOTE: original value = 174.66, greet value = 1043.87, CO2 emissions for DRI-EAF Steel production w/ 83% H2 and 0% scrap (kg CO2e/metric tonne annual steel lab production)
        # steel_NG_supply_EI = ((greet2['Steel']['B260'].value +                                                      # NOTE: original value = 13, greet value = 12.67, these are upstream emissions of NG, not exact emissions of NG used in steel process, Upstream Natural Gas emissions for DRI-EAF Steel production (g CO2e/MJ)
        #                         (greet2['Steel']['B250'].value * VOC_to_CO2e) +                  
        #                         (greet2['Steel']['B251'].value * CO_to_CO2e) + 
        #                         (greet2['Steel']['B258'].value * CH4_gwp_to_CO2e) + 
        #                         (greet2['Steel']['B259'].value * N2O_gwp_to_CO2e)
        #                         ) / mmbtu_to_MJ
        #                         )
        # steel_iron_ore_EI = ((greet2['Steel']['B92'].value +                                                        # NOTE: original value = 0.048, greet value = 0.045 Iron ore production emissions for use in DRI-EAF Steel production (kg CO2e/kg iron ore)
        #                         (greet2['Steel']['B82'].value * VOC_to_CO2e) +                     
        #                         (greet2['Steel']['B83'].value * CO_to_CO2e) +
        #                         (greet2['Steel']['B90'].value * CH4_gwp_to_CO2e) + 
        #                         (greet2['Steel']['B91'].value * N2O_gwp_to_CO2e)
        #                         ) * (g_to_kg / ton_to_kg)
        #                     )
        # # TODO: confirm / validate we can pull these values from Greet with Masha (compare original values with GREET values)
        # steel_H2O_EI = (MMBTU_NG_to_kg_CO2e / greet2['Steel']['B249'].value)                                        # TODO: check conversion, original value = 0.00013, greet value = 15.33, Water consumption emissions for use in DRI-EAF Steel production (kg CO2e/gal H20)      
        # steel_H2O_consume = (greet2['Steel']['Y113'].value * (gal_H2O_to_MT/ton_to_MT))                             # NOTE: original value = 0.8037, greet value = 6.296, H2O consumption for DRI-EAF Steel production w/ 83% H2 and 0% scrap (metric tonne H2O/metric tonne steel production)
        # steel_H2_consume = (greet2['Steel']['AK66'].value * (mmbtu_to_MJ/MJLHV_per_kg_h2) * (kg_to_MT/ton_to_MT))   # NOTE: original value = 0.06596, greet value = 0.08184, Hydrogen consumption for DRI-EAF Steel production w/ 83% H2 regardless of scrap (metric tonnes H2/metric tonne steel production)
        # steel_NG_consume = ((greet2['Steel']['AE63'].value +                                                        # NOTE: original value = 0.71657, greet value = 4.5729, Natural gas consumption for DRI-EAF Steel production (GJ/ton steel)
        #                         greet2['Steel']['AG63'].value +                                    
        #                         greet2['Steel']['AK63'].value + 
        #                         greet2['Steel']['AM63'].value
        #                         ) * (mmbtu_to_GJ / ton_to_MT)
        #                     )
        # steel_lime_consume = (greet2['Steel']['AM68'].value)                                                        # NOTE: original value = 0.01812, greet value = 0.01269, Lime consumption for DRI-EAF Steel production (metric tonne lime/metric tonne steel production)
        # steel_iron_ore_consume = (greet2['Steel']['AM69'].value)                                                    # NOTE: original value = 1.629, greet value = 1.82333, Iron ore consumption for DRI-EAF Steel production (metric tonne iron ore/metric tonne steel production)
        # steel_PO_consume = (greet2['Steel']['Y108'].value * (mmbtu_to_MWh/ton_to_MT))                               # NOTE: original value = 0.5502, greet value = 12.259, Total Energy consumption for DRI-EAF Steel production w/ 83% H2 and 0% scrap (MWh/metric tonne steel production)
        # greet2.close()

        data_dict = {
                     # Hardcoded values
                     'NH3_boiler_EI':NH3_boiler_EI,
                     'smr_NG_combust':smr_NG_combust,
                     'smr_HEX_eff':smr_HEX_eff,
                     'smr_NG_supply':smr_NG_supply,
                     'ccs_PO_consume':ccs_PO_consume,
                     'atr_steam_prod':atr_steam_prod,
                     'atr_ccs_steam_prod':atr_ccs_steam_prod,
                     'grid_trans_losses':grid_trans_losses,
                     'fuel_to_grid_curr':fuel_to_grid_curr,
                     'fuel_to_grid_futu':fuel_to_grid_futu,
                     # Renewable infrastructure embedded EI and h2 production via water electrolysis
                     'wind_capex_EI':wind_capex_EI,
                     'solar_pv_capex_EI':solar_pv_capex_EI,
                     'nuclear_PWR_capex_EI':nuclear_PWR_capex_EI,
                     'nuclear_BWR_capex_EI':nuclear_BWR_capex_EI,
                     'coal_capex_EI':coal_capex_EI,
                     'gas_capex_EI':gas_capex_EI,
                     'hydro_capex_EI':hydro_capex_EI,
                     'bio_capex_EI':bio_capex_EI,
                     'geothermal_EGS_capex_EI':geothermal_EGS_capex_EI,
                     'geothermal_flash_capex_EI':geothermal_flash_capex_EI,
                     'geothermal_binary_capex_EI':geothermal_binary_capex_EI,
                     'pem_ely_PO_consume':pem_ely_PO_consume,
                     'pem_ely_stack_capex_EI':pem_ely_stack_capex_EI,
                     'pem_ely_stack_and_BoP_capex_EI':pem_ely_stack_and_BoP_capex_EI,
                     'alk_ely_stack_capex_EI':alk_ely_stack_capex_EI,
                     'alk_ely_stack_and_BoP_capex_EI':alk_ely_stack_and_BoP_capex_EI,
                     'soec_ely_stack_capex_EI':soec_ely_stack_capex_EI,
                     'soec_ely_stack_and_BoP_capex_EI':soec_ely_stack_and_BoP_capex_EI,
                     'battery_LFP_residential_EI':battery_LFP_residential_EI,
                     'battery_LFP_commercial_EI':battery_LFP_commercial_EI,
                     # Steam methane reforming (SMR) and Autothermal Reforming (ATR)
                     'smr_ccs_NG_consume':smr_ccs_NG_consume,
                     'smr_ccs_RNG_consume':smr_ccs_RNG_consume,
                     'smr_ccs_NG_PO_consume':smr_ccs_NG_PO_consume,
                     'smr_ccs_RNG_PO_consume':smr_ccs_RNG_PO_consume,
                     'smr_ccs_steam_prod':smr_ccs_steam_prod,
                     'smr_ccs_perc_capture':smr_ccs_perc_capture,
                     'atr_ccs_NG_consume':atr_ccs_NG_consume,
                     'atr_ccs_RNG_consume':atr_ccs_RNG_consume,
                     'atr_ccs_NG_PO_consume':atr_ccs_NG_PO_consume,
                     'atr_ccs_RNG_PO_consume':atr_ccs_RNG_PO_consume,
                     'atr_ccs_perc_capture':atr_ccs_perc_capture,
                     'smr_NG_consume':smr_NG_consume,
                     'smr_RNG_consume':smr_RNG_consume,
                     'smr_NG_PO_consume':smr_NG_PO_consume,
                     'smr_RNG_PO_consume':smr_RNG_PO_consume,
                     'smr_steam_prod':smr_steam_prod,
                     'atr_NG_consume':atr_NG_consume,
                     'atr_RNG_consume':atr_RNG_consume,
                     'atr_NG_PO_consume':atr_NG_PO_consume,
                     'atr_RNG_PO_consume':atr_RNG_PO_consume,
                     # Ammonia (NH3)
                     'green_NH3_H2_consume':green_NH3_H2_consume,
                     'green_NH3_PO_consume':green_NH3_PO_consume,
                     'blue_NH3_H2_consume':blue_NH3_H2_consume,
                     'blue_NH3_PO_consume':blue_NH3_PO_consume,
                     'conventional_NH3_H2_consume':conventional_NH3_H2_consume,
                     'conventional_NH3_PO_consume':conventional_NH3_PO_consume,
                     # Steel
                     'steel_lime_EI':steel_lime_EI,
                     'steel_CH4_prod':steel_CH4_prod,
                     'steel_CO2_prod':steel_CO2_prod,
                     'steel_NG_supply_EI':steel_NG_supply_EI,
                     'steel_iron_ore_EI':steel_iron_ore_EI,
                     'steel_H2O_EI':steel_H2O_EI,
                     'steel_H2O_consume':steel_H2O_consume,
                     'steel_H2_consume':steel_H2_consume,
                     'steel_NG_consume':steel_NG_consume,
                     'steel_lime_consume':steel_lime_consume,
                     'steel_iron_ore_consume':steel_iron_ore_consume,
                     'steel_PO_consume':steel_PO_consume,
                    }
        
        # Dump data to yaml file
        yaml_file = open(self.filename, mode="w+")
        yaml.dump(data_dict, yaml_file, default_flow_style=False)
        yaml_file.close()

    def format_data(self):
        """
        """
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"{self.filename} does not exist. Try `download_resource` first.")

        yaml_file = open(self.filename, mode='r')
        self.data = yaml.load(yaml_file, Loader=yaml.SafeLoader)
        yaml_file.close()

#Adhoc testing
if __name__ == '__main__':
    test = GREETData()
#     print(test.data)

#NOTE: Runtime ~ 1m16s to fully parse greet
#NOTE: Runtime <2s with greet preprocessed (load yaml to dict)
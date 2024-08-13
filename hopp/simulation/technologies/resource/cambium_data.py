import os
import csv
from pathlib import Path
from typing import Union
import numpy as np

from hopp.utilities.keys import get_developer_nrel_gov_key, get_developer_nrel_gov_email
from hopp.utilities.log import hybrid_logger as logger
from hopp.simulation.technologies.resource.resource import Resource
from hopp import ROOT_DIR

CAMBIUM_BASE_URL = "https://scenarioviewer.nrel.gov/api/get-data-cache/"

class CambiumData(Resource):
    """
    Class to manage Cambium emissions and grid mix generation data

    Args:
        lat: (float) latitude
        lon: (float) longitude
        year: (int) year
        project_uuid: (string) uuid of the cambium project (Cambium 2022 vs 2023), can be found in the scenarioviewer.nrel.gov URL after selecting the project in the viewer
            Default: Cambium 2023 == '0f92fe57-3365-428a-8fe8-0afc326b3b43'
            Available: 
                Cambium 2022 == '82460f06-548c-4954-b2d9-b84ba92d63e2'
                Cambium 2023 == '0f92fe57-3365-428a-8fe8-0afc326b3b43'
        scenario: (string) scenario name to query as it appears in the Cambium Scenario Viewer
            Default: 'Mid-case with 100% decarbonization by 2035'
            Available: 
                'High demand growth', 'High natural gas prices', 'High renewable energy cost', 'Low natural gas prices', 'Low renewable energy cost', 'Mid-case', 'Mid-case with 100% decarbonization by 2035', 'Mid-case with 95% decarbonization by 2050'
    #       Documentation: See pdf for additional information on each scenario https://www.nrel.gov/docs/fy24osti/88507.pdf
        location_type: (string) geographic resolution of cambium emissions and grid generation data
            Default: 'GEA Regions 2023'
            Available:
                'GEA Regions 2023', 'Nations'
        time_type: (string) time resolution of data
            Default: 'hourly'
            Available:
                'hourly' == 8760 array for entire queried year
                'annual' == single annual value
        path_resource: directory where to save downloaded files
        filepath: file path of resource file to load
        use_api: Make an API call even if there is an existing data file. Default == False
        kwargs: additiona keyword arguments

    """
    filename_map = {
        'project_uuid': {'82460f06-548c-4954-b2d9-b84ba92d63e2':'Cambium22',
                         '0f92fe57-3365-428a-8fe8-0afc326b3b43':'Cambium23'},
        'scenario': {'High demand growth':'HighDemandGrowth',
                     'High natural gas prices':'HighNGPrices',
                     'High renewable energy cost':'HighRenewableCost',
                     'Low natural gas prices':'LowNGPrices',
                     'Low renewable energy cost':'LowRenewableCost',
                     'Mid-case':'MidCase',
                     'Mid-case with 100% decarbonization by 2035':'MidCase100by2035', 
                     'Mid-case with 95% decarbonization by 2050':'MidCase95by2050'},
        'location_type': {'GEA Regions 2023':'GEA',
                          'Nations':'Nation'},
    }
    def __init__(
        self,
        lat: float,
        lon: float,
        year: int,
        project_uuid: str = '0f92fe57-3365-428a-8fe8-0afc326b3b43',
        scenario: str = 'Mid-case with 100% decarbonization by 2035',
        location_type: str = 'GEA Regions 2023',
        time_type: str = 'hourly',
        path_resource: Union[str, Path] = ROOT_DIR / "simulation" / "resource_files",
        filepath: Union[str, Path] ="",
        use_api: bool = False,
        **kwargs
    ):
        # Run init of Resource super class
        super().__init__(lat,lon,year)

        # Check if path_resource is a directory, if yes define as self.path_resource attribute
        if os.path.isdir(path_resource):
            self.path_resource = path_resource
        
        # update attribute with cambium directory
        self.path_resource = os.path.join(self.path_resource, 'cambium')

        # Force override internal definitions if passed in
        self.__dict__.update(kwargs)

        # Define the filepath and file name for the resource file
        if filepath == "":
            filepath = os.path.join(self.path_resource, 
                                    str(lat) + "_" + str(lon) + "_" + str(self.filename_map['project_uuid'][project_uuid]) + "_" +
                                    str(self.filename_map['scenario'][scenario]) + "_" + str(time_type) + "_" +
                                    str(self.filename_map['location_type'][location_type]) + "_" + str(year) + ".csv")
        self.filename = filepath

        # Check if the download directory exists (hopp/simulation/resource_files/cambium), if not make the directory
        self.check_download_dir()

        # If the resource file does not exist in directory or use_api flag == True, download the data
        if not os.path.isfile(self.filename) or use_api:
            self.download_resource()

        # TODO: verify purpose / functionality in wind_resource.py and solar_resource.py
        self.format_data()

        logger.info("CambiumData: {}".format(self.filename))

    def download_resource(self):
        pass
        #TODO: update logic with API calls

    def format_data(self):
        """
        """
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"{self.filename} does not exist. Try `download_resource` first.")

        self.data = self.filename
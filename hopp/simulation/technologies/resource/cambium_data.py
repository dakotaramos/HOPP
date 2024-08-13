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
        lat: latitude
        lon: longitude
        year: year
        path_resource: directory where to save downloaded files
        filepath: file path of resource file to load
        use_api: Make an API call even if there is an existing data file. Default == False
        kwargs: additiona keyword arguments
        
    """
import sys 
import os

from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pyreadr ## to read R objcets




os.chdir("/home/cwijes1/projects/jgi_rna_seq/jgi_jupyter_notebooks/jgi_full_atlas_11_22_2021/jgi_full_atlas_with_sp_nacl/data/wgcna_analysis/con_tre_multi_ion_removed_wgcna/")
sys.path.insert(0, '/home/cwijes1/projects/jgi_rna_seq/scripts/tau')
sys.path.insert(0, '/home/cwijes1/projects/jgi_rna_seq/scripts/network_analysis_utils/')

# import custom functions
import calculate_tau as filter_functions
import network_utils
import wgcna_analysis_helper_functions as wgcna_functions
#from helper_functions import run_pca
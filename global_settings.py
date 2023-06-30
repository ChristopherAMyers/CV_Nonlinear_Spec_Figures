import matplotlib.pyplot as plt
import os
import platform

script_dir = os.path.abspath(os.path.dirname(__file__))
plt.style.use(script_dir + '/style.mplstyle')

data_root_dir = '/Users/cmyers/Library/CloudStorage/OneDrive-UniversityofCaliforniaMerced/crysel_violet/2DES/shao-yu_replication/'
if platform.system().lower() == 'windows':
    data_root_dir = 'C:\\Users\\chris\\OneDrive - University of California Merced\crysel_violet\\2DES\\shao-yu_replication\\'

os.makedirs('png', exist_ok=True)
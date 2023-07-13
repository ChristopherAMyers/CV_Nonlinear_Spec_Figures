import matplotlib.pyplot as plt
import os
from os.path import join
import platform

#   load the plotting style in the same directory as this script
script_dir = os.path.abspath(os.path.dirname(__file__))
plt.style.use(script_dir + '/style.mplstyle')

#   folder to hold plot images
os.makedirs('png', exist_ok=True)

#   main directory that holds spectra calculations
data_root_dir = '/Users/cmyers/Library/CloudStorage/OneDrive-UniversityofCaliforniaMerced/crysel_violet/2DES/shao-yu_replication/'
if platform.system().lower() == 'windows':
    data_root_dir = 'C:\\Users\\chris\\OneDrive - University of California Merced\crysel_violet\\2DES\\shao-yu_replication\\'



#   Each 'job_dir' must contain the vee_traj1.dat file and 
#   vee__2nd_order_cumulant_transient_absorption_spec.txt file.
#
#   If Only using the SE decay and not the difference between SE and GSB,
#   then set 'job_dir' equal to 'job_dir_SE'.
job_info = {
    'stripped':
    {
        'job_dir': join(data_root_dir, 'gs-aimd-longer', 'stripped'),
        'job_dir_SE': join(data_root_dir, 'gs-aimd-longer', 'stripped_SE'),
        'plt_args': {
            'label': 'Stripped', 
            'color': '#377E22'
            }
    },
    'mm_4hb':
    {
        'job_dir': join(data_root_dir, 'gs-aimd-longer', 'mm_4hb'),
        'job_dir_SE': join(data_root_dir, 'gs-aimd-longer', 'mm_4hb_SE'),
        'plt_args': {
            'label': '4-Peripheral', 
            'color': '#21ADEF'
            }
    },
    'mm_C4':
    {
        'job_dir': join(data_root_dir, 'gs-aimd-longer', 'mm_C4'),
        'job_dir_SE': join(data_root_dir, 'gs-aimd-longer', 'mm_C4_SE'),
        'plt_args': {
            'label': '4-Axial', 
            'color': '#D321FF'
            }
    },
    'mm_qm1':
    {
        'job_dir': join(data_root_dir, 'gs-aimd-longer', 'mm_qm1'),
        'job_dir_SE': join(data_root_dir, 'gs-aimd-longer', 'mm_qm1_SE'),
        'plt_args': {
            'label': 'Full MM', 
            'color': 'red'
            }
    },
    'qm2':
    {
        'job_dir': join(data_root_dir, 'gs-aimd-longer', 'qm2'),
        'job_dir_SE': join(data_root_dir, 'gs-aimd-longer', 'qm2_SE'),
        'plt_args': {
            'label': 'QM + MM', 
            'color': 'blue'
            }
    },
}

from voidfindertk import io
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from voidfindertk.analysis_tools import analysis_voids,join_box_void
filename='./clean_popvds.dat'
void = io.read_popcorn_output(filename)
box = io.read_table('halos_fede.dat') 

void_out = join_box_void(box=box,voids=void)
output = analysis_voids(**void_out,n_items=100)
print(output.keys())

sns.histplot(data=pd.DataFrame(output), x='velocity_norm', bins=50, kde=True).set(title=f'nvoids:{len(output['velocity_norm'])}')
#sns.histplot(data=velocity_voids_df, x='n_tracers', bins=50, kde=True).set(title=f'nvoids:{len(velocity_voids)}')

# #plt.savefig('popcorn_velocity_hist.jpg')
# #plt.savefig('popcorn_tracers_hist.jpg')
plt.show()

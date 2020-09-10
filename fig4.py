


"""

Figure 4: Geometry of the AIS in RGC.

"""

from brian2 import *
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy import stats
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import gridspec
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import seaborn as sns

### Figure parameters
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Laoding the results of analyses
df_cells = pd.read_excel('RGC_electrical_properties.xlsx')

dates = array(df_cells['Date'])[:-3]
retinas = array(df_cells['Retina'])[:-3]
cells = array(df_cells['Cell'])[:-3]
ages = array(df_cells['Age'])[:-3]

ais_starts = array(df_cells['AIS start (um)'])[:-3]
ais_lengths = array(df_cells['AIS length (um)'])[:-3]

# remove nans
dates = dates[~isnan(ais_lengths)]
ais_starts = ais_starts[~isnan(ais_lengths)]
ais_lengths = ais_lengths[~isnan(ais_lengths)]

### Load the images of an example cell
# path_images = 'data/RGC data/images/fig4/' 
path_images = '/Users/sarah/Documents/data/Martijn Sierksma/images/20191114 P10/Retina C/cell 1/' 


### Statistics
print ('Mean Delta:', mean(ais_starts), '+-', std(ais_starts))
print ('Mean L:', mean(ais_lengths), '+-', std(ais_lengths))
print ('N cells:', len(ais_starts))

### Figure
fig = figure('AIS position and axonal currents', figsize=(11, 4))
gs = gridspec.GridSpec(1, 2, width_ratios=[2,1]) 
gs1 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[0], wspace=0)

# Panel A: labellings in an exmaple cell

# biocytin
ax1 = fig.add_subplot(gs1[0])
img1 = mpimg.imread(path_images + 'fig_biocytin.jpg')
ax1.imshow(img1)
ax1.axis('off')
ax1.set_title('biocytin', color='m')
ax1.annotate("A", xy=(-0.2,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=16, weight='bold')

# ankyrin G
ax2 = fig.add_subplot(gs1[1], sharey=ax1)
img2 = mpimg.imread(path_images + 'fig_ankG.jpg')
ax2.imshow(img2, cmap='gray')
ax2.axis('off')
ax2.set_title('ankyrin G', color='g')

# masked anyrin G
ax3 = fig.add_subplot(gs1[2], sharey=ax1)
img3 = mpimg.imread(path_images + 'fig_masked_ankG.jpg')
ax3.imshow(img3, cmap='gray')
ax3.axis('off')
ax3.set_title('ankyrin G', color='g')

# merge
ax4 = fig.add_subplot(gs1[3], sharey=ax1)
img4 = mpimg.imread(path_images + 'fig_merge.jpg')
ax4.imshow(img4, cmap='gray')
ax4.axis('off')
ax4.set_title('merge')

# Panel B: AIS start position vs length
gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[1], wspace=0.5)
ax5 = fig.add_subplot(gs2[0])

# Regression
slope_dl, intercept_dl, r_dl, p_dl, _ = stats.linregress(ais_lengths, ais_starts)
print ('Regression Delta-L', 'r=', r_dl, 'p=', p_dl)
print ('Pearson Delta-L', stats.pearsonr(ais_lengths, ais_starts))

sns.scatterplot(x=ais_lengths, y=ais_starts, color='k', ax=ax5)
ax5.plot(linspace(20, 45, 50), intercept_dl + slope_dl * linspace(20, 45, 50), 'k-')
ax5.set_ylim(0,20)
ax5.set_xlim(10,50)
ax5.set_xlabel('$L$ ($\mu$m)')
ax5.set_ylabel('$\Delta$ ($\mu$m)')
ax5.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=16, weight='bold')

tight_layout()

### Saving the figure
save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'
# fig.savefig(save_path + "fig4.pdf", bbox_inches='tight')

# fig.savefig(save_path + "fig4.png", dpi=300)












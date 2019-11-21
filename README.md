# Nipah_phylogenetics
Collection of scripts used for "Inference of Nipah virus Evolution, 1999-2015" 


## Enrichment Probes

How does one make enrichment probes?  It's surprisingly simple.  The user generates a fasta file containing sequences of interest.  A multiline or single line fasta can be used as input.  These sequences can represent only a limited subset of the available sequences, or perhaps representative sequences (maybe for each viral clade?).  Some limited data suggests that an enrichment probe set works better if there are many diverse sequences in your fasta.

The [enrichment probe script](/python_multi_seq_fasta_test_hybrid_oligos_V2.py) will slice a fasta sequence to generate enrichment probe sequences.  The probe size is controled by "j" and the spacing between the probes is controled by "i".  Currently the script will generate probes of 80 bp long with a spacing of 200 bp between the probes.  If you want to change the probe size (100 or 120 bp length?) just change the value of j.  If you want the probes to have no gaps between them, then change the value of i to j+1.

Also, this was written in Python 2.  If you want to use it with Python 3, you probably just need to change the "print" syntax.

Script usage:
```
python python_multi_seq_fasta_test_hybrid_oligos_V2.py > enrichment_probes.txt
```

## Enrichment probes used in "Inference of Nipah virus Evolution, 1999-2015"

File containing probe sequences can be found [here](/Nipah_oligos.xls)
Probes can be purchase from [Twist Biosciences](https://twistbioscience.com).  Please note that probes need to be 5' biotinylated.

My last conversation with the folks at Twist Biosciences indicated that they could make approximately 6,500 probes at a total probe concentration of 30ug.

## Assembly script

File for our in-house assembly pipeline include the [mapper](/map_and_dedup_Bangladesh_NEBNextV1.sh), [wrapper script](array_script_NEBNext_V1.sh), and [input_file](file_names.txt).  The script that does the heavy lifting for genome assembly is [here](/map_and_dedup_Bangladesh_NEBNextV1.sh).  The other scrips are wrapped scripts for input to CDC's High performance computing cluster (uses Sun Grid Engine).

Usage:
```
sh array_script_NEBNext_V1.sh
```

## Read Strandedness

This script was originally from from Jason [Ladner](https://github.com/jtladner/Scripts/blob/master/strand_specific_coverage/strandspec_covplot/strandspec_covplot_v1.1.py), but I made some edits and used this [verison](/strandspec_covplot_v1.2.py).  

Usage:
Please see Jason' [readme](https://github.com/jtladner/Scripts/tree/master/strand_specific_coverage/strandspec_covplot).  The usage is the same.

## Chimeric Reads

This script is from Jason Ladner.  Check out his github [page](https://github.com/jtladner/Scripts/tree/master/chimeric_reads) for original [script](https://github.com/jtladner/Scripts/blob/master/chimeric_reads/chimeric_reads_v3.6.2.py).

## Visualizing Geographic Spread from a BEAST Continuous Diffusion Analysis

Hold on to your hats, folks.  This script has a lot of moving parts.  The original geographic spread script was made by Gytis Dudas and can be found [here](https://github.com/evogytis/baltic/blob/master/curonia.ipynb).  Kudos to him for this excellent piece of work!  This script was used to make the really cool Ebola virus spread video found [here](https://www.youtube.com/watch?v=j4Ut4krp8GQ).

######How does my modified curionia script work?
I had to break it down into several steps to understand and get the program working.
1. Load a beast tree and geojson for your region of interest.  Are the geographic polygons correctly loaded?

```
#This code is from curonia.py from Gytis Dudas but modified by Shannon Whitmer.
#module load Python/3.6.1
#jupyter notebook

#Point to baltic.
import imp
bt = imp.load_source('baltic','/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/baltic.py') ## point to where baltic repo was cloned

#import libraries
%matplotlib inline
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon ## for polygons
from matplotlib.collections import PatchCollection ## for polygons too
from matplotlib.colors import LinearSegmentedColormap ## for colour maps
from matplotlib import gridspec ## for composite figures
import matplotlib.patheffects as path_effects ## for elegant text
from IPython.display import clear_output
from IPython.display import HTML

import datetime
import math
import time
import sys

import unicodedata
# import unidecode ## for removing diacritics from example geoJSON

import numpy as np
from scipy.interpolate import UnivariateSpline ## used to smooth counts of lineages in each location at any given time
from scipy.interpolate import interp1d ## used to linearly interpolate between data points used in colouring polygons
from sklearn.decomposition import IncrementalPCA ## used to identify PCA1 when automatically producing a colour map

#import bezier ## custom arbitrary order Bezier curves
import requests ## used to fetch examples from internet
import json ## used for importing JSONs
try:
    from StringIO import StringIO as sio
    from cStringIO import StringIO as csio
except ImportError:
    from io import StringIO as sio
    from io import BytesIO as csio
    
def removeDiacritics(string):
    """
    Removes diacritic marks from unicode.
    """
#    output=None
#    if isinstance(string, str):
#        output=string
#    elif isinstance(string, unicode):
#        output=string.encode('utf-8')
#        output=unidecode.unidecode(string)
#        nkfd_form = unicodedata.normalize('NFKD', unicode(string))
#        output= u"".join([c for c in nkfd_form if not unicodedata.combining(c)])
#        output = ''.join((c for c in unicodedata.normalize('NFD', string) if unicodedata.category(c) != 'Mn'))
#        output=unicodedata.normalize('NFKD', string).encode('ASCII', 'ignore')
#    return output
#    return unicodedata.normalize('NFKD', string).encode('ASCII', 'ignore')
    return unicodedata.normalize('NFKD', string)

def calendarTimeline(start_date, end_date, infmt='%Y-%m-%d',outfmt='%Y-%b',optfmt=None,month_step=12):
    """
    Given two calendar dates returns a list of calendar dates at monthly (by default) intervals.
    """
    current_date = datetime.datetime.strptime(start_date,infmt)
    ending_date = datetime.datetime.strptime(end_date,infmt)
    
    timeline=[]
    while current_date <= ending_date:
        if optfmt and current_date.month!=1:
            d=datetime.datetime.strftime(current_date,optfmt)
        else:
            d=datetime.datetime.strftime(current_date,outfmt)
        #print('Current Date:%s'%d)
        timeline.append(d)
        carry, new_month = divmod(current_date.month - 1 + month_step, 12)
        new_month += 1
        current_date = current_date.replace(year=current_date.year + carry,month=new_month)
    return timeline

typeface='Helvetica Neue' ## set default matplotlib font and font size
mpl.rcParams['font.weight']=300
mpl.rcParams['axes.labelweight']=300
mpl.rcParams['font.family']=typeface
mpl.rcParams['font.size']=22

#frame='<iframe style="border: 0; width: 400px; height: 472px;" src="https://bandcamp.com/EmbeddedPlayer/album=29809561/size=large/bgcol=333333/linkcol=e99708/artwork=small/transparent=true/" seamless><a href="http://romowerikoito.bandcamp.com/album/nawam-r">NAWAMAR by Romowe Rikoito</a></iframe>'

print('Done!')
#HTML(frame)

start='1971-09-01' ## start date of animation
end='2015-03-18' ## end date of animation

xtimeline=calendarTimeline(start,end,'%Y-%m-%d','%Y-%m-%d') ## create timeline from start to end of animation delimited with months


#import geoJSON:
json_map=json.load(open('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/location_variable_testing_2/GEOJSON/bgd_admin2.geojson','r')) ## read from (hopefully saved) local copy

print('Done!')

#Convert geoJSON into matplotlib ploygons:
features=json_map['features']
location_points={} ## location points will be stored here
polygons={} ## polygons will be stored here

locName='ADM2_EN' ## key name for each feature

for loc in features: ## iterate through features (locations)
    poly = np.asarray(loc['geometry']['coordinates']) ## get coordinates
    location=removeDiacritics(loc['properties'][locName]) ## standardised location name (remove diacritics)
#     print(location.encode().decode('utf-8'))
    if location not in ['Isla Sala y Gomez'] and 'Gal' not in location: ## ignore Isla Sala y Gomez
        polygons[location]=[]
        location_points[location]=[]
        if loc['geometry']['type']=='MultiPolygon': ## multiple parts detected
            for part in np.asarray(poly): ## iterate over each component polygon
                for coords in np.asarray(part): ## iterate over coordinates
                    coords=np.array(coords)
                    xs=coords[:,0] ## longitudes
                    ys=coords[:,1] ## latitudes

                    location_points[location].append(np.vstack(zip(xs,ys))) ## append coordinates to location's list of coordinates
        if loc['geometry']['type']=='Polygon': ## location is single part
            for coords in np.asarray(poly): ## iterate over coordinates
                coords=np.array(coords)
                xs=coords[:,0] ## longitudes
                ys=coords[:,1] ## latitudes

                location_points[location].append(np.vstack(zip(xs,ys))) ## append coordinates to location's list of coordinates

        complete_location=[]
        for part in location_points[location]: ## iterate over each component of a location
            complete_location.append(Polygon(part,True)) ## create a polygon for each component of a location

        polygons[location]=complete_location ## assign list of polygons to a location
#     elif location=='Isla Sala y Gomez': ## if location is Isla Sala y Gomez - print a geoJSON entry example
#         print('example geoJSON entry:\n%s\n\nnote that only the coordinate field is called\n'%(loc))
        
print('polygons loaded:\n%s'%(polygons.keys()))
```

2. Can i import my Cauchy MCC tree and pull coordiate data out of it?
```
#Can I import my Cauchy MCC tree and pull coordinate data out of it??

#Point to baltic.
import imp
bt = imp.load_source('baltic','/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/baltic.py') ## point to where baltic repo was cloned

#Test code to import my Cauchy MCC tree:
tree_path='/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/Nipah_Cauchy_new_HKYG_UCLN_skygrid_1234_log.time.MCC.trees'

ll=bt.loadNexus(tree_path)

#Import geojson
df=gpd.read_file('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/location_variable_testing_2/GEOJSON/SEAsia_forBEAST_V2.geojson')


# locTrait='location' ## name of locations in the tree
latTrait='location1'
longTrait='location2'
#print([x.traits for x in ll.Objects])
print([x.traits[longTrait] for x in ll.Objects if longTrait in x.traits])
print([x.traits[latTrait] for x in ll.Objects if latTrait in x.traits])
print([x.index for x in ll.Objects if latTrait in x.traits])

for k in bt_tree.Objects: ## iterate over branches
    locA=''
    if longTrait in k.parent.traits:
        long_parent=k.parent.traits[longTrait] ## get location of parent
        lat_parent=k.parent.traits[latTrait] ## get location of parent
        print('Parent long and lat%.5f %.5f'%(long_parent, lat_parent))
    long_current=k.traits[longTrait] ## get location of branch
    lat_current=k.traits[latTrait] ## get location of branch
    print('Current long and lat%.5f %.5f'%(long_current, lat_current))
    
    
    
    gdf.plot(ax=ax, color='red')
    
df.plot()

print('Done!')
```
3. 





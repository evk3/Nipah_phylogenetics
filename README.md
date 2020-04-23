# Nipah phylogenetics
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

File containing probe sequences can be found [here](/Nipah_oligos.xlsx).
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
Please see Jason's [readme](https://github.com/jtladner/Scripts/tree/master/strand_specific_coverage/strandspec_covplot).  The usage is the same.

## Chimeric Reads

This script is from Jason Ladner.  Check out his github [page](https://github.com/jtladner/Scripts/tree/master/chimeric_reads) for original [script](https://github.com/jtladner/Scripts/blob/master/chimeric_reads/chimeric_reads_v3.6.2.py).

## Genome Coverage Plots and Minor Variant Plots

A question was raised during review of the manuscript as to whether the use of enrichment probes biased the enriched reads.  We directly tested this by sequencing a paired subset of our viral isolates using enrichment and non-enrichment based methods.

![Supplementary_FigureX1](/Supplementary_FigureX1_V1_crop_flat?raw=true "Title")

<b>Supplementary Figure X1:</b>  Read coverage and minor variant plots for paired Nipah virus isolates sequenced using enrichment-based and non-enrichment based NGS methods.  <b>A)</b> Read coverage plots for a subset of paired viral isolates sequenced using enrichment-based and non-enrichment based NGS methods.  Plots show mean coverage versus genome position (bp) and loess-smoothed coverage estimates (dark, solid line).  <b>B)</b> Minor variant frequency versus genome position.  Subplots on y- and x-axes indicate density of minor variant locations across the genome and minor variant frequencies. 

![Supplementary_FigureX2](/Minor_variants_reviewers_8_to_share?raw=true "Title")

<b>Supplementary Figure X2:</b> Read coverage and minor variant plots for all new genomes sequenced using enrichment and non-enrichment based NGS methods, similar to Supplementary Figure X1.  <b>A)</b> Read coverage plots.  Plots show mean coverage versus genome position (bp) and loess-smoothed coverage estimates (dark, solid line).  <b>B)</b> Minor variant frequency versus genome position.  Subplots on y- and x-axes indicate density of minor variant locations across the genome and minor variant frequencies.   <b>C)</b> Minor frequency versus genome position.  Points are colored according to specimen number.  Left – minor variants for all new genomes.  Right – minor variants for all new genomes with variants for specimen 201206119 removed.

Link to raw variant data found [here](/Minor_variants_reviewers_8_to_share.xlsx).  Please note that genomes were aligned to maintain position consistency.  Position #1 corresponds to base #1 in JN808863.  At base 6558 there is a 0, 6, or 12 bp gap in the alignment.  All bases are re-aligned at base 6570.

## Visualizing Geographic Spread from a BEAST Continuous Diffusion Analysis

Hold on to your hats, folks.  This script has a lot of moving parts.  The original geographic spread script was made by Gytis Dudas and can be found [here](https://github.com/evogytis/baltic/blob/master/curonia.ipynb).  Kudos to him for this excellent piece of work!  This script was used to make the really cool Ebola virus spread video found [here](https://www.youtube.com/watch?v=j4Ut4krp8GQ).

My final plotting ipython notebook can be found [here](/curonia_nipah_V1.ipynb)

I had to break the curionia script down into several steps to understand and get the program working:

###### 1. Load a beast tree and geojson for your region of interest.  Are the geographic polygons correctly loaded?

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

###### 2. Can I import my Cauchy MCC tree and pull coordiate data out of it?
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
###### 3. Can I import a shapefile and add a matplot polygon to it?
```
#Can I import a shapefile and add a matplotlib polygon to it?
import descarteslabs as dl
import numpy as np
import geopandas as gpd

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


shape_file=gpd.read_file('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/location_variable_testing_2/GEOJSON/bangladesh-latest-free/gis_osm_water_a_free_1.shp')



#Rename rows by fclass column values:
shape_file=shape_file.set_index('fclass')
rivers_shape_file=shape_file.loc[['river', 'water']]
shape_file.info()

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
        
#print('polygons loaded:\n%s'%(polygons.keys()))

plt.figure(figsize=(15,15),facecolor='w') ## start figure
gs = gridspec.GridSpec(2, 1,height_ratios=[4,1]) ## define subplots
ax1 = plt.subplot(gs[0]) ## map here

ax1.plot()

for loc in polygons.keys():
    ax1.add_collection(PatchCollection(polygons[loc],facecolor='lightgrey',edgecolor='w',zorder=1)) ## plot location polygons

#for index, row in shape_file.iterrows():
#    if index == 'river':
#        print(index)
#        print(row['geometry'])
#        ax1.add_collection(PatchCollection(row['geometry'],facecolor='lightblue',edgecolor='',zorder=2))
        
        #internal_polygons[k.index]=([float(i) for i in k.traits[longPolygon]],[float(l) for l in k.traits[latPolygon]])
        #print(internal_polygons[k.index][0])
        #print(internal_polygons[k.index][1])
        #location_points_internal[k.index]=(np.vstack(zip(internal_polygons[k.index][0],internal_polygons[k.index][1])))
        #uncertainty_polygons[k.index].append(Polygon(location_points_internal[k.index], True))



ax1.set_aspect(1) ## equal aspect ratio
ax1.set_xlim(88.0, 93)
ax1.set_ylim(21.5, 26.75)

rivers_shape_file.plot(ax=ax1, facecolor='lightblue')

#plt.show()
```

###### 4. Can I animate by phylogentic tree by time and add it next to my map?
```
#Can I import a shapefile and add a matplotlib polygon to it?
import descarteslabs as dl
import numpy as np
import geopandas as gpd

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


shape_file=gpd.read_file('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/location_variable_testing_2/GEOJSON/bangladesh-latest-free/gis_osm_water_a_free_1.shp')



#Rename rows by fclass column values:
shape_file=shape_file.set_index('fclass')
rivers_shape_file=shape_file.loc[['river', 'water']]
shape_file.info()

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
        
#print('polygons loaded:\n%s'%(polygons.keys()))

plt.figure(figsize=(15,15),facecolor='w') ## start figure
gs = gridspec.GridSpec(2, 1,height_ratios=[4,1]) ## define subplots
ax1 = plt.subplot(gs[0]) ## map here

ax1.plot()

for loc in polygons.keys():
    ax1.add_collection(PatchCollection(polygons[loc],facecolor='lightgrey',edgecolor='w',zorder=1)) ## plot location polygons

#for index, row in shape_file.iterrows():
#    if index == 'river':
#        print(index)
#        print(row['geometry'])
#        ax1.add_collection(PatchCollection(row['geometry'],facecolor='lightblue',edgecolor='',zorder=2))
        
        #internal_polygons[k.index]=([float(i) for i in k.traits[longPolygon]],[float(l) for l in k.traits[latPolygon]])
        #print(internal_polygons[k.index][0])
        #print(internal_polygons[k.index][1])
        #location_points_internal[k.index]=(np.vstack(zip(internal_polygons[k.index][0],internal_polygons[k.index][1])))
        #uncertainty_polygons[k.index].append(Polygon(location_points_internal[k.index], True))



ax1.set_aspect(1) ## equal aspect ratio
ax1.set_xlim(88.0, 93)
ax1.set_ylim(21.5, 26.75)

rivers_shape_file.plot(ax=ax1, facecolor='lightblue')

#plt.show()
```

###### 5. Can I import my Cauchy MCC tree and pull coordinate and uncertainty polygons out of it?
```
#Can I import my Cauchy MCC tree and pull coordinate data and polygons out of it??

#Point to baltic.
import imp
bt = imp.load_source('baltic','/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/baltic.py') ## point to where baltic repo was cloned

#Test code to import my Cauchy MCC tree:
tree_path='/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/Nipah_Cauchy_new_HKYG_UCLN_skygrid_1234_log.time.MCC.trees'

bt_tree=bt.loadNexus(tree_path)

#Import geojson
#df=gpd.read_file('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/location_variable_testing_2/GEOJSON/SEAsia_forBEAST_V2.geojson')


popCentres={} ## dictionary with point coordinates
internal_polygons={} ##dictionary with internal node polygon estimates

longTrait='location2'
latTrait='location1'

longPolygon='location1_80%HPD_1'
latPolygon='location2_80%HPD_1'


for k in bt_tree.Objects: ## iterate over branches
    popCentres[k.index]=(float(k.traits[longTrait]),float(k.traits[latTrait] )) ## assign longitude and latitude to location
    
    if k.branchType == 'node' and k.index != 'Root':
        internal_polygons[k.index]=([float(i) for i in k.traits[longPolygon]],[float(l) for l in k.traits[latPolygon]])
        
        
print('migration centres:\n%s\n'%(popCentres.keys()))
print('dictionary format for migration centre coordinates:\n%s\n'%(popCentres))
    
print('internal polygons:\n%s\n'%(internal_polygons.keys()))
print('dictionary format for internal polygon coordinates:\n%s\n'%(internal_polygons))


print('Done!')
```

###### 6. Try plotting Bezier curves onto a map with a single time frame
```
#Objective: Try plotting Bezier curves with points from tree.

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
import imp

import unicodedata
# import unidecode ## for removing diacritics from example geoJSON

import numpy as np
from scipy.interpolate import UnivariateSpline ## used to smooth counts of lineages in each location at any given time
from scipy.interpolate import interp1d ## used to linearly interpolate between data points used in colouring polygons
#from sklearn.decomposition import IncrementalPCA ## used to identify PCA1 when automatically producing a colour map

import bezier ## custom arbitrary order Bezier curves
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

def Bezier_control(pointA,pointB,height,frac):
    """ 
    Given a line defined by 2 points A & B, 
    find a third point at a given distance (height) that defines a line perpendicular to line AB which intercepts AB 
    at fraction (frac) along AB.
    Equation derived by Luiz Max Fagundes de Carvalho (University of Edinburgh).
    """
    x1,y1=pointA
    x2,y2=pointB

    sign=1
    if x1>x2:
        sign=-1

    slope = (y2-y1) / (x2-x1)
    d=np.sqrt((y2-y1)**2 + (x2-x1)**2) ## distance between points
    
    h=np.sqrt(height**2+(d*frac)**2) ## distance between desired height and point along line

    n1=x1+h*np.cos(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign ## magic
    n2=y1+h*np.sin(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign

    return (n1,n2) ## return third point's coordinate
#**********************************************end of Bezier control function*****************************************

plt.figure(figsize=(15,18),facecolor='w') ## start figure
gs = gridspec.GridSpec(2, 1,height_ratios=[4,1]) ## define subplots
ax1 = plt.subplot(gs[0]) ## map here
ax2 = plt.subplot(gs[1]) ## plot to show relationship between distances between migration centres and the distance of the control point to the straight line connecting migration centres

#Point to baltic.
bt = imp.load_source('baltic','/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/baltic.py') ## point to where baltic repo was cloned

#Test code to import my Cauchy MCC tree:
tree_path='/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/Nipah_Cauchy_new_HKYG_UCLN_skygrid_1234_log.time.MCC.trees'

bt_tree=bt.loadNexus(tree_path)

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
        
#print('polygons loaded:\n%s'%(polygons.keys()))
#********************************end of pulling in map**************************************************

popCentres={} ## dictionary with point coordinates

#Pull in latitude and longitude values from continuous diffusion analysis.
#IE - just parse tree and make a dictionary with node index as key and latitude and longitude values.
longTrait='location1'
latTrait='location2'

for k in bt_tree.Objects: ## iterate over branches
    popCentres[k.index]=(float(k.traits[longTrait]),float(k.traits[latTrait] )) ## assign longitude and latitude to location
        

migration_function={i:{j:None for j in popCentres.keys() if i!=j} for i in popCentres.keys()} ## matrix of pairs of locations that will contain the Bezier function
all_distances=[] ## keep track of all distances
#control=lambda d:1/(np.abs(np.log(d))) ## function that will convert the distance between migration points into a distance 
                                #to be used for finding the control point for a Bezier function
control=lambda d:-1+0.1/float(d)**0.15+0.5
#control=lambda d:30*(1/d)**0.5

for i in popCentres.keys(): ## iterate over locations
    A=popCentres[i] ## fetch origin
    lon,lat=A
    #groupA=aggregation[i] ## get origin's high order group
    #c=group_colours[groupA](normalised_coordinates[i]) ## fetch colour map of the group, colour of location is determined by index of location's coordinates along PCA1
    
    #Set c to black for initial testing.
    c='k'
    
    ax1.scatter(lon,lat,s=30,facecolor=c,edgecolor='w',lw=2,zorder=100) ## plot migration centre
    
    for j in popCentres.keys(): ## iterate over locations again
        if i!=j: ## if not self
            B=popCentres[j] ## fetch destination
            d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between location A and location B
            
            all_distances.append(d) ## remember actual distince
            
            bez_points=[A] ## Bezier curve will start at point A
            #bez_points.append(Bezier_control(A,B,1/d ,0.1))
            #bez_points.append(Bezier_control(A,B,0.01     ,0.2))
            #bez_points.append(Bezier_control(A,B,-0.02    ,0.3))
            #bez_points.append(Bezier_control(A,B,-d**0.5,0.4))
            bez_points.append(Bezier_control(A,B,control(d),0.01)) ## head towards the first control point, which is perpendicular to line AB at distance control(d), and 0.01 of the way along line AB
            bez_points.append(Bezier_control(A,B,0.0,0.05)) ## head towards second control point, directly on the line AB, 0.1 of the way along line AB
            bez_points.append(B) ## Bezier curve will finish at point B
            bez_points=np.array(bez_points).transpose()
            curve = bezier.Curve(np.asfortranarray(bez_points),degree=len(bez_points)) ## Bezier curve object
            migration_function[i][j]=curve.evaluate_multi ## only interested in calling the evaluate_multi function of the curve, which will return a list of coordinates, given a numpy array of fractions along the line

for loc in polygons.keys():
    ax1.add_collection(PatchCollection(polygons[loc],facecolor='lightgrey',edgecolor='w',zorder=1)) ## plot location polygons

dict_locs=list(popCentres.keys())
print(dict_locs)
for locA in dict_locs:
    #print('In outer for loop')
    for locB in dict_locs[:dict_locs.index(locA)]: ## iterate over pairs of locations
        A=popCentres[locA]
        #print('A')
        #print(locA)
        #print(A)
        B=popCentres[locB]
        #print(locB)
        #print(B)
        d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between locations A and B
        sD=sorted(all_distances)
        
        if locA!=locB and (d in sD[:30] or d in sD[-7:]): ## if locations aren't the same and are some of the longest  or shortest distances - plot away
            print('d')
            print(d)
            Bezier_smooth=35 ## number of coordinates at which to compute the Bezier curve (more = smoother curve)
            eval_Bezier=np.linspace(0.0,1.0,Bezier_smooth) ## numpy array of 35 values going from 0.0 to 1.0
            migration=migration_function[locA][locB](eval_Bezier) ## compute Bezier curve coordinates along the path
            xs,ys=migration ## unpack Bezier coordinates
            print(A)
            print(B)
            print(xs)
            print(ys)
            for q in range(len(xs)-1): ## iterate through Bezier line segments with fading alpha and reducing width
                x1,y1=xs[q],ys[q] ## coordinates of current point
                x2,y2=xs[q+1],ys[q+1] ## coordinates of next point
                
                segL=(q+1)/float(len(xs)) ## fraction along length of Bezier line
                
                #Commenting out code to just plot fc as black. Haven't initialized aggregation yet.
                #if aggregation[locA]!=aggregation[locB]: ## locations in different high-order groups
                #    fc=group_colours[aggregation[locA]](normalised_coordinates[locA]) ## colour by origin colour
                #else: ## locations in same high-order group
                fc='r' ## colour black
                    
                ax1.plot([x1,x2],[y1,y2],lw=7*segL,alpha=1,color=fc,zorder=99,solid_capstyle='round') ## plot actual lineage with width proportional to position along Bezier curve
                ax1.plot([x1,x2],[y1,y2],lw=10*segL,alpha=1,color='w',zorder=98,solid_capstyle='round') ## plot white outline underneath

#Don't want to worry about this yet:
#for stretch in international_border: ## plot international border
#    xs,ys=zip(*stretch)
#    ax1.plot(xs,ys,color='k',zorder=11,lw=2)
#    ax1.plot(xs,ys,color='w',zorder=10,lw=5)
    
ax1.set_aspect(1) ## equal aspect ratio
ax1.set_xlim(88.0, 93)
ax1.set_ylim(21.5, 26.75)

ax2.plot(sorted(all_distances),list(map(control,sorted(all_distances))),ls='--',color='k') ## plot distances between points against their control point distances
ax2.set_ylabel('control point distance')
ax2.set_xlabel('distance between migration centres')

plt.show()
```

###### 7. OK, let's combine the animated tree, map and Bezier curves:
```
#Objective: Try plotting Bezier curves with points from tree.

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
import imp

import unicodedata
# import unidecode ## for removing diacritics from example geoJSON

import numpy as np
from scipy.interpolate import UnivariateSpline ## used to smooth counts of lineages in each location at any given time
from scipy.interpolate import interp1d ## used to linearly interpolate between data points used in colouring polygons
#from sklearn.decomposition import IncrementalPCA ## used to identify PCA1 when automatically producing a colour map

import bezier ## custom arbitrary order Bezier curves
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

def Bezier_control(pointA,pointB,height,frac):
    """ 
    Given a line defined by 2 points A & B, 
    find a third point at a given distance (height) that defines a line perpendicular to line AB which intercepts AB 
    at fraction (frac) along AB.
    Equation derived by Luiz Max Fagundes de Carvalho (University of Edinburgh).
    """
    x1,y1=pointA
    x2,y2=pointB

    sign=1
    if x1>x2:
        sign=-1

    slope = (y2-y1) / (x2-x1)
    d=np.sqrt((y2-y1)**2 + (x2-x1)**2) ## distance between points
    
    h=np.sqrt(height**2+(d*frac)**2) ## distance between desired height and point along line

    n1=x1+h*np.cos(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign ## magic
    n2=y1+h*np.sin(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign

    return (n1,n2) ## return third point's coordinate
#**********************************************end of Bezier control function*****************************************

plt.figure(figsize=(15,18),facecolor='w') ## start figure
gs = gridspec.GridSpec(2, 1,height_ratios=[4,1]) ## define subplots
ax1 = plt.subplot(gs[0]) ## map here
ax2 = plt.subplot(gs[1]) ## plot to show relationship between distances between migration centres and the distance of the control point to the straight line connecting migration centres

#Point to baltic.
bt = imp.load_source('baltic','/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/baltic.py') ## point to where baltic repo was cloned

#Test code to import my Cauchy MCC tree:
tree_path='/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/Nipah_Cauchy_new_HKYG_UCLN_skygrid_1234_log.time.MCC.trees'

bt_tree=bt.loadNexus(tree_path)

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
        
#print('polygons loaded:\n%s'%(polygons.keys()))
#********************************end of pulling in map**************************************************

popCentres={} ## dictionary with point coordinates

#Pull in latitude and longitude values from continuous diffusion analysis.
#IE - just parse tree and make a dictionary with node index as key and latitude and longitude values.
longTrait='location1'
latTrait='location2'

for k in bt_tree.Objects: ## iterate over branches
    popCentres[k.index]=(float(k.traits[longTrait]),float(k.traits[latTrait] )) ## assign longitude and latitude to location
        

migration_function={i:{j:None for j in popCentres.keys() if i!=j} for i in popCentres.keys()} ## matrix of pairs of locations that will contain the Bezier function
all_distances=[] ## keep track of all distances
#control=lambda d:1/(np.abs(np.log(d))) ## function that will convert the distance between migration points into a distance 
                                #to be used for finding the control point for a Bezier function
control=lambda d:-1+0.1/float(d)**0.15+0.5
#control=lambda d:30*(1/d)**0.5

for i in popCentres.keys(): ## iterate over locations
    A=popCentres[i] ## fetch origin
    lon,lat=A
    #groupA=aggregation[i] ## get origin's high order group
    #c=group_colours[groupA](normalised_coordinates[i]) ## fetch colour map of the group, colour of location is determined by index of location's coordinates along PCA1
    
    #Set c to black for initial testing.
    c='k'
    
    ax1.scatter(lon,lat,s=30,facecolor=c,edgecolor='w',lw=2,zorder=100) ## plot migration centre
    
    for j in popCentres.keys(): ## iterate over locations again
        if i!=j: ## if not self
            B=popCentres[j] ## fetch destination
            d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between location A and location B
            
            all_distances.append(d) ## remember actual distince
            
            bez_points=[A] ## Bezier curve will start at point A
            #bez_points.append(Bezier_control(A,B,1/d ,0.1))
            #bez_points.append(Bezier_control(A,B,0.01     ,0.2))
            #bez_points.append(Bezier_control(A,B,-0.02    ,0.3))
            #bez_points.append(Bezier_control(A,B,-d**0.5,0.4))
            bez_points.append(Bezier_control(A,B,control(d),0.01)) ## head towards the first control point, which is perpendicular to line AB at distance control(d), and 0.01 of the way along line AB
            bez_points.append(Bezier_control(A,B,0.0,0.05)) ## head towards second control point, directly on the line AB, 0.1 of the way along line AB
            bez_points.append(B) ## Bezier curve will finish at point B
            bez_points=np.array(bez_points).transpose()
            curve = bezier.Curve(np.asfortranarray(bez_points),degree=len(bez_points)) ## Bezier curve object
            migration_function[i][j]=curve.evaluate_multi ## only interested in calling the evaluate_multi function of the curve, which will return a list of coordinates, given a numpy array of fractions along the line

for loc in polygons.keys():
    ax1.add_collection(PatchCollection(polygons[loc],facecolor='lightgrey',edgecolor='w',zorder=1)) ## plot location polygons

dict_locs=list(popCentres.keys())
print(dict_locs)
for locA in dict_locs:
    #print('In outer for loop')
    for locB in dict_locs[:dict_locs.index(locA)]: ## iterate over pairs of locations
        A=popCentres[locA]
        #print('A')
        #print(locA)
        #print(A)
        B=popCentres[locB]
        #print(locB)
        #print(B)
        d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between locations A and B
        sD=sorted(all_distances)
        
        if locA!=locB and (d in sD[:30] or d in sD[-7:]): ## if locations aren't the same and are some of the longest  or shortest distances - plot away
            print('d')
            print(d)
            Bezier_smooth=35 ## number of coordinates at which to compute the Bezier curve (more = smoother curve)
            eval_Bezier=np.linspace(0.0,1.0,Bezier_smooth) ## numpy array of 35 values going from 0.0 to 1.0
            migration=migration_function[locA][locB](eval_Bezier) ## compute Bezier curve coordinates along the path
            xs,ys=migration ## unpack Bezier coordinates
            print(A)
            print(B)
            print(xs)
            print(ys)
            for q in range(len(xs)-1): ## iterate through Bezier line segments with fading alpha and reducing width
                x1,y1=xs[q],ys[q] ## coordinates of current point
                x2,y2=xs[q+1],ys[q+1] ## coordinates of next point
                
                segL=(q+1)/float(len(xs)) ## fraction along length of Bezier line
                
                #Commenting out code to just plot fc as black. Haven't initialized aggregation yet.
                #if aggregation[locA]!=aggregation[locB]: ## locations in different high-order groups
                #    fc=group_colours[aggregation[locA]](normalised_coordinates[locA]) ## colour by origin colour
                #else: ## locations in same high-order group
                fc='r' ## colour black
                    
                ax1.plot([x1,x2],[y1,y2],lw=7*segL,alpha=1,color=fc,zorder=99,solid_capstyle='round') ## plot actual lineage with width proportional to position along Bezier curve
                ax1.plot([x1,x2],[y1,y2],lw=10*segL,alpha=1,color='w',zorder=98,solid_capstyle='round') ## plot white outline underneath

#Don't want to worry about this yet:
#for stretch in international_border: ## plot international border
#    xs,ys=zip(*stretch)
#    ax1.plot(xs,ys,color='k',zorder=11,lw=2)
#    ax1.plot(xs,ys,color='w',zorder=10,lw=5)
    
ax1.set_aspect(1) ## equal aspect ratio
ax1.set_xlim(88.0, 93)
ax1.set_ylim(21.5, 26.75)

ax2.plot(sorted(all_distances),list(map(control,sorted(all_distances))),ls='--',color='k') ## plot distances between points against their control point distances
ax2.set_ylabel('control point distance')
ax2.set_xlabel('distance between migration centres')

plt.show()
```
###### 8. Final code that worked:
```
#OK, tree and map are animated, yeay!
#Let's add in uncertainty polygons and change the color of the clades.

#This code is from curonia.py from Gytis Dudas.
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

import descarteslabs as dl
import geopandas as gpd

import unicodedata
# import unidecode ## for removing diacritics from example geoJSON

import numpy as np
from scipy.interpolate import UnivariateSpline ## used to smooth counts of lineages in each location at any given time
from scipy.interpolate import interp1d ## used to linearly interpolate between data points used in colouring polygons
from sklearn.decomposition import IncrementalPCA ## used to identify PCA1 when automatically producing a colour map

import bezier ## custom arbitrary order Bezier curves
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
        timeline.append(d)
        carry, new_month = divmod(current_date.month - 1 + month_step, 12)
        new_month += 1
        current_date = current_date.replace(year=current_date.year + carry,month=new_month)
    return timeline

def Bezier_control(pointA,pointB,height,frac):
    """ 
    Given a line defined by 2 points A & B, 
    find a third point at a given distance (height) that defines a line perpendicular to line AB which intercepts AB 
    at fraction (frac) along AB.
    Equation derived by Luiz Max Fagundes de Carvalho (University of Edinburgh).
    """
    x1,y1=pointA
    x2,y2=pointB

    sign=1
    if x1>x2:
        sign=-1

    slope = (y2-y1) / (x2-x1)
    d=np.sqrt((y2-y1)**2 + (x2-x1)**2) ## distance between points
    
    h=np.sqrt(height**2+(d*frac)**2) ## distance between desired height and point along line

    n1=x1+h*np.cos(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign ## magic
    n2=y1+h*np.sin(np.arctan(height/float(d)/frac)+np.arctan(slope))*sign

    return (n1,n2) ## return third point's coordinate
#**********************************************end of Bezier control function*****************************************


def animate(frame):
    #### Primary plotting (map)
    ax1.lines=[line for line in ax1.lines if 'borders' in line.get_label()] ## reset lines (except borders) and texts in the plot
    #print(ax1.collections)
    ax1.collections=[patch for patch in ax1.collections if patch.get_label() in ('borders', 'leaf_node', 'internal_node')] ##reset patches (except map) 
    #ax1.collections=[]
    #print(ax1.collections)
    #ax1.points=[] ## reset points
    
    ax1.texts=[] ## remove all text from plot

    if len(animation_grid)-1==frame: ## if at last frame
        next_time=animation_grid[frame] ## current frame is next frame
    else:
        next_time=animation_grid[frame+1] ## next frame
        
    current_time=animation_grid[frame]
    
    
    effects=[path_effects.Stroke(linewidth=4, foreground='white'),
                 path_effects.Stroke(linewidth=0.5, foreground='k')] ## black text, white outline
    ax1.text(0.60,0.95,'Time: %.3f'%(current_time),size=20,transform=ax1.transAxes,zorder=1000,path_effects=effects) ## add text to indicate current time point
    

    #Find the branches in tree that we want to trace.
    exists=[k for k in bt_tree.Objects if k.parent!=bt_tree.root and k.parent.absoluteTime!=None and k.parent.absoluteTime<=current_time<=k.absoluteTime] ## identify lineages that exist at current timeslice

    #lineage_locations=[c.traits[locTrait] for c in exists if c.traits[locTrait]!='Not Available'] ## identify locations where lineages are present
    #presence=set(lineage_locations) ## all locations where lineages currently are in the tree
    
    #Changed above code to have presence trace all of the tree branches.
    presence=set(exists)
    
    update=10 ## update progress bar every X frames
    
    #### Secondary plotting (tree)
    Ls2=[x for x in ax2.lines if 'Colour' not in str(x.get_label())] ## fetch all lines in plot
    partials=[x for x in ax2.lines if 'partial' in str(x.get_label())] ## fetch all tree branches in progress
    finished_lines=[x for x in ax2.lines if 'finished' in str(x.get_label())] ## fetch all tree branches that are finished
    finished_points=[x for x in ax2.collections if 'finished' in str(x.get_label())] ## fetch all tip circles that are finished
    
    finished_labels=[str(x.get_label()) for x in finished_lines]+[str(x.get_label()) for x in finished_points] ## combine everything that's finished (branches + tip circles)
    partial_labels=[str(x.get_label()) for x in partials] ## partially plotted branches
    
    if frame>0 and frame%update==0: ## progress bar
        clear_output()
        timeElapsed=(time.time() - t0)/60.0
        progress=int((frame*(50/float(len(animation_grid)))))
        percentage=frame/float(len(animation_grid))*100
        rate=timeElapsed/float(frame)
        ETA=rate*(len(animation_grid)-frame)
        sys.stdout.write("[%-50s] %6.2f%%  frame: %5d %10s  time: %5.2f min  ETA: %5.2f min (%6.5f s/operation) %s %s %s" % ('='*progress,percentage,frame,animation_grid[frame],timeElapsed,ETA,rate,len(partials),len(finished_lines),len(finished_points)))
        sys.stdout.flush()

        
    #***********************Animate map!******************************************************************
    smooth=20 ## how many segments Bezier lines will have

    #tracking_length=365.0 ## number of days over which to plot the lineage
    depth=1
    departure_condition = lambda f:f-f ## determines how far away (in Bezier fraction) the tail of the migrating lineage is
    transition_point=1.0 ## determines the fixed time point along a branch at which migration happens
    locA=''
    #***********************Animate map!******************************************************************
    
    ####
    ## COMMENT this bit out if you don't want the tree to appear out of the time arrow
    ####
    
    for ap in bt_tree.Objects: ## iterate over branches
        idx='%s'%(ap.index) ## get unique id of branch
        xp=ap.parent.absoluteTime ## get parent's time

        x=ap.absoluteTime ## get branch's time
        y=ap.y ## get branch's y coordinate

        #location=ap.traits[locTrait] ## get branch's location
        
                     
        ## aggregate individual locations into higher-order groups (e.g. subdivisions into a country)

        group=aggregation[str(ap.index)] ## get location's group
        c=group_colours[group] ## get group colour map
        #c=cmap(normalised_coordinates[location]) ## get colour
        
        #Setting C to black for initial testing.
        #c='k'
        
        if xp!=None and xp<=current_time<x: ## branch is intersected by time arrow
            if 'partial_%s'%(idx) in partial_labels: ## if branch was partially drawn before
                l=[w for w in partials if 'partial_%s'%(idx)==str(w.get_label())][-1] ## get branch line
                l.set_data([xp,current_time],[y,y]) ## adjust its end coordinate to be time arrow
            else: ## branch is intersected, but not drawn before
                ax2.plot([xp,current_time],[y,y],lw=branchWidth,color=c,zorder=99,label='partial_%s'%(ap.index)) ## draw branch ending at time arrow, label as partially drawn
                
        if x<=current_time: ## time arrow passed branch - add it to finished class
            if 'partial_%s'%(idx) in partial_labels: ## if branch has been partially drawn before
                l=[w for w in partials if 'partial_%s'%(idx)==str(w.get_label())][-1] ## get branch
                l.set_data([xp,x],[y,y]) ## set end coordinate to be actual end coordinate
                l.set_label('finished_%s'%(idx)) ## set its label to finished
                
            if 'finished_%s'%(idx) not in finished_labels: ## branch has not been drawn before at all
                ax2.plot([xp,x],[y,y],lw=branchWidth,color=c,zorder=99,label='finished_%s'%(ap.index)) ## draw branch, add to finished class
                
            if 'partial_%s'%(idx) in partial_labels or 'finished_%s'%(idx) not in finished_labels: 
                if ap.branchType=='leaf': ## if leaf
                    ax2.scatter(x,y,s=tipSize,facecolor=c,edgecolor='none',zorder=102,label='finished_%s'%(ap.index)) ## add tip circle
                    ax2.scatter(x,y,s=tipSize*2,facecolor='k',edgecolor='none',zorder=101,label='finished_%s'%(ap.index)) ## add tip circle outline underneath
                elif ap.branchType=='node': ## if node
                    yl=ap.children[0].y ## get y coordinates of first and last child
                    yr=ap.children[-1].y
                    ax2.plot([x,x],[yl,yr],lw=branchWidth,color=c,zorder=99,label='finished_%s'%(ap.index)) ## plot vertical bar for node
    ####
    ## COMMENT this bit out if you don't want the tree to appear out of the time arrow
    ####
                
    #****************************end of animate tree***************************************************************  
    
    #****************************beginning of animate map***************************************************************
        #print(k.parent.absoluteTime)
        #print(k.absoluteTime)
        #print(time_point)
        locA=ap.parent.index ## get index of parent
        locB=ap.index ## get index of branch
        #print(popCentres[locA])
        #print(popCentres[locB])
        #print('Current time:')
        #print(current_time)
        #print("branch absolute time")
        #print(ap.absoluteTime)
        
        groupA=aggregation[str(ap.index)] ## get origin's high order group
        c=group_colours[groupA] ## fetch colour map of the group, colour of location is determined by index of location's coordinates along PCA1
    
        #Set c to black for initial testing.
        #c='k'
    
        A=popCentres[ap.index] ## fetch origin
        #print(ap.index)
        #print(A)
        lon,lat=A
        
        #Plot all points on the map with the code below:
        #if ap.branchType == 'leaf':
        #    ax1.scatter(lon,lat,s=30,facecolor=c,label='leaf_node', edgecolor='w',lw=2,zorder=99) ## plot migration centre
        #elif ap.branchType == 'node':
        #    ax1.scatter(lon,lat,s=30,facecolor='grey',label='internal_node', edgecolor='w',lw=2,zorder=98) ## plot migration centre
        
        #Plot root point:
        if ap.index == 0:
            ax1.scatter(lon,lat,s=30,facecolor='grey',label='internal_node', edgecolor='w',lw=2,zorder=98) ## plot migration centre
       
    
        if locA!='' and locA!=locB and ap.parent.absoluteTime!=None and ap.parent.absoluteTime<current_time<=ap.absoluteTime: ## if parent location is present, locations don't match, and branch is intersected
            #print('Current time:')
            #print(current_time)
            #print("branch absolute time")
            #print(ap.absoluteTime)
            
            transition=ap.parent.absoluteTime+ap.length*transition_point ## migration time
            tracking_length=ap.length
            #print(current_time)
            #print(ap.absoluteTime)
            #print('In if loop')
            #print((transition-time_point)/float(depth))
            begin_bezier=((current_time-ap.parent.absoluteTime)/(ap.absoluteTime-ap.parent.absoluteTime)) 
            #begin_bezier=(current_time/transition) ## migrating lineage's head
            end_bezier=departure_condition(begin_bezier) ## migrating lineage's tail
            #print(begin_bezier)
            #print(end_bezier)
        
            #begin_bezier,end_bezier = np.clip([begin_bezier,end_bezier],0.0,1.0) ## restrict to interval [0,1]
        
            #print(begin_bezier)
            #print(end_bezier)

        
            if end_bezier<=1.0 and begin_bezier>=0.0 and begin_bezier<=1.0: ## if lineage should still be visible (head departed and tail hasn't arrived)
                points=np.linspace(end_bezier,begin_bezier,smooth) ## get a bunch of points between the start and end of migrating lineage
                #print(points)
                #print(popCentres[locA])
                #print(popCentres[locB])
            
                #If an an internal node, plot uncertanty estimates for internal nodes.
                if ap.branchType == 'node' or ap.index == 0:
                    ax1.add_collection(PatchCollection(uncertainty_polygons[ap.index],facecolor=c,alpha=0.1,label='uncertainty',edgecolor='w',zorder=80)) ## plot location polygons

                
                #Plot points and have them appear and disappear by time.
                #Hard to view the data this way, though...
                if ap.branchType == 'leaf':
                    ax1.scatter(lon,lat,s=120,facecolor=c,label='leaf_node', edgecolor='w',lw=2,zorder=100) ## plot migration centre
                elif ap.branchType == 'node':
                    ax1.scatter(lon,lat,s=30,facecolor=c,label='internal_node', edgecolor='w',lw=2,zorder=100) ## plot migration centre
            
                bezier_line=migration_function[locA][locB](points) ## get coordinates of the migrating lineage
                #print(bezier_line)
                xs,ys=bezier_line
                for q in range(len(xs)-1): ## iterate through Bezier line segments with fading alpha and reducing width
                    x1,y1=xs[q],ys[q] ## get coordinates for current segment's start
                    x2,y2=xs[q+1],ys[q+1] ## get coordinates for current segment's end
                    segL=(q+1)/float(len(xs)) ## fraction along length of Bezier line
                
                    #Commenting the color stuff out, just plot one color.
                    #if aggregation[locA]!=aggregation[locB]: ## locations in different high-order groups
                    fc=group_colours[aggregation[str(ap.index)]] ## colour by origin colour
                    #else: ## locations in same high-order group
                    #    fc='k' ## colour black
                    #fc='r'
                
                    ax1.plot([x1,x2],[y1,y2],lw=7*segL,alpha=1,color=fc,zorder=101,solid_capstyle='round') ## plot actual lineage
                    ax1.plot([x1,x2],[y1,y2],lw=10*segL,alpha=1,color='w',zorder=100,solid_capstyle='round') ## plot underlying white background to help lineages stand out
    #***********************END of Animate map!******************************************************************
    
    for l in Ls2: ## iterate over lines in tree
        if 'time' in l.get_label(): ## if line is time arrow
            l.set_data([current_time,current_time],[0,1]) ## adjust time arrow
    
    #rivers_shape_file.plot(ax=ax1, facecolor='lightblue')

#*************************************end of animate definition*******************************************************


                    
#*************************************************************beginning of program main******************************

typeface='Helvetica Neue' ## set default matplotlib font and font size
mpl.rcParams['font.weight']=300
mpl.rcParams['axes.labelweight']=300
mpl.rcParams['font.family']=typeface
mpl.rcParams['font.size']=22

figWidth=10 ## map figure width
dpi=90 ## dots per inch for each .png (90 used in the final version)
every=12 ## put months labels every number of months

#Test code to import my Cauchy MCC tree:
tree_path='/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/Nipah_Cauchy_new_HKYG_UCLN_skygrid_1234_log.time.MCC.trees'

bt_tree=bt.loadNexus(tree_path)

popCentres={} ## dictionary with point coordinates
internal_polygons={} ## dictionary with coordinates of internal node uncertainty polygons
location_points_internal={}
uncertainty_polygons={}

#Pull in latitude and longitude values from continuous diffusion analysis.
#IE - just parse tree and make a dictionary with node index as key and latitude and longitude values.
longTrait='location1'
latTrait='location2'

longPolygon='location1_80%HPD_1'
latPolygon='location2_80%HPD_1'

#Initialize dictionaries with keys of tree indexes, but no values.
for r in bt_tree.Objects: ## iterate over branches
    uncertainty_polygons[r.index]=[]
    location_points_internal[r.index]=[]

#print('uncertainty polygons')
#print(uncertainty_polygons)
#print('locaiton points internal')
#print(location_points_internal)

for k in bt_tree.Objects: ## iterate over branches
    popCentres[k.index]=(float(k.traits[longTrait]),float(k.traits[latTrait] )) ## assign longitude and latitude to location
    
    #Grab internal polygons and uncertanty estimates for internal nodes.
    #Make internal_ploygon with Polygon from internal HPD_80%_estimates.
    if k.branchType == 'node' and k.index != 'Root':
        internal_polygons[k.index]=([float(i) for i in k.traits[longPolygon]],[float(l) for l in k.traits[latPolygon]])
        #print(internal_polygons[k.index][0])
        #print(internal_polygons[k.index][1])
        location_points_internal[k.index]=(np.vstack(zip(internal_polygons[k.index][0],internal_polygons[k.index][1])))
        uncertainty_polygons[k.index].append(Polygon(location_points_internal[k.index], True))
    
    
#Find the right coordinates for the map and tree.
#Width and height ratios are the lat and lon values from the map.
map_width=(93.0-88.0)
map_height=(26.75-21.5)
ratio=map_width/float(map_height) ## aspect ratio of map

plt.clf() 
plt.cla()
plt.figure(figsize=(figWidth*2,figWidth*ratio),facecolor='w') ## start figure
gs = gridspec.GridSpec(ncols=2, nrows=1,wspace=0.0,hspace=0.0) ## define subplots
ax1 = plt.subplot(gs[0]) ## map here
ax2 = plt.subplot(gs[1]) ## Tree

#frame='<iframe style="border: 0; width: 400px; height: 472px;" src="https://bandcamp.com/EmbeddedPlayer/album=29809561/size=large/bgcol=333333/linkcol=e99708/artwork=small/transparent=true/" seamless><a href="http://romowerikoito.bandcamp.com/album/nawam-r">NAWAMAR by Romowe Rikoito</a></iframe>'

print('Done!')
#HTML(frame)


#import geoJSON:
json_map=json.load(open('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/location_variable_testing_2/GEOJSON/bgd_admin2.geojson','r')) ## read from (hopefully saved) local copy

print('Done!')

#Pull in shapefile to include rivers and water on map.
shape_file=gpd.read_file('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/location_variable_testing_2/GEOJSON/bangladesh-latest-free/gis_osm_water_a_free_1.shp')

#Rename rows by fclass column values:
shape_file=shape_file.set_index('fclass')
rivers_shape_file=shape_file.loc[['river', 'water']]
#shape_file.info()

aggregation={'0':'group 0',
             '1':'group 2',
             '2':'group 2',
             '3':'group 2', 
             '4':'group 2',
             '773':'group 2',
             '3906':'group 2',
             '6319':'group 2',
             '6320':'group 2',
             '6321':'group 2',
             '7082':'group 2', 
             '7083':'group 2', 
             '7084':'group 2',
             '7085':'group 2',
             '7850':'group 2', 
             '7851':'group 2', 
             '7852':'group 2', 
             '8621':'group 2', 
             '8622':'group 2', 
             '12769':'group 2', 
             '19547':'group 2',
             '19548':'group 2',
             '19549':'group 2',
             '23140':'group 2', 
             '28708':'group 2', 
             '36218':'group 2', 
             '36219':'group 2', 
             '36987':'group 2', 
             '42968':'group 2', 
             '42969':'group 2',
             '42970':'group 2',
             '43733':'group 2',
             '50299':'group 2',
             '60809':'group 2',
             '67347':'group 1',
             '67348':'group 1',
             '68116':'group 1',
             '68117':'group 1',
             '68118':'group 1',
             '68119':'group 1',
             '68891':'group 1',
             '71353':'group 1',
             '71354':'group 1',
             '72122':'group 1',
             '72122':'group 1',
             '77072':'group 1',
             '77073':'group 1',
             '77834':'group 1',
             '77835':'group 1',
             '77836':'group 1',
             '77837':'group 1',
             '78593':'group 1',
             '78594':'group 1',
             '82243':'group 1',
             '89994':'group 1',
             '89995':'group 1',
             '90767':'group 1',
             '95074':'group 1',
             '95075':'group 1',
             '98927':'group 1',
             '98928':'group 1',
             '102855':'group 1',
             '102856':'group 1',
             '106741':'group 1',
             '106742':'group 1',
             '110295':'group 1'}

group_colours={'group 0': 'k',
              'group 1': 'b',
              'group 2': 'r'}

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
        
#print('polygons loaded:\n%s'%(polygons.keys()))
#********************************end of pulling in map**************************************************

#***************************************#Make Matrix of location distances*******************************************
control=lambda d:1/((np.sqrt(np.abs(np.log(d))))*2) ## function that will convert the distance between migration points into a distance 
                                #to be used for finding the control point for a Bezier function

all_distances=[] ## keep track of all distances
migration_function={i:{j:None for j in popCentres.keys() if i!=j} for i in popCentres.keys()} ## matrix of pairs of locations that will contain the Bezier function

for i in popCentres.keys(): ## iterate over locations
    A=popCentres[i] ## fetch origin
    lon,lat=A
    groupA=aggregation[str(i)] ## get origin's high order group
    c=group_colours[groupA] ## fetch colour map of the group, colour of location is determined by index of location's coordinates along PCA1
    
    #Set c to black for initial testing.
    #c='k'
    
    for j in popCentres.keys(): ## iterate over locations again
        if i!=j: ## if not self
            B=popCentres[j] ## fetch destination
            
            #Code below is historic from Gytis Dudas.  Not an appropriate way to calculate geographic distance.
            #Recommend using geopy to calculate geodesic distance instead!
            #Leaving code as is, because it is not used for plotting virus spread in map.
            d=np.sqrt((B[1]-A[1])**2 + (B[0]-A[0])**2) ## calculate distance between location A and location B
            
            all_distances.append(d) ## remember actual distince
            
            bez_points=[A] ## Bezier curve will start at point A
            #bez_points.append(Bezier_control(A,B,1/d ,0.1))
            #bez_points.append(Bezier_control(A,B,0.01     ,0.2))
            #bez_points.append(Bezier_control(A,B,-0.02    ,0.3))
            #bez_points.append(Bezier_control(A,B,-d**0.5,0.4))
            #bez_points.append(Bezier_control(A,B,control(d),0.01)) ## head towards the first control point, which is perpendicular to line AB at distance control(d), and 0.01 of the way along line AB
            bez_points.append(Bezier_control(A,B,0.01,0.01)) ## head towards the first control point, which is perpendicular to line AB at distance control(d), and 0.01 of the way along line AB
            bez_points.append(Bezier_control(A,B,0.0,0.1)) ## head towards second control point, directly on the line AB, 0.1 of the way along line AB
            bez_points.append(B) ## Bezier curve will finish at point B
            bez_points=np.array(bez_points).transpose()
            curve = bezier.Curve(np.asfortranarray(bez_points),degree=len(bez_points)) ## Bezier curve object
            migration_function[i][j]=curve.evaluate_multi ## only interested in calling the evaluate_multi function of the curve, which will return a list of coordinates, given a numpy array of fractions along the line
#***************************************Make Matrix of location distances*******************************************

# locTrait='location' ## name of locations in the tree
locTrait='location2_median'

t0 = time.time() ## time how long animation takes


start='1971-01-01' ## start date of animation
end='2017-01-01' ## end date of animation

animation_grid=list(np.linspace(bt.decimalDate(start),bt.decimalDate('1984-01-01'),20)) ## few frames from 2011 to 2015
animation_grid+=list(np.linspace(bt.decimalDate('1984-01-02'),bt.decimalDate(end),140)) ## many frames from 2015 to 2016
#animation_grid+=list(np.linspace(bt.decimalDate('2018-06-01'),bt.decimalDate(end),20)) ## medium frames from 2016 to 2018

tipSize=20 ## size of tip circles in tree
branchWidth=2 ## line width of tree branches

animation_duration=120.0 ## seconds

print('Start of animation: %.2f\nEnd: %.2f'%(min(animation_grid),max(animation_grid)))

print(len([x for x in bt_tree.Objects if locTrait not in x.traits]))

print('Number of frames to animate: %d'%(len(animation_grid)))
            
global travelers ## the animation will need to have information to traveling lineages
travelers=[x for x in bt_tree.Objects if x!=bt_tree.root] ## find lineages that have travelled - they're what's going to be animated
print('\nNumber of travelling lineages: %d (%.3f%% of all lineages)'%(len(travelers),len(travelers)/float(len(bt_tree.Objects))*100))


#ax2 = plt.axes()## ax2 is tree


xtimeline=calendarTimeline(start,end,'%Y-%m-%d','%Y-%m-%d') ## create timeline from start to end of animation delimited with months
xpos=[bt.decimalDate(b) for b in xtimeline] ## convert calendar dates to decimal time
xlabels=[(bt.convertDate(b,'%Y-%m-%d','%b\n%Y') if '-01-' in b else bt.convertDate(b,'%Y-%m-%d','%b')) if (int(b.split('-')[1])+every-1)%every==0 else '' for b in xtimeline] ## month or year-month (if January) tick labels every given number of months


################
## Tertiary plot begins - TREE
################

####
## UNCOMMENT if you'd like the tree to be plotted in grey initially and get coloured over time
####
# iterate over objects in tree
for k in bt_tree.Objects:
    #loc=k.traits[locTrait] ## get branch location
    group=aggregation[str(k.index)]
    c=group_colours[group] ## get colour map
    #c=cmap(normalised_coordinates[loc]) ## fixed colour for location
    
    #set c = black for now.
    c='k'
    
    mask_cmap=mpl.cm.Greys ## colour map that goes on top of lines
    #grey_colour=mask_cmap((sorted_group.index(group)+1)/float(len(sorted_group)+2)) ## get colour based on index of location in list, avoiding white colour
    grey_colour="grey"
    
    y=k.y
    yp=k.parent.y
    
    x=k.absoluteTime
    xp=k.parent.absoluteTime
    
    if k.branchType=='leaf':
        ax2.scatter(x,y,s=tipSize,facecolor=grey_colour,edgecolor='none',zorder=102,label='LeafBW_%s'%(k.index)) ## plot black and white tip circle on top
        ax2.scatter(x,y,s=tipSize,facecolor=c,edgecolor='none',zorder=101,label='LeafColour_%s'%(k.index)) ## plot colour tip circle underneath black and white tip circle
        ax2.scatter(x,y,s=tipSize*2,facecolor='k',edgecolor='none',zorder=100,label='Colour') ## black outline underneath every tip
        
    elif k.branchType=='node':
        yl=k.children[0].y
        yr=k.children[-1].y
       
        if xp==0.0:
            xp=x

        ax2.plot([x,x],[yl,yr],color=grey_colour,lw=branchWidth,zorder=99,label='NodeHbarBW_%s'%(k.index))
        ax2.plot([x,x],[yl,yr],color=c,lw=branchWidth,zorder=98,label='NodeHbarColour_%s'%(k.index))
        
    ax2.plot([xp,x],[y,y],color=grey_colour,lw=branchWidth,zorder=99,label='BranchBW_%s'%(k.index)) ## plot black and white branch on top
    ax2.plot([xp,x],[y,y],color=c,lw=branchWidth,zorder=98,label='BranchColour_%s'%(k.index)) ## plot colour branch underneath black and white branch (note the zorder attribute)
    
####
## UNCOMMENT if you'd like the tree to be plotted in grey initially and get coloured over time
####

#plot all of the points right away.
#for k in bt_tree.Objects:
#    A=popCentres[k.index] ## fetch origin
#    lon,lat=A
    #ax1.scatter(lon,lat,s=30,facecolor=c,edgecolor='w',lw=2,zorder=100) ## plot migration centre

    
for loc in polygons.keys():
    ax1.add_collection(PatchCollection(polygons[loc],facecolor='lightgrey',label='borders', edgecolor='w',zorder=1)) ## plot location polygons
    
#ax1.plot()## ax1 is map
ax1.set_aspect(1) ## equal aspect ratio
ax1.set_xlim(88.0, 93)
ax1.set_ylim(21.5, 26.75)

#******************************************************************************
#Add rivers_shape_file plotting here:
rivers_shape_file.plot(ax=ax1, facecolor='lightblue', label='borders')

ax2.axvline(xpos[0],color='k',lw=3,label='time',zorder=200) ## add time arrow to indicate current time

ax2.set_xticks([x+1/24.0 for x in xpos]) ## add ticks, tick labels and month markers
ax2.set_xticklabels(xlabels) ## set x axis labels
[ax2.axvspan(xpos[x],xpos[x]+1,facecolor='k',edgecolor='none',alpha=0.04) for x in range(0,len(xpos),2)] ## set grey vertical bars for x axis

ax2.xaxis.tick_bottom() ## make tree plot pretty
ax2.yaxis.tick_left()
[ax2.spines[loc] for loc in ax2.spines]

ax2.tick_params(axis='x',size=0) ## no ticks
ax2.tick_params(axis='y',size=0)
ax2.set_xticklabels([]) ## no tick labels
ax2.set_yticklabels([])

ax2.set_xlim(min(xpos),max(xpos)) ## axis limits
ax2.set_ylim(-bt_tree.ySpan*0.01,bt_tree.ySpan*1.01)
################
## Tertiary plot ends - TREE
################
            
for frame in range(len(animation_grid)): ## iterate through each frame
    animate(frame) ## animate will modify the map, tree and cases
    plt.savefig('/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/animate_tree_V2/ani_frame_%05d.png'%(frame), format='png',bbox_inches='tight',dpi=dpi) ## save individual frames for stitching up using 3rd party software (e.g. FFMpeg)
    
print('\n\nDONE!')

print('\nTime taken: %.2f minutes'%((time.time() - t0)/60.0))

print('Suggested frame rate for an animation lasting %s seconds using %s frames: %s'%(animation_duration,len(animation_grid),int(len(animation_grid)/animation_duration))) ## frame rate
plt.show()
```

###### 9. Convert individual png images into a movie:
```
%%bash

frames=/scicomp/home/evk3/Diagnostics/Nipah/GLM_2018/alignment_without_mojiang_cedar_bat/bang_geo_only/animate_tree_V2/
ffmpeg_path=/scicomp/home/evk3/setup/ffmpeg-4.1

cd  $frames; $ffmpeg_path/ffmpeg -framerate 2 -start_number 0 -i ani_frame_%05d.png -pix_fmt yuv420p -b:a 64k -vf scale="2160:trunc(ow/a/2)*2" nipah_25Apr2019_animation.HD.264.mp4
```




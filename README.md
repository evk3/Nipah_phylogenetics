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

My final potting ipython notebook can be found [here](/curonia_nipah_V1.p)

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

###### 2. Can i import my Cauchy MCC tree and pull coordiate data out of it?
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

###### 5. Can I import by Cauchy MCC tree and pull coordinate and uncertainty polygons out of it?
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

###### 6. Try plotting Bexier curves onto a map with a single time frame
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
test
```





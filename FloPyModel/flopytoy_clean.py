#----------------------------------------------------------------------------
#
# Created on Sun Mar 29 15:06:36 2020
# @author: dbabrams
# FloPy Toy Model code, runs MODFLOW 2005
# Unless otherwise stated, all units are in feet and days.
#
#----------------------------------------------------------------------------

#%% 

'''Import packages.'''
#----------------------------------------------------------------------------
import flopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
#----------------------------------------------------------------------------

#%% 

'''Create a MODFLOW model object and run with MODFLOW 2005.'''
#----------------------------------------------------------------------------
modelname = "my_model2"
m = flopy.modflow.Modflow(modelname, exe_name = 'mf2005')
#----------------------------------------------------------------------------

#%% 

'''Define model grid, river cells, and top and bottom elevations.'''
#----------------------------------------------------------------------------
# Define Model Grid
xlo = 2727053
xhi = 3165901
ylo = 2528219
yhi = 2800979
Lx = xhi - xlo # Width of the model domain
Ly = yhi - ylo # Height of the model domain
nlay = 2 # Number of model layers
nrow = 50 # Number of rows
ncol = 75 # Number of columns
dx = Lx/ncol # grid spacing (x-direction)
dy = Ly/nrow # grid spacing (y-direction)

# Define river cells
df = pd.read_csv('https://raw.githubusercontent.com/dbabrams/G572_Mahomet_Example/develop_abrams/RiverElevationData/majorriverelevations.csv')
df['lay'] = 0
df['row'] = nrow- np.floor((df['lamy']-ylo)/dy)-1
df = df[df['row']>=0]
df = df[df['row']<nrow]
df['col'] = np.floor((df['lamx']-xlo)/dx)
df = df[df['col']>=0]
df = df[df['col']<ncol]
df['stage'] = df['VALUE']
df['cond'] = 50000.
df['rbot'] = df['stage'] - 5.
df = df.drop(['wkt_geom','VALUE','lamx','lamy'], axis=1)
df = df.drop_duplicates(subset=['row', 'col'])
df_riv = df

# Top elevation array
df = pd.read_csv('https://raw.githubusercontent.com/dbabrams/G572_Mahomet_Example/develop_abrams/elevations/l1_top.csv')
df['row'] = nrow- np.floor((df['lamy']-ylo)/dy)-1
df = df[df['row']>=0]
df = df[df['row']<nrow]
df['col'] = np.floor((df['lamx']-xlo)/dx)
df = df[df['col']>=0]
df = df[df['col']<ncol]
df['elev'] = df['VALUE']
df = df.drop(['wkt_geom','VALUE','lamx','lamy'], axis=1)
df = df.drop_duplicates(subset=['row', 'col'])
ztop = np.zeros([nrow,ncol])
for index, values in df.iterrows():
    ztop[np.int(values['row']),np.int(values['col'])]=values['elev']

# burning the river elevations into the top elevation grid    
for index, values in df_riv.iterrows():
    ztop[np.int(values['row']),np.int(values['col'])]=values['stage']
    
# Bottom elevation array
df = pd.read_csv('https://raw.githubusercontent.com/dbabrams/G572_Mahomet_Example/develop_abrams/elevations/l1_bot.csv')
df['row'] = nrow- np.floor((df['lamy']-ylo)/dy)-1
df = df[df['row']>=0]
df = df[df['row']<nrow]
df['col'] = np.floor((df['lamx']-xlo)/dx)
df = df[df['col']>=0]
df = df[df['col']<ncol]
df['elev'] = df['VALUE']
df = df.drop(['wkt_geom','VALUE','lamx','lamy'], axis=1)
df = df.drop_duplicates(subset=['row', 'col'])
zbot = np.zeros([nrow,ncol])
for index, values in df.iterrows():
    zbot[np.int(values['row']),np.int(values['col'])]=values['elev']

# Ensuring that the model is 10 ft thick everywhere
for idx, x in np.ndenumerate(ztop):
    if zbot[idx]+10>=x:
        zbot[idx]= x-1.
        
zbot1 = zbot + 2/3*(ztop-zbot)
zbot2 = zbot        
zbot = [zbot1, zbot2]

#%%         
'''Create the Discretization package'''
#----------------------------------------------------------------------------       
#for index, values in df_riv.iterrows():
#    zbot[np.int(values['row']),np.int(values['col'])]=values['stage']-50        
        
nper = 1 #specify number of stress periods
steady = [True] #specify if stress period is transient or steady-state

# create flopy discretization object
# length and time are feet (1) and days (4).
# See https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/index.html?dis.htm 
dis = flopy.modflow.ModflowDis(model=m, nlay=nlay, nrow=nrow, ncol=ncol, 
                               delr=dx, delc=dy, top=ztop, botm=zbot, 
                               itmuni = 4, lenuni = 1, 
                               nper=nper, steady=steady)
#----------------------------------------------------------------------------

#%%


'''Create the Basic Package, which contains ibound and starting heads'''
#----------------------------------------------------------------------------
# Create ibound as array of ints (1), indicating all cells are active
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
#ibound[:, :, 0] = -1 # Designate left boundary cells as constant head
#ibound[:, :, -1] = -1 # Designate right boundary cells as constant head

# Create starting head array, must be floats.
#strt = 800*np.ones((nlay, nrow, ncol), dtype=np.float32) #set every cell to 5.0
#strt[:, :, 0] = 10. #set left side head to 10 ft
#strt[:, :, -1] = 0. #set right side head to 0 ft

#Create flopy bas object

bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=10000.)
#----------------------------------------------------------------------------
#%%

'''Create a river package'''
#----------------------------------------------------------------------------
# https://flopy.readthedocs.io/en/latest/source/flopy.modflow.mfriv.html
rivs = {0: df_riv.to_numpy()}
riv = flopy.modflow.ModflowRiv(model=m, stress_period_data=rivs)
#----------------------------------------------------------------------------
#%%

'''Create the Layer Property Flow Package, which contains information about
hydraulic conductivity and other information about how to calculate flow'''
#----------------------------------------------------------------------------
hi = 250.
# point 2: increased low k
lo = 10.
thresh = 25.

#layer 1
df = pd.read_csv('https://raw.githubusercontent.com/dbabrams/G572_Mahomet_Example/develop_abrams/hydraulicconductivity/l1_k.csv')
df['row'] = nrow- np.floor((df['lamy']-ylo)/dy)-1
df = df[df['row']>=0]
df = df[df['row']<nrow]
df['col'] = np.floor((df['lamx']-xlo)/dx)
df = df[df['col']>=0]
df = df[df['col']<ncol]
df['hk'] = df['VALUE']
df = df.drop(['wkt_geom','VALUE','lamx','lamy'], axis=1)
df_hk1 = df.drop_duplicates(subset=['row', 'col'])
hk1 = hi*np.ones([nrow,ncol])
for index, values in df_hk1.iterrows():
    if values['hk']<=thresh:
        hk1[np.int(values['row']),np.int(values['col'])]=lo
        
# Point 3: add high k where rivers are present
for index, values in df_riv.iterrows():
    hk1[np.int(values['row']),np.int(values['col'])]=hi   

# layer 2
df = pd.read_csv('https://raw.githubusercontent.com/dbabrams/G572_Mahomet_Example/develop_abrams/hydraulicconductivity/l9_k.csv')
df['row'] = nrow- np.floor((df['lamy']-ylo)/dy)-1
df = df[df['row']>=0]
df = df[df['row']<nrow]
df['col'] = np.floor((df['lamx']-xlo)/dx)
df = df[df['col']>=0]
df = df[df['col']<ncol]
df['hk'] = df['VALUE']
df = df.drop(['wkt_geom','VALUE','lamx','lamy'], axis=1)
df = df.drop_duplicates(subset=['row', 'col'])
hk2 = hi*np.ones([nrow,ncol])
for index, values in df.iterrows():
    if values['hk']<=thresh:
        hk2[np.int(values['row']),np.int(values['col'])]=lo
        
# add high k where rivers are present
#for index, values in df_riv.iterrows():
#    hk2[np.int(values['row']),np.int(values['col'])]=hi   

hk = [hk1, hk2]
vk = [hk1/10, hk2/10] #define vertical hydraulic conductivity

#define layer type as convertible (1), must be an integer
#for more information, see https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/index.html?dis.htm
laytyp = np.ones((nlay,), dtype=np.int32)

# create the LPF object
lpf = flopy.modflow.ModflowLpf(model=m, hk=hk, vka=vk, laytyp=laytyp, ipakcb=1)
#----------------------------------------------------------------------------
#%%

'''Create a drain package'''

#----------------------------------------------------------------------------
newdata= {}
i = 0
for index, values in df_hk1.iterrows():
    if values['hk']<=thresh:
        rn = np.int(values['row'])
        cn = np.int(values['col'])
        newdata[i] = [0,rn,cn,ztop[rn,cn],100.]
        i = i+1
        
df_drn = pd.DataFrame.from_dict(newdata, orient='index')
        
drns = {0: df_drn.to_numpy()}
drn = flopy.modflow.ModflowDrn(model=m, stress_period_data=drns)
#----------------------------------------------------------------------------
#%%

'''Create a recharge package'''
#----------------------------------------------------------------------------
hi = 0.003
lo = 0.003
df = pd.read_csv('https://raw.githubusercontent.com/dbabrams/G572_Mahomet_Example/develop_abrams/hydraulicconductivity/l1_k.csv')
df['row'] = nrow- np.floor((df['lamy']-ylo)/dy)-1
df = df[df['row']>=0]
df = df[df['row']<nrow]
df['col'] = np.floor((df['lamx']-xlo)/dx)
df = df[df['col']>=0]
df = df[df['col']<ncol]
df['rchg'] = df['VALUE']
df = df.drop(['wkt_geom','VALUE','lamx','lamy'], axis=1)
df = df.drop_duplicates(subset=['row', 'col'])
rchg = hi*np.ones([nrow,ncol])
for index, values in df.iterrows():
    if values['rchg']<=thresh:
        rchg[np.int(values['row']),np.int(values['col'])]=lo

rech = flopy.modflow.ModflowRch(model=m,rech = rchg)
#----------------------------------------------------------------------------
#%%

'''Create the Output Control Package'''
#----------------------------------------------------------------------------
#create oc stress period data. 
spd = {(0, 0): ['print head', 'print budget', 'save head', 'save budget']}
#create output control object
oc = flopy.modflow.ModflowOc(model=m, stress_period_data=spd, compact=True)
#----------------------------------------------------------------------------

#%%

'''Create the PCG Solver Object'''
#----------------------------------------------------------------------------
# for the time being, we will use default settings with the solver
#nwt = flopy.modflow.ModflowNwt(model=m)
pcg = flopy.modflow.ModflowPcg(model=m)
#----------------------------------------------------------------------------

#%%

'''Write MODFLOW input files.'''
#----------------------------------------------------------------------------
m.write_input()
#----------------------------------------------------------------------------

#%%

'''Run the model'''
#----------------------------------------------------------------------------
# Executve the model run
success, mfoutput = m.run_model(pause=False, report=True)
# Report back if the model did not successfully complete
if not success:
    raise Exception('MODFLOW did not terminate normally.')
#----------------------------------------------------------------------------
    
#%%   
    
'''Extract binary data from head and flow files'''
#----------------------------------------------------------------------------
#extract binary data from head file as flopy head object
headobj = flopy.utils.binaryfile.HeadFile(modelname+'.hds')
#extract head data from head object
head = headobj.get_data(totim=1.0)

#extract binary data from budget file as flopy budget object
budgobj = flopy.utils.binaryfile.CellBudgetFile(modelname+'.cbc')
#extract flow data from budget object, define face over which flow is reported
frf = budgobj.get_data(text='flow right face', totim=1.0)
fff = budgobj.get_data(text='flow front face', totim=1.0)
#----------------------------------------------------------------------------

#%%

'''Plot grid and boundary conditions'''
#----------------------------------------------------------------------------
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.PlotMapView(model=m, layer=0)
#grid = modelmap.plot_grid()
ib = modelmap.plot_ibound()
bc = modelmap.plot_bc('drn')
bc = modelmap.plot_bc('riv')
#add labels and legend
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Ibound', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='blue',label='Const. Head',ec='black'),
                   mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                   mp.patches.Patch(color='black',label='Inactive Cell',ec='black'),
                   mp.patches.Patch(color='green',label='River',ec='black'),
                   mp.patches.Patch(color='yellow',label='Drain',ec='black')],
                   bbox_to_anchor=(1.5,1.0))
#----------------------------------------------------------------------------

#%%

'''Plot results'''
#----------------------------------------------------------------------------
layviz = 1

plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=m, layer=layviz) #use plotmapview to attach plot to model
#grid = modelmap.plot_grid() #plot model grid
bc = modelmap.plot_bc('riv')
contour_levels = np.linspace(head[layviz].min(),head[layviz].max(),11) #set contour levels for contouring head
head_contours = modelmap.contour_array(head[layviz], levels=contour_levels) #create head contours
#flows = modelmap.plot_discharge(frf[0], fff[0], head=head) #create discharge arrows

#display parameters
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Steady-State Model, Flow(ft^3/d) and Head(ft) Results', fontsize = 15, fontweight = 'bold')
plt.colorbar(head_contours,aspect=5)
plt.show(modelmap)
#----------------------------------------------------------------------------
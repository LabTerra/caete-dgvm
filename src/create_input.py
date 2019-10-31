import pickle as pkl
from create_input import data_in

# Location of the input data downloaded from ISIMIP-PIK server via rsync
files = "/d/c1/homes/amazonfaceme/jpdarela/dlds/dlds_daily" 

# INstantiate our input data object
data_in = data_in(files)

# TODO parrallelization of gridcell input creation/update

## Use with a mask (creates a bunch og gridcells)
#data_in = data_in.data_dict(mask=ct.mask)

# OR

## Update a existent instance of data_in
with open("data_in-instance.pkl", 'rb') as fh:
    data_in = pkl.load(fh)

# Apply data_dict creation method:

grid_points = [(239, 183),
               (219, 184),
               (404, 186),
               (440, 61),
               (572, 55)]


for x, y in grid_points:
    for var in data_in.varnames:
        data_in.data_dict(varname=var, nx=x, ny=y)

with open("data_in-instance.pkl", 'wb') as fh:
    pkl.dump(data_in, fh)


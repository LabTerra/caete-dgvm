
# Functionalities to write model outputs 

import numpy as np
from netCDF4 import Dataset as dt
from caete_module import global_pars as gp

nt = gp.ntimes
npls = gp.npls

mask_fpath = './mask12.npy'
NO_DATA = [-9999.0, -9999.0]
lsmk = np.load(mask_fpath)

monthly_out = ['aresp',
               'clit',
               'csoil',
               'emaxm',
               'evapm',
               'hresp',
               'lai',
               'npp',
               'photo',
               'rcm',
               'rg',
               'rm',
               'runom',
               'tsoil',
               'wsoil',
               'wue',
               'cue'] 

npls_out = ['area',
            'area0',
            'clin',
            'cfin',
            'cwin',
            'cleaf',
            'cawood',
            'cmass',
            'cfroot']

def mask_gen(nlayers):
    mask1 = lsmk[0]
    if nlayers == 1:
        return mask1
    z = np.zeros(shape=(nlayers, mask1.shape[0],mask1.shape[1]),dtype=np.bool)
    for i in range(nlayers):
        z[i,:,:] = mask1
    return z

def flt_attrs():
    return {'header'  : ['long_name',                 'unit',           'standart_name', 'ldim'],
            
            'rsds'    : ['short_wav_rad_down',        'W m-2',                   'rsds',    nt], 
            'wind'    : ['wind_velocity',             'm s-1',                   'wind',    nt],
            'ps'      : ['sur_pressure',              'Pa',                        'ps',    nt],
            'tas'     : ['sur_temperature_2m',        'celcius',                  'tas',    nt],
            'tsoil'   : ['soil_temperature',          'celcius',            'soil_temp',    nt],
            'pr'      : ['precipitation',             'Kg m-2 month-1',            'pr',    nt],      
            'wsoil'   : ['soil_water_content-wsoil',  'kg m-2',                  'mrso',    nt],
            'evapm'   : ['evapotranpiration',         'kg m-2 day-1',              'et',    nt],
            'emaxm'   : ['potent. evapotranpiration', 'kg m-2 day-1',           'etpot',    nt],
            'runom'   : ['total_runoff',              'kg m-2 day-1',            'mrro',    nt],
            'aresp'   : ['autothrophic respiration',  'kg m-2 year-1',             'ar',    nt],
            'photo'   : ['photosynthesis',            'kg m-2 year-1',             'ph',    nt],
            'npp'     : ['net primary productivity',  'kg m-2 year-1',            'npp',    nt],
            'lai'     : ['Leaf Area Index',           'm2 m-2',                   'LAI',    nt],
            'rcm'     : ['stomatal resistence',       's m-1',                    'rcm',    nt],
            'hresp'   : ['heterothrophic respiration','kg m-2 year-1',             'hr',    nt],
            'clit'    : ['Litter Carbon',             'Kg m-2',                  'clit',    nt],
            'csoil'   : ['Soil Carbon',               'Kg m-2',                 'csoil',    nt],
            'rm'      : ['maintenance respiration',   'kg m-2 year-1',             'rm',    nt],
            'rg'      : ['growth respiration',        'kg m-2 year-1',             'rg',    nt],
            'wue'     : ['water use efficiency',      'ad',                       'wue',    nt],
            'cue'     : ['Carbon use efficiency',     'ad',                       'cue',    nt],
            'cawood'  : ['C in abovewgrownd wood',    'kg m-2',                'cawood',  npls],
            'cfroot'  : ['C in fine roots',           'kg m-2',                'cfroot',  npls],
            'cleaf'   : ['C in leaves',               'kg m-2',                 'cleaf',  npls],
            'area'    : ['occupation coefficient',    '%',                       'area',  npls],
            'area0'   : ['occupation coefficient',    '%',                  'init_area',  npls],
            'cwin'    : ['init C in abovewgrownd wood','kg m-2',               'cawood',  npls],
            'cfin'    : ['init C in fine roots',       'kg m-2',               'cfroot',  npls],
            'cmass'   : ['total Carbon -Biomass',      'kg m-2',                'cmass',  npls],
            'clin'    : ['init C in leaves',           'kg m-2 ',               'cleaf',  npls]}


def write_CAETE_output(nc_filename, arr, var, pls_mode=False):

    ldim = flt_attrs()[var][3]
    if ldim > 1:
        t, la, lo,  = arr.shape
        one_layer = False
    else:
        la,lo = arr.shape
        one_layer = True

    lsmk_internal = mask_gen(ldim)
    # create netcdf file
    rootgrp = dt(nc_filename, mode='w', format='NETCDF3_CLASSIC')
    #dimensions  & variables
        
    rootgrp.createDimension("latitude", la)
    rootgrp.createDimension("longitude", lo)

    if var in monthly_out  :
        rootgrp.createDimension("time", None)
        time = rootgrp.createVariable(varname="time", datatype=np.float32,
                                      dimensions=("time",))
    elif var in npls_out:
        rootgrp.createDimension("pls", npls)
        pls = rootgrp.createVariable(varname="pls", datatype=np.int32,
                                      dimensions=("pls",))
    
    latitude  = rootgrp.createVariable(varname="latitude", datatype=np.float32,dimensions=("latitude",))
    longitude = rootgrp.createVariable(varname="longitude", datatype=np.float32, dimensions=("longitude",))

    if var in monthly_out:
        var_ = rootgrp.createVariable(varname = 'annual_cycle_mean_of_'+str(flt_attrs()[var][2]), datatype=np.float32,
                                       dimensions=("time","latitude","longitude",),
                                       fill_value=NO_DATA[0])
    elif var in npls_out:
        var_ = rootgrp.createVariable(varname = str(flt_attrs()[var][2]), datatype=np.float32,
                                       dimensions=("pls","latitude","longitude",),
                                       fill_value=NO_DATA[0])
    elif one_layer:
        var_ = rootgrp.createVariable(varname = str(flt_attrs()[var][2]), datatype=np.float32,
                                       dimensions=("latitude","longitude",),
                                       fill_value=NO_DATA[0])

   
    #attributes
    ## rootgrp
    rootgrp.description = flt_attrs()[var][0] + " caete-v1.0 OUTPUT"
    rootgrp.source = "CAETE model outputs"
    ## time
    if var in monthly_out:
        time.units = "days since 1850-01-01 00:00:00.0"
        time.calendar = "noleap"
        time.axis='T'

    if var in npls_out:
        pls.units = '1'
        pls.axis= u'T'
        
    ## lat
    latitude.units = u"degrees_north"
    latitude.long_name=u"latitude"
    latitude.standart_name =u"latitude"
    latitude.axis = u'Y'
    ## lon
    longitude.units = "degrees_east"
    longitude.long_name = "longitude"
    longitude.standart_name = "longitude"
    longitude.axis = u'X'
    ## var
    var_.long_name=flt_attrs()[var][0]
    var_.units = flt_attrs()[var][1]
    var_.standard_name=flt_attrs()[var][2]
    var_.missing_value=NO_DATA[0]

    ## WRITING DATA
    if var in monthly_out:
        time[:] = np.array([15.5, 45., 74.5, 105., 135.5, 166.,
                               196.5, 227.5, 258., 288.5, 319., 349.5])
    if var in npls_out:
        pls[:] = np.arange(1,npls+1)
    
    longitude[:] = np.arange(-179.75, 180, 0.5)
    latitude[:] =  np.arange(-89.75, 90, 0.5)
     
    if ldim > 1:
        var_[:,:,:] = np.fliplr(np.ma.masked_array(arr, lsmk_internal))
    else:
        var_[:,:] = np.flipud(np.ma.masked_array(arr, lsmk_internal))
    rootgrp.close()

# ids

file_attrs_dict = flt_attrs()



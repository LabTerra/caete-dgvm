#-*-coding:utf-8-*-
# "CAETÊ"

"""
Copyright 2017- LabTerra

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
# import glob
# import time
# from math import ceil
# import pickle
# import multiprocessing as mp
# import psutil


from netCDF4 import num2date, MFDataset
import numpy as np
import plsgen as pls
from caete_module import global_par as gp
from caete_module import caete as model
from caete_module import photo as model_funcs

npls = gp.npls

# Mask for model execution
mask = np.load('random_amazon_mask_315.npy')

# Create the semi-random table of Plant Life Strategies
d_at = pls.table_gen(npls)

class data_in:

    """Program to extract CAETÊ input data from ISIMIP netCDF4 files"""


    def __init__(self, inputs_folder):

        self.root_dir = os.getcwd()

        try:
            self.files = sorted(glob.glob1(inputs_folder, '*.nc4'))
            self.NotWork = False

        except:
            self.files = None
            self.NotWork = True

        self.inputs_folder = inputs_folder
        self.metadata = {}
        self.celldata = {}


        self.varnames = ['tas', 'ps', 'pr', 'rsds', 'hurs']

#        self.varmask = [True, True, True, True, True, True, True, False, False]

        return None


    def _open_dts(self, varname):

        """Use de MFdataset class for read the array in a netCDF file"""

        #TODO Implement this mechanism without changing directory
        #     e.g. use os.path module

        assert varname in self.varnames, 'Unknown Variable: %s' % varname

        # Insert this block in a try/except to catch IO errors and return to the correct dir
        os.chdir(self.inputs_folder)

        fpath = varname + '_' + '*' + '.nc4'
        dt = MFDataset(fpath, 'r')

        os.chdir(self.root_dir)

        return dt


    def data_dict(self, varname=None, nx=None, ny=None, mask=None):

        """ An object that stores input data for caete

           Create a data dictionary (self.cell_data) with the input data for CAETÊ

           varname: str variable name (look self.varnames)
           nx; ny: int (x, y) point for the gridcell
           mask: np.array(dtype = bool) mask for creation of several gridcells 
        """

        if mask is None and varname is not None:

            assert varname in self.varnames, 'Unknown Variable: %s' % varname
            assert nx is not None and ny is not None, "In the '1 gricdell' mode ny and nx need to be present"

            k = varname + '_' + str(ny) + '-' + str(nx)

            cells_done = self.celldata.keys()

            if k in cells_done:
                return None

            print("Extracting var %s - grd %d, %d" % (varname, ny, nx))

            with self._open_dts(varname) as fh:
                self.celldata[k] = {'varname' : varname,
                                    'calendar' : fh.variables['time'].calendar,
                                    'time_unit': fh.variables['time'].units,
                                    'var_unit' : fh.variables[varname].units,
                                    'init_date': num2date(fh.variables['time'][0],
                                                             fh.variables['time'].units,
                                                             fh.variables['time'].calendar),
                                    'end_date' : num2date(fh.variables['time'][-1],
                                                             fh.variables['time'].units,
                                                             fh.variables['time'].calendar),
                                    'var_data' : fh.variables[varname][:, ny, nx],
                                    'time_data' : fh.variables['time'][:],
                                    'lat_unit' : fh.variables['lat'].units,
                                    'latitude' : fh.variables['lat'][:],
                                    'lon_unit' : fh.variables['lon'].units,
                                    'longitude' : fh.variables['lon'][:],
                                    'ny' : fh.dimensions['lat'].size,
                                    'nx' : fh.dimensions['lon'].size,
                                    'len' : fh.variables['time'][:].size}

        else:
            assert len(mask.shape) == 2
            assert ny is None and ny is None and varname is None

            dim1 = mask.shape[0]
            dim2 = mask.shape[1]

            #variables = [self.varnames[i] for i in self.varmask if i]
            print(self.varnames)

            for va in self.varnames:
                print("Extracting var %s" % va)
                with self._open_dts(va) as fh:
                    for Y in range(dim1):
                        for X in range(dim2):
                            if not mask[Y, X]:
                                print("...gridcell(%d, %d)" % (Y, X))
                                k = va + '_' + str(Y) + '-' + str(X)
                                self.celldata[k] = {'varname' : va,
                                                    'calendar' : fh.variables['time'].calendar,
                                                    'time_unit': fh.variables['time'].units,
                                                    'var_unit' : fh.variables[va].units,
                                                    'init_date': num2date(fh.variables['time'][0],
                                                                             fh.variables['time'].units,
                                                                             fh.variables['time'].calendar),
                                                    'end_date' : num2date(fh.variables['time'][-1],
                                                                             fh.variables['time'].units,
                                                                             fh.variables['time'].calendar),
                                                    'var_data' : fh.variables[va][:, Y, X],
                                                    'time_data' : fh.variables['time'][:],
                                                    'lat_unit' : fh.variables['lat'].units,
                                                    'latitude' : fh.variables['lat'][:],
                                                    'lon_unit' : fh.variables['lon'].units,
                                                    'longitude' : fh.variables['lon'][:],
                                                    'ny' : fh.dimensions['lat'].size,
                                                    'nx' : fh.dimensions['lon'].size,
                                                    'len' : fh.variables['time'][:].size}


class gridcell_dyn:
    
    """
    Defines the gridcell object - This object stores all the input data,
    the data comming from model runs for each grid point and all the metadata
    describing the life cycle of the gridcell
    """
    
    def __init__(self, x, y):
    
        """Construct the gridcell object"""
    
        # CELL Identifiers
        self.x = x                            # Grid point x coordinate 
        self.y = y                            # Grid point y coordinate
        self.name = str(y) + '-' + str(x)     # Name of the gridcell in the form of 'y-x'
        self.pos = (int(self.x), int(self.y)) # Tuple of (x, y) gridcell position
        
        #TODO
        self.neighbours = []                  # Neighbours pos of the gridcell 
        
        # Time related attributes and execution control
        #TODO
        self.complete = False   # Indicates if the gridcell already have runs
        self.run_id = None      # String defining the stage of execution (spinup, hist, proj)
        self.run_index = []
        self.index_dicts = {}   # Dict containing the metadata about executed runs
        self.ndays = 0          # Number of days for the current execution
        self.last_index = 0     # index of the last simulated day

        # This is the version that we want
        # Time attributes
        self.init_date = None
        self.time_index = None  # Array with the time stamps
        self.calendar = None    # Calendar name
        self.time_unit = None   # Time unit 
        self.start_date = None  
        self.end_date = None

        # Input data
        self.filled = False     # Indicates when the gridcell is filled with input data
        self.pr = None
        self.ps = None
        self.rsds = None
        self.tas = None
        self.rhs = None

        # Spinup data
        self.clin = None
        self.cfin = None
        self.cwin = None

        # OUTPUTS
        self.emaxm = None
        self.tsoil = None
        self.photo = None
        self.aresp = None
        self.npp = None
        self.lai = None
        self.clit = None  # 2D
        self.csoil = None # 3D
        self.hresp = None

        self.rcm = None
        self.f5 = None
        self.runom = None
        self.evapm = None
        self.wsoil = None
        self.rm = None
        self.rg = None

        self.cleaf = None
        self.cawood = None
        self.cfroot = None
        self.area = None
        self.wue = None
        self.cue = None
        self.cdef = None

        self.wfim = None
        self.gfim = None
        self.sfim = None
        self.dl_final = None
        self.dw_final = None
        self.dr_final = None
        self.clf = None
        self.caf = None
        self.cff = None

        self.nitro_min = None
        self.phop_lab = None
        self.vcmax = None
        self.specific_la = None
        self.nupt = None
        self.pupt = None
        self.litter_l = None
        self.cwd = None
        self.litter_fr = None
        self.lnr = None
        self.storage_pool = None
        
        # TODO 
        # Define the time control system: Storage, Save outputs, flush data from grid - cell - Try to mantain only metadata
        def update_neighbors(self):
            raise NotImplementedError

        def get_ndays(self):
            """return the number of days to be simulated"""
            pass


        def set_last_index(self):

            """Update the last_index attribute"""

            if len(self.run_index) == 0:
                return 1
            else:
                self.last_index = run_index[-1]
                return 0


# defining functions that drive model throghout the grid data & save outputs


def chunks(lst, chunck_size):
    """Yield successive n-sized chunks from lst."""

    for i in range(0, len(lst), chunck_size):
        yield lst[i:i + chunck_size]


## Prototype for a function that reads the input data and using a data_in instance 

# def init_gridcell(grd, mode='test'):
#     """ Organize model input """

#     # __future__ =
#     # implement the bigdict creation here! Only need to run one time. Save data as a .pkl
#     # implement Concurrent computation to acelerate this process - ONLY IN SOMBRERO
#     # not to run in your home computer
#     if mode == 'test':
#         dt1 = data_in('/home/jpdarela/Desktop/caete_jp/dev/CAETE_JP/model_inputs/dlds_daily/')

#         # fill dt1.metadata
#         dt1.data_dict('pr', grd.x, grd.y)
#         dt1.data_dict('ps', grd.x, grd.y)
#         dt1.data_dict('rsds', grd.x, grd.y)
#         dt1.data_dict('tas', grd.x, grd.y)
#         dt1.data_dict('hurs', grd.x, grd.y)
#         return dt1
#     return None


def init_caete_dyn(grd, dt1, name, init_date=None, end_date=None):
    
    """ Prepare the gridcell instance for run.
        TODO: Adapt to the dynamics structure

        grd: gridcell_dyn instance
        dt1: data_in instance
        name: str a name for the gridcell
        init_date; end_date: datetime objects"""

    # TODO

    # Check for the array size of caete_module.global_par.nt1
    # Control Time execution

    # Implement Mechanism to catch and execute model in a specified time span
    # Use module netcdf functions

    assert dt1 is not None, 'INVALID INPUT DATA'

    nx = str(grd.x)
    ny = str(grd.y)
    grd.name = name

    grd.pr = dt1.celldata['pr_' + ny + '-' + nx]['var_data'][:]
    grd.ps = dt1.celldata['ps_' + ny + '-' + nx]['var_data'][:]
    grd.rsds = dt1.celldata['rsds_' + ny + '-' + nx]['var_data'][:]
    grd.tas = dt1.celldata['tas_' + ny + '-' + nx]['var_data'][:]
    grd.rhs = dt1.celldata['hurs_' + ny + '-' + nx]['var_data'][:]

    assert grd.pr.size == grd.ps.size, 'ps is different from pr'
    assert grd.pr.size == grd.rsds.size, 'rsds is different from pr'
    assert grd.pr.size == grd.tas.size, 'tas is different from pr'
    assert grd.pr.size == grd.rhs.size, 'rhs is different from pr'

    grd.calendar = dt1.celldata['ps_' + ny + '-' + nx]['calendar']
    grd.time_index = dt1.celldata['ps_' + ny + '-' + nx]['time_data']
    grd.time_unit = dt1.celldata['ps_' + ny + '-' + nx]['time_unit']

    grd.filled = True


def spin(grd):
    global RUN
    RUN += int(gp.nt1)
    print("run %s" % RUN//gp.nt1)


def run_dyn(grd, at=np.copy(d_at, order='F')):
    """apply CAETÊ in dynamical mode in gridcell grd for nt1 days"""

    # starting variables
    # execute the model according to grd time variables

    ndays = 1

    RUN = 0
    w0 = np.zeros(shape=npls,) + 0.01
    g0 = np.zeros(shape=npls,)
    s0 = np.zeros(shape=npls,)

    dcl = np.zeros(npls,)
    dca = np.zeros(npls,)
    dcf = np.zeros(npls,)

    if grd.run_id == 'spinup':
        grd.clin, grd.cfin, grd.cwin = model_funcs.spinup2(0.5, d_at)
    
    # --------------------
    ##  Some assertions
    assert int(gp.nt1) > ndays, "different shapes in T i m e "
    # assert some at attributes

    # inputs to the subroutine caete_dyn
    # run,dt,w0,g0,s0, dcl,dca,dcf,prec,temp,p0,par,rhs&
    #    &,cleaf_ini,cawood_ini,cfroot_ini

    # subroutine caete_dyn(run,dt,w0,g0,s0,dcl,dca,dcf,prec,temp,p0,par,rhs&
    #  &,cleaf_ini,cawood_ini,cfroot_ini,emaxm,tsoil,photo_pft,aresp_pft&
    #  &,npp_pft,lai_pft,c_litter,c_soil,het_resp,rcm_pft,f51,runom_pft&
    #  &,evapm_pft,wsoil_pft,rm_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft&
    #  &,grid_area,wue,cue,cdef,wfim,gfim,sfim,dl_final,dw_final,dr_final&
    #  &,clf,caf,cff,nitro_min,phop_lab,vcmax,specific_la,nupt,pupt&
    #  &,litter_l,cwd,litter_fr,lnr,storage_pool)

    outputs = model.caete_dyn(grd.ndays, grd.x, grd.y, RUN, at, w0, g0, s0, dcl, dca, dcf, grd.pr,
                           grd.tas, grd.ps, grd.rsds,
                           grd.rhs, grd.clin, grd.cwin, grd.cfin)

    # TODO
    # atualizar attributos de tempo de grd
    # Cut off no_data tails
    # Include the cut of the time variables
    def cut_tail(in_array):
        return in_array[0:ndays-1]
    
    grd.emaxm = outputs[0]
    grd.tsoil = outputs[1]
    grd.photo = outputs[2]
    grd.aresp = outputs[3]
    grd.npp = outputs[4]
    grd.lai = outputs[5]
    grd.clit = outputs[6]  # 2D
    grd.csoil = outputs[7] # 3D
    grd.hresp = outputs[8]

    grd.rcm = outputs[9]
    grd.f5  = outputs[10]
    grd.runom = outputs[11]
    grd.evapm = outputs[12]
    grd.wsoil = outputs[13]
    grd.rm = outputs[14]
    grd.rg = outputs[15]

    grd.cleaf = outputs[16]
    grd.cawood = outputs[17]
    grd.cfroot = outputs[18]
    grd.area = outputs[19]
    grd.wue = outputs[20]
    grd.cue = outputs[21]
    grd.cdef = outputs[22]
        
    grd.wfim = outputs[23]
    grd.gfim = outputs[24]
    grd.sfim = outputs[25]
    grd.dl_final = outputs[26]
    grd.dw_final = outputs[27]
    grd.dr_final = outputs[28]
    grd.clf = outputs[29]
    grd.caf = outputs[30]
    grd.cff = outputs[31]

    grd.nitro_min = outputs[32]
    grd.phop_lab = outputs[33]
    grd.vcmax = outputs[34]
    grd.specific_la = outputs[35]
    grd.nupt = outputs[36]
    grd.pupt = outputs[37]
    grd.litter_l = outputs[38]
    grd.cwd = outputs[39]
    grd.litter_fr = outputs[40]
    grd.lnr = outputs[41]
    grd.storage_pool = outputs[42]


    def spin(grd1, r_number):
        """ Continuation Runs """

        outputs = model.caete_dyn(grd.x, grd.y, r_number, at, grd1.wfim, grd1.gfim, grd1.sfim,
                               grd1.dl_final, grd1.dw_final, grd1.dr_final, grd1.pr,
                               grd1.tas, grd1.ps, grd1.rsds, grd1.rhs,
                               grd1.clf, grd1.caf, grd1.cff)

        grd1.emaxm = np.hstack((grd1.emaxm, outputs[0]))
        grd1.tsoil = np.hstack((grd1.tsoil, outputs[1]))
        grd1.photo = np.hstack((grd1.photo, outputs[2]))
        grd1.aresp = np.hstack((grd1.aresp, outputs[3]))
        grd1.npp = np.hstack((grd1.npp, outputs[4]))
        grd1.lai = np.hstack((grd1.lai, outputs[5]))
        grd1.clit = np.hstack((grd1.clit, outputs[6]))
        grd1.csoil = np.hstack((grd1.csoil, outputs[7]))
        grd1.hresp = np.hstack((grd1.hresp, outputs[8]))

        grd1.rcm = np.hstack((grd1.rcm, outputs[9]))
        grd1.f5 = np.hstack((grd1.f5, outputs[10]))
        grd1.runom = np.hstack((grd1.runom, outputs[11]))
        grd1.evapm = np.hstack((grd1.evapm, outputs[12]))
        grd1.wsoil = np.hstack((grd1.wsoil, outputs[13]))
        grd1.rm = np.hstack((grd1.rm, outputs[14]))
        grd1.rg = np.hstack((grd1.rg, outputs[15]))

        grd1.cleaf = np.hstack((grd1.cleaf, outputs[16]))
        grd1.cawood = np.hstack((grd1.cawood, outputs[17]))
        grd1.cfroot = np.hstack((grd1.cfroot, outputs[18]))
        grd1.area = np.vstack((grd1.area, outputs[19]))
        grd1.wue = np.hstack((grd1.wue, outputs[20]))
        grd1.cue = np.hstack((grd1.cue, outputs[21]))
        grd1.cdef = np.hstack((grd1.cdef, outputs[22]))
        
        grd1.wfim = outputs[23]
        grd1.gfim = outputs[24]
        grd1.sfim = outputs[25]
        grd1.dl_final = outputs[26]
        grd1.dw_final = outputs[27]
        grd1.dr_final = outputs[28]
        grd1.clf = outputs[29]
        grd1.caf = outputs[30]
        grd1.cff = outputs[31]

        grd1.nitro_min = np.hstack((grd1.nitro_min, outputs[32]))
        grd1.phop_lab = np.hstack((grd1.phop_lab, outputs[33]))
        grd1.vcmax = np.hstack((grd1.vcmax, outputs[34]))
        grd1.specific_la = np.hstack((grd1.specific_la, outputs[35]))
        grd1.nupt = np.hstack((grd1.nupt, outputs[36]))
        grd1.pupt = np.hstack((grd1.pupt, outputs[37]))
        grd1.litter_l = np.hstack((grd1.litter_l, outputs[38]))
        grd1.cwd = np.hstack((grd1.cwd, outputs[39]))
        grd1.litter_fr = np.hstack((grd1.litter_fr, outputs[40]))
        grd1.lnr = np.hstack((grd1.lnr, outputs[41]))
        grd1.storage_pool = np.hstack((grd1.storage_pool, outputs[42]))

# GAMBIARRA NERVOSA

    print('\n--run 1\n')
    RUN += int(gp.nt1)
    spin(grd, RUN)

    print('\n--run 2\n')
    RUN += int(gp.nt1)
    spin(grd, RUN)

    print('\n--run 3\n')
    RUN += int(gp.nt1)
    spin(grd, RUN)

    print('\n--run 4\n')
    RUN += int(gp.nt1)
    spin(grd, RUN)

    print('\n--run 5\n')
    RUN += int(gp.nt1)
    spin(grd, RUN)


    # print('run 6\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 7\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 8\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 9\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 10\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 11\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 12\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 13\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 14\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 15\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 16\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 17\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 18\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 19\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 20\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

   # print('run 6\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 7\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 8\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 9\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)

    # print('run 10\n')
    # RUN += int(gp.nt1)
    # spin(grd, RUN)
        # cont_run(grd,RUN)
# #---------------           

def model_flush(grd):
    """collect garbage grd"""
    pass

def rm_apply(gridcell_obj):
    """apply the model--- to be employed in multiprocessing"""
    pass

if __name__ == "__main__":
    pass

    # Abaixo a implementação do caete estacionário
    
    # if os.path.exists('/home/ecologia/results7'):
    #     pass
    # else:
    #    os.mkdir('/home/ecologia/results7')

    # if os.path.exists('/home-ext/ecologia/outputs_nc'):
    #     pass
    # else:
    #    os.mkdir('/home-ext/ecologia/outputs_nc')

    # if os.path.exists('/home-ext/ecologia/pkl'):
    #     pass
    # else:
    #     os.mkdir('/home-ext/ecologia/pkl')

    # # iniciando o modelo para todas (land) as celulas do grid
    # # este bloco de código cria uma lista que contém n instâncias da
    # # classe gridcell (L-579) e depois aplica o modelo a cada uma delas.
    # # Preferi deixar esta classe o mais limpa possível, pois as suans instầncias são modificada em
    # # paralelo e para isso precisam ser pickladas (veja o modulo pickle e também o modulo multiprocessing)
    
    # # Por isso as funções que iniciam e aplicam o modelo a cada instancia de gridcell
    # # não são métodos propriamente ditos e sim, funções que manipulam estas instâncias.
    
    # land_data = []
    # id_n = 0
    # print('init caete --- %d = npls' % npls, end='---> ')
    # print(time.ctime())
    # for Y in range(ny):
    #     for X in range(nx):
    #         if not mask[Y][X]:
    #             id_n += 1
    #             grd_cell = gridcell(X, Y)
    #             init_caete(grd_cell)
    #             land_data.append(grd_cell)

    # # del input arrays
    # del(global_pr)
    # del(global_ps)
    # del(global_tas)
    # del(global_rhs)
    # del(global_rsds)
    # del(npp_init)
    # x = dir()
    # print('dir() before model run')
    # print(x)

    # # divide land_data when data is too big - npls > i
    # if npls < 5:
        
    #     with mp.Pool(processes=24,maxtasksperchild=20) as p:
    #         result = p.map(rm_apply, land_data)
            
    #     print('\nModelo aplicado a %d localidades-- %d PLSs' % (id_n,npls), end='--->')
    #     print(time.ctime())

    #     print('\nCALCULANDO SOMAS', end='--->')
    #     print(time.ctime())
   
    #     #for grd in result:
    #     #    grd_dict(grd)

    #     print('\nSalvando resultados')
    #     for v in varlist:
    #         print(v, end='-')
    #         assemble(result, v)

    #     print('terminado', end='---: ')
    #     print(time.ctime())

    # else:
    #     id = 0
    #     files = []
    #     f_result = []
    #     for lst in chunks(land_data, ceil(len(land_data)/10)):
    #         with mp.Pool(processes=24, maxtasksperchild=20) as p:
    #             result = p.map(rm_apply, lst)

    #         pkl_name = '/home-ext/ecologia/pkl/result' + '_' + str(id)
    #         with open(pkl_name, 'wb') as fh:
    #             pickle.dump(result,fh)
    #         files.append(pkl_name)
    #         id += 1

    #     print('\nModelo aplicado a %d localidades-- %d PLSs' % (len(land_data),npls), end='--->')
    #     del(land_data)
    #     print(time.ctime())

    #     for fh in files:
    #         with open(fh,'rb') as fh2:
    #             tmp = pickle.load(fh2)
    #         f_result += tmp[:]

    #     del(tmp)
    #     x = dir()

    #     print('dir() in grd_result + assemble ')
    #     print(x)

    #     print('\nSalvando resultados')

    #     for v in varlist:
    #         print(v, end='-')
    #         assemble(f_result, v)

    #     print('terminado', end='---: ')
    #     print(time.ctime())

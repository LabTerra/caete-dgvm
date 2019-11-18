import os
import glob
from netCDF4 import num2date, MFDataset


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

        # TODO Implement this mechanism without changing directory
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

           Create a data dictionary (self.celldata) with the input data for CAETÊ

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

            # TODO (jp) as variáveis que são extraidas aqui nesta seção são redundantes. Elas são comuns á todas
            # as células de grid. Acho que a melhor solução é criar un atributo para armazrnas as variáveis que são redundantes
            with self._open_dts(varname) as fh:
                self.celldata[k] = {
                    'varname':
                    varname,
                    'calendar':
                    fh.variables['time'].calendar,
                    'time_unit':
                    fh.variables['time'].units,
                    'var_unit':
                    fh.variables[varname].units,
                    'init_date':
                    num2date(fh.variables['time'][0],
                             fh.variables['time'].units,
                             fh.variables['time'].calendar),
                    'end_date':
                    num2date(fh.variables['time'][-1],
                             fh.variables['time'].units,
                             fh.variables['time'].calendar),
                    'var_data':
                    fh.variables[varname][:, ny, nx],
                    'time_data':
                    fh.variables['time'][:],
                    'lat_unit':
                    fh.variables['lat'].units,
                    'latitude':
                    fh.variables['lat'][:],
                    'lon_unit':
                    fh.variables['lon'].units,
                    'longitude':
                    fh.variables['lon'][:],
                    'ny':
                    fh.dimensions['lat'].size,
                    'nx':
                    fh.dimensions['lon'].size,
                    'len':
                    fh.variables['time'][:].size
                }

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
                                self.celldata[k] = {
                                    'varname':
                                    va,
                                    'calendar':
                                    fh.variables['time'].calendar,
                                    'time_unit':
                                    fh.variables['time'].units,
                                    'var_unit':
                                    fh.variables[va].units,
                                    'init_date':
                                    num2date(fh.variables['time'][0],
                                             fh.variables['time'].units,
                                             fh.variables['time'].calendar),
                                    'end_date':
                                    num2date(fh.variables['time'][-1],
                                             fh.variables['time'].units,
                                             fh.variables['time'].calendar),
                                    'var_data':
                                    fh.variables[va][:, Y, X],
                                    'time_data':
                                    fh.variables['time'][:],
                                    'lat_unit':
                                    fh.variables['lat'].units,
                                    'latitude':
                                    fh.variables['lat'][:],
                                    'lon_unit':
                                    fh.variables['lon'].units,
                                    'longitude':
                                    fh.variables['lon'][:],
                                    'ny':
                                    fh.dimensions['lat'].size,
                                    'nx':
                                    fh.dimensions['lon'].size,
                                    'len':
                                    fh.variables['time'][:].size
                                }

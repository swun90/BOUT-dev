# Provides a class BoutData which makes access to code
# inputs and outputs easier
#
# Intended to make access to BOUT++ data sets easier, 
# with one eye on the API needed for PyXPAD graphical interface.
# 

import os
import glob

from collect import collect

try:
    from boututils import DataFile
except ImportError:
    print("ERROR: boututils.DataFile couldn't be loaded")
    raise

class BoutDim(object):
    """
    Represents a dimension
    
    """
    name = ""
    label = ""
    units = ""
    data = None
    
class BoutVar(object):
    """
    Data representing an output variable
    
    Data members
    ============
    
    data   NumPy array of the data

    dim     A list of dimensions, each of which contains:
      - name 
      - label  Short axis label (e.g. "Time (sec)")
      - units  (e.g. "s")
      - data   Axis values (NumPy array)
    order   Index of time dimension
    time    A shortcut to the time data (dim[order].data). May be None
    """
    def __init__(self):
        self.data = None
        
        

        self.dim = [t]
        self.order = 0
        self.time = t.data
    
class BoutData(object):
    """
    Data members
    ============

    dim     A list of dimensions, each of which contains:
      - name 
      - label  Short axis label (e.g. "Time (sec)")
      - units  (e.g. "s")
      - data   Axis values (NumPy array)

    
    options   A dictionary of settings
    """
    def __init__(self, path=".", prefix="BOUT.dmp"):
        """
        Initialise BoutData object
        """
        self._path = path
        self._prefix = prefix
        
        # Label for this data
        self.label = path

        # Check that the path contains some data

        # Dimensions
        t = BoutDim()
        t.name  = "Time"
        t.label = "Time"
        t.data  = collect("t_array", path=path, prefix=prefix)
        
        self.dimensions = [t]
        
        # Options
        self.options = {}
        #with open( os.path.join(self._path, "BOUT.inp"), "r") as f:
        
        # Available variables
        self.varNames = []
        
        file_list = glob.glob(os.path.join(path, prefix+"*.nc"))
        file_list.sort()
        if file_list == []:
            raise ValueError("ERROR: No data files found")
        
        with DataFile(file_list[0]) as f:
            # Get variable names
            self.varNames = f.keys()
            
            # Dimensions
            
        
    def keys(self):
        """
        Return a list of available variable names
        """
        return self.varNames
    
    def read(self, name):
        """
        Read a variable of given name
        
        Returns a BoutVar object
        """
        try:
            var = BoutVar()
            
            # Collect the data from the repository
            var.data = collect(name, path=self._path, prefix=self._prefix)
            
            var.name = name
            var.source = self._path

            return var
        except:
            return None
    
    def __getitem__(self, name):
        """
        A shortcut, which just returns the data
        
        BoutData["name"] is equivalent to BoutData.read(name).data

        since this is what is required most of the time
        
        """
        try:
            # Collect the data from the repository
            data = collect(name, path=self._path, prefix=self._prefix)

            return data
        except:
            return None

        

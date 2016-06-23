# save data into NetCDF file.
# © 2016 SUN Chuankui <sunck@pku.edu.cn>
from boututils import DataFile

def save2nc(file_name,**var):
  """save2nc(file_name, varName=var, ...)
    save variables to netCDF file. 

    'file_name' is the exported file name, and the suffix ".nc" will be append automatically.
    'varName' is the name saved in exported files,
    'var' is the variable saved. varName and var do NOT need the quotation mark.

    E.g. save2nc('file',t=t,x=psi,P0=P0_origin) will save variables t as t, psi as x, P0_origin as P0 into a netcdf file named 'file.nc'.
  """
  if file_name[-3:] != '.nc':
    file_name += '.nc'
  f = DataFile(file_name, write=True, create=True)
  for v in var:
    try:
      varg = var[v]
      f.write(v,varg)
    except:
      print 'check the source code "save2nv.py" for more infomation'
  f.close()

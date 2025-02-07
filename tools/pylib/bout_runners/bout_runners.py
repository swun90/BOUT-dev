#!/usr/bin/env python

"""Classes for running one or several mpi-runs with BOUT++ at once.
   Read the docstring of 'basic_runner', or refer to the user manual of
   BOUT++ for more info. Examples can be found in
   BOUT/examples/bout_runners_example."""

# NOTE: This document uses folding. A hash-symbol followed by three {'s
# denotes the start of a fold, and a hash-symbol followed by three }'s
# denotes the end of a fold
__authors__ = 'Michael Loeiten'
__email__   = 'mmag@fysik.dtu.dk'
__version__ = '1.012'
__date__    = '2016.02.17'

import os
import re
import itertools
import glob
import timeit
import datetime
import shutil
from numbers import Number
import numpy as np
from boututils.run_wrapper import shell, launch, getmpirun
from boututils.options import BOUTOptions
from boututils.datafile import DataFile

#{{{class basic_runner
# As a child class uses the super function, the class must allow an
# object as input
class basic_runner(object):
#{{{docstring
    """Class for mpi running one or several runs with BOUT++.
    Calling self.execute_runs() will run your BOUT++ program with the possible
    combinations given in the member data using the mpi runner.

    Before each run, a folder system, based on the member data, rooted
    in self._directory, will be created. The BOUT.inp of self._directory
    is then copied to the execution folder.

    A log-file for the run is stored in self._directory

    By default self._directory = 'data', self._nproc = 1 and
    self._allow_size_modification = False

    self._program_name is by default set to the same name as any .o files in the
    folder where an instance of the object is created. If none is found
    the creator tries to run make.

    All other data members are set to None by default.

    The data members will override the corresponding options given in
    self._directory/BOUT.inp.

    See the doctring of the constructor (__int__) for options.
    See BOUT/examples/bout_runners_example for examples."""
#}}}

#{{{__init__
    def __init__(self,\
                 nproc        = 1,\
                 directory    = 'data',\
                 solver       = None,\
                 mms          = None,\
                 atol         = None,\
                 rtol         = None,\
                 mxstep       = None,\
                 grid_file    = None,\
                 nx           = None,\
                 ny           = None,\
                 nz           = None,\
                 zperiod      = None,\
                 zmin         = None,\
                 zmax         = None,\
                 dx           = None,\
                 dy           = None,\
                 dz           = None,\
                 MXG          = None,\
                 MYG          = None,\
                 NXPE         = None,\
                 ixseps1      = None,\
                 ixseps2      = None,\
                 jyseps1_1    = None,\
                 jyseps1_2    = None,\
                 jyseps2_1    = None,\
                 jyseps2_2    = None,\
                 symGlobX     = None,\
                 symGlobY     = None,\
                 ddx_first    = None,\
                 ddx_second   = None,\
                 ddx_upwind   = None,\
                 ddx_flux     = None,\
                 ddy_first    = None,\
                 ddy_second   = None,\
                 ddy_upwind   = None,\
                 ddy_flux     = None,\
                 ddz_first    = None,\
                 ddz_second   = None,\
                 ddz_upwind   = None,\
                 ddz_flux     = None,\
                 nout         = None,\
                 timestep     = None,\
                 additional   = None,\
                 series_add   = None,\
                 restart      = None,\
                 restart_from = None,\
                 cpy_source   = None,\
                 cpy_grid     = None,\
                 sort_by      = None,\
                 make         = None,\
                 allow_size_modification = False):
        #{{{docstring
        """The constructor of the basic_runner.

        All the member data is set to None by default. If the
        data members are not set, the values from BOUT.inp will be used.
        The exception is nproc (default = 1), directory (default =
        'data') and allow_size_modification (default = False), which
        always needs to be set.

        Input:
        nproc        -    The number of processors to use in the mpirun (int)
        directory    -    The directory of the BOUT.inp file (str)
        solver       -    The solver to be used in the runs (str or
                          iterable)
        mms          -    Whether or not mms should be run (bool)
        atol         -    Absolute tolerance (number or iterable)
        rtol         -    Relative tolerance (number or iterable)
        mxstep       -    Max internal step pr output step (int or
                          iterable)
        grid_file    -    The grid file (str or iterable)
        nx           -    Number of nx in the run (int or iterable)
        ny           -    Number of ny in the run (int or iterable)
        nz           -    Number of nz in the run (int or iterable)
        zperiod      -    Domain size in  multiple of fractions of 2*pi
                          (int or iterable)
        zmin         -    Minimum range of the z domain
        zmax         -    Maximum range of the z domain
        dx           -    Grid size in the x direction (Number or iterable)
        dy           -    Grid size in the x direction (Number or iterable)
        dz           -    Grid size in the x direction (Number or iterable)
        MXG          -    The number of guard cells in the x direction
                          (int)
        MYG          -    The number of guard cells in the y direction
                          (int)
        NXPE         -    Numbers of processors in the x direction
        ixseps1      -    Separatrix location for 'upper' divertor (int
                          or iterable)
        ixseps2      -    Separatrix location for 'lower' divertor (int
                          or iterable)
        jyseps1_1    -    Branch cut location 1_1 [see user's manual for
                          details] (int or iterable)
        jyseps1_2    -    Branch cut location 1_2 [see user's manual for
                          details] (int or iterable)
        jyseps2_1    -    Branch cut location 2_1 [see user's manual for
                          details] (int or iterable)
        jyseps2_2    -    Branch cut location 2_2 [see user's manual for
                          details] (int or iterable)
        symGlobX     -    symmetricGLobalX: x defined symmetrically between
                          0 and 1 (bool)
        symGlobY     -    symmetricGLobalX: y defined symmetrically (bool)
        ddx_first    -    Method used for for first ddx terms (str or iterable)
        ddx_second   -    Method used for for second ddx terms (str or iterable)
        ddx_upwind   -    Method used for for upwind ddx terms (str or iterable)
        ddx_flux     -    Method used for for flux ddx terms (str or iterable)
        ddy_first    -    Method used for for first ddy terms (str or iterable)
        ddy_second   -    Method used for for second ddy terms (str or iterable)
        ddy_upwind   -    Method used for for upwind ddy terms (str or iterable)
        ddy_flux     -    Method used for for flux ddy terms (str or iterable)
        ddz_first    -    Method used for for first ddz terms (str or iterable)
        ddz_second   -    Method used for for second ddz terms (str or iterable)
        ddz_upwind   -    Method used for for upwind ddz terms (str or iterable)
        ddz_flux     -    Method used for for flux ddz terms (str or iterable)
        nout         -    Number of outputs stored in the *.dmp.* files
                          (int or iterable)
        timestep     -    The time between each output stored in the
                          *.dmp.* files (int or iterable)
        additional   -    Additional option for the run given on the form
                          ('section_name','variable name', values) or as
                          iterable on the same form, where values can be
                          any value or string or an iterable of those
        series_add   -    The same as above, with the exception that
                          no combination will be performed between the elements
                          during a run
        restart      -    Wheter or not to use the restart files
                          ('overwrite' or 'append')
        restart_from -    Path to restart from
        cpy_source   -    Wheter or not to copy the source files to the
                          folder of the *.dmp.* files (bool)
        cpy_grid     -    Wheter or not to copy the grid files to the
                          folder of the *.dmp.* files (bool)
        sort_by      -    Defining what will be the fastest running
                          variable in the run, which can be useful if one
                          for example would like to 'group' the runs before
                          sending it to a post processing function (see
                          the docstring of the run function for more
                          info). The possibilities are
                         'spatial_domain',
                         'temporal_domain',
                         'solver',
                         'ddx_first',
                         'ddx_second',
                         'ddx_upwind',
                         'ddx_flux',
                         'ddy_first',
                         'ddy_second',
                         'ddy_upwind',
                         'ddy_flux',
                         'ddz_first',
                         'ddz_second',
                         'ddz_upwind',
                         'ddz_flux',
                         any 'variable_name' from additional or series_add
                         an iterable consisting of several of these. If
                         an iterable is given, then the first element is
                         going to be the fastest varying variable, the
                         second element is going to be the second fastest
                         varying variable and so on.
        make         -   Whether or not to make the program (bool)

        allow_size_modification - Whether or not to allow bout_runners
                                  modify nx and ny in order to find a
                                  valid split of the domain (bool)
        """
        #}}}

        # Setting the member data
        self._nproc           = nproc
        self._directory       = directory
        self._solver          = self._set_member_data(solver)
        self._mms             = mms
        self._atol            = self._set_member_data(atol)
        self._rtol            = self._set_member_data(rtol)
        self._mxstep          = self._set_member_data(mxstep)
        self._grid_file       = self._set_member_data(grid_file)
        self._nx              = self._set_member_data(nx)
        self._ny              = self._set_member_data(ny)
        self._nz              = self._set_member_data(nz)
        self._zperiod         = self._set_member_data(zperiod)
        self._zmin            = self._set_member_data(zmin)
        self._zmax            = self._set_member_data(zmax)
        self._dx              = self._set_member_data(dx)
        self._dy              = self._set_member_data(dy)
        self._dz              = self._set_member_data(dz)
        self._MXG             = MXG
        self._MYG             = MYG
        self._NXPE            = self._set_member_data(NXPE)
        self._ixseps1         = self._set_member_data(ixseps1)
        self._ixseps2         = self._set_member_data(ixseps2)
        self._jyseps1_1       = self._set_member_data(jyseps1_1)
        self._jyseps1_2       = self._set_member_data(jyseps1_2)
        self._jyseps2_1       = self._set_member_data(jyseps2_1)
        self._jyseps2_2       = self._set_member_data(jyseps2_2)
        self._symGlobX        = symGlobX
        self._symGlobY        = symGlobY
        self._ddx_first       = self._set_member_data(ddx_first)
        self._ddx_second      = self._set_member_data(ddx_second)
        self._ddx_upwind      = self._set_member_data(ddx_upwind)
        self._ddx_flux        = self._set_member_data(ddx_flux)
        self._ddy_first       = self._set_member_data(ddy_first)
        self._ddy_second      = self._set_member_data(ddy_second)
        self._ddy_upwind      = self._set_member_data(ddy_upwind)
        self._ddy_flux        = self._set_member_data(ddy_flux)
        self._ddz_first       = self._set_member_data(ddz_first)
        self._ddz_second      = self._set_member_data(ddz_second)
        self._ddz_upwind      = self._set_member_data(ddz_upwind)
        self._ddz_flux        = self._set_member_data(ddz_flux)
        self._nout            = self._set_member_data(nout)
        self._timestep        = self._set_member_data(timestep)
        self._additional      = additional
        self._series_add      = series_add
        self._restart         = restart
        self._restart_from    = restart_from
        self._cpy_source      = cpy_source
        self._cpy_grid        = cpy_grid
        self._sort_by         = self._set_member_data(sort_by)
        self._make            = make
        self._allow_size_modification = allow_size_modification

        # Make some space to distinguish from the rest of the terminal
        print("\n")

        # Initializing self._warnings and self._error
        # self._warnings will be filled with warnings
        # self._errors will be filled with errors
        # The warnings and errors will be printed when the destructor is called
        self._warnings = []
        self._errors   = []

        # Check if make is a boolean
        if self._make is not None:
            if type(self._make) != bool:
                self._errors.append("TypeError")
                raise TypeError("make must be boolean if set")

        # Set self._program_name from the *.o file. Make the program if
        # the *.o file is not found
        self._set_program_name()

        # Make the file if make is True
        if self._make:
            self._run_make()

        # Obtain the MPIRUN
        self._MPIRUN = getmpirun()

        # The run type is going to be written in the run.log file
        self._run_type = 'basic'

        #{{{ Set self._additional and self._series_add correctly
        # self._additional must be on a special form (see
        # basic_error_checker).
        if self._additional is not None:
            if not(hasattr(self._additional, "__iter__")) or\
               (type(self._additional) == str) or\
               (type(self._additional) == dict):
                # Put additional as a double iterable
                self._additional = [(self._additional)]
            else:
                if not(hasattr(self._additional[0], "__iter__")) or\
                   (type(self._additional[0]) == str) or\
                   (type(self._additional) == dict):
                    # Put self._additional as an iterable
                    self._additional = [self._additional]
        # Do the same for series_add
        if self._series_add is not None:
            if not(hasattr(self._series_add, "__iter__")) or\
               (type(self._series_add) == str) or\
               (type(self._series_add) == dict):
                # Put series_add as a double iterable
                self._series_add = [(self._series_add)]
            else:
                if not(hasattr(self._series_add[0], "__iter__")) or\
                   (type(self._series_add[0]) == str) or\
                   (type(self._series_add) == dict):
                    # Put self._series_add as an iterable
                    self._series_add = [self._series_add]
        #}}}

        # Check that nproc is given correctly
        if type(self._nproc) != int:
            message  = "nproc is of wrong type\n"+\
                       "nproc must be given as an int"
            self._errors.append("TypeError")
            raise TypeError(message)

        #{{{ Set NYPE from NXPE and nproc
        if self._NXPE is not None:
            # Make self._NYPE as an appendable list
            self._NYPE = []

            # Check that NXPE is of correct type
            check_if_int = [\
                    (self._NXPE,  'NXPE'),\
                   ]
            self._check_for_correct_type(var = check_if_int,\
                                          the_type = int,\
                                          allow_iterable = True)

            #Check that NXPE and nproc is consistent
            for cur_NXPE in self._NXPE:
                if (self._nproc % cur_NXPE) != 0:
                    self._errors.append("RuntimeError")
                    message = "nproc =" + str(self._nproc) + " not divisible by"+\
                              " NXPE = " + str(cur_NXPE) +\
                              " (the number of processors in the x direction)"
                    raise RuntimeError(message)

                # Append NYPE
                self._NYPE.append(int(self._nproc/cur_NXPE))
        else:
            self._NYPE = None
        #}}}

        # Check if the instance is set correctly
        self._check_for_basic_instance_error()
#}}}

#{{{__del__
    def __del__(self):
        """The destructor will print all the warning and error messages"""

        # Switch to see if error occured
        error_occured = False

        # If errors occured
        if len(self._errors) > 0:
            message = "! A " + self._errors[0] + " occurred. !"
            # Find the boarder length
            len_boarder = len(message)
            # Print the message
            print("\n"*2 + "!"*len_boarder)
            print(message)
            print('!'*len_boarder + "\n"*2)
            error_occured = True
        if len(self._warnings) > 0:
            print('\n'*3 + 'The following WARNINGS were detected:')
            print('-'*80)
            for warning in self._warnings:
                print(warning + '\n')
            print('-'*80 + '\n'*3)
        elif len(self._warnings) > 0 and not(error_occured):
            print('\n'*3 + ' ' + '~'*69)
            print("| No WARNINGS detected before instance destruction in"+\
                  " 'bout_runners'. |")
#}}}

#{{{execute_runs
    def execute_runs(self,\
                     remove_old = False,\
                     post_processing_function = None,\
                     post_process_after_every_run = True,\
                     **kwargs):
        #{{{docstring
        """
        Makes a run for each of the combination given by the member data.

        Input
        remove_old                  - boolean telling whether old run
                                      files should be deleted or not
        post_processing_function    - a function to be called after
                                      one or several run. This function
                                      must accept the string of
                                      self._dmp_folder if
                                      post_process_after_each_run is
                                      True, and a list of dmp folders if
                                      post_process_after_each_run is
                                      False
        post_process_after_each_run - boolean telling whether
                                      post_processing_function
                                      should be called after each run
                                      (if True), or after the number of
                                      runs decided by self._sort_by
                                      (see the constructor of
                                      basic_runner for more info)
        **kwargs                    - parameters to be passed to the
                                      post_processing_function
        """
        #}}}

        # Check for errors in the run function
        self._error_check_for_run_input(remove_old,\
                                         post_processing_function,\
                                         post_process_after_every_run)

        # Create the run log
        self._create_run_log()

        # We check that the given combination of nx and ny is
        # possible to perform with the given nproc
        if (self._nx is not None) and (self._ny is not None):
            self._get_correct_domain_split()

        # Get the combinations of the member functions
        possibilities = self._get_possibilities()
        combinations = self._get_combinations(possibilities)

        # If we are not running the post processing function after every
        # run, make an appendable list over all the runs which will be
        # passed as an input parameter to the post processing function
        if not(post_process_after_every_run):
            list_of_dmp_folders = []

        # Print either 'now running' or 'now submitting'
        self._print_run_or_submit()

        # Set self._len_group if post_processing_function is set, but
        # self._sort_by is None
        if (post_processing_function is not None) and\
           (not(post_process_after_every_run)) and\
           (self._len_group == None):
               # self._len_group is to a number by _get_swapped_input_list
               # (which is called if self._sort_by is not None)
               # If there are no sorting, self._len_group will be None
               # We will make self._len_group the length of the
               # number of runs here
               self._len_group = len(combinations)

        # The run
        for run_no, combination in enumerate(combinations):

            # Get the folder to store the data
            skip_run = self._prepare_dmp_folder(combination)
            if skip_run:
                # Skip this run
                continue

            if remove_old:
                # Remove old data
               self._remove_data()

            # Copy the grid (if any) if cpy_grid files is True
            if (self._cpy_grid) and (self._grid_file is not None):
                combination_list = combination.split()
                # Loop through all the combinations
                for elem in combination_list:
                    # Find the grid
                    if elem[0:4] == 'grid':
                        # Remove grid=, so that only the path remains
                        cur_grid = elem.replace('grid=', '')
                        # Copy the grid file
                        shutil.copy2(cur_grid, self._dmp_folder)

            # Check if the run has been performed previously
            do_run = self._check_if_run_already_performed()
            # Do the actual runs
            if do_run:
                # Call the driver for a run
                self._run_driver(combination, run_no)

            # If we would like to call a post_processing function
            if post_processing_function is not None:
                if post_process_after_every_run:
                    # Call the post processing function
                    self._call_post_processing_function(\
                            function = post_processing_function,\
                            folders  = self._dmp_folder,\
                            **kwargs)
                else:
                    # Append the dmp folder to the list of dmp folders
                    list_of_dmp_folders.append(self._dmp_folder)
                    # If the run_no+1 is divisible by self._len_group
                    if ((run_no+1) % self._len_group == 0):
                        # Call the post processing function
                        self._call_post_processing_function(\
                                function = post_processing_function,\
                                folders  = list_of_dmp_folders,\
                                **kwargs)
                        # Reset the list_of_dmp_folders
                        list_of_dmp_folders = []
#}}}

#{{{_run_driver
    def _run_driver(self, combination, run_no):
        """The machinery which actually performs the run and eventual
        post_processing"""

        # Get the time when the run starts
        start = datetime.datetime.now()
        # Do the run
        output, run_time = self._single_run(combination)
        # Print info to the log file for the runs
        self._append_run_log(start, run_no, run_time)
        print('\n')

        # Perform eventual post-processing
#}}}

#{{{ Functions called by the constructor
#{{{_set_member_data
    def _set_member_data(self, input_parameter):
        """Returns the input_parameter as a list if it is different than None,
        and if it is not iterable"""

       # If the input_data is not set, the value in BOUT.inp will
       # be used
        if input_parameter is not None:
            # If the input_data is not an iterable, or if it is a
            # string: Put it to a list
            if not(hasattr(input_parameter, "__iter__")) or\
               (type(input_parameter)) == str:
                input_parameter = [input_parameter]

        return input_parameter
#}}}

#{{{_set_program_name
    def _set_program_name(self):
        """Set self._program_name from the *.o file. Make the program if
        the *.o file is not found"""

        # Find the *.o file
        o_files = glob.glob("*.o")
        if len(o_files) > 0:
            # Pick the first instance as the name
            self._program_name = o_files[0].replace('.o', '')
        else:
            # Check if there exists a make
            make_file = glob.glob("*make*")
            if len(make_file) > 0:
                # Run make
                self._run_make()
                # Set the make flag to False, so it is not made again
                self._make = False
                # Search for the .o file again
                o_files = glob.glob("*.o")
                if len(o_files) > 0:
                    self._program_name = o_files[0].replace('.o', '')
                else:
                    self._program_name = False
                    message = 'The constructor could not make your'+\
                              ' program'
                    self._errors.append("RuntimeError")
                    raise RuntimeError(message)
            else:
                self._errors.append("RuntimeError")
                raise RuntimeError("No make file found in current" +\
                                   " directory")
#}}}

#{{{_check_for_basic_instance_error
    def _check_for_basic_instance_error(self):
        """Check if there are any type errors when creating the object"""

        #{{{Check if directory has the correct type
        if type(self._nproc) != int:
            message  = "nproc is of wrong type\n"+\
                       "nproc must be given as an int"
            self._errors.append("TypeError")
            raise TypeError(message)
        if type(self._directory) != str:
            message  = "directory is of wrong type\n"+\
                       "directory must be given as a str"
            self._errors.append("TypeError")
            raise TypeError(message)
        #}}}

        #{{{Check if MXG and MYG has the correct type
        # Check if MXG and MYG is given as a single int
        # One should not be allowed to specify MXG and MYG as an
        # iterable, as MXG is used to find the correct split, and
        # because it in principle could be incompatible with the method
        # (first, second etc.) used
        check_if_int = [\
                        (self._MXG,  'MXG'),\
                        (self._MYG,  'MYG'),\
                       ]
        self._check_for_correct_type(var = check_if_int,\
                                      the_type = int,\
                                      allow_iterable = False)
        #}}}

        #{{{Check if BOUT.inp exsists in the self._directory
        # Check if there are any BOUT.inp files in the self._directory
        inp_file = glob.glob(self._directory + "/BOUT.inp")
        if len(inp_file) == 0:
            self._errors.append("RuntimeError")
            raise RuntimeError("No BOUT.inp files found in '" +\
                                self._directory + "'")
        #}}}

        #{{{Check grid_file are strings, that they exsist, and one can sort
        if self._grid_file is not None:
            # Set a variable which is has length over one if the test fails
            not_found = []
            if type(self._grid_file) == str:
                # See if the grid_file can be found
                grid_file = glob.glob(self._grid_file)
                # The grid_file cannot be found
                if len(grid_file) == 0:
                    not_found.append(self._grid_file)
            # If several grid files are given
            elif hasattr(self._grid_file, "__iter__"):
                for elem in self._grid_file:
                    # See if the grid_file can be found
                    grid_file = glob.glob(elem)
                    # The grid_file cannot be found
                    if len(grid_file) == 0:
                        not_found.append(elem)
            if len(not_found) > 0:
                message =  "The following grid files were not found\n"
                message += "\n".join(not_found)
                self._errors.append("RuntimeError")
                raise RuntimeError(message)
            if (self._sort_by is not None) and ('grid_file' in self._sort_by):
                # Set a success flag
                success = True
                # The start name of the files
                start_name = 'grid_file'
                # Check if grid file is iterable
                if hasattr(self._grid_file, "__iter__"):
                    for grid in grid_file:
                        if grid[0:len(start_name)] != start_name:
                            success = False
                else:
                    # Only one grid file
                    if self._grid_file[0:len(start_name)] != start_name:
                        success = False
                if not(success):
                    message =  "The name of the grid file must start with"+\
                               " 'grid_file' in order to sort by them."
                    self._errors.append("RuntimeError")
                    raise RuntimeError(message)

        #}}}

        #{{{Check nx, ny, nz, zperiod, nout, mxstep, separatrix are int/iterable
        check_if_int = [\
            (self._nx       , 'nx')       ,\
            (self._ny       , 'ny')       ,\
            (self._nz       , 'nz')       ,\
            (self._zperiod  , 'zperiod')  ,\
            (self._nout     , 'nout')     ,\
            (self._mxstep   , 'mxstep')   ,\
            (self._ixseps1  , 'ixseps1')  ,\
            (self._ixseps2  , 'ixseps2')  ,\
            (self._jyseps1_1, 'jyseps1_1'),\
            (self._jyseps1_2, 'jyseps1_2'),\
            (self._jyseps2_1, 'jyseps2_1'),\
            (self._jyseps2_2, 'jyseps2_2'),\
            ]

        self._check_for_correct_type(var = check_if_int,\
                                      the_type = int,\
                                      allow_iterable = True)
        #}}}

        #{{{Check timestep, atol, rtol, zmin/max, dx, dy, dz is Number/iterable
        # Check if the following is a number
        check_if_number = [\
            (self._timestep, 'timestep'),\
            (self._zmin    , 'zmin')    ,\
            (self._zmax    , 'zmax')    ,\
            (self._dx      , 'dx')      ,\
            (self._dy      , 'dy')      ,\
            (self._dz      , 'dz')      ,\
            (self._atol    , 'atol')    ,\
            (self._rtol    , 'rtol')     \
            ]

        self._check_for_correct_type(var = check_if_number,\
                                    the_type = Number,\
                                    allow_iterable = True)
        #}}}

        #{{{Check if solver, grid_file, methods and sort_by is str/list of str
        # Check if instance is string, or an iterable containing strings
        check_if_string = [\
            (self._solver    , 'solver')    ,\
            (self._grid_file , 'grid_file') ,\
            (self._ddx_first , 'ddx_first') ,\
            (self._ddx_second, 'ddx_second'),\
            (self._ddx_upwind, 'ddx_upwind'),\
            (self._ddx_flux  , 'ddx_flux')  ,\
            (self._ddy_first , 'ddy_first') ,\
            (self._ddy_second, 'ddy_second'),\
            (self._ddy_upwind, 'ddy_upwind'),\
            (self._ddy_flux  , 'ddy_flux')  ,\
            (self._ddz_first , 'ddz_first') ,\
            (self._ddz_second, 'ddz_second'),\
            (self._ddz_upwind, 'ddz_upwind'),\
            (self._ddz_flux  , 'ddz_flux')  ,\
            (self._sort_by   , 'sort_by')    \
            ]

        self._check_for_correct_type(var = check_if_string,\
                                    the_type = str,\
                                    allow_iterable = True)
        #}}}

        #{{{Check if solver is set to the correct possibility
        # Check if the solver is possible
        # From /include/bout/solver.hxx
        possible_solvers = [\
            'cvode',\
            'pvode',\
            'ida',\
            'petsc',\
            'karniadakis',\
            'rk4',\
            'euler',\
            'rk3ssp',\
            'power',\
            'arkode'\
            ]

        # Do the check if the solver is set
        if self._solver is not None:
            self._check_if_set_correctly(var = (self._solver, 'solver'),\
                                          possibilities = possible_solvers)
        #}}}

        #{{{Check if the methods is set to the correct possibility
        # Check if ddx or ddy is possible
        possible_method = [\
            'C2',\
            'C4',\
            'W2',\
            'W3'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self._ddx_first , 'ddx_first') ,\
            (self._ddx_second, 'ddx_second'),\
            (self._ddy_first , 'ddy_first') ,\
            (self._ddy_second, 'ddy_second')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] is not None:
                self._check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)

        # Check if ddz is possible
        possible_method.append('FFT')

        # Make a list of the variables
        the_vars = [\
            (self._ddz_first , 'ddz_first') ,\
            (self._ddz_second, 'ddz_second') \
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] is not None:
                self._check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)

        # Check for upwind terms
        possible_method = [\
            'U1',\
            'U4',\
            'W3'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self._ddx_upwind, 'ddx_upwind'),\
            (self._ddy_upwind, 'ddy_upwind'),\
            (self._ddz_upwind, 'ddz_upwind')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] is not None:
                self._check_if_set_correctly(var          = var,\
                                              possibilities = possible_method)

        # Check for flux terms
        possible_method = [\
            'SPLIT',\
            'NND'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self._ddx_flux  , 'ddx_flux'),\
            (self._ddy_flux  , 'ddy_flux'),\
            (self._ddz_flux  , 'ddz_flux')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] is not None:
                self._check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)
        #}}}

        #{{{Check if sort_by is set to the correct possibility
        # Appendable list
        possible_sort_by = []

        # Append the 1st element of sort_checks if the 0th elements of
        # sort_checks != None
        sort_checks = [\
            (self._nx,         "spatial_domain")  ,\
            (self._ny,         "spatial_domain")  ,\
            (self._nz,         "spatial_domain")  ,\
            (self._dx,         "spatial_domain")  ,\
            (self._dy,         "spatial_domain")  ,\
            (self._dz,         "spatial_domain")  ,\
            (self._ixseps1,    "spatial_domain")  ,\
            (self._ixseps2,    "spatial_domain")  ,\
            (self._jyseps1_1,  "spatial_domain")  ,\
            (self._jyseps1_2,  "spatial_domain")  ,\
            (self._jyseps2_1,  "spatial_domain")  ,\
            (self._jyseps2_2,  "spatial_domain")  ,\
            (self._symGlobX,   "spatial_domain")  ,\
            (self._symGlobY,   "spatial_domain")  ,\
            (self._timestep,   "temporal_domain") ,\
            (self._nout,       "temporal_domain") ,\
            (self._solver,     "solver")          ,\
            (self._mms,        "solver")          ,\
            (self._atol,       "solver")          ,\
            (self._rtol,       "solver")          ,\
            (self._mxstep,     "solver")          ,\
            (self._ddx_first,  "ddx_first")       ,\
            (self._ddx_second, "ddx_second")      ,\
            (self._ddx_upwind, "ddx_upwind")      ,\
            (self._ddx_flux,   "ddx_flux")        ,\
            (self._ddy_first,  "ddy_first")       ,\
            (self._ddy_second, "ddy_second")      ,\
            (self._ddy_upwind, "ddy_upwind")      ,\
            (self._ddy_flux,   "ddy_flux")        ,\
            (self._ddz_first,  "ddz_first")       ,\
            (self._ddz_second, "ddz_second")      ,\
            (self._ddz_upwind, "ddz_upwind")      ,\
            (self._ddz_flux,   "ddz_flux")        ,\
            (self._grid_file,  "grid_file")        \
            ]

        for sort_check in sort_checks:
            if sort_check[0] is not None:
                if not(sort_check[1] in possible_sort_by):
                    possible_sort_by.append(sort_check[1])

        # Append the additionals and series_add
        # If additional is set
        if self._additional is not None:
            for additional in self._additional:
                # The additional now contains a tuple of three elements
                # We would like to extract the section (if any) and variable
                # and append them to the possibilities list
                # If the section is empty
                if additional[0] == '':
                    section = ''
                else:
                    section = additional[0] + ':'
                possible_sort_by.append(section + additional[1])
        # Do the same for series_add
        if self._series_add is not None:
            for series_add in self._series_add:
                if series_add[0] == '':
                    section = ''
                else:
                    section = series_add[0] + ':'
                possible_sort_by.append(section + series_add[1])

        # Make a list of the variables
        the_vars = [\
            (self._sort_by, 'sort_by')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] is not None:
                self._check_if_set_correctly(var           = var,\
                                              possibilities = possible_sort_by)
        #}}}

        #{{{Check if restart is set correctly
        if self._restart is not None:
            if type(self._restart) != str:
                self._errors.append("TypeError")
                raise TypeError ("restart must be set as a string when set")

        possible_method = [\
            'overwrite',\
            'append'\
            ]

        # Make a list of the variables
        the_vars = [\
            (self._restart, 'restart')\
            ]

        for var in the_vars:
            # Do the check if the method is set
            if var[0] is not None:
                self._check_if_set_correctly(var           = var,\
                                              possibilities = possible_method)
        #}}}

        #{{{Check if restart_from is set correctly
        if self._restart_from is not None:
            # Throw warning if restart is None
            if self._restart is None:
                message = "restart_from will be ignored as restart = None"
                self._warning_printer(message)
                self._warnings.append(message)

            if type(self._restart_from) != str:
                self._errors.append("TypeError")
                raise TypeError ("restart_from must be set as a string when set")

            # Check if any restart files are present
            if len(glob.glob(os.path.join(self._restart_from,'*restart*')))== 0:
                self._errors.append("FileNotFoundError")
                raise FileNotFoundError("No restart files found in " +\
                                 self._restart_from)
        #}}}

        #{{{Check for options set in both member data and in the grid file
        if self._grid_file is not None:
            # Check if the following variables are found in the grid
            # file
            check_if_in_grid =[\
                    (self._nx       , "nx")               ,\
                    (self._ny       , "ny")               ,\
                    (self._nz       , "nz")               ,\
                    (self._dx       , "dx")               ,\
                    (self._dy       , "dy")               ,\
                    (self._dz       , "dz")               ,\
                    (self._MXG      , "MXG")              ,\
                    (self._MYG      , "MYG")              ,\
                    (self._NXPE     , "NXPE")              ,\
                    (self._NYPE     , "NYPE")              ,\
                    (self._ixseps1  , "ixseps1")          ,\
                    (self._ixseps2  , "ixseps2")          ,\
                    (self._jyseps1_1, "jyseps1_1")        ,\
                    (self._jyseps1_2, "jyseps1_2")        ,\
                    (self._jyseps2_1, "jyseps2_1")        ,\
                    (self._jyseps2_2, "jyseps2_2")        ,\
                    (self._symGlobX , "symmmetricGlobalX"),\
                    (self._symGlobY , "symmmetricGlobalY") \
                    ]
            for var in check_if_in_grid:
                # If the variable is set
                if var[0] is not None:
                    # Loop through the grid files
                    for grid_file in self._grid_file:
                        # Open (and automatically close) the grid files
                        f = DataFile(grid_file)
                        # Search for mesh data in the grid file
                        grid_variable = f.read(var[1])
                        # If the variable is found
                        if grid_variable is not None:
                            self._errors.append("TypeError")
                            message  = var[1] + " was specified both in the "
                            message += "driver and in the grid file.\n"
                            message += "Please remove " + var[1]
                            message += " from the driver if you would "
                            message += "like to run with a grid file."
                            raise TypeError(message)
        #}}}

        #{{{If grid files are set: Use nx, ny and nz values in the grid file
        if self._grid_file is not None:
            # Make a dict of appendable lists
            spatial_domain = {'nx':[], 'ny':[], 'nz':[]}
            for grid_file in self._grid_file:
                # Open (and automatically close) the grid files
                f = DataFile(grid_file)
                # Search for nx, ny and nz in the grid file
                mesh_types = ["nx", "ny", "nz"]
                for mesh_type in mesh_types:
                    grid_variable = f.read(mesh_type)
                    # If the variable is found
                    if grid_variable is not None:
                        spatial_domain[mesh_type].append(grid_variable)
            # Check that the lengths of nx, ny and nz are the same
            # unless they are not found
            len_nx = len(spatial_domain['nx'])
            len_ny = len(spatial_domain['ny'])
            len_nz = len(spatial_domain['nz'])
            if len_nx != 0:
                self._nx = spatial_domain['nx']
            if len_ny != 0:
                self._ny = spatial_domain['ny']
            if len_nz != 0:
                self._nz = spatial_domain['nz']
        #}}}

        #{{{Check that nx, ny and nz are of the same length
        if self._nx is not None and self._ny is not None:
            self._check_if_same_len((self._nx, 'nx'), (self._ny, 'ny'))
        if self._nx is not None and self._nz is not None:
            self._check_if_same_len((self._nx, 'nx'), (self._nz, 'nz'))
        if self._ny is not None and self._nz is not None:
            self._check_if_same_len((self._ny, 'ny'), (self._nz, 'nz'))
        #}}}

        #{{{Check that NXPE and NYPE are of the same length as nx, ny, nz
        if self._nx is not None and self._NXPE is not None:
            self._check_if_same_len((self._nx, 'nx'), (self._NXPE, 'NXPE'))
        if self._ny is not None and self._NXPE is not None:
            self._check_if_same_len((self._ny, 'ny'), (self._NXPE, 'NXPE'))
        if self._nz is not None and self._NXPE is not None:
            self._check_if_same_len((self._nz, 'nz'), (self._NXPE, 'NXPE'))

        if self._nx is not None and self._NYPE is not None:
            self._check_if_same_len((self._nx, 'nx'), (self._NYPE, 'NYPE'))
        if self._ny is not None and self._NYPE is not None:
            self._check_if_same_len((self._ny, 'ny'), (self._NYPE, 'NYPE'))
        if self._nz is not None and self._NYPE is not None:
            self._check_if_same_len((self._nz, 'nz'), (self._NYPE, 'NYPE'))
        #}}}

        #{{{Check (zperiod), (zmin, zmax) and (dz) is not set simultaneously
        if (self._zperiod is not None and\
           (self._zmin != None or self._zmax != None)):
            self._errors.append("TypeError")
            message = "zperiod and zmin or zmax cannot be set simultaneously."
            raise TypeError(message)
        elif (self._dz is not None and\
             (self._zmin != None or self._zmax != None)):
            self._errors.append("TypeError")
            message = "dz and zmin or zmax cannot be set simultaneously."
            raise TypeError(message)
        elif (self._zperiod is not None and self._dz):
            self._errors.append("TypeError")
            message = "dz and zperiod cannot be set simultaneously."
            raise TypeError(message)
        #}}}

        #{{{Check that dz is not set
        # dz is currently set throught zmin and zmax
        if self._dz is not None:
            self._errors.append("TypeError")
            message  = "dz can currently just be set through zmin and zmax\n"
            message += "dz = 2*pi*(zmax-zmin)/(MZ-1)"
            raise TypeError (message)
        #}}}

        #{{{Check that dx, dy and dz are of the same length
        if self._dx is not None and self._dy is not None:
            self._check_if_same_len((self._dx, 'dx'), (self._dy, 'dy'))
        if self._dx is not None and self._dz is not None:
            self._check_if_same_len((self._dx, 'dx'), (self._dz, 'dz'))
        if self._dy is not None and self._dz is not None:
            self._check_if_same_len((self._dy, 'dy'), (self._dz, 'dz'))
        #}}}

        #{{{Check that (dx, nx), (dy, ny) and (dz,nz) are of the same length
        if self._dx is not None and self._nx is not None:
            self._check_if_same_len((self._dx, 'dx'), (self._nx, 'nx'))
        if self._dy is not None and self._ny is not None:
            self._check_if_same_len((self._dy, 'dy'), (self._ny, 'ny'))
        if self._nz is not None and self._dz is not None:
            self._check_if_same_len((self._dz, 'dz'), (self._nz, 'nz'))
        #}}}

        #{{{ Check that timestep and nout have the same len
        if self._timestep is not None and self._nout is not None:
            self._check_if_same_len((self._timestep, 'timestep'),\
                                   (self._nout, 'nout'))
        #}}}

        #{{{Check that additional and series_add are on the correct form
        self._error_check_additional((self._additional, 'additional'))
        self._error_check_additional((self._series_add, 'series_add'))
        #}}}

        #{{{Check that self._series_add[:][2] have the same length
        if self._series_add is not None:
            # Make the second indices iterable if they are not already
            for index in range(len(self._series_add)):
                if not(hasattr(self._series_add[index][2], "__iter__")) or\
                   (type(self._series_add[index][2]) == str) or\
                   (type(self._series_add[index][2]) == dict):
                    # Check if the type is a tuple
                    if (type(self._series_add[index]) != tuple):
                        self._series_add[index][2]=[self._series_add[index][2]]
                    else:
                        # We are dealing with tuples
                        # Cast to list
                        self._series_add[index] = list(self._series_add[index])
                        self._series_add[index][2]=[self._series_add[index][2]]
                        # Recast to tuple
                        self._series_add[index] = tuple(self._series_add[index])

            # Collect all second indices
            second_indices = [elems[2] for elems in self._series_add]
            # Find the length of the second indices
            lengths = [len(elem) for elem in second_indices\
                       if (type(elem)!=str and type(elem)!=dict)]
            # Check if any string or dicts were given
            if len(second_indices) != len(lengths):
                message  = "series_add is on the wrong form.\n"
                message += "series_add should be on the form\n"
                message += "series_add=\ \n"
                message +=\
                        "     [(section1, name1, [value1-1, value1-2,... ,"+\
                        "value1-n]),\ \n"
                message +=\
                        "      (section2, name2, [value2-1, value2-2,... ,"+\
                        "value2-n]),\ \n"
                message +=\
                        "       ...])\n"
                self._errors.append("TypeError")
                raise TypeError(message)

            # Check that the length of the second indices are the same
            # L.count(value) -> integer -- return number of occurrences
            # of value
            # stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical
            if not(lengths.count(lengths[0]) == len(lengths)):
                message = "The length of the second index of the elements"+\
                          " of series_add must be the same"
                self._errors.append("TypeError")
                raise TypeError(message)
        #}}}

        #{{{Check mms, symGlobX, symGlobY, cpy_src/grid, allow_size_mod is bool
        check_if_bool = [\
            (self._mms                    , 'mms')                    ,\
            (self._symGlobX               , 'symGlobX')               ,\
            (self._symGlobY               , 'symGlobY')               ,\
            (self._cpy_source             , 'cpy_source')             ,\
            (self._cpy_grid               , 'cpy_grid')               ,\
            (self._allow_size_modification, 'allow_size_modification') \
            ]

        self._check_for_correct_type(var = check_if_bool,\
                                    the_type = bool)
        #}}}

        #{{{Check grid_file == None if cpy_grid==True
        if (self._grid_file is None) and (self._cpy_grid == True):
            # Raise error
            self._errors.append("TypeError")
            message = "Cannot copy the grid files if none exists in "+\
                      " 'grid_file'"
            raise TypeError(message)
        #}}}

        #{{{Check that zmin and zmax has the same length
        if (self._zmin is not None) and (self._zmax is not None):
            self._check_if_same_len((self._zmin, 'zmin'),\
                                     (self._zmax, 'zmax'))

        #}}}
#}}}
#}}}

#{{{Functions called by _check_for_basic_instance_error
    #{{{_error_check_additional
    def _error_check_additional(self, input_member):
        #{{{docstring
        """
        Checks that the input_member is on the following form:

        input_member = [(section1, name1, [value1-1, value1-2, ...]),
                        (section2, name2, [value2-1, value2-2, ...]),
                        ...]

        Input:
        Either self._additional or self._series_add
        input_member[0]     -   the input data
        input_member[1]     -   the name of the input data
        """
        #}}}

        # If input_member is set
        if input_member[0] is not None:
            # Set a success variable that will fail if anything goes
            # wrong
            success = True
            # If we need to change the elements, make sure that the
            # input element is not a tuple
            if type(input_member[0]) == tuple:
                success = False
            # Loop through all elements in input_member
            for elem in input_member[0]:
                # Check if self._addition is iterable, but not a string
                # or dict
                if (hasattr(elem, "__iter__")) and\
                   (type(elem) != str) and\
                   (type(elem) != dict):
                    if type(elem[0]) == str:
                        # Check that the second element (the name) is a
                        # string
                        if type(elem[1]) != str:
                            success = False
                        # If more than three elements are given
                        if len(elem) != 3:
                            success = False
                    # elem[0] is not a string
                    else:
                        success = False
                # elem is not iterable or is a dict or a string
                else:
                    success = False
            if not(success):
                message  = input_member[1]+" is on the wrong form.\n"
                message += input_member[1]+" should be on the form\n"
                message += input_member[1]+"=\ \n"
                message +=\
                        "     [(section1, name1, [value1-1, value1-2,...]),\ \n"
                message +=\
                        "      (section2, name2, [value2-1, value2-2,...]),\ \n"
                message +=\
                        "       ...])\n"
                self._errors.append("TypeError")
                raise TypeError(message)
        #}}}
#}}}

#{{{ Functions called by the execute_runs function
#{{{_error_check_for_run_input
    def _error_check_for_run_input(self                        ,\
                                   remove_old                  ,\
                                   post_processing_function    ,\
                                   post_process_after_every_run \
                                   ):
        """Check if there are any type errors in input for the run
        function"""

        #{{{Check if remove_old is of the correct type
        check_if_bool = [\
            (remove_old , 'remove_old'),\
            ]

        self._check_for_correct_type(var = check_if_bool,\
                                    the_type = bool)
        #}}}

        #{{{Check if remove_old and restart is set on the same time
        if remove_old == True and self._restart is not None:
            self._errors.append("RuntimeError")
            raise RuntimeError("You should not remove old data if you"\
                               " want a restart run")
        #}}}

        #{{{Check that the post_processing_function is a fuction
        if (post_processing_function is not None) and\
           (not(hasattr(post_processing_function, '__call__'))):
            self._errors.append("RuntimeError")
            message = "post_process_after_every_run must be a"+\
                      " function"
            raise RuntimeError(message)
        #}}}

        #{{{Check that the post_process_after_every_run is not set alone
        if (post_process_after_every_run is not None) and\
           (type(post_processing_function) == None):
            self._errors.append("RuntimeError")
            message = "post_process_after_every_run can only be set if"+\
                      " post_processing_function is given"
            raise RuntimeError(message)
        #}}}

        #{{{Check that the post_process_after_every_run is a boolean
        if (post_process_after_every_run is not None) and\
           (type(post_process_after_every_run) != bool):
            self._errors.append("RuntimeError")
            message = "post_process_after_every_run must be set to"+\
                      " a boolean when set"
            raise RuntimeError(message)
        #}}}

        # Check for errors in a child class
        self._check_for_child_class_errors(
                                   remove_old                  ,\
                                   post_processing_function    ,\
                                   post_process_after_every_run \
                                          )
#}}}

#{{{Functions called by _error_check_for_run_input
    #{{{_check_for_child_class_errors
    def _check_for_child_class_errors(
                                   self                        ,\
                                   remove_old                  ,\
                                   post_processing_function    ,\
                                   post_process_after_every_run \
                                   ):
        """Function which check for errors in a child class.

        Here a virtual function"""
        pass
    #}}}
#}}}

#{{{_create_run_log
    def _create_run_log(self):
        """Makes a run_log file if it doesn't exists"""

        # Checks if run_log exists
        self._run_log = self._directory + "/run_log.txt"
        if os.path.isfile(self._run_log) == False:
            # The header
            header = ['start_time', 'run_type', 'run_no', 'run_time_H:M:S', 'dump_folder']
            header_format = '{:<19}   {:<9}   {:<6}   {:<17}   {:<}'
            # Create the log file, and print the header
            with open(self._run_log , "w") as f:
                f.write(header_format.format(*header) + '\n')

        # Preparation of the run
        print("\nRunning with inputs from '" + self._directory + "'")
#}}}

#{{{_get_correct_domain_split
    def _get_correct_domain_split(self):
        """Checks that the grid can be split in the correct number of
        processors. If not, vary the number of points until value is found."""

        #{{{First we check if self._MXG is given
        if self._MXG is None:
            # We need to find MXG in BOUT.inp. We use BOUTOption for
            # this
            # Object initialization
            myOpts = BOUTOptions(self._directory)
            # Get MXG as string
            local_MXG = myOpts.root['MXG']
            # Evaluate the string
            local_MXG = eval(local_MXG)
        else:
            local_MXG = self._MXG
        #}}}

        # If NXPE is not set, we will try to find a optimal grid size
        # Flag to determine if a warning should be printed
        produce_warning = False
        print("\nChecking the grid split for the meshes\n")
        if self._NXPE is None:
            #{{{ If NXPE is not set
            for size_nr in range(len(self._nx)):
                print("Checking nx=" + str(self._nx[size_nr]) +\
                      " and ny=" + str(self._ny[size_nr]))
                # Check to see if succeeded
                init_split_found = False
                cur_split_found  = False
                add_number = 1
                # Counter to see how many times the while loop has been
                # called
                count = 0

                #{{{While cur_split_found == False
                while cur_split_found == False:
                    # The same check as below is performed internally in
                    # BOUT++ (see boutmesh.cxx under
                    # if(options->isSet("NXPE")))
                    for i in range(1, self._nproc+1, 1):
                        MX = self._nx[size_nr] - 2*local_MXG
                        # self._nproc is called NPES in boutmesh
                        if (self._nproc % i == 0) and \
                           (MX % i == 0) and \
                           (self._ny[size_nr] % (self._nproc/i) == 0):
                            # If the test passes
                            cur_split_found = True

                    # Check if cur_split_found is true, eventually
                    # update the add_number
                    add_number, produce_warning = self._check_cur_split_found(\
                                                             cur_split_found,\
                                                             produce_warning,\
                                                             add_number,\
                                                             size_nr,\
                                                             using_nx = True,\
                                                             using_ny = True)


                    #{{{ Check if the split was found the first go.
                    # This will be used if self_allow_size_modification is
                    # off, or if we are using a grid file
                    if count == 0 and cur_split_found:
                        init_split_found = True
                    #}}}

                    # Add one to the counter
                    count += 1
                #}}}

                # Check if initial split succeeded
                self._check_init_split_found(init_split_found, size_nr,\
                                             test_nx = True, test_ny = True,
                                             produce_warning = produce_warning)
            #}}}
        else:
            #{{{ If NXPE is set
            # Check if NXPE and NYPE is set consistently with nproc
            self._check_NXPE_or_NYPE(type_txt = 'NXPE', local_MXG = local_MXG)
            self._check_NXPE_or_NYPE(type_txt = 'NYPE')
            #}}}
#}}}

#{{{_get_possibilities
    def _get_possibilities(self):
        """ Returns the list of the possibilities. In get_combinations
        the elements of this list is going to be put together to a list
        of strings which will be used when making a run."""

        #{{{Set combination of nx, ny and nz (if not set in grid_file)
        # Appendable list
        spatial_grid_possibilities = []
        if (self._grid_file is None):
            # Dictionary where
            # - the first element is the variable itself
            # - the second element is the section of the variable
            # - the third element is an appendable list
            spatial_grid_str = {\
                                'nx'     :[self._nx,      'mesh:', []],\
                                'ny'     :[self._ny,      'mesh:', []],\
                                'nz'     :[self._nz,      'mesh:', []],\
                                'dx'     :[self._dx,      'mesh:', []],\
                                'dy'     :[self._dy,      'mesh:', []],\
                                'dz'     :[self._dz,      'mesh:', []],\
                                'zperiod':[self._zperiod, '',      []],\
                                'zmin'   :[self._zmin,    '',      []],\
                                'zmax'   :[self._zmax,    '',      []],\
                               }
            # Store the keys as an own variable
            keys = list(spatial_grid_str.keys())
            # Append the different dimension to the list of strings
            for key in keys:
                # If the variable is not empty
                if spatial_grid_str[key][0] is not None:
                    # Fill the appendable list with the elements from
                    # the variable
                    for elem in spatial_grid_str[key][0]:
                        spatial_grid_str[key][2].append(\
                            spatial_grid_str[key][1] + key + '=' + str(elem)\
                                )

            # The goal is to combine the these strings to one string
            # Find the largest length
            lengths = [len(spatial_grid_str[key][2]) for key in keys]
            max_len = np.max(lengths)
            # Make the strings the same length
            for key in keys:
                # We do this by filling it with empty strings
                while len(spatial_grid_str[key][2]) <= max_len:
                    spatial_grid_str[key][2].append('')

            # Append this to the spatial grid possibilities as a string
            for number in range(max_len):
                # Make a list
                current_grid = [spatial_grid_str[key][2][number] for key in\
                                keys]
                # Join the strings in the list and append
                spatial_grid_possibilities.append(' '.join(current_grid))
        #}}}

        #{{{Set the combination of timestep and nout if is not None
        # Appendable lists
        temporal_grid_possibilities = []
        timestep_str = []
        nout_str     = []
        # Append the different time options to the list of strings
        if self._timestep is not None:
            for timestep in self._timestep:
                timestep_str.append('timestep=' + str(timestep))
        if self._nout is not None:
            for nout in self._nout:
                nout_str.append('nout=' + str(nout))
        # Combine the strings to one string
        # Find the largest length
        max_len = np.max([len(timestep_str), len(nout_str)])
        # Make the strings the same length
        if len(timestep_str) < max_len:
            timestep_str.append('')
        if len(nout_str) < max_len:
            nout_str.append('')
        # Append the temporal grid possibilities as a string
        for number in range(max_len):
            # Make a list
            current_times = [timestep_str[number],\
                             nout_str[number]\
                            ]
            # Join the strings in the list and append
            temporal_grid_possibilities.append(' '.join(current_times))
        #}}}

        #{{{Set the combination of the series_add option if is not None
        # Appendable list
        series_add_possibilities = []
        if self._series_add is not None:
            # Dictionary to handle the data, where the key is going to
            # be the element number in self._series_add, and the values
            # are going to be the sub dictionary defined below
            all_info = {}
            # Loop through all elements and fill the dictionary
            for nr, elem in enumerate(self._series_add):
                # Put in the sub dictionary
                all_info[nr] = {'values':None,\
                                'section_and_var':None,\
                                'sec_var_vals':[]}
                # Fill the values
                all_info[nr]['values'] = elem[2]
                # Fill the section and variable key
                all_info[nr]['section_and_var'] = elem[0] + ":" + elem[1] + "="
                # Fill in the combinations
                for val in all_info[nr]['values']:
                    all_info[nr]['sec_var_vals'].append(\
                        all_info[nr]['section_and_var'] + str(val)\
                            )

            # Make an appendable list
            all_sec_var_vals = []
            for key in all_info.keys():
                all_sec_var_vals.append(all_info[key]['sec_var_vals'])

            # Zip the sec_var_vals together (* unpacks), join them with
            # a space, and append them to series_add_possibilities
            for one_possibility in zip(*all_sec_var_vals):
                series_add_possibilities.append(' '.join(one_possibility))
        #}}}

        #{{{Put non-iterable variables into a list if they are not set to None
        # This makes the member data iterable, and usable in
        # generate_possibilities
        if self._MXG is not None:
            self._MXG = [self._MXG]
        if self._MYG is not None:
            self._MYG = [self._MYG]
        if self._mms is not None:
            self._mms = [self._mms]
        if self._symGlobX is not None:
            self._symGlobX = [self._symGlobX]
        if self._symGlobY is not None:
            self._symGlobY = [self._symGlobY]
        #}}}

        #{{{List of tuple of variables to generate possibilities from
        tuple_of_variables = [\
            (self._solver,     "solver", "type")             ,\
            (self._mms,        "solver", "mms")              ,\
            (self._atol,       "solver", "atol")             ,\
            (self._rtol,       "solver", "rtol")             ,\
            (self._mxstep,     "solver", "mxstep")           ,\
            (self._MXG,        "",       "MXG")              ,\
            (self._MYG,        "",       "MYG")              ,\
            (self._NXPE,       "",       "NXPE")             ,\
            (self._NYPE,       "",       "NYPE")             ,\
            (self._grid_file,  "",       "grid")             ,\
            (self._ddx_first,  "ddx",    "first")            ,\
            (self._ddx_second, "ddx",    "second")           ,\
            (self._ddx_upwind, "ddx",    "upwind")           ,\
            (self._ddx_flux,   "ddx",    "flux")             ,\
            (self._ddy_first,  "ddy",    "first")            ,\
            (self._ddy_second, "ddy",    "second")           ,\
            (self._ddy_upwind, "ddy",    "upwind")           ,\
            (self._ddy_flux,   "ddy",    "flux")             ,\
            (self._ddz_first,  "ddz",    "first")            ,\
            (self._ddz_second, "ddz",    "second")           ,\
            (self._ddz_upwind, "ddz",    "upwind")           ,\
            (self._ddz_flux,   "ddz",    "flux")             ,\
            (self._ixseps1,    "mesh",   "ixseps1")          ,\
            (self._ixseps2,    "mesh",   "ixseps2")          ,\
            (self._jyseps1_1,  "mesh",   "jyseps1_1")        ,\
            (self._jyseps1_2,  "mesh",   "jyseps1_2")        ,\
            (self._jyseps2_1,  "mesh",   "jyseps2_1")        ,\
            (self._jyseps2_2,  "mesh",   "jyseps2_2")        ,\
            (self._symGlobX,   "mesh",   "symmetricGlobalX") ,\
            (self._symGlobY,   "mesh",   "symmetricGlobalY")  \
            ]
        #}}}

        #{{{Append the additional option to tuple of variables if set
        if self._additional is not None:
            for additional in self._additional:
                # If the last element of additional is not iterable we need
                # put them into a list to make them iterable (in order to
                # use them in generate_possibilities)
                if (not(hasattr(additional[2], "__iter__"))) or\
                   (type(additional[2]) == str):
                    # We have to specify the whole additional, as this can
                    # be given as a tuple, and tuples does not support item
                    # assignments
                    additional = (additional[0],\
                                  additional[1],\
                                  [additional[2]])
                # Append the additional to tuple of variables
                tuple_of_variables.append(\
                                (additional[2],\
                                additional[0],\
                                additional[1])\
                                )
        #}}}

        #{{{List of the possibilities of the variables
        # Start out with the already generated
        # spatial_grid_possibilities and temporal_grid_possibilities
        list_of_possibilities = [spatial_grid_possibilities,\
                                 temporal_grid_possibilities,\
                                 series_add_possibilities]

        # Append the possibilities to the list of possibilities
        for var in tuple_of_variables:
            list_of_possibilities.append(\
                    self._generate_possibilities(var[0], var[1], var[2])\
                    )
        #}}}

        # Return the list_of possibilities
        return list_of_possibilities
#}}}

#{{{_get_combinations
    def _get_combinations(self, input_list):
        """The input_list is a list with lists as element.
        Returns a list of all combinations between the elements of the
        input_list."""

        # Remove empty elements in input_list in order for
        # itertools.product to work
        input_list = [elem for elem in input_list if elem != []]

        # If we would like to sort the input list (choose which variable
        # to be the fastest varying)
        if self._sort_by is not None:
            # Swap the list corresponding to the sort_by statement so
            # that that list will be the last. The itertools.product
            # will then make that list the fastest varying in the list
            input_list = self._get_swapped_input_list(input_list)
        else:
            # Initialize this member data to None
            self._len_group = None

        # The last element in the input_list will be the fastest varying
        # element
        all_combinations_as_tuple = list(itertools.product(*input_list))

        # all_combination_as_tuple is a list with tuples as elements
        # We would like to combine the elements in these tuples to one
        # string
        # Make an appendable list
        all_combinations_as_strings = []

        # Loop over the elements in the list containing tuples
        for a_tuple in all_combinations_as_tuple:
            # Join the elements in a tuple and store it
            all_combinations_as_strings.append(' '.join(a_tuple))

        return all_combinations_as_strings
#}}}

#{{{_print_run_or_submit
    def _print_run_or_submit(self):
        """Prints 'Now running'"""
        print("\nNow running:")
#}}}

#{{{_prepare_dmp_folder
    def _prepare_dmp_folder(self, combination):
        """
        Set the folder to dump data in based on the input from the
        combination.

        - Copy the input file to the final folder.
        - Copy restart files if restart_from is set (can set skip_run=True)
        - Copy files if restart is set to overwrite
        - Copy the source files to the final folder is cpy_source is True.

        Returns skip_run = True if there are any troubles with the copying
        """
        # Obtain folder names
        folder_name = self._get_folder_name(combination)
        self._dmp_folder = os.path.join(self._directory, folder_name)
        # If the last character is '/', then remove it
        if self._dmp_folder[-1] == '/':
            self._dmp_folder = self._dmp_folder[:-1]

        # Create folder if it doesn't exists
        self._create_folder(self._dmp_folder)
        # If self._dmp_folder contains anything other than
        # self._directory
        if self._dmp_folder != self._directory:
            # Copy the input file into this folder
            src = os.path.join(self._directory, 'BOUT.inp')
            shutil.copy2(src, self._dmp_folder)

        # Copy restart files if restart_from is set
        # skip_run is set to False by default
        skip_run = False
        if self._restart and self._restart_from:
            # Copy the files to restart
            skip_run = self._copy_restart_files()

        # Save files if restart is set to "overwrite"
        if self._restart == 'overwrite':
            self._move_old_runs()

        # Copy the source files if cpy_source is True
        if self._cpy_source:
            # This will copy all C++ files to the dmp_folder
            cpp_extension= ['.cc', '.cpp', '.cxx', '.C', '.c++',\
                            '.h',  '.hpp', '.hxx', '.h++']
            # Copy for all files in the extension
            for extension in cpp_extension:
                file_names = glob.glob('*' + extension)
                for a_file in file_names:
                    shutil.copy2(a_file, self._dmp_folder)
        return skip_run
#}}}

#{{{_remove_data
    def _remove_data(self):
        """Removes *.nc and *.log files from the dump directory"""

        print("Removing old data")
        remove_extensions = ['dmp.*', 'fail.*', 'restart.*', 'log.*', 'cpy']
        files_to_rm = []
        for extension in remove_extensions:
            files_to_rm.extend(\
                    glob.glob(os.path.join(self._dmp_folder, "*." + extension)))

        # Cast to set (unique values)
        files_to_rm = set(files_to_rm)
        for f in files_to_rm:
            os.remove(f)

        # Remove dirs
        folder_to_rm = glob.glob(os.path.join(self._dmp_folder, "run*"))
        # Filter to only inlcude folders
        folder_to_rm = [f for f in folder_to_rm if os.path.isdir(f)]
        for f in folder_to_rm:
            os.removedirs(f)
#}}}

#{{{_check_if_run_already_performed
    def _check_if_run_already_performed(self):
        """
        Checks if the run has been run previously.

        If restart is set, and no files are found, a warning will be
        printed.

        Returns
        True    - The run will be performed
        False   - The run will NOT be performed
        """

        dmp_files = glob.glob(os.path.join(self._dmp_folder, '*.dmp.*'))
        # If no BOUT.inp files are found or if self._restart is not set
        # (meaning that the run will be done even if files are found)
        if len(dmp_files) != 0 and self._restart is None:
            print('Skipping the run as *.dmp.* files was found in '\
                  + self._dmp_folder)
            print('To overwrite old files, run with'+\
                  ' self.execute_runs(remove_old=True)\n')
            return False
        # Either no files are found, or restart is set
        else:
            if len(dmp_files) == 0 and\
               self._restart is not None and\
               self._restart_from is None:
                message = "'restart' was set to " +self._restart+\
                          ", but no dmp files found."+\
                          " Setting 'restart' to None"
                self._restart = None
                self._warning_printer(message)
                self._warnings.append(message)
            return True
#}}}

#{{{_call_post_processing_function
    def _call_post_processing_function(\
                    self           ,\
                    function = None,\
                    folders  = None,\
                    **kwargs):
        """Function which calls the post_processing_function"""

        function(folders, **kwargs)

#}}}
#}}}

#{{{Function called by _set_program_name
#{{{_run_make
    def _run_make(self):
        """Make cleans and makes the .cxx program"""

        print("Make clean eventually previously compiled\n")
        command = "make clean"
        status, output = shell(command, pipe=True)
        print("Making the .cxx program\n")
        command = "make"
        status, output = shell(command, pipe=True)
        print(output)
        # Check if any errors occurred
        if status != 0:
            self._errors.append("RuntimeError")
            raise RuntimeError("Error encountered during make.")
#}}}
#}}}

#{{{ Functions called by the basic_error_checker
#{{{_check_for_correct_type
    def _check_for_correct_type(self,\
                                 var            = None,\
                                 the_type       = None,\
                                 allow_iterable = None):
        """Checks if a variable has the correct type

        Input:
        var            - a tuple consisting of
                         var[0] - the variable (a data member)
                         var[1] - the name of the variable given as a string
        the_type       - the data type
        allow_iterable - if an iterable with the element as type is
                         allowed
        """

        # Set a variable which is False if the test fails
        success = True
        for cur_var in var:
            # There is an option that the variable could be set to None,
            # and that the default value from BOUT.inp will be used
            if cur_var[0] is not None:
                # Check for the correct type
                if isinstance(cur_var[0], the_type) == False:
                    # Check if it is an iterable if iterables are
                    # allowed
                    if allow_iterable and\
                       hasattr(cur_var[0], "__iter__") and\
                       type(cur_var[0]) != dict:
                        for elem in cur_var[0]:
                            # Check for the correct type
                            if isinstance(elem, the_type) == False:
                                success = False
                    else:
                        # Neither correct type, nor iterable
                        success = False
                if not(success):
                    message  = cur_var[1] + " is of wrong type\n"+\
                               cur_var[1] + " must be " + the_type.__name__
                    if allow_iterable:
                        # If iterable is allowed, then add this
                        message += " or an iterable with " + the_type.__name__
                        message += " as elements."
                    self._errors.append("TypeError")
                    raise TypeError(message)
#}}}

#{{{_check_if_set_correctly
    def _check_if_set_correctly(self,\
                                 var           = None,\
                                 possibilities = None):
        """Check if a variable is set to a possible variable. Called by
        the error checkers"""

        # Set a variable which is False if the test fails
        success = True

        # Due to the check done in check_for_correct_type: If the
        # variable is not a string it will be an iterable
        if type(var[0]) != str:
            for elem in var[0]:
                # Check if the element is contained in the possibilities
                if not(elem in possibilities):
                    success = False
        else:
            # The variable was a string
            if not(var[0] in possibilities):
                success = False

        if not(success):
            message = var[1] + " was not set to a possible option.\n"+\
                      "The possibilities are \n" + "\n".join(possibilities)
            self._errors.append("TypeError")
            raise TypeError(message)
#}}}

#{{{_check_if_same_len
    def _check_if_same_len(self, object1 = None, object2 = None):
        """Checks if object1 and object2 has the same length

        Input:
        object1 - a tuple of the object [0] and its name [1]
        object2 - a tuple an object [0] different than object1 together with
                  its name [1]
        """

        try:
            len_dim1 = len(object1[0])
        # If object1 does not have length
        except TypeError:
            len_dim1 = 1
        try:
            len_dim2 = len(object2[0])
        # If object2 does not have length
        except TypeError:
            len_dim2 = 1

        if len_dim1 != len_dim2:
            message = object1[1] + " and " + object2[1] + " must have the same"
            message += " length when specified"
            self._errors.append("RuntimeError")
            raise RuntimeError (message)
#}}}
#}}}

#{{{ Functions called by _get_correct_domain_split
    #{{{_check_cur_split_found
    def _check_cur_split_found(self, cur_split_found,\
                               produce_warning,\
                               add_number, size_nr,\
                               using_nx = None, using_ny = None):
        #{{{docstring
        """
        Checks if the current split is found.

        Will add a number if not found.

        Input:
        cur_split_found     -   whether or not the current split was
                                found
        produce_warning     -   if a warning should be produced
        add_number          -   the number added to nx and/or ny
        size_nr             -   index of the current nx and/or ny
        using_nx            -   if add_number should be added to nx
        using_ny            -   if add_number should be added to ny

        Output:
        add_number          -   the number to eventually be added the
                                next time
        produce_warning     -   whether or not a warning should be
                                produced
        """
        #}}}

        # If the value tried is not a good value
        if cur_split_found == False:
            # Produce a warning
            produce_warning = True
            if using_nx:
                self._nx[size_nr] += add_number
            if using_ny:
                self._ny[size_nr] += add_number

            print("Mismatch, trying "+ str(self._nx[size_nr]) +\
                  "*" + str(self._ny[size_nr]))
            add_number = (-1)**(abs(add_number))\
                         *(abs(add_number) + 1)
        else:
            # If no warnings has been produced so far
            if not(produce_warning):
                produce_warning = False

        return add_number, produce_warning
    #}}}

    #{{{_check_init_split_found
    def _check_init_split_found(self, init_split_found, size_nr,\
                                test_nx = None, test_ny = None,\
                                produce_warning = None):
        #{{{docstring
        """
        Check if the initial split was a good choice when checking the grids.

        Will raise eventual errors.
        Input:
        init_split_found    -   boolean revealing whether or not a good
                                split was found on the first trial
        size_nr             -   the index of the current nx, ny or NXPE
                                under consideration
        test_nx             -   whether or not the test was run on nx
        test_ny             -   whether or not the test was run on ny
        produce_warning     -   whether or not a warning should be
                                produced
        """
        #}}}

        #{{{ If the initial split did not succeed
        if not(init_split_found):
            # If modification is allowed
            if not(self._allow_size_modification) or\
                  (self._grid_file != None):
                # If the split fails and the a grid file is given
                if self._grid_file is not None:
                    self._errors.append("RuntimeError")
                    message = "The grid can not be split using the"+\
                              " current number of nproc.\n"+\
                              "Suggest using "
                    if test_nx:
                        message += "nx = " + str(self._nx[size_nr]) + " "
                    if test_ny:
                        message += "ny = " + str(self._ny[size_nr]) + " "
                    message += " with the current nproc"
                    raise RuntimeError(message)
                # If the split fails and no grid file is given
                else:
                    self._errors.append("RuntimeError")
                    message  = "The grid can not be split using the"+\
                               " current number of nproc.\n"
                    message += "Setting allow_size_modification = True"+\
                               " will allow modification of the grid"+\
                               " so that it can be split with the"+\
                               " current number of nproc"
                    raise RuntimeError(message)
        #}}}

        #{{{ When the good value is found
        print("Successfully found the following good values for the mesh:")
        message = ''
        if test_nx:
            message += "nx = " + str(self._nx[size_nr]) + " "
        if test_ny:
            message += "ny = " + str(self._ny[size_nr])

        print(message + "\n")
        #}}}

        #{{{ Make the warning if produced
        if produce_warning:
            message = "The mesh was changed to allow the split given by nproc"
            self._warning_printer(message)
            self._warnings.append(message)
        #}}}
    #}}}

    #{{{_check_NXPE_or_NYPE
    def _check_NXPE_or_NYPE(self, type_txt  = None,\
                            local_MXG       = None,\
                            produce_warning = None\
                            ):
        """
        Check if NXPE or NYPE is consistent with nproc

        Input
        type_txt        -   can be either 'NXPE' or 'NYPE' and is specifying
                            whether NXPE or NYPE should be checked
        local_MXG       -   the current MXG
        produce_warning -   whether or not a warning should be produced
        """

        for size_nr in range(len(self._nx)):
            # Check the type
            if type_txt == 'NXPE':
                print("Checking nx = " + str(self._nx[size_nr]) +\
                      " with NXPE = " + str(self._NXPE[size_nr]))
            elif type_txt == 'NYPE':
                print("Checking ny = " + str(self._ny[size_nr]) +\
                      " with NYPE = " + str(self._NYPE[size_nr]))
            # Check to see if succeeded
            init_split_found = False
            cur_split_found  = False
            add_number = 1
            # Counter to see how many times the while loop has been
            # called
            count = 0

            #{{{While cur_split_found == False
            while cur_split_found == False:
                # The same check as below is performed internally in
                # BOUT++ (see boutmesh.cxx under
                # if((MX % NXPE) != 0)
                # and
                # if((MY % NYPE) != 0)
                if type_txt == 'NXPE':
                    MX = self._nx[size_nr] - 2*local_MXG
                    # self._nproc is called NPES in boutmesh
                    if (MX % self._NXPE[size_nr]) == 0:
                        # If the test passes
                        cur_split_found = True
                    # Check if cur_split_found is true, eventually
                    # update the add_number
                    add_number, produce_warning = self._check_cur_split_found(\
                                                             cur_split_found,\
                                                             produce_warning,\
                                                             add_number,\
                                                             size_nr,\
                                                             using_nx = True,\
                                                             using_ny = False)
                elif type_txt == 'NYPE':
                    MY = self._ny[size_nr]
                    # self._nproc is called NPES in boutmesh
                    if (MY % self._NYPE[size_nr]) == 0:
                        # If the test passes
                        cur_split_found = True
                    # Check if cur_split_found is true, eventually
                    # update the add_number
                    add_number, produce_warning = self._check_cur_split_found(\
                                                             cur_split_found,\
                                                             produce_warning,\
                                                             add_number,\
                                                             size_nr,\
                                                             using_nx = True,\
                                                             using_ny = False)

                #{{{ Check if the split was found the first go.
                # This will be used if self_allow_size_modification is
                # off, or if we are using a grid file
                if count == 0 and cur_split_found:
                    init_split_found = True
                #}}}

                # Add one to the counter
                count += 1
            #}}}

            # Check if initial split succeeded
            if type_txt == 'NXPE':
                self._check_init_split_found(init_split_found, size_nr,\
                                             test_nx = True, test_ny = False,
                                             produce_warning = produce_warning)
            elif type_txt == 'NYPE':
                self._check_init_split_found(init_split_found, size_nr,\
                                             test_nx = False, test_ny = True,
                                             produce_warning = produce_warning)
    #}}}
#}}}

#{{{Function called by _prepare_dmp_folder
#{{{_get_folder_name
    def _get_folder_name(self, combination):
        """Returning the folder name where the data will be stored.

        If all options are given the folder structure should be on the
        form solver/method/nout_timestep/mesh/additional/grid"""

        # Combination is one of the combination of the data members
        # which is used as the command line arguments in the run
        combination = combination.split()

        #{{{Append from eventual grid file
        # If there is a grid file, we will extract the values from the
        # file, and put it into this local combination variable, so that
        # a proper dmp folder can be made on basis on the variables
        # A flag to see whether or not the grid file was found
        grid_file_found = False
        # Check if grid is in element, and extract its path
        for elem in combination:
            if elem[0:5] == 'grid=':
                cur_grid = elem.replace('grid=','')
                grid_file_found = True

        # If the grid file is found, open it
        if grid_file_found:
            # Open (and automatically close) the grid files
            f = DataFile(cur_grid)
            # Search for mesh types in the grid file
            mesh_types = [\
                          ("mesh:", "nx"       ) ,\
                          ("mesh:", "ny"       ) ,\
                          ("mesh:", "nz"       ) ,\
                          ("mesh:", "zperiod"  ) ,\
                          ("mesh:", "zmin"     ) ,\
                          ("mesh:", "zmax"     ) ,\
                          ("mesh:", "dx"       ) ,\
                          ("mesh:", "dy"       ) ,\
                          ("mesh:", "dz"       ) ,\
                          ("mesh:", "ixseps1"  ) ,\
                          ("mesh:", "ixseps2"  ) ,\
                          ("mesh:", "jyseps1_1") ,\
                          ("mesh:", "jyseps1_2") ,\
                          ("mesh:", "jyseps2_1") ,\
                          ("mesh:", "jyseps2_2") ,\
                          ("",      "MXG")       ,\
                          ("",      "MYG")       ,\
                          ]
            for mesh_type in mesh_types:
                grid_variable = f.read(mesh_type[1])
                # If the variable is found
                if grid_variable is not None:
                    # Append it to the combinations list
                    combination.append(mesh_type[0] +\
                                       mesh_type[1] +\
                                       "=" + str(grid_variable))
        #}}}

        # Make lists for the folder-type, so that we can append the
        # elements in the combination folders if it is found
        solver        = []
        method        = []
        nout_timestep = []
        mesh          = []
        additional    = []
        grid_file     = []

        # We will loop over the names describing the methods used
        # Possible directional derivatives
        dir_derivatives = ['ddx', 'ddy', 'ddz']

        # Check trough all the elements of combination
        for elem in combination:

            # If 'solver' is in the element
            if 'solver' in elem:
                # Remove 'solver:' and append it to the mesh folder
                cur_solver = elem.replace('solver:','')
                cur_solver = cur_solver.replace('=','_')
                # Append it to the solver folder
                solver.append(cur_solver)

            # If nout or timestep is in the element
            elif ('nout' in elem) or\
                 ('timestep' in elem):
                # Remove '=', and append it to the
                # nout_timestep folder
                nout_timestep.append(elem.replace('=','_'))

            # If any quantity related to mesh is in the combination
            elif ('mesh' in elem) or\
                 ('MXG' in elem) or\
                 ('MYG' in elem) or\
                 ('NXPE' in elem) or\
                 ('NYPE' in elem) or\
                 ('zperiod' in elem) or\
                 ('zmin' in elem) or\
                 ('zmax' in elem) or\
                 (('dx' in elem) and not('ddx' in elem)) or\
                 (('dy' in elem) and not('ddy' in elem)) or\
                 (('dz' in elem) and not('ddz' in elem)):
                # Remove 'mesh:', and append it to the mesh folder
                cur_mesh = elem.replace('mesh:','')
                cur_mesh = cur_mesh.replace('=','_')
                mesh.append(cur_mesh)

            # If a grid file is in the combination
            elif (elem[0:4] == 'grid'):
                # Remove .grd .nc and =
                cur_grid = elem.replace('.grd','')
                cur_grid = cur_grid.replace('.nc','')
                cur_grid = cur_grid.replace('=','_')
                grid_file.append(cur_grid)

            # If the element is none of the above
            else:
                # It could either be a dir derivative
                # Set a flag to state if any of the dir derivative was
                # found in the combination
                dir_derivative_set = False
                # If any of the methods are in combination
                for dir_derivative in dir_derivatives:
                    if dir_derivative in elem:
                        # Remove ':', and append it to the
                        # method folder
                        cur_method = elem.replace(':','_')
                        cur_method = cur_method.replace('=','_')
                        method.append(cur_method)
                        dir_derivative_set = True

                # If the dir_derivative_set was not set, the only
                # possibility left is that the element is an
                # 'additional' option
                if not(dir_derivative_set):
                    # Replace ':' and '=' and append it to the
                    # additional folder
                    cur_additional = elem.replace(':','_')
                    cur_additional = cur_additional.replace('=','_')
                    cur_additional = cur_additional.replace('"','-')
                    cur_additional = cur_additional.replace("'",'-')
                    cur_additional = cur_additional.replace('(',',')
                    cur_additional = cur_additional.replace(')',',')
                    additional.append(cur_additional)

        # We sort the elements in the various folders alphabetically,
        # to ensure that the naming convention is always the same, no
        # matter how the full combination string looks like
        # Sort alphabetically
        solver.sort()
        #{{{ Manual sort solver
        # We want 'type' to be first, and 'atol' and 'rtol' to be last
        sort_these = [\
                      ('type',0) ,\
                      ('atol',-1),\
                      ('rtol',-1) \
                     ]
        # Loop through everything we want to sort
        for sort_this in sort_these:
            # Flag to check if found
            found_string = False
            for elem_nr, elem in enumerate(solver):
                if sort_this[0] in elem:
                    swap_nr = elem_nr
                    # Set the flag that the string is found
                    found_string = True
            # If type was found
            if found_string != False:
                # Swap the elements in the solver
                solver[sort_this[1]], solver[swap_nr] =\
                        solver[swap_nr], solver[sort_this[1]]
        #}}}
        method.sort()
        nout_timestep.sort()
        mesh.sort()
        additional.sort()
        grid_file.sort()

        # Combine the elements in the various folders
        solver              = ['_'.join(solver)]
        method              = ['_'.join(method)]
        nout_timestep       = ['_'.join(nout_timestep)]
        mesh                = ['_'.join(mesh)]
        additional          = ['_'.join(additional)]
        grid_file           = ['_'.join(grid_file)]

        # Put all the folders into the combination_folder
        combination_folder = [\
                              solver       ,\
                              method       ,\
                              nout_timestep,\
                              mesh         ,\
                              additional   ,\
                              grid_file     \
                             ]
        # We access the zeroth element (if given) as the folders are
        # given as a list
        combination_folder = [folder[0] for folder in combination_folder\
                              if (folder != []) and (folder !=[''])]

        # Make the combination folder as a string
        combination_folder = '/'.join(combination_folder)

        return combination_folder
#}}}

#{{{_create_folder
    def _create_folder(self, folder):
        """Creates a folder if it doesn't exists"""

        if not os.path.exists(folder):
            os.makedirs(folder)
            print(folder + " created\n")
#}}}

#{{{_copy_restart_files
    def _copy_restart_files(self):
        """
        Function which copies restart files from self._restart_from
        """
        # Check for files in dmp_folder
        if len(glob.glob(os.path.join(self._dmp_folder,'*restart*'))) !=0 or\
           len(glob.glob(os.path.join(self._dmp_folder,'*dmp*'))) !=0:
            message = "Restart or dmp files was found in " + self._dmp_folder +\
                      " when restart_from was set. Run skipped."
            self._warning_printer(message)
            self._warnings.append(message)
            skip_run = True
        else:
            skip_run = False

        if not skip_run:
            print("\nCopying files from {0} to {1}\n".\
                  format(self._restart_from, self._dmp_folder))

            # Files with these extension will be given the
            # additional extension .cpy when copied to the destination
            # folder
            extensions_w_cpy = ['inp', 'log.*']

            if self._cpy_source:
                extensions_w_cpy.extend(['cc' , 'cpp'  , 'cxx', 'C'  , 'c++',\
                                        'h'  , 'hpp'  , 'hxx', 'h++'])

            # Additional files that will be copied to the destination
            # folder
            extensions = [*extensions_w_cpy, 'restart.*']

            if self._restart == "append":
                extensions.append("dmp.*")

            # Copy for all files in the extension
            for extension in extensions:
                file_names =\
                    glob.glob(os.path.join(self._restart_from, '*.'+extension))
                for cur_file in file_names:
                    # Check if any of the extensions matches the current
                    # string (must strip the '.' as "in" does not accept
                    # wildcards
                    if any([ewc.split('.')[0] in cur_file
                            for ewc in extensions_w_cpy]):
                        # Add ".cpy" to the file name (without the path)
                        name = os.path.split(cur_file)[-1] + '.cpy'
                        shutil.copy2(cur_file, os.path.join(self._dmp_folder, name))
                    else:
                        shutil.copy2(cur_file, self._dmp_folder)

        return skip_run
#}}}

#{{{Save _move_old_runs
    def _move_old_runs(self):
        """Move old runs if restart is set to 'overwrite'"""
        print("Moving old runs\n")
        # Check for folders in the dmp directory
        directories = [\
                       name for name in\
                       os.listdir(self._dmp_folder) if\
                       os.path.isdir(os.path.join(\
                                    self._dmp_folder, name))\
                      ]
        # Find occurrences of 'run' in these folders
        prev_runs = [name for name in directories if 'run' in name]
        # Check that the list is not empty
        if len(prev_runs) != 0:
            # Sort the folders alphabetically
            prev_runs.sort()
            # Pick the last of prev_runs
            prev_runs = prev_runs[-1]
            # Pick the number from the last run
            # First split the string
            overwrite_nr = prev_runs.split('_')
            # Pick the last element of overwrite_nr, and cast it
            # to an integer
            overwrite_nr = int(overwrite_nr[-1])
            # Add one to the overwrite_nr, as we want to create
            # a new directory
            overwrite_nr += overwrite_nr
        else:
            # Set the overwrite_nr
            overwrite_nr = 1
        # Create the folder for the previous runs
        self._create_folder(\
                os.path.join(self._dmp_folder, 'run_' +\
                             str(overwrite_nr)))

        extensions_to_move = ['cpy', 'log.*', 'dmp.*',\
                              'cc' , 'cpp'  , 'cxx'  , 'C'  , 'c++',\
                              'h'  , 'hpp'  , 'hxx'  , 'h++']

        for extension in extensions_to_move:
            file_names =\
                glob.glob(os.path.join(self._dmp_folder, '*.'+extension))

            # Cast to unique file_names
            file_names = set(file_names)

            # Move the files
            for cur_file in file_names:
                dst = os.path.join(self._dmp_folder,\
                                   "run_" + str(overwrite_nr))
                shutil.move(cur_file, dst)
#}}}
#}}}

#{{{Function called by _run_driver
#{{{_single_run
    def _single_run(self, combination):
        """Makes a single MPIRUN of the program"""

        # Get the command to be used
        command = self._get_command_to_run(combination)

        # Time how long the time took
        tic = timeit.default_timer()

        # Launch the command
        status, out = launch(command,\
                             runcmd = self._MPIRUN,\
                             nproc = self._nproc,\
                             pipe = True,\
                             verbose = True)

        # If the run returns an exit code other than 0
        if status != 0:
            message = "! An error occurred. Printing the output to stdout !"
            print("\n" + "!"*len(message))
            print(message)
            print("!"*len(message) + "\n")
            print(out)
            self._errors.append("RuntimeError")
            message =  "An error occurred the run."
            message += " Please see the output above for details."
            # Search if parantheses are present, but without ' or "
            if     ('(' in combination and\
                   not(    re.search(r'\"(.*)\(', combination)\
                        or re.search(r"\'(.*)\(", combination)))\
                or (')' in combination and\
                   not(   re.search(r'\)(.*)\"', combination)
                       or re.search(r"\)(.*)\'", combination))):
                message = 'A "(" and/or ")" symbol seem to have appeared in the'
                message += " command line.\nIf this true, you can avoid"
                message += " this problem by adding an extra set of"
                message += " quotation marks. For example\n\n"
                message += "additional=('variable', 'bndry_xin',"
                message += " '\"dirichlet_o4(0.0)\")'\n"
                message += "rather than\n"
                message += "additional=('variable', 'bndry_xin',"
                message += " 'dirichlet_o4(0.0))'"
            else:
                message =  "An error occurred the run."
                message += " Please see the output above for details."
            raise RuntimeError(message)

        # Estimate elapsed time
        toc = timeit.default_timer()
        elapsed_time = toc - tic

        return out, elapsed_time
#}}}

#{{{_append_run_log
    def _append_run_log(self, start, run_no, run_time):
        """Appends the run_log"""

        # Convert seconds to H:M:S
        run_time = str(datetime.timedelta(seconds=run_time))

        start_time = (str(start.year) + '-' + str(start.month) + '-' +\
                      str(start.day) + '.' + str(start.hour) + ":" +\
                      str(start.minute) + ':' + str(start.second))

        # If the run is restarted with initial values from the last run
        if self._restart:
            dmp_line = self._dmp_folder + '-restart-'+self._restart
            if self._restart_from:
                dmp_line += " from " + self._restart_from
        else:
            dmp_line = self._dmp_folder

        # Line to write
        line = [start_time, self._run_type, run_no, run_time, dmp_line]
        # Opens for appending
        log_format = '{:<19}   {:^9}   {:^6}   {:<17}   {:<}'
        with open(self._run_log , "a") as f:
            f.write(log_format.format(*line) + '\n')
#}}}
#}}}

#{{{Function called by _get_possibilities
#{{{_generate_possibilities
    def _generate_possibilities(self, variables=None, section=None, name=None):
        """Generate the list of strings of possibilities"""

        if variables is not None:
            # Set the section name correctly
            if section != "":
                section = section + ":"
            else:
                section = ""
            # Set the combination of the variable
            var_possibilities = []
            # Find the number of different dimensions
            for var in variables:
                var_possibilities.append(section + name + '=' + str(var))
        else:
            var_possibilities = []

        return var_possibilities
#}}}
#}}}

#{{{Functions called by _get_combinations
#{{{_get_swapped_input_list
    def _get_swapped_input_list(self, input_list):
        """Finds the element in the input list, which corresponds to the
        self._sort_by criterion. The element is swapped with the last
        index, so that itertools.product will make this the fastest
        varying variable"""

        # We make a sort list containing the string to find in the
        # input_list
        sort_list = []

        # We loop over the elements in self._sort_by to find what
        # string we need to be looking for in the elements of the lists
        # in input_list
        for sort_by in self._sort_by:
            # Find what list in the input_list which contains what we
            # would sort by

            #{{{ If we would like to sort by the spatial domain
            if sort_by == 'spatial_domain':
                # nx, ny and nz are all under the section 'mesh'
                find_in_list = 'mesh'
            #}}}

            #{{{ If we would like to sort by the temporal domain
            elif sort_by == 'temporal_domain':
                # If we are sorting by the temporal domain, we can either
                # search for timestep or nout
                if self._timestep is not None:
                    find_in_list = 'timestep'
                elif self._nout is not None:
                    find_in_list = 'nout'
            #}}}

            #{{{ If we would like to sort by the method
            elif (sort_by == 'ddx_first') or\
                 (sort_by == 'ddx_second') or\
                 (sort_by == 'ddx_upwind') or\
                 (sort_by == 'ddx_flux') or\
                 (sort_by == 'ddy_first') or\
                 (sort_by == 'ddy_second') or\
                 (sort_by == 'ddy_upwind') or\
                 (sort_by == 'ddy_flux') or\
                 (sort_by == 'ddz_first') or\
                 (sort_by == 'ddz_second') or\
                 (sort_by == 'ddz_upwind') or\
                 (sort_by == 'ddz_flux'):
                find_in_list = sort_by.replace('_',':')
            #}}}

            #{{{ If we would like to sort by the solver
            elif sort_by == 'solver':
                find_in_list = sort_by
            #}}}

            #{{{ If we would like to sort by anything else
            else:
                find_in_list = sort_by
            #}}}

            # Append what to be found in the input_list
            sort_list.append(find_in_list)

        # For all the sort_list, we would like check if the match
        # can be found in any of the elements in input_list
        # Appendable list
        lengths = []
        for sort_nr, sort_by_txt in enumerate(sort_list):
            # Make a flag to break the outermost loop if find_in_list is
            # found
            break_outer = False
            # Loop over the lists in the input_list to find the match
            for elem_nr, elem in enumerate(input_list):
                # Each of the elements in this list is a string
                for string in elem:
                    # Check if fins_in_list is in the string
                    if sort_by_txt in string:
                        # If there is a match, store the element number
                        swap_from_index = elem_nr
                        # Check the length of the element (as this is
                        # the number of times the run is repeated, only
                        # changing the values of sort_by [defining a
                        # group])
                        lengths.append(len(elem))
                        # Break the loop to save time
                        break_outer = True
                        break
                # Break the outer loop if find_in_list_is_found
                if break_outer:
                    break

            # As it is the last index which changes the fastest, we swap the
            # element where the find_in_list was found with the last element
            input_list[swap_from_index], input_list[-(sort_nr + 1)] =\
                    input_list[-(sort_nr + 1)], input_list[swap_from_index]

        # The number of runs in one 'group'
        # Initialize self._len_group with one as we are going to
        # multiply it with all the elements in lengths
        self._len_group = 1
        for elem in lengths:
            self._len_group *= elem

        return input_list
#}}}
#}}}

#{{{Function called by _single_run
#{{{_get_command_to_run
    def _get_command_to_run(self, combination):
        """ Returns a string of the command which will run the BOUT++
        program"""

        # Creating the arguments
        arg = " -d " + self._dmp_folder + " " + combination

        # If the run is set to overwrite
        if self._restart == 'overwrite':
            arg += ' restart'
        elif self._restart == 'append':
            arg += ' restart append'

        # Replace excessive spaces with a single space
        arg = ' '.join(arg.split())
        command = "./" + self._program_name + " " + arg

        return command
#}}}
#}}}

#{{{Functions called from several places in the code
    #{{{_warning_printer
    def _warning_printer(self, message):
        """Function for printing warnings"""

        print('\n'*3 + '*'*37 + 'WARNING' + '*'*36)
        # Makes sure that no more than 80 characters are printed out at
        # the same time
        message_chunks=[]
        for chunk in self._message_chunker(message):
            rigth_padding = ' '*(76 - len(chunk))
            print('* ' + chunk + rigth_padding + ' *')
        print('*'*80 + '\n'*3)
    #}}}

    #{{{_message_chunker
    def _message_chunker(self, message, chunk=76):
        """Generator used to chop a message so it doesn't exceed some
        width"""

        for start in range(0, len(message), chunk):
            yield message[start:start + chunk]
    #}}}
#}}}
#}}}



#{{{class PBS_runner
class PBS_runner(basic_runner):
#{{{docstring
    """Class for mpi running one or several runs with BOUT++.
    Works like the basic_runner, but submits the jobs to a Portable
    Batch System (PBS).

    For the additional member data, see the docstring of __init__.

    For more info check the docstring of bout_runners.
    """
#}}}

# The constructor
#{{{__init__
    def __init__(self,\
                 BOUT_nodes            = 1         ,\
                 BOUT_ppn              = 1         ,\
                 BOUT_walltime         = None      ,\
                 BOUT_queue            = None      ,\
                 BOUT_mail             = None      ,\
                 BOUT_run_name         = None      ,\
                 post_process_nproc    = None      ,\
                 post_process_nodes    = None      ,\
                 post_process_ppn      = None      ,\
                 post_process_walltime = None      ,\
                 post_process_queue    = None      ,\
                 post_process_mail     = None      ,\
                 post_process_run_name = None      ,\
                 **kwargs):
        #{{{docstring
        """The constructor of the PBS_runner.

        All the member data is set to None by default, with the
        exception of BOUT_nodes (default=1) and BOUT_ppn (default = 4).

        Input:
        BOUT_nodes              -    Number of nodes for one submitted
                                     BOUT job
        BOUT_ppn                -    Processors per node for one
                                     submitted BOUT job
        BOUT_walltime           -    Maximum wall time for one submitted
                                     BOUT job
        BOUT_queue              -    The queue to submit the BOUT jobs
        BOUT_mail               -    Mail address to notify when a BOUT job
                                     has finished
        BOUT_run_name           -    Name of the BOUT run on the cluster
                                     (optional)
        post_process_nproc      -    Total number of processors for one
                                     submitted post processing job
        post_process_nodes      -    Number of nodes for one submitted
                                     post processing job
        post_process_ppn        -    Processors per node for one
                                     submitted BOUT job
        post_process_walltime   -    Maximum wall time for one
                                     submitting post processing job
        post_process_queue      -    The queue to submit the post
                                     processing jobs
        post_process_mail       -    Mail address to notify when a post
                                     processing job has finished
        post_process_run_name   -    Name of the post processing run on the
                                     cluster (optional)
        **kwargs                -    As the constructor of bout_runners
                                     is called, this additional keyword
                                     makes it possible to specify the
                                     member data of bout_runners in the
                                     constructor of PBS_runner (i.e.
                                     nprocs = 1
                                     is an allowed keyword argument in
                                     the constructor of PBS_runner).
                                     For a full list of possible
                                     keywords, see the docstring of the
                                     bout_runners constructor.
        """
        #}}}

        # Note that the constructor accepts additional keyword
        # arguments (**kwargs). These must match the keywords of the
        # parent class 'basic_runner', which is called by the 'super'
        # function below

        # Call the constructor of the superclass
        super(PBS_runner, self).__init__(**kwargs)

        # Options set for the BOUT runs
        self._BOUT_nodes            = BOUT_nodes
        self._BOUT_ppn              = BOUT_ppn
        self._BOUT_walltime         = BOUT_walltime
        self._BOUT_mail             = BOUT_mail
        self._BOUT_queue            = BOUT_queue
        self._BOUT_run_name         = BOUT_run_name
        # Options set for the post_processing runs
        self._post_process_nproc    = post_process_nproc
        self._post_process_nodes    = post_process_nodes
        self._post_process_ppn      = post_process_ppn
        self._post_process_walltime = post_process_walltime
        self._post_process_mail     = post_process_mail
        self._post_process_queue    = post_process_queue
        self._post_process_run_name = post_process_run_name

        # Options set for all runs
        self._run_type      = 'basic_PBS'

        # Error check the input data
        self._check_for_PBS_instance_error()

        # Initialize the jobid returned from the PBS
        self._PBS_id = []
#}}}

# The run_driver
#{{{_run_driver
    def _run_driver(self, combination, run_no):
        """The machinery which actually performs the run"""

        # Submit the job to the queue
        self._single_submit(combination, run_no, append_to_run_log = True)
#}}}

#{{{Functions called by the constructor
    #{{{_check_for_PBS_instance_error
    def _check_for_PBS_instance_error(self):
        """Check if there are any type errors when creating the object"""

        #{{{Check if BOUT_ppn and BOUT_nodes have the correct type
        # BOUT_ppn and BOUT_nodes are set by default, however, we must check that
        # the user has not given them as wrong input
        if type(self._BOUT_ppn) != int:
            message  = "BOUT_ppn is of wrong type\n"+\
                       "BOUT_ppn must be given as a int"
            self._errors.append("TypeError")
            raise TypeError(message)
        if type(self._BOUT_nodes) != int:
            message  = "BOUT_nodes is of wrong type\n"+\
                       "BOUT_nodes must be given as a int"
            self._errors.append("TypeError")
            raise TypeError(message)
        #}}}

        #{{{Check that nprocs, BOUT_nodes and BOUT_ppn is consistent
        if self._nproc > (self._BOUT_nodes * self._BOUT_ppn):
            message = 'Must have nproc <= BOUT_nodes * BOUT_ppn'
            self._errors.append("TypeError")
            raise TypeError(message)
        #}}}

        #{{{Check all the proper post_process data is set if any is set
        check_if_set = [\
                        self._post_process_nproc,\
                        self._post_process_nodes,\
                        self._post_process_ppn,\
                       ]
        # All elements of check_if_set must be set if any is set
        not_None = 0
        for check in check_if_set:
            if check is not None:
                not_None += 1

        if (not_None != 0) and (not_None != len(check_if_set)):
            message = "If any of post_process_nproc, post_process_nodes,"+\
                       " post_process_ppn and post_process_walltime is"+\
                       " set, all others must be set as well."
            self._errors.append("TypeError")
            raise TypeError(message)
        #}}}

        #{{{Check if post_process_ppn and post_process_nodes is int if set
        check_if_int = [\
                        (self._post_process_nodes, 'post_process_nodes') ,\
                        (self._post_process_ppn,   'post_process_ppn') \
                       ]
        self._check_for_correct_type(var = check_if_int,\
                                      the_type = int,\
                                      allow_iterable = False)
        #}}}

        #{{{Check that post_process_nprocs,nodes,ppn is consistent if set
        if self._post_process_nproc is not None:
            if self._post_process_nproc > \
                    (self._post_process_nodes * self._post_process_ppn):
                message = 'Must have post_process_nproc <= post_process_nodes'+\
                          '* post_process_ppn'
                self._errors.append("TypeError")
                raise TypeError(message)
        #}}}

        #{{{Check if walltime, mail and queue is a string if set
        check_if_str = [\
                        (self._BOUT_walltime,         'BOUT_walltime')        ,\
                        (self._BOUT_mail,             'BOUT_mail')            ,\
                        (self._BOUT_queue,            'BOUT_queue')           ,\
                        (self._BOUT_run_name,         'BOUT_run_name')        ,\
                        (self._post_process_walltime, 'BOUT_walltime')        ,\
                        (self._post_process_mail,     'post_process_mail')    ,\
                        (self._post_process_queue,    'post_process_queue')   ,\
                        (self._post_process_run_name, 'post_process_run_name') \
                       ]
        self._check_for_correct_type(var = check_if_str,\
                                      the_type = str,\
                                      allow_iterable = False)
        #}}}

        #{{{Check that walltime is on correct format
        # A list to loop over
        walltimes = []
        # Append the walltimes if set
        if self._BOUT_walltime is not None:
            walltimes.append((self._BOUT_walltime,\
                              'BOUT_walltime'))
        if self._post_process_walltime is not None:
            walltimes.append((self._post_process_walltime,\
                              'post_process_walltime'))

        # Loop over the walltimes
        for walltime in walltimes:
            # Set a flag which states whether or not the check was
            # successful
            success = True
            # Split the walltime string
            walltime_list = walltime[0].split(':')
            # Check that the list has three elements
            if len(walltime_list) == 3:

                # Check that seconds is on the format SS
                if len(walltime_list[2]) == 2:
                    # Check that the last element (seconds) is a digit (int)
                    if walltime_list[2].isdigit():
                        # Check that the element is less than 59
                        if int(walltime_list[2]) > 59:
                            success = False
                    # Seconds is not a digit
                    else:
                        success = False
                # Seconds is not on the format SS
                else:
                    success = False


                # Do the same for the second last element (minutes)
                if len(walltime_list[1]) == 2:
                    # Check that the last element (seconds) is a digit (int)
                    if walltime_list[1].isdigit():
                        if int(walltime_list[1]) > 59:
                            success = False
                    # Minutes is not a digit
                    else:
                        success = False
                # Seconds is not on the format SS
                else:
                    success = False

                # Check that the first element (hours) is a digit
                if not(walltime_list[0].isdigit()):
                        success = False

            # walltime_list does not have three elements
            else:
                success = False

            if not(success):
                message = walltime[1] + " must be on the form H...H:MM:SS"
                self._errors.append("TypeError")
                raise TypeError(message)
        #}}}
    #}}}
#}}}

#{{{Functions called by _error_check_for_run_input
    #{{{_check_for_child_class_errors
    def _check_for_child_class_errors(
                                   self                        ,\
                                   remove_old                  ,\
                                   post_processing_function    ,\
                                   post_process_after_every_run \
                                   ):
        """Function which check for errors in a child class."""

        # Check member data is set if post_processing_function is not None
        if post_processing_function is not None:
            check_if_set = [\
                            self._post_process_nproc,\
                            self._post_process_nodes,\
                            self._post_process_ppn,\
                           ]
            # All elements of check_if_set must be set if any is set
            not_None = 0
            for check in check_if_set:
                if check is not None:
                    not_None += 1

            if (not_None != 0) and (not_None != len(check_if_set)):
                message = "post_process_nproc, post_process_nodes,"+\
                          " and post_process_ppn and must"+\
                          " be set if post_processing_function is set."
                self._errors.append("TypeError")
                raise TypeError(message)
    #}}}
#}}}

#{{{Functions called by the execute_runs
    #{{{ _print_run_or_submit
    def _print_run_or_submit(self):
        """Prints 'Submitting'"""
        print("\nSubmitting:")
    #}}}
#}}}

#{{{Functions called by _run_driver
    #{{{_single_submit
    def _single_submit(self, combination, run_no, append_to_run_log = None):
        """Submit a single BOUT job and submit the jobid to self._PBS_id"""

        # Get the script (as a string) which is going to be
        # submitted
        job_string = self._get_job_string(run_no,\
                                          combination,\
                                          append_to_run_log)

        # The submission
        PBS_id = self._submit_to_PBS(job_string)
        self._PBS_id.append(PBS_id)
    #}}}

    #{{{_call_post_processing_function
    def _call_post_processing_function(\
                                       self            ,\
                                       function = None ,\
                                       folders  = None ,\
                                       **kwargs         \
                                       ):
        """Function which submits the post processing to the PBS

        This is done by making a self deleting temporary python file
        that will be called by a PBS script."""

        #{{{ Create a python script, calling the post-processing function
        # Get the start_time (to be used in the name of the file)
        start_time = self._get_start_time()

        # The name of the file
        python_name = "tmp_" + function.__name__ + '_'+start_time+'.py'

        # Make the script
        python_tmp  = '#!/usr/bin/env python\n'
        python_tmp += 'import os\n'
        # Import the post processing function
        python_tmp += 'from ' + function.__module__ +\
                      ' import ' + function.__name__ + '\n'
        # Convert the keyword args to proper arguments
        # Appendable list
        arguments = []
        for key in kwargs.keys():
            if type(kwargs[key]) != str:
                # If the value is not a string, we can append it directly
                arguments.append(str(key) + '=' + str(kwargs[key]))
            else:
                # If the value is a string, we need to put quotes around
                arguments.append(str(key) + '="' + str(kwargs[key]) +'"')

        # Put a comma in between the arguments
        arguments = ', '.join(arguments)
        # Call the post processing function
        if type(folders) == list:
            python_tmp+=function.__name__+"("+str(folders)+","+arguments+")\n"
        elif type(folders) == str:
            python_tmp+=function.__name__+"('"+str(folders)+"',"+arguments+")\n"
        # When the script has run, it will delete itself
        python_tmp += "os.remove('" + python_name + "')\n"

        # Write the python script
        with open(python_name, "w") as f:
            f.write(python_tmp)
        #}}}

        #{{{Create and submit the shell script
        # Creating the job string
        if self._post_process_run_name is None:
            job_name = 'post_process_' + function.__name__ + '_'+ start_time
        else:
            job_name = self._post_process_run_name

        # Get core of the job string
        job_string = self._create_PBS_core_string(\
                                job_name         = job_name                   ,\
                                nodes            = self._post_process_nodes   ,\
                                ppn              = self._post_process_ppn     ,\
                                walltime         = self._post_process_walltime,\
                                mail             = self._post_process_mail    ,\
                                queue            = self._post_process_queue    \
                                )
        # Call the python script in the submission

        job_string += 'python ' + python_name + '\n'
        job_string += 'exit'

        # Create the dependencies
        dependencies = ':'.join(self._PBS_id)

        # Submit the job
        print('\nSubmitting the post processing function "' +\
                function.__name__ + '"\n')
        PBS_id = self._submit_to_PBS(job_string, dependent_job = dependencies)
        #}}}
    #}}}
#}}}

#{{{ Functions called by _single_submit
    #{{{_get_job_string
    def _get_job_string(self, run_no, combination, append_to_run_log):
        """Make a string which will saved as a shell script before being sent to
        the PBS queue."""

        #{{{Make the job name based on the combination
        # Split the name to a list
        combination_name = combination.split(' ')
        # Remove whitespace
        combination_name = [element for element in combination_name\
                            if element != '']
        # Collect the elements
        combination_name = '_'.join(combination_name)
        # Replace bad characters
        combination_name = combination_name.replace(':','')
        combination_name = combination_name.replace('=','-')

        # Name of job
        if self._BOUT_run_name is None:
            job_name =\
                combination_name + '_' + self._directory + '_' + str(run_no)
        else:
            job_name = self._BOUT_run_name
        #}}}

        #{{{Make the main command that will be used in the PBS script
        command = self._get_command_to_run( combination )
        command = 'mpirun -np ' + str(self._nproc) + ' ' + command

        # Print the command
        print(command + '\n')
        #}}}

        #{{{ Creating the core job string
        job_string = self._create_PBS_core_string(\
                                job_name         = job_name           ,\
                                nodes            = self._BOUT_nodes   ,\
                                ppn              = self._BOUT_ppn     ,\
                                walltime         = self._BOUT_walltime,\
                                mail             = self._BOUT_mail    ,\
                                queue            = self._BOUT_queue    \
                                )
        #}}}

        if append_to_run_log:
            #{{{ Get the time for start of the submission
            start = datetime.datetime.now()
            start_time = (str(start.year) + '-' + str(start.month) + '-' +\
                          str(start.day) + '.' + str(start.hour) + ":" +\
                          str(start.minute) + ':' + str(start.second))
            #}}}

            #{{{ Start the timer
            job_string += 'start=`date +%s`\n'
            # Run the bout program
            job_string += command + '\n'
            # end the timer
            job_string += 'end=`date +%s`\n'
            # Find the elapsed time
            job_string += 'time=$((end-start))\n'
            # The string is now in seconds
            # The following procedure will convert it to H:M:S
            job_string += 'h=$((time/3600))\n'
            job_string += 'm=$((($time%3600)/60))\n'
            job_string += 's=$((time%60))\n'
            #}}}

            #{{{ Append to the run log
            # Ideally we would check if any process were writing to
            # run_log.txt
            # This could be done with lsof command as described in
            # http://askubuntu.com/questions/14252/how-in-a-script-can-i-determine-if-a-file-is-currently-being-written-to-by-ano
            # However, lsof is not available on all clusters

            # Using the same formatting as in _append_run_log, we are going
            # to echo the following to the run_log when the run is finished
            job_string += "echo '" +\
                          '{:<19}'.format(start_time)     + " "*3  +\
                          '{:^9}'.format(self._run_type)  + " "*3  +\
                          '{:^6}'.format(str(run_no))     + " "*3  +\
                          "'$h':'$m':'$s"                 + " "*10 +\
                          '{:<}'.format(self._dmp_folder) + " "*3  +\
                          " >> $PBS_O_WORKDIR/" + self._directory  +\
                          "/run_log.txt\n"
            #}}}

        # Exit the qsub
        job_string += 'exit'

        return job_string
    #}}}
#}}}

#{{{Functions called by _submit_to_PBS
#{{{_get_start_time
    def _get_start_time(self):
        """Returns a string of the current time down to micro precision"""

        # The time is going to be appended to the  job name and python name
        time_now = datetime.datetime.now()
        start_time = str(getattr(time_now, 'hour')) + '-' +\
                     str(getattr(time_now, 'minute'))+ '-' +\
                     str(getattr(time_now, 'second'))
        # In case the process is really fast, so that more than one job
        # is submitted per second, we add a microsecond in the
        # names for safety
        start_time += '-' + str(getattr(time_now,'microsecond'))
        return start_time
#}}}
#}}}

#{{{Functions called by several functions
#{{{_create_PBS_core_string
    def _create_PBS_core_string(\
                                self           ,\
                                job_name = None,\
                                nodes    = None,\
                                ppn      = None,\
                                walltime = None,\
                                mail     = None,\
                                queue    = None \
                                ):
        """Creates the core of a PBS script as a string"""

        # Shebang line
        job_string = '#!/bin/bash\n'
        # The job name
        job_string += '#PBS -N ' + job_name + '\n'
        job_string += '#PBS -l nodes=' + str(nodes) + ':ppn=' + str(ppn)  + '\n'
        # If walltime is set
        if walltime is not None:
            # Wall time, must be in format HOURS:MINUTES:SECONDS
            job_string += '#PBS -l walltime=' + walltime + '\n'
        # If submitting to a specific queue
        if queue is not None:
            job_string += '#PBS -q ' + queue + '\n'
        job_string += '#PBS -o ' + os.path.join(self._dmp_folder, job_name) +\
                      '.log' + '\n'
        job_string += '#PBS -e ' + os.path.join(self._dmp_folder, job_name) +\
                      '.err' + '\n'
        # If we want to be notified by mail
        if mail is not None:
            job_string += '#PBS -M e' + mail + '\n'
        # #PBS -m abe
        # a=aborted b=begin e=ended
        job_string += '#PBS -m e ' + '\n'
        # cd to the folder you are sending the qsub from
        job_string += 'cd $PBS_O_WORKDIR ' + '\n'

        return job_string
#}}}

#{{{_submit_to_PBS
    def _submit_to_PBS(self, job_string, dependent_job=None):
        """Saves the job_string as a shell script, submits it and
        deletes it. Returns the output from PBS as a string"""

        # Create the name of the temporary shell script
        # Get the start_time used for the name of the script
        start_time = self._get_start_time()
        script_name = 'tmp_'+start_time+'.sh'

        # Save the string as a script
        with open(script_name, "w") as shell_script:
                shell_script.write(job_string)

        # Submit the jobs
        if dependent_job==None:
            # Without dependencies
            command = "qsub ./"+script_name
            status, output = shell(command, pipe=True)
        else:
            # If the length of the depend job is 0, then all the jobs
            # have completed, and we can carry on as usual without
            # dependencies
            if len(dependent_job) == 0:
                command = "qsub ./"+script_name
                status, output = shell(command, pipe=True)
            else:
                # With dependencies
                command = "qsub -W depend=afterok:" + dependent_job + " ./"+script_name
                status, output = shell(command, pipe=True)

        # Check for success
        if status != 0:
            if status == 208:
                message = "Runs finished before submission of the post"+\
                          " processing function. When the runs are done:"+\
                          " Run again with 'remove_old = False' to submit"+\
                          " the function."
                self._warnings.append(message)
            else:
                print("\nSubmission failed, printing output\n")
                print(output)
                self._errors.append("RuntimeError")
                message = "The submission failed with exit code " + str(status) +\
                          ", see the output above"
                raise RuntimeError(message)

        # Trims the end of the output string
        output = output.strip(' \t\n\r')

        # Delete the shell script
        try:
            os.remove(script_name)
        except FileNotFoundError:
            # Do not raise an error
            pass

        return output
#}}}
#}}}
#}}}



#{{{if __name__ == '__main__':
if __name__ == '__main__':
    """If bout_runners is run as a script, it will just call the demo
    function"""

    print("\n\nTo find out about the bout_runners, please read the user's"+\
          " manual, or have a look at 'BOUT/examples/bout_runners_example'")
#}}}

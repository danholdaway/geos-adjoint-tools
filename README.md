# GEOS Adjoint Monitoring Tools

### Use of ops_adjoint_monitoring.py:

This tool will plot the RMS and horizontal contour plots of the output of the GEOS adjoint as used for observation impacts in operations.

Several arguments are required when calling this tool:

- -g (or -\-geos_ver). This is the version of GEOS, default 522
- -v (or -\-fp_or_fpp). Whether the experiemnt is fp or fpp, default fp
- -d (or -\-date). The analysis date in question, default 20191001
- -f (or -\-fields). List of space seperated fields to be ploted, default u v tv sphu
- -l (or -\-level_to_plot). Level to plot. Can also be maxrms to plot the level coresponsing to maximum RMS.

Example:

    module load other/python/GEOSpyD/Ana2019.03_py3.7
    python ops_adjoint_monitoring.py -g 525 -v fpp -d 20190826 -f u v -l 50

This software requires netcdf, cartopy and matplotlib

README

    PyGS (Python Galaxy Survey) is designed to be an interactive class to deal with large-scale redshift surveys.
    
    Its ideal is to make the run-of-the-mill operations as simple as possible, but to provide enough flexibility
    in the background to be of real usefulness. While intended as an interactive class, for use with the python
    shell or IPython, it can of course also be scripted.
    
    Many of the more intensive routines are written in Fortran and made available through F2PY.
    
    Instantiating the class requires a filename, or folder, which contains the survey information (redshifts,
    angles, magnitudes etc.). Once the import of the data has been successful, a directory structure will be
    created which is intended to allow speed-of-access and computation, and ease-of-reuse for future 
    calculations. See the documentation for create_sample for detailed explanations.
    
    Most calculations have default options which will be used in the absence of user-specified options, and these
    are deemed to be the most common usages of the method. However, there is always the opportunity to tweak these
    as needed. Generally there are 'convenience' functions for performing a series of calculations, where it is
    common to do these together.
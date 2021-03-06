'''
Created on May 22, 2012

@author: jeven

    
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
    
'''
#!/usr/bin/python


#============================================
# IMPORT ALL NECESSARY MODULES
#============================================
import sys
import os
import errno
import time
from fort.DFT import dft

import cosmolopy.distance as cd
from cosmolopy import *

from scipy.interpolate import griddata

#import mayavi.mlab as mm
try:
    import matplotlib
    matplotlib.use("TkAgg")
    #matplotlib.interactive(True)
    import matplotlib.pyplot as plt
    from matplotlib.image import *
except:
    sys.exit("Please install Matplotlib")
try:
    import numpy as np
except:
    sys.exit("Please install Numpy")
try:
    import numpy.lib.recfunctions as nlrf
except:
    sys.exit("Please install Numpy")
try:
    import asciitable as atab
except:
    sys.exit("Please install asciitable module")
    
    
    
class PyGS(object):
    '''
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
    '''


    def __init__(self,filename,radians=True, cosmology = fidcosmo):
        '''
        Imports the survey and creates a project directory
        
        Instantiating the class requires a filename or directory in which is kept (solely) the survey information.
        This step is designed to be as seamless as possible, though survey information is kept in a variety of formats.
        To this end, note the following:
        
        The default expected input is one ascii file with columns designated (in order) as:
        
        [redshift, RA, DEC, magnitude, comoving_distance, luminosity]
        
         NOTE: only these quantities are used in any calculations in the class as yet, and therefore no other columns
               will be read. In future versions this may change.
               
        If a folder is specified, a list of variables and the corresponding files must be submitted.
        
        If a single file is specified in a different order than that given, the order must be specified.
        
        If the survey information is in a completely different format, then the best thing to to do is reformat it 
        yourself and try again. And email me the problem and I'll try to update the class to deal with it.
        
        
        Once the information is read in, a project directory is automatically created. This serves to speed up future
        calculations and simplify the layout of results from your survey. The project folder is by default named
        <survey_filename>_Project/ and included is a hidden file ".PyLSSproject" which serves to tell this class that
        the folder is in fact formatted in the default way. Within this folder the directory structure looks like:
        
        Project_Folder/
            .PyLSSproject
            correctly_formatted_survey.dat
            Sub_Sample_1/
                .PyLSSsample
                sample_description.txt
                Data/
                    correctly_formatted_sample.dat
                    Mocks/
                        correctly_formatted_mock_1.dat
                        ...
                    Randoms/
                        correctly_formatted_random_1.dat
                        ...
                Calculation_Type_1/
                    results.dat
                    ...
                ...
            ...
            
        Each data file has a default header which lets the class know how the calculation was performed and some
        other information which is useful for archiving purposes. This serves to let the class know whether or not
        to perform a calculation again, or merely read in an existing calculation.
        
        A note on cosmologies: changing the fiducial cosmology will potentially change a good many of the calculations
        for a given survey (especially those relying on distances). This is reflected in the "Sub_Sample" directory
        level. That is, two identical subsamples of the main survey will be called different subsamples if their
        cosmology is specified differently.
                
        '''
        self.is_project_file = False
        self.is_sample_file = False
        
        if "_PyLSSProject_" in os.listdir(os.path.split(filename)[0]):
            self.is_project_file = True
        
        if "_PyLSSSample_" in os.listdir(os.path.split(filename)[0]):
            self.is_sample_file = True
          
        self.cosmology = cosmology
          
        if not self.is_project_file and not self.is_sample_file:
            self.initial_file_import(filename, radians)
            self.create_project_dir(filename)
        elif self.is_project_file: 
            self.survey_file_import(filename)
            self.project_dir = os.path.split(filename)[0]
         
            c_dist = cd.comoving_distance(self.survey_data["redshift"],**self.cosmology)
            if 'mag' in self.survey_data.dtype.names:
                lum = self.LuminosityFromAppMag(c_dist)
                self.survey_data = nlrf.append_fields(self.survey_data, names=("c_dist","lum"), data=(c_dist,lum), 
                                                      dtypes=(float,float),usemask=False,asrecarray=True)
                self.survey_data.dtype.names = ("redshift",'ra','dec','mag','c_dist','lum')
            else:
                self.survey_data = nlrf.append_fields(self.survey_data, names="c_dist", data=c_dist, 
                                                      dtypes=float,usemask=False,asrecarray=True)
                self.survey_data.dtype.names = ("redshift",'ra','dec','mag','c_dist')
            
            self.create_sample_dir()
        else:
            self.survey_file_import(filename)
            self.project_dir = os.path.split(os.path.split(filename)[0])[0]
            self.sample_dir = os.path.split(filename)[0]
            self.sample_name = os.path.split(filename)[1]
        
        # Define some overall parameters of the sample.
        self.N = self.survey_data['redshift'].shape[0]
        self.survey_name = os.path.split(self.project_dir)[1]
        
        # Make a Plots Directory
        self.plots_dir = self.sample_dir+'/SurveyPlots'    
        try:
            os.makedirs(self.plots_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        
        # Make a 1-D Power Spectrum Directory
        self.power1d_dir = self.sample_dir+'/PowerSpecs_1D'    
        try:
            os.makedirs(self.power1d_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
    
        # Make a 2-D Power Spectrum Directory
        self.power2d_dir = self.sample_dir+'/PowerSpecs_2D'    
        try:
            os.makedirs(self.power2d_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
            
        # Make a Filtered surver Directory
        self.filtered_dir = self.filtered_dir+'/FilteredData'    
        try:
            os.makedirs(self.filtered_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
            
    def initial_file_import(self, filename,radians=True):
        """
        This method imports the file for the first time (the hard part...)
        
        
        """
        # If the filename given is actually a directory, we assume that all
        # files within it are components of the survey and will joined together
        # into one survey, with the title of the folder name.
        if os.path.isdir(filename):
            splits = os.listdir(filename)
            datafile = open(filename+'/'+splits[0],"r")
            self.survey_data = atab.read(datafile.read())
            datafile.close()
            for split_file in splits[1::]:
                datafile = open(filename+'/'+split_file,"r")
                self.survey_data = np.hstack(self.survey_data,atab.read(datafile.read()))
                datafile.close()

        # Otherwise if the filename is actually a file, just read it in.
        elif os.path.isfile(filename):
            datafile = open(filename,"r")
            self.survey_data = atab.read(datafile.read())
            datafile.close()
         
        # Print out the fields read in.
        print "The field names of the data just read are: "
        print self.survey_data.dtype.names   
        
        if 'redshift' in self.survey_data.dtype.names and 'ra' in self.survey_data.dtype.names and 'dec' in self.survey_data.dtype.names and 'mag' in self.survey_data.dtype.names:
            self.survey_data = self.survey_data[["redshift","ra","dec","mag"]]
            return

        
        # Attempt to automatically get correctly formatted names from fields
        name_list = list(self.survey_data.dtype.names)
        
        print "Attempting automatic conversion to PyLSS format"
        for i,fname in enumerate(name_list):
            if fname is 'z':
                name_list[i] = 'redshift'
            elif fname is 'Redshift':
                name_list[i] = 'redshift'
            elif fname is 'Z':
                name_list[i] = 'redshift'
            elif fname is "RA":
                name_list[i] = 'ra'
            elif fname is "Ra":
                name_list[i] = 'ra'
            elif fname is "DEC":
                name_list[i] = 'dec'
            elif fname is "Dec":
                name_list[i] = 'dec'
            elif fname is "Mag":
                name_list[i] = 'mag'
            elif fname is "MAG":
                name_list[i] = 'mag'
            elif fname is "appmag":
                name_list[i] = 'mag'
            elif fname is "app_mag":
                name_list[i] = 'mag'
            elif fname is "mag_r":
                name_list[i] = 'mag'
        
        new_tuple = tuple(name_list)
        self.survey_data.dtype.names = new_tuple    
        if 'redshift' in self.survey_data.dtype.names and 'ra' in self.survey_data.dtype.names and 'dec' in self.survey_data.dtype.names and 'mag' in self.survey_data.dtype.names:
            self.survey_data = self.survey_data[["redshift","ra","dec","mag"]]
            return
    
        # Automatic conversion must have failed, so get the user to input the names
        print "Automatic conversion has failed. Please input column names manually - "
        
        redshift_column = 'ttt'
        ra_column = 'xxx'
        dec_column = 'vvv'
        mag_column = 'ggg'
        
        while redshift_column not in name_list:
            redshift_column = input("Which column name contains redshift info?")
        while ra_column not in name_list:
            ra_column = input("Which column name contains Right Ascension info?")
        while dec_column not in name_list:
            dec_column = input("Which column name contains Declination info?")
            
        while mag_column not in name_list:
            mag_column = input("Which column name contains Apparent Magnitude info? (choose r_band if more than one, or write None if none)")
            if not mag_column:
                break


        for i,fname in enumerate(name_list):
            if fname == redshift_column:
                name_list[i] = 'redshift'
            if fname == ra_column:
                name_list[i] = 'ra'
            if fname == dec_column:
                name_list[i] = 'dec'
            if fname == mag_column:
                name_list[i] = 'mag'

        new_tuple = tuple(name_list)
        self.survey_data.dtype.names = new_tuple    
        if 'redshift' in self.survey_data.dtype.names and 'ra' in self.survey_data.dtype.names and 'dec' in self.survey_data.dtype.names and 'mag' in self.survey_data.dtype.names:
            self.survey_data = self.survey_data[["redshift","ra","dec","mag"]]
            return
        elif 'redshift' in self.survey_data.dtype.names and 'ra' in self.survey_data.dtype.names and 'dec' in self.survey_data.dtype.names:
            self.survey_data = self.survey_data[["redshift","ra","dec"]]
            return
                

        if not radians:
            self.survey_data['ra'] = self.survey_data['ra']*np.pi/180.0
            self.survey_data['dec'] = self.survey_data['dec']*np.pi/180.0
        
            
    def create_project_dir(self,filename):
        """
        Creates the directory structure into which all calculations will be placed.
        """
        # Create the Project Directory itself (check to see whether it has been made previously)
        self.project_dir = filename.partition('.')[0]+"_Project"
        print self.project_dir
        try:
            os.makedirs(self.project_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
            else:
                print "A project folder for this survey has already been created, project contents will be overwritten."
                do_continue = input("Continue? (y/n)")
                
                if do_continue is "n":
                    print "Will use the corresponding file from inside the project"
                    self.is_project_file = True
                    return
                

        # Change the focus inside the project directory
        os.chdir(self.project_dir)
        
        # Create the project file identifier
        project_file = open("_PyLSSProject_","w")
        project_file.close()
        
        # Rewrite the survey information into the project format
        survey_file = open(os.path.split(filename)[1],"w")
        survey_file.write("# File: "+os.path.split(filename)[1]+"\n")
        survey_file.write("# Created on "+time.strftime("%d/%m/%Y  %H:%M:%s")+"\n")
        #survey_file.write("# Columns: Redshift, Right Ascension, Declination, Magnitude \n")
        survey_file.write("\n")
        
        atab.write(self.survey_data,survey_file)
        
        survey_file.close()
        
    def survey_file_import(self,filename):
        """
        A simple function to import the survey data (not a subsample)
        """
        datafile = open(filename,"r")
        self.survey_data = atab.read(datafile.read())
        datafile.close()
        
        
    def create_sample_dir(self):
        """
        Creates a directory containing the subsample
        """
        
        self.sample_name = input("Please enter a name for the subsample")
        self.sample_dir = self.project_dir+'/'+self.sample_name
        try:
            os.makedirs(self.sample_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
            else:
                print "This name has already been taken. If you continue, you will overwrite the contents."
                do_continue = input("Continue? (y/n)")
                
                if do_continue is "n":
                    print "Will use the corresponding file from inside the sample"
                    self.is_sample_file = True
                    return
                
        os.chdir(self.sample_dir)
        
        # Create the sample file identifier
        sample_file = open("_PyLSSSample_","w")
        sample_file.close()
        
        # Rewrite the survey information into the sample format
        survey_file = open(self.sample_name+'.dat',"w")
        survey_file.write("# File: "+self.sample_name+".dat\n")
        survey_file.write("# Created on "+time.strftime("%d/%m/%Y  %H:%M:%s")+"\n")
        
        # Write the properties of the sample
        survey_file.write("\n")
        survey_file.write("# Sample conditions: \n")
        
        entered = False
        while not entered:
            self.get_condition = input("Would you like to specify sample conditions? ('y'/'n')")
            if self.get_condition != 'y' and self.get_condition != 'n':
                print "You entered an invalid command, please type 'y' or 'n'"
            else:
                entered = True
                
        if self.get_condition == 'y':
            try:
                z_min = input("What is the minimum redshift preferred?")
            except SyntaxError:
                z_min =  None
            try:
                z_max = input("What is the maximum redshift preferred?")
            except SyntaxError:
                z_max =  None
            try:
                ra_min = input("What is the minimum Right Ascension preferred (in radians!)?")
            except SyntaxError:
                ra_min =  None
                
            try:
                ra_max = input("What is the maximum Right Ascension preferred (in radians!)?")
            except SyntaxError:
                ra_max =  None
            
            try:
                dec_min = input("What is the minimum Declination preferred (in radians!)?")
            except SyntaxError:
                dec_min =  None
            try:
                dec_max = input("What is the maximum Declination preferred (in radians!)?")
            except SyntaxError:
                dec_max =  None
            try:
                mag_min = input("What is the minimum apparent magnitude preferred?")
            except SyntaxError:
                mag_min =  None
            try:
                mag_max = input("What is the maximum apparent magnitude preferred?")
            except SyntaxError:
                mag_max =  None
            try:
                dist_min = input("What is the minimum comoving distance preferred (in Mpc/h)?")
            except SyntaxError:
                dist_min =  None
            try:
                dist_max = input("What is the maximum comoving distance preferred (in Mpc/h)?")
            except SyntaxError:
                dist_max =  None
            try:
                lum_min = input("What is the minimum luminosity preferred (in units of log(M_sun))?")
            except SyntaxError:
                lum_min =  None
            try:
                lum_max = input("What is the maximum luminosity preferred (in units of log(M_sun))?")
            except SyntaxError:
                lum_max =  None
            
        else:
            z_min = None
            z_max = None
            ra_min = None
            ra_max = None
            dec_min = None
            dec_max = None
            mag_min = None
            mag_max = None
            dist_min = None
            dist_max = None
            lum_min = None
            lum_max = None

        if z_min:
            self.survey_data = self.survey_data[self.survey_data['redshift']>z_min]
            survey_file.write("# Minimum Redshift: "+str(z_min)+'\n')
        if z_max:
            self.survey_data = self.survey_data[self.survey_data['redshift']<z_max]
            survey_file.write("# Maximum Redshift: "+str(z_max)+'\n')
        if ra_min:
            self.survey_data = self.survey_data[self.survey_data['ra']>ra_min]
            survey_file.write("# Minimum Right Asc.: "+str(ra_min)+'\n')
        if ra_max:
            self.survey_data = self.survey_data[self.survey_data['ra']<ra_max]
            survey_file.write("# Maximum Right Asc.: "+str(ra_max)+'\n')
        if dec_min:
            self.survey_data = self.survey_data[self.survey_data['dec']>dec_min]
            survey_file.write("# Minimum Declination: "+str(dec_min)+'\n')
        if dec_max:
            self.survey_data = self.survey_data[self.survey_data['dec']<dec_max]
            survey_file.write("# Maximum Declination: "+str(dec_max))
        if mag_min:
            self.survey_data = self.survey_data[self.survey_data['mag']>mag_min]
            survey_file.write("# Minimum App. Mag.: "+str(mag_min)+'\n')
        if mag_max:
            self.survey_data = self.survey_data[self.survey_data['mag']<mag_max]
            survey_file.write("# Maximum App. Mag.: "+str(mag_max)+'\n')
        if dist_min:
            self.survey_data = self.survey_data[self.survey_data['c_dist']>dist_min]
            survey_file.write("# Minimum Com. Distance: "+str(dist_min)+'\n')
        if dist_max:
            self.survey_data = self.survey_data[self.survey_data['c_dist']<dist_max]
            survey_file.write("# Maximum Com. Distance: "+str(dist_max)+'\n')
        if lum_min:
            self.survey_data = self.survey_data[self.survey_data['lum']>lum_min]
            survey_file.write("# Minimum Luminosity: "+str(lum_min)+'\n')
        if lum_max:
            self.survey_data = self.survey_data[self.survey_data['lum']<lum_max]
            survey_file.write("# Maximum Luminosity: "+str(lum_max)+'\n')
            
        
        # Cosmological Parameters.
        survey_file.write("# ------------------------------\n")
        survey_file.write("# Cosmological Parameters of Sample\n")
        survey_file.write("# Omega_Lambda: "+str(self.cosmology["omega_lambda_0"])+'\n')
        survey_file.write("# Omega_M: "+str(self.cosmology["omega_M_0"])+'\n')
        survey_file.write("# Omega_k: "+str(self.cosmology["omega_k_0"])+'\n')
        survey_file.write("# h: "+str(self.cosmology["h"])+'\n')
        survey_file.write("# =================================\n\n")
        
        atab.write(self.survey_data,survey_file)
        
        survey_file.close()
        
    def LuminosityFromAppMag(self,c_dist):
        """
        Calculates the luminosity of apparent magnitudes based on distance. VERY SIMPLE.
        
        Units of log10(L_Sun)
        """
        
        lum = np.log10((0.0813/((9.461E15)**2))*c_dist**2) - 0.4*self.survey_data["mag"]
        
        return lum
        
    def PlotHist(self,quantity="redshift",bins=None):
        """
        Plots a Histogram of the data for any quantity
        """
        
        if not bins:
            bins = self.N/10000
            if bins < 30:
                bins = 30
        
        plt.clf()
        plt.title("Histogram of "+quantity+" values for Sample "+self.sample_name+" from "+self.survey_name)
        plt.hist(self.survey_data[quantity], bins)
        plt.xlabel(quantity)
        plt.ylabel("Number")
        plt.savefig(self.plots_dir+'/Histogram_'+quantity+'_'+str(bins)+"bins.eps")
        plt.show()
        
    def PlotPolarDot(self,distance=False):
        """
        Makes a polar plot of the data, TODO: with an optional filter condition
        
        If distance is true, plots against comoving distance rather than redshift.
        """
        
        plt.clf()
        plt.title("Polar Plot of sample "+self.sample_name+" from "+self.survey_name)
        
        if not distance:
            plt.polar(self.survey_data['ra'],self.survey_data['redshift'],'.',markersize=1)
            plt.savefig(self.plots_dir+'/PolarPlot_Redshift.eps')
        else:
            plt.polar(self.survey_data['ra'],self.survey_data['c_dist'],'.',markersize=1)
            plt.savefig(self.plots_dir+'/PolarPlot_Distance.eps')
        plt.show()    
    def PlotPolarSmooth(self,bins=None):
        """
        Makes a polar plot of the data smoothed by a filter., TODO: with an optional filter condition
        
        If distance is true, plots against comoving distance rather than redshift.
        """
        x_p= self.survey_data['c_dist']*np.cos(self.survey_data['ra'])
        y_p = self.survey_data['c_dist']*np.sin(self.survey_data['ra'])
        
        if not bins:
            bins = np.ceil(np.sqrt(self.N/10000))
            if bins < 30:
                bins = 30
        
            
            
        smooth = np.histogram2d(x_p,y_p, ((bins,bins)))[0]
    
        plt.clf()
        plt.title("Smoothed Polar Plot of sample "+self.sample_name+" from "+self.survey_name)
        plt.imshow(smooth,interpolation="Gaussian",aspect='equal')
        plt.tick_params(axis='both', labelbottom=False,labelleft=False)
        plt.savefig(self.plots_dir+'/PolarPlot_Smooth.eps')
        plt.show()
        
    def PlotSkyDot(self):
        """
        Makes a 2-D flat sky-plot of the data points.
        """
        plt.clf()
        plt.title("Sky Plot of sample "+self.sample_name+" from "+self.survey_name)
        plt.plot(self.survey_data['ra'],self.survey_data['dec'],'.',markersize=1)
        plt.xlabel("Right Ascension")
        plt.ylabel("Declination")
        plt.savefig(self.plots_dir+'/SkyPlot.eps')
        plt.show()
        
    def PlotSkySmooth(self,bins=None):
        """
        Makes a 2-D flat sky-plot of the data points, smoothed by Gaussian filter.
        """
        plt.clf()
        plt.title("Smoothed Sky Plot of sample "+self.sample_name+" from "+self.survey_name)
        
        if not bins:
            bins = np.ceil(np.sqrt(self.N/10000))
            if bins < 30:
                bins = 30
            
        smooth = np.histogram2d(self.survey_data['ra'],self.survey_data['dec'], ((bins,bins)))[0]
    
        plt.clf()
        plt.imshow(smooth,interpolation=None,aspect='equal',extent=(np.min(self.survey_data['ra']),
                    np.max(self.survey_data['ra']),np.min(self.survey_data['dec']),np.max(self.survey_data['dec'])))
        plt.tick_params(axis='both')
        plt.savefig(self.plots_dir+'/SkyPlot_Smooth.eps')
        plt.show()
        
        
    def PowerSpec_1d_fast(self,quantity="redshift",bins=None,min_wave=0.0001,max_wave=1000.0,pad=1,in_wavelength=True):
        """
        Computes the 1D power spectrum and plots it.
        
        This method uses the numpy fft algorithm to calculate the real to real DFT of any quantity in the GS data,
        then converts this to a power spectrum (merely squares it). As such, it needs what I call 'density values'
        as its input, which here are merely histogram values. It will only calculate as many DFT points as there are
        bins in the histogram, which is fine for frequency-based plots, but makes wavelength (inverse frequency) plots
        very angular. To increase the number of points you can pad the original data with zeroes on either end, but
        even this does not ensure a smooth wavelength plot. To get a very nice looking wavelength plot, use
        PowerSpec_1d_slow().
        
        
        INPUT PARAMETERS:
        quantity    [redshift]    : the quantity type used for the dft. Generally 'redshift', 'c_dist', 'ra' or 'dec'
        bins        [None]        : the number of bins to use in the histogram to get 'densities'. 
                                    If none, an 'optimal' bin number is used.
        min_wave    [0.0001]      : The minimum wavelength to plot.
        max_wave    [1000.0]      : The maximum wavelength to plot
        pad         [1]           : How much the original data length is multiplied by (zeroes padded in the extra space).
                                    The higher then padding, the longer the calculation, but the better the resolution.
        wavelength  [True]        : Whether to plot the power against wavelength. If False, plots against frequency.
        
        OUTPUT PARAMETERS:
        ps                        : array of power spectrum values
        freq                      : array of frequencies
        wvl                       : array of wavelengths
        """
        
        # First bin the quantity (higher number of bins tends towards phase calculation)
        if not bins:
            bins = self.N/1
            if bins < 30:
                bins = 30
                
        hist = np.histogram(self.survey_data[quantity],bins)
        n = pad*hist[0].size
        ps = np.abs(np.fft.rfft(hist[0],n=n)/np.sqrt(n))**2

        step = hist[1][2]-hist[1][1]
        freq = np.fft.fftfreq(n, step)
        
        wavelength = 1.0/freq[:np.floor(n/2+1)]

        plt.clf()
        plt.title("1D P.S. for "+quantity+" of sample "+self.sample_name.partition('.')[0]+" from "+self.survey_name)
        if in_wavelength:
            plt.plot(wavelength[(wavelength<max_wave)&(wavelength>min_wave)],ps[(wavelength<max_wave)&(wavelength>min_wave)],'b-')
            plt.xlabel("Wavelength")
        else:
            plt.plot(freq[(wavelength<max_wave)&(wavelength>min_wave)]*2*np.pi,ps[(wavelength<max_wave)&(wavelength>min_wave)],'b-')
            plt.xlabel("Frequency (radians)")
        plt.ylabel("Power")
        if in_wavelength:
            plt.savefig(self.power1d_dir+'/'+quantity+'_bins'+str(bins)+'_pad'+str(pad)+'_wavelength.eps')
        else:
            plt.savefig(self.power1d_dir+'/'+quantity+'_bins'+str(bins)+'_pad'+str(pad)+'_frequency.eps')
        plt.show()
        
        return ps,freq*2.0*np.pi,wavelength

    def PowerSpec_1d_slow(self,quantity="redshift",bins=None,min_wave=0.0001,max_wave=1000.0,n_waves=1000):
        """
        Computes the 1D power spectrum and plots it.
        
        This method uses a slower DFT method, rather than an FFT. However, the bonus is that the increments in
        wavelength or regular and completely customizable (producing smooth plots). 
        
        INPUT PARAMETERS
        quantity    ['redshift']    : the quantity to analyse.
        bins        [None]          : the number of bins to use in the histogram , if None, performs a 'phase'
                                        DFT, which is to say each data point counts for 1.0 at its precise location.
        min_wave    [0.0001]        :the minimum wavelength to calculate
        max_wave    [1000.0]        :the maximum wavelength ot calculate
        n_waves     [1000]          :the number of wavelengths to calculate
        
        OUTPUT PARAMETERS
        ps    : the power spectrum
        wvl   : the wavelengths
        """
        wavelengths = np.asfortranarray(np.linspace(min_wave,max_wave,n_waves)) 
        
        if bins:        
            hist,edges = np.histogram(self.survey_data[quantity],bins)
            centres = []
            for i,edge in enumerate(edges[:-1]):
                centres = centres + [(edge+edges[i+1])/2]
            centres = np.asfortranarray(centres)
            ps = dft.dft_one(np.asfortranarray(hist),centres,wavelengths)
        else:
            ps = dft.phasedft_one(np.asfortranarray(self.survey_data[quantity]),wavelengths)
        

        plt.clf()
        plt.title("1D P.S. for "+quantity+" of sample "+self.sample_name.partition('.')[0]+" from "+self.survey_name)
        plt.plot(wavelengths,ps,'b-')
        plt.xlabel("Wavelength")
        plt.ylabel("Power")
        plt.show()
        plt.savefig(self.power1d_dir+'/'+quantity+'_bins'+str(bins)+'_wavelength_smooth.eps')
        
        return ps,wavelengths
        
    def PowerSpec_2d_fast(self,quantity=("c_dist","ra"),bins=None,min_wave=(0.0001,0.0001),max_wave=(1000.0,1000.0),
                     pad=(1,1),wavelength=True):
        """
        Returns and plots the 2d power spectra of the given quantities
        
        This is the equivalent procedure to PowerSpec_1d_fast for 2-dimensions. However, in this case,
        the plotting becomes more difficult. The default plots are density images, with colours showing
        peaks. The irregular nature of the wavelength plots is pronounced in 2D however. 
        
        To produce the wavelength plots, a regular grid is created and the original power spectrum
        is interpolated onto it. This produces very striated plots however. To get nice plots in 
        wavelength, use PowerSpec_2d_slow().
        
        INPUT PARAMETERS
        quantity (tuple) [('c_dist','ra')]    :the quantities to use on (x,y) axes.
        bins (int)        [None]              : the number of bins to use in both directions
        min_wave (tuple)    [(0.0001,0.0001)] : the minimum wavelengths to use
        max_wave (tuple)    [(1000.0,1000.0)] :the maximum wavelengths to use.
        pad      (tuple)    [(1,1)]            :the factor use for padding with zeroes.
        wavelength (bool)    [True]            :whether to plot against wavelength
        
        OUTPUT PARAMETERS
        cps    : the 2D power spectrum in the range given
        cfreq1 : the frequencies in the first quantity
        cfreq2 : the frequencies in the second quantity
        wvl1   : the wavelengths in the first quantity
        wvl2   : the wavlengths in the second quantity
        
        """
        
        
        # First bin the quantity (higher number of bins tends towards phase calculation)
        if not bins:
            bins = [np.ceil(np.sqrt(self.N/100)),np.ceil(np.sqrt(self.N/100))]
            if bins[0] < 30:
                bins = [30,30]
                
        hist = np.histogram2d(self.survey_data[quantity[0]],self.survey_data[quantity[1]], bins)
        
        n1 = pad[0]*(hist[1].size -1)
        n2 = pad[1]*(hist[2].size -1)

        ps = np.abs(np.fft.rfft2(hist[0],s=(n1,n2))/np.sqrt(n1*n2))**2
        
        ps = ps[np.floor(n1/2+1):,1:]
        
        step1 = hist[1][2]-hist[1][1]
        step2 = hist[2][2]-hist[2][1]
        
        freq1 = np.fft.fftfreq(n1, step1)
        freq2 = np.fft.fftfreq(n2, step2)
        freq1 = freq1[1:np.floor(n1/2)]
        freq2 = freq2[1:np.floor(n2/2)]
        
        wavelength1 = 1.0/freq1
        wavelength2 = 1.0/freq2
        
        condition_1 = (wavelength1<max_wave[0])&(wavelength1>min_wave[0])
        condition_2 = (wavelength2<max_wave[1])&(wavelength2>min_wave[1])

        
        wvl1 = wavelength1[condition_1]
        wvl2 = wavelength2[condition_2]
        cfreq1 = 2.0*np.pi*freq1[condition_1]
        cfreq2 = 2.0*np.pi*freq2[condition_2]

        cps = ps[condition_1,:]
        cps = cps[:,condition_2]
        
        if not wavelength:
            plt.clf()
            plt.title("2D P.S. for "+quantity[0]+"and "+quantity[1]+" of sample "+self.sample_name.partition('.')[0]+" from "+self.survey_name)
            plt.imshow(cps,interpolation="gaussian",aspect=np.max(cfreq1)/np.max(cfreq2),extent=(np.min(cfreq1),np.max(cfreq1),np.min(cfreq2),np.max(cfreq2)))
            plt.xlabel(quantity[0]+" frequency")
            plt.ylabel(quantity[1]+ " frequency")
            plt.show()
            plt.savefig(self.power2d_dir+'/'+quantity[0]+'AND'+quantity[1]+'_bins'+str(bins[0])+','+str(bins[1])+'_pad'+str(pad[0])+','+str(pad[1])+'_frequency.eps')
        
        else: 
            grid_1, grid_2 = np.mgrid[np.min(wvl1):np.max(wvl1):300j,np.min(wvl2):np.max(wvl2):300j]
            points = np.zeros((wvl2.size*wvl1.size,2))
            for i,wv in enumerate(wvl2):
                points[i*wvl1.size:(i+1)*wvl1.size,0] = wvl1
            for i,wv in enumerate(wvl2):
                points[i*wvl1.size:(i+1)*wvl1.size,1] = wv
             
            print wvl2.size*wvl1.size
            print cps.size    
            grid_z0 = griddata(points, np.reshape(cps,(wvl1.size*wvl2.size)), (grid_1, grid_2), method='linear')
            
            plt.title("2D P.S. for "+quantity[0]+"and "+quantity[1]+" of sample "+self.sample_name.partition('.')[0]+" from "+self.survey_name)
            plt.imshow(grid_z0.T,interpolation="gaussian",aspect=np.max(wvl1)/np.max(wvl2),origin='lower',extent=(np.min(wvl1),np.max(wvl1),np.min(wvl2),np.max(wvl2)))
            plt.xlabel(quantity[0]+" wavelength")
            plt.ylabel(quantity[1]+ " wavelength")
            plt.show()
            plt.savefig(self.power2d_dir+'/'+quantity[0]+'AND'+quantity[1]+'_bins'+str(bins[0])+','+str(bins[1])+'_pad'+str(pad[0])+','+str(pad[1])+'_frequency.eps')

        return cps, cfreq1, cfreq2, wvl1, wvl2 
    
    def PowerSpec_2d_slow(self,quantity=("c_dist","ra"),bins=None,min_wave=(0.0001,0.0001),max_wave=(1000.0,1000.0),n_waves=(1000,1000)):
        """
        Returns and plots the 2d power spectra of the given quantities
        
        This is the equivalent procedure to PowerSpec_1d_fast for 2-dimensions. However, in this case,
        the plotting becomes more difficult. The default plots are density images, with colours showing
        peaks. The irregular nature of the wavelength plots is pronounced in 2D however. 
        
        To produce the wavelength plots, a regular grid is created and the original power spectrum
        is interpolated onto it. This produces very striated plots however. To get nice plots in 
        wavelength, use PowerSpec_2d_slow().
        
        INPUT PARAMETERS
        quantity (tuple) [('c_dist','ra')]    :the quantities to use on (x,y) axes.
        bins (int)        [None]              : the number of bins to use in both directions, if none, do phase.
        min_wave (tuple)    [(0.0001,0.0001)] : the minimum wavelengths to use
        max_wave (tuple)    [(1000.0,1000.0)] :the maximum wavelengths to use.
        n_waves  (tuple)    [(1000,1000)]     : the number of wavlengths to use.
        wavelength (bool)    [True]            :whether to plot against wavelength
        
        OUTPUT PARAMETERS
        ps    : the 2D power spectrum in the range given
        wvl1   : the wavelengths in the first quantity
        wvl2   : the wavlengths in the second quantity
        
        """
        
        wvl1 = np.asfortranarray(np.linspace(min_wave[0],max_wave[0],n_waves[0]))
        wvl2 = np.asfortranarray(np.linspace(min_wave[1],max_wave[1],n_waves[1]))
        
        # First bin the quantity (higher number of bins tends towards phase calculation)
        if bins:
            hist,x_edges,y_edges = np.histogram2d(self.survey_data[quantity[0]],self.survey_data[quantity[1]], bins)       
            x_centres = []
            for i,edge in enumerate(x_edges[:-1]):
                x_centres = x_centres + [(edge+x_edges[i+1])/2]
            y_centres = []
            for i,edge in enumerate(y_edges[:-1]):
                y_centres = y_centres + [(edge+y_edges[i+1])/2]
            x_centres = np.asfortranarray(x_centres)
            y_centres = np.asfortranarray(y_centres)
            
            ps = dft.dft_two(np.asfortranarray(hist),x_centres,y_centres,wvl1,wvl2)
        else:
            ps = dft.phasedft_two(np.asfortranarray(self.survey_data[quantity[0]]),np.asfortranarray(self.survey_data[quantity[1]]),
                                  wvl1,wvl2)

            
        plt.title("2D P.S. for "+quantity[0]+"and "+quantity[1]+" of sample "+self.sample_name.partition('.')[0]+" from "+self.survey_name)
        plt.imshow(ps.T,interpolation="gaussian",aspect=np.max(wvl1)/np.max(wvl2),origin='lower',extent=(np.min(wvl1),np.max(wvl1),np.min(wvl2),np.max(wvl2)))
        plt.xlabel(quantity[0]+" wavelength")
        plt.ylabel(quantity[1]+ " wavelength")
        plt.show()
        plt.savefig(self.power2d_dir+'/'+quantity[0]+'AND'+quantity[1]+'_bins'+str(bins[0])+','+str(bins[1])+'_frequency_smooth.eps')

        return ps, wvl1, wvl2  
    
    def PowerSpec_3d_fast(self,quantity=("c_dist","ra","dec"),bins=None,min_wave=(0.0001,0.0001,0.0001),
                          max_wave=(1000.0,1000.0,1000.0), pad=(1,1,1),wavelength=False):
        """
        Returns and plots the 2d power spectra of the given quantities
        
        At the moment this does not work correctly. This will be fixed soon.
        
        This is the equivalent procedure to PowerSpec_1d_fast for 2-dimensions. However, in this case,
        the plotting becomes more difficult. The default plots are density images, with colours showing
        peaks. The irregular nature of the wavelength plots is pronounced in 2D however. 
        
        To produce the wavelength plots, a regular grid is created and the original power spectrum
        is interpolated onto it. This produces very striated plots however. To get nice plots in 
        wavelength, use PowerSpec_2d_slow().
        
        INPUT PARAMETERS
        quantity (tuple) [('c_dist','ra')]    :the quantities to use on (x,y) axes.
        bins (int)        [None]              : the number of bins to use in both directions
        min_wave (tuple)    [(0.0001,0.0001)] : the minimum wavelengths to use
        max_wave (tuple)    [(1000.0,1000.0)] :the maximum wavelengths to use.
        pad      (tuple)    [(1,1)]            :the factor use for padding with zeroes.
        wavelength (bool)    [True]            :whether to plot against wavelength
        
        OUTPUT PARAMETERS
        cps    : the 2D power spectrum in the range given
        cfreq1 : the frequencies in the first quantity
        cfreq2 : the frequencies in the second quantity
        wvl1   : the wavelengths in the first quantity
        wvl2   : the wavlengths in the second quantity
        
        """
        
        
        # First bin the quantity (higher number of bins tends towards phase calculation)
        if not bins:
            bins = np.array([np.ceil(np.sqrt(self.N/100)),np.ceil(np.sqrt(self.N/100)),np.ceil(np.sqrt(self.N/100))])
            if bins[0] < 30:
                bins = np.array([30,30,30])
        print bins.shape
        print self.survey_data[quantity[0]]       
        hist,x_edges,y_edges,z_edges = np.histogramdd([self.survey_data[quantity[0]],self.survey_data[quantity[1]],self.survey_data[quantity[2]]], bins)
        
        n1 = pad[0]*(x_edges.size -1)
        n2 = pad[1]*(y_edges.size -1)
        n3 = pad[2]*(z_edges.size -1)

        ps = np.abs(np.fft.rfftn(hist,s=(n1,n2,n3))/np.sqrt(n1*n2*n3))**2
        
        ps = ps[np.floor(n1/2+1):,np.floor(n2/2+1):,1:]
        
        step1 = x_edges[2]-x_edges[1]
        step2 = y_edges[2]-y_edges[1]
        step3 = z_edges[2]-z_edges[1]
        
        freq1 = np.fft.fftfreq(n1, step1)
        freq2 = np.fft.fftfreq(n2, step2)
        freq3 = np.fft.fftfreq(n3, step2)
        freq1 = freq1[1:np.floor(n1/2)]
        freq2 = freq2[1:np.floor(n2/2)]
        freq3 = freq3[1:np.floor(n3/2)]
        
        wavelength1 = 1.0/freq1
        wavelength2 = 1.0/freq2
        wavelength3 = 1.0/freq3
        
        condition_1 = (wavelength1<max_wave[0])&(wavelength1>min_wave[0])
        condition_2 = (wavelength2<max_wave[1])&(wavelength2>min_wave[1])
        condition_3 = (wavelength3<max_wave[2])&(wavelength3>min_wave[2])

        
        wvl1 = wavelength1[condition_1]
        wvl2 = wavelength2[condition_2]
        wvl3 = wavelength3[condition_3]
        cfreq1 = 2.0*np.pi*freq1[condition_1]
        cfreq2 = 2.0*np.pi*freq2[condition_2]
        cfreq3 = 2.0*np.pi*freq2[condition_3]

        cps = ps[condition_1,:,:]
        cps = cps[:,condition_2,:]
        cps = cps[:,:,condition_3]
        
        return cps, cfreq1, cfreq2, cfreq3, wvl1, wvl2 ,wvl3
    
    def PowerSpec_3d_slow(self,quantity=("c_dist","ra","dec"),bins=None,min_wave=(0.0001,0.0001,0.0001),
                          max_wave=(1000.0,1000.0,1000.0),n_waves=(1000,1000,1000)):
        """
        Returns and plots the 2d power spectra of the given quantities
        
        This is the equivalent procedure to PowerSpec_1d_fast for 2-dimensions. However, in this case,
        the plotting becomes more difficult. The default plots are density images, with colours showing
        peaks. The irregular nature of the wavelength plots is pronounced in 2D however. 
        
        To produce the wavelength plots, a regular grid is created and the original power spectrum
        is interpolated onto it. This produces very striated plots however. To get nice plots in 
        wavelength, use PowerSpec_2d_slow().
        
        INPUT PARAMETERS
        quantity (tuple) [('c_dist','ra')]    :the quantities to use on (x,y) axes.
        bins (int)        [None]              : the number of bins to use in both directions, if none, do phase.
        min_wave (tuple)    [(0.0001,0.0001)] : the minimum wavelengths to use
    
            
        """
        
        wvl1 = np.asfortranarray(np.linspace(min_wave[0],max_wave[0],n_waves[0]))
        wvl2 = np.asfortranarray(np.linspace(min_wave[1],max_wave[1],n_waves[1]))
        wvl3 = np.asfortranarray(np.linspace(min_wave[2],max_wave[2],n_waves[2]))
        
        # First bin the quantity (higher number of bins tends towards phase calculation)
        if bins:
            hist,x_edges,y_edges,z_edges = np.histogramdd(self.survey_data[quantity[0]],self.survey_data[quantity[1]],self.survey_data[quantity[2]], bins)       
            x_centres = []
            for i,edge in enumerate(x_edges[:-1]):
                x_centres = x_centres + [(edge+x_edges[i+1])/2]
            y_centres = []
            for i,edge in enumerate(y_edges[:-1]):
                y_centres = y_centres + [(edge+y_edges[i+1])/2]
            z_centres = []
            for i,edge in enumerate(z_edges[:-1]):
                z_centres = z_centres + [(edge+z_edges[i+1])/2]
            x_centres = np.asfortranarray(x_centres)
            y_centres = np.asfortranarray(y_centres)
            z_centres = np.asfortranarray(z_centres)
            
            ps = dft.dft_three(np.asfortranarray(hist),x_centres,y_centres,z_centres,wvl1,wvl2,wvl3)
        else:
            ps = dft.phasedft_two(np.asfortranarray(self.survey_data[quantity[0]]),np.asfortranarray(self.survey_data[quantity[1]]),
                                  np.asfortranarray(self.survey_data[quantity[2]]),wvl1,wvl2,wvl3)


        return ps, wvl1, wvl2,wvl3      
   
    #def filtering(self,quantity=['c_dist'],kmin=[35.0],kmax=[45.0],n_k=[10],min_phase=0.5):
    #    """
    #    Filters the galaxies based on whether they are in phase with a selected Fourier peak.
    #    
    #    SHOULD BE ABLE TO DO ALL ANALYSIS WITH THE FILTERED GALAXIES BASED ON REMAINS CONDITION
    #    IS IT BETTER TO DO THIS BY WRITING OUT THE 'FILTERED SURVEY'?
    #    """
    #   
    #    from PyGS.fort.filter import filters
    #    if type(quantity) == type([1,2,3]):
    #        dim = len(quantity)
    #        if dim > 3:
    #            sys.exit("Too many quantities in filter (must be 3 or less)")
    #        if (len(kmin) != dim) or (len(kmax) != dim) or (len(n_k) !=dim):
    #            sys.exit("All input values must have same dimension")
    #                                                        
    #    if dim == 1:
    #        filters.filter_1d(self.survey_data[quantity[0]],kmin[0],kmax[0],n_k[0],min_phase)
    #        self.remains = np.asfortranarray(filters.remains)
    #
    #        filter_file = open(self.filtered_dir+'/1D_filter_'+kmin[0]+'to'+kmax[0])    
    #        filter_file.write("# Fourier filter with following characteristics:\n")
    #        filter_file.write("# kmin: "+kmin[0]+'\n')
    #        filter_file.write("# kmax: "+ kmax[0]+'\n')
    #        filter_file.write("# n_k: "+ n_k[0]+'\n')
    #        filter_file.write("#======================")
    #    
    #        atab.write(self.survey_data[self.remains],filter_file)
    #    
    #        filter_file.close()
    #        
    #    if dim == 3:
    #        filters.filter_3d(self.survey_data[quantity[0],quantity[1],quantity[2]],kmin,kmax,n_k,min_phase)
    #        self.remains = np.asfortranarray(filters.remains)
    #    
    #    filter_file = open(self.filtered_dir+'/'+dim+'D_filter_'+kmin[0])    
    #    survey_file.write("# Omega_k: "+str(self.cosmology["omega_k_0"])+'\n')
    #    survey_file.write("# h: "+str(self.cosmology["h"])+'\n')
    #    survey_file.write("# =================================\n\n")
        
        atab.write(self.survey_data,survey_file)
        
        survey_file.close()

























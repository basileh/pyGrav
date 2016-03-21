#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 06 17:51:43 2014

@author: Basile HECTOR

##########################
        pyGrav 
##########################

Module which contains the main classes of the program:

Base class is an objct of type ChannelList, which basically contains channel
lists such as found in CG5 ascii output files (grav, tilts, temp, sd, time...)

Derived classes follow a logical hierarchy, where each 'subclass' are 
instanciated as an element of a dictionnary from the parent class:
class Campaign
has a dictionnary of 'Survey'-types objects
    class Survey
    has a dictionnary of 'Loop'-types objects
        class Loop
        has a dictionnary of 'Station'-types objects
            class Station
        
Each of these objects is derived from a ChannelList object. The instance of
the base class 'Campaign' therefore contains all the data set (and several 
times)
            
These class also have specific properties and populating, writing, handling and
processing functions, called from the main program.

"""

#import os, sys, glob, numpy
from time import time
import glob, os
from copy import *
from PyQt4 import QtGui,QtCore
from datetime import *
from synthetic_tides import *
from scipy import sqrt
from scipy.interpolate import interp1d
from collections import OrderedDict
import numpy as np

from matplotlib.dates import date2num,num2date

###############################################################################
class ChannelList(object):     
    """ChannelList-type object
    
    Properties:
    
    list objects (time series):
    line             gravity line
    station          station name
    alt              altitude (m)
    grav             gravity value (mgal)
    sd               Standard Deviation (mgal)
    tiltx            tilt x value (arcsec)
    tilty            tilt y value (arcsec)
    temp             temperature variations (k)
    etc              cg5 earth tide corrections (mgal)
    dur              reading duration (sec)
    rej              number of rejected values
    t                datetime object of central second of the reading
    corr_g           corrected g channel
    keepdata         1 if value is to be kept, 0 if not
    
    status objects:
    keepitem        1 if item is to be kept (for inversion for instance)
                    0 otherwise
                    
    functions:
    - extract_subset
    - read_c_file
    - read_modif_c_file
    - populateFromSubDictionnaries
                
    """
    line=[]
    station=[] 
    alt=[]
    grav=[]  
    sd=[]
    tiltx=[]
    tilty=[]
    temp=[]
    etc=[]
    dur=[]
    rej=[]
    t=[]
    corr_g=[]
    keepdata=[]
    keepitem=int
  
    def __init__(self):
        """
        """
        self.line=[]
        self.station=[]
        self.alt=[]
        self.grav=[]
        self.sd=[]
        self.tiltx=[]
        self.tilty=[]
        self.temp=[]
        self.etc=[]
        self.dur=[]
        self.rej=[]
        self.t=[]
        self.corr_g=[]
        self.keepdata=[]
        self.keepitem=1
        
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print 'length of field t: %d'%(len(self.t))       
        print 'first station: %s'%(self.station[0])   
        print 'last station: %s'%(self.station[len(self.station)-1]) 
        print 'station availables: %s'%(set(self.station))
        return "you're welcome"
        
    def extract_subset(self,start_date,end_date):
        """Function for extracting shorter time series        
        Provide a starting and ending dates and the function will return an
        other ChannelList object with a subset of the original        
        """
        templist=ChannelList()
        tsub=[t for t in self.t if date2num(t)>=date2num(start_date) and date2num(t)<=date2num(end_date)]        
        indexes=[self.t.index(t) for t in tsub]
        indexes.sort()
        templist.line=[self.line[ind] for ind in indexes]
        templist.station=[self.station[ind] for ind in indexes]
        templist.alt=[self.alt[ind] for ind in indexes]
        templist.grav=[self.grav[ind] for ind in indexes]
        templist.sd=[self.sd[ind] for ind in indexes]
        templist.tiltx=[self.tiltx[ind] for ind in indexes]
        templist.tilty=[self.tilty[ind] for ind in indexes]
        templist.temp=[self.temp[ind] for ind in indexes]
        templist.etc=[self.etc[ind] for ind in indexes]
        templist.dur=[self.dur[ind] for ind in indexes]
        templist.rej=[self.rej[ind] for ind in indexes]
        templist.t=[self.t[ind] for ind in indexes]
        templist.keepdata=[self.keepdata[ind] for ind in indexes]
        return templist
     
    def read_c_file(self,filename):
        """
        read a generic c file and populate channels
        """                 
        try:
            #essaye d'ouvrir le fichier
            fh = open(filename, 'r')
            i=0

            PBAR = ProgressBar(total=len([1 for line in  open(filename, 'r')]),textmess='Load c data')   
            
            print "number of lines: %d"%len([1 for line in  open(filename, 'r')])
            PBAR.show()

            for line in fh:    
                PBAR.progressbar.setValue(i)
                i+=1
                # Clean line
                line = line.strip()
                # Skip blank and comment lines
                if (not line) or (line[0] == '/') or (line[0] == 'L'): continue

                vals=line.split()
                # fill object properties:
                self.line.append(0)
                self.station.append(int(vals[0]))
                self.alt.append(0)
                self.grav.append(float(vals[15]))
                self.sd.append(float(vals[2]))
                self.dur.append(int(vals[3]))
                self.rej.append(int(vals[4]))
                self.tiltx.append(float(vals[5]))
                self.tilty.append(float(vals[6]))
                self.temp.append(float(vals[7]))
                self.etc.append(float(vals[8]))                
                self.t.append(datetime(int(vals[11][4:6])+2000,int(vals[11][2:4]),\
                int(vals[11][0:2]),int(vals[12][0:2]),int(vals[12][2:4]),\
                int(vals[12][4:6])))             
                self.keepdata.append(1)     
                                                                                    
        except IOError:
            #si ça ne marche pas, affiche ce message et continue le prog
            print 'No file : %s' %(filename)            
        except ValueError:
            print 'pb at line %d : check raw data file'%(i)
        except IndexError:
            print 'pb at line %d : check raw data file: possibly last line?'%(i)         
            
    def read_modif_c_file(self,filename):
        """
        read a slighlty modified generic CGxTool 'c' file and populate channels
        One more column is added to the 'c' file containing 1 and 0 whether 
        data line is kept or not (for drift adjustment purpose)
        
        To be modified: only available for dates > 2000
        """         
        
        try:
            #essaye d'ouvrir le fichier
            fh = open(filename, 'r')
            i=0
            #PBAR = ProgressBar(total=len([1 for line in  open(filename, 'r')]),textmess='Load c data')               
            print "number of lines: %d"%len([1 for line in  open(filename, 'r')])
            #PBAR.show()
            for line in fh:    
                #PBAR.progressbar.setValue(i)
                i+=1
                # Clean line
                line = line.strip()
                # Skip blank and comment lines
                if (not line) or (line[0] == '/') or (line[0] == 'L'): continue

                vals=line.split()
                # fill object properties:
                self.line.append(0)
                self.station.append(int(vals[0]))
                self.alt.append(0)
                #self.grav.append(float(vals[15]))
                self.grav.append(float(vals[1]))                
                self.sd.append(float(vals[2]))
                self.dur.append(int(vals[3]))
                self.rej.append(int(vals[4]))
                self.tiltx.append(float(vals[5]))
                self.tilty.append(float(vals[6]))
                self.temp.append(float(vals[7]))
                self.etc.append(float(vals[8]))                
                self.t.append(datetime(int(vals[11][4:6])+2000,int(vals[11][2:4]),\
                int(vals[11][0:2]),int(vals[12][0:2]),int(vals[12][2:4]),\
                int(vals[12][4:6])))             
                self.keepdata.append(int(vals[16]))
                                                                                    
        except IOError:
            #si ça ne marche pas, affiche ce message et continue le prog
            print 'No file : %s' %(filename)            
        except ValueError:
            print 'pb at line %d : check raw data file'%(i)
        except IndexError:
            print 'pb at line %d : check raw data file: possibly last line?'%(i)            
                        
    def populateFromSubDictionnaries(self,subdict):
        """
        When subdictionnaries (stations, loops, surveys) are filled (their 
        time series) prior to the upstream campaign structure, this should be
        filled, because its time series are used for some functions (such as
        tide and atmosphere corrections)
        
        Arguments:
        - subdict:          subdictionaries to fill from
        """                
        for keyssubobj,subobj in sorted(subdict.iteritems(), key=lambda x: x[1].t[1]):    
            for i in range(len(subobj.t)):              
                self.line.append(subobj.line[i])
                self.station.append(subobj.station[i])
                self.alt.append(subobj.alt[i])
                self.grav.append(subobj.grav[i])
                self.sd.append(subobj.sd[i])
                self.tiltx.append(subobj.tiltx[i])
                self.tilty.append(subobj.tilty[i])
                self.temp.append(subobj.temp[i])
                self.etc.append(subobj.etc[i])
                self.dur.append(subobj.dur[i])
                self.rej.append(subobj.rej[i])
                self.t.append(subobj.t[i])
                if subobj.corr_g:
                    self.corr_g.append(subobj.corr_g[i])            
###############################################################################
class Campaign(ChannelList):
    """Campaign-type object

    Derived (inheritance) from a ChannelList-type object 
    all properties and functions from ChannelList apply
    
    Properties:

    survey_dic		dictionary of survey objects. Keys are start dates as 
                       given by datetime.toordinal()


    Functions:
    - readRawDataFile
    - populateSurveyDic
    - readModifCDataFilesAndPopulateSurveyDic
    - saveProcessedDoubledifferencesData
    - writeStationLocationFile
    - saveProcessedData
    - corrHeightBouguer
    - driftInversion
    """

    survey_dic = {}
  
    def __init__(self):
        ChannelList.__init__(self)
        """
        Initializes to an empty dictionary.
    
        """
        self.survey_dic = {}
        self.name="Campaign"
    
    def readRawDataFile(self,filename):  
        """ read a raw ascii data text file extracted from CG5
     
        sometimes bad values are written by the CG5 soft, if such case happen,
        an error is raised and the user is asked for checking the data file 
        manually
        """    
   
        try:
            #essaye d'ouvrir le fichier
            fh = open(filename, 'r')
            i=0
            PBAR = ProgressBar(total=len([1 for line in  open(filename, 'r')]),textmess='Load raw data')   
            print "number of lines: %d"%len([1 for line in  open(filename, 'r')])
            PBAR.show()
            for line in fh:    
                PBAR.progressbar.setValue(i)
                i+=1
                # Clean line
                line = line.strip()
                # Skip blank and comment lines
                if (not line) or (line[0] == '/') or (line[0] == 'L'): continue
        	     #parse string line first with respect to '/' caracters (used in the date format), 
        	     #then with ':' (used for the time display), eventually with the classic ' '
                vals_temp1=line.split('/')
                vals_temp2=vals_temp1[0].split(':')
                vals_temp3=vals_temp2[0].split()                
                vals_temp4=vals_temp2[2].split()

                # fill object properties:
                self.line.append(float(vals_temp3[0]))
                self.station.append(float(vals_temp3[1]))
                self.alt.append(float(vals_temp3[2]))
                self.grav.append(float(vals_temp3[3]))
                self.sd.append(float(vals_temp3[4]))
                self.tiltx.append(float(vals_temp3[5]))
                self.tilty.append(float(vals_temp3[6]))
                self.temp.append(float(vals_temp3[7]))
                self.etc.append(float(vals_temp3[8]))
                self.dur.append(int(vals_temp3[9]))
                self.rej.append(int(vals_temp3[10]))
                self.t.append(datetime(int(vals_temp4[3]),int(vals_temp1[1]),\
                int(vals_temp1[2]),int(vals_temp3[11]),int(vals_temp2[1]),\
                int(vals_temp4[0])))             
                self.keepdata.append(1)                                                                                         
        except IOError:
            #si ça ne marche pas, affiche ce message et continue le prog
            print 'No file : %s' %(filename)            
        except ValueError:
            print 'pb at line %d : check raw data file'%(i)
        except IndexError:
            print 'pb at line %d : check raw data file: possibly last line?'%(i)  
            
    def populateSurveyDic(self,options,base_station,start_end_dates=None,time_threshold=0,base_cycling_station=0):
        """populate survey dictionary following user-defined rules

        fill as many survey objects as wanted, each survey is then made of
        different loops, defined following user's rules. Each loop is made of
        different stations
        
        Arguments:
        options:                1: provide a list of starting and ending dates
                                   for survey identification
                                2: automatic survey identification
                                
        base station:           base station name (number)                        

        start_end_dates:        list of tuples with starting and ending dates 
                                if options = 1

        time_threshold:         time interval threshold that separate two 
                                surveys if options = 2
 
        base_cycling_station:   1 if base station = cycling station, 0 otherwise   

        if options = 2, and base station = cycling station, a number of data to
        keep after and before the entrance in a loop is defined (currently 15)
        and could need to be passed as an argument.
        """            
    
        if options==1:
            for survey_tuple in start_end_dates:
                # store into survey dictionary: keys are start dates of survey
                key_id=survey_tuple[0].toordinal()
                temp_camp=deepcopy(self)
                temp_surv=temp_camp.extract_subset(survey_tuple[0],survey_tuple[1])  
                # create a new Survey-type object:                                
                temp_surv2=Survey(temp_surv,str(datetime.fromordinal(key_id).date()))
                temp_surv2.populateLoopDic(1,base_station)
                self.survey_dic[key_id]=temp_surv2
                
                
        if options==2:
            # non zero difference indices
            diff_sta=np.diff(self.station)            
            ind=[i for i in range(len(self.station)-1) if diff_sta[i] != 0]
            ind.append(len(self.station)-1)
            # diff t (object of type datetime.timedelta)
            t_temp=[self.t[i] for i in ind]
            # add first and last value
            t_temp.insert(0,self.t[0])
            t_temp.append(self.t[len(self.t)-1])
            diff_t=np.diff(t_temp)
            print diff_t
            # define number of data to retain before the start and after the 
            # end of a survey (base station): could be passed as an argument
            n_to_keep=15            
            
            #initiate test value (1=survey found, 0=no survey)
            test=0
            k=0
            for i in ind:
                print i
                print diff_t[k].days*24
                #print diff_t[i].seconds/60/60
                diff_hr=diff_t[k].days*24 + diff_t[k].seconds/60/60
                print diff_hr                
                if diff_hr >= time_threshold :
                    print 'ok'
                    # if we are quitting a survey
                    if test==1:
                        if base_cycling_station==1:
                            if ind[k-1]+n_to_keep <= len(self.t)-1:
                                tstop=self.t[ind[k-1]+n_to_keep]
                            else :
                                tstop=self.t[len(self.t)-1]
                                
                        elif base_cycling_station==0:
                            tstop=self.t[ind[k-1]]
                        key_id=tstart.toordinal()  
                        print key_id
                        print str(datetime.fromordinal(key_id).date())
                        print tstart
                        print tstop                        
                        temp_camp=deepcopy(self)
                        temp_surv=temp_camp.extract_subset(tstart,tstop)  
                        # create a new Survey-type object:                                
                        temp_surv2=Survey(temp_surv,str(datetime.fromordinal(key_id).date()))
                        temp_surv2.populateLoopDic(1,base_station)
                        test=0
                        self.survey_dic[key_id]=temp_surv2                                
                    # if we have entered a survey                
                    if test==0:
                        if base_cycling_station==1:
                            tstart=self.t[i-n_to_keep]
                        elif base_cycling_station==0:
                            tstart=self.t[i]
                        test=1                        
                k=k+1     

    def readModifCDataFilesAndPopulateSurveyDic(self,filename):
        """
        Populate a survey dictionnary for all survey identified in the 
        hierarchical file. The format of the input file is the following:
        Directory 'root directory where all data is stored'
        Survey: 'surveyname' nloops: 'number of loops' directory: 'survey directory'
        Loop: 'loop name' filename: 'loop modified c file name'
        Loop: 'loop name' filename: 'loop modified c file name'
        Loop: 'loop name' filename: 'loop modified c file name'
        ...
        Survey: 'surveyname' nloops: 'number of loops' directory: 'survey directory'
        Loop: 'loop name' filename: 'loop modified c file name'
        ...

        where survey name is survey date as string yyyy-mm-dd, and number of 
        loops is the number of loops that follow

        """        
        # open hierarchical file
        f=open(filename,'r')  
        for line in f:
            # Clean line
            line = line.strip() 
            vals=line.split()
            
            if vals[0]=='Directory':
                rootdir=vals[1]
            
            if vals[0]=='Survey:':
                #instanciate a survey object
                survtemp=Survey(ChannelList(),vals[1])
                surveydate=datetime.strptime(vals[1],'%Y-%m-%d')
                survtemp.keepitem=1
                keysurv=int(surveydate.toordinal())
                print "process survey: %d (%s)"%(keysurv,vals[1])
                #keysurv=str(surveydate.toordinal())
                #print "process survey: %s (%s)"%(keysurv,vals[1])                
                nloops=int(vals[3])
                survdir=vals[5]
                k=0
               
            if vals[0]=='Loop:':
                k=k+1
                channeltemp=ChannelList()
                loopdir=rootdir+os.sep+survdir+os.sep+vals[3]
                channeltemp.read_modif_c_file(loopdir)
                temp_loop=Loop(channeltemp,vals[1])
                temp_loop.keepitem=1
                temp_loop.populateStationDic()
                survtemp.loop_dic[int(vals[1])]=temp_loop 
                if k==nloops:
                    #the current survey is filled, store it in the self
                    self.survey_dic[keysurv]=survtemp   
        

    def calculateDoubleDifferences(self,referencesurvey,option):
        """
        Compute double difference from a reference survey
        the reference survey is passed as an argument
        std is calculated as sqrt(std_currenstation² + std_referencestation²)
        
        the function returns a Campaign object, in which each output_dic of 
        each Survey contains tuples of gravity double differences (see 
        output_dic format under Survey objects)
        
        Parameters:
        - reference_survey:                 explicit
        - option:                           1 = classic double differences, 
                                            2 = dd from network mean      
        """
        
        #instanciate a Campaign object to store double differences
        DD=Campaign()
        
        if option==2:
            #only for option 2:
            mean_ref_surv=np.mean([sta[1] for staid,sta in referencesurvey.output_dic.iteritems()])        
        
        for keysurv,surv in self.survey_dic.iteritems():
            if surv.keepitem==1:
                #empty survey
                survtmp=Survey(ChannelList(),surv.name)             
                DD.survey_dic[keysurv]=survtmp
        
                if option==1: #classic double differences
                    for staid,sta in referencesurvey.output_dic.iteritems():
                        try:
                            currdata=self.survey_dic[keysurv].output_dic[staid]                    
                            if currdata:
                                grav=currdata[1]-sta[1]                                
                                std=sqrt(currdata[2]*currdata[2]+sta[2]*sta[2])                            
                                DD.survey_dic[keysurv].output_dic[staid]=(sta[0],grav,std)
                        except KeyError:
                            pass
                elif option==2: #double difference from network mean:
                    mean_curr_surv=[]
                    for staid,sta in referencesurvey.output_dic.iteritems():
                        try:
                            currdata=surv.output_dic[staid]                    
                            if currdata:
                                mean_curr_surv.append(currdata[1])
                        except KeyError:
                            # means station is not available: do not append vector for mean calculation
                            pass
                    mean_curr_surv=np.mean(mean_curr_surv)         
                    
                    for staid,sta in referencesurvey.output_dic.iteritems():
                        try:
                            currdata=surv.output_dic[staid]                    
                            if currdata:
                                grav=currdata[1]-mean_curr_surv-(sta[1]-mean_ref_surv)
                                std=sqrt(currdata[2]*currdata[2]+sta[2]*sta[2])
                                DD.survey_dic[keysurv].output_dic[staid]=(sta[0],grav,std)
                        except KeyError:
                            # do nothing
                            pass
                        
        return DD                        
        
    def saveProcessedDoubledifferencesData(self,filenameGrav,filenameSTD,DD,referencesurvey):
        """
        Save Double difference files:
        write a doubledifferencesGrav.dat file and a
        doubledifferencesSTD.dat file containing a 1 line header, 
        a survey date column, and grav or std values for each station (µgal)
        
        Another (transposed) output file is also proposed
        
        Parameters:
        - filenameGrav:              filename for gravity double differences
        - filenameSTD:               filename for STD double differences
        - DD:                        a Campaign object, in which each 
                                     output_dic of each Survey contains tuples 
                                     of gravity double differences (see 
                                     output_dic format under Survey objects)                                   
        """
        
        
        ###########
        #FORMAT 1:#
        ###########
        survids=[survid for survid,surv in DD.survey_dic.iteritems()]
        survids=np.sort(survids)
        file1=open(filenameGrav+"1.dat",'w')
        file2=open(filenameSTD+"1.dat",'w')
        # print header        
        file1.write("date        ")
        file2.write("date        ")               
        for staid,sta in referencesurvey.output_dic.iteritems():
            file1.write("%6d "%sta[0])
            file2.write("%6d "%sta[0])
        file1.write("\n")    
        file2.write("\n")            
        # print data    
        for surv in survids:
            # print date:
            datetmp=str(datetime.fromordinal(surv))
            file1.write("%11s "%datetmp[0:10])
            file2.write("%11s "%datetmp[0:10])
            # print data
            for staid,sta in referencesurvey.output_dic.iteritems():
                try:
                    currdata=DD.survey_dic[surv].output_dic[staid]
                    if currdata:
                        grav=currdata[1]*1000
                        std=currdata[2]*1000
                        file1.write("%6.1f "%grav)
                        file2.write("%6.1f "%std)   
                except KeyError:
                    # means station is not available
                    file1.write("   NaN ")
                    file2.write("   NaN ")   
            file1.write("\n")    
            file2.write("\n")                      
        file1.close()
        file2.close()
        
        ###########
        #FORMAT 2:#
        ###########
        file1=open(filenameGrav+"2.dat",'w')
        file2=open(filenameSTD+"2.dat",'w')
        file1.write("     ")
        file2.write("     ")        
        for surv in survids:
            datetmp=str(datetime.fromordinal(surv))
            file1.write(" %11s"%datetmp[0:10])
            file2.write(" %11s"%datetmp[0:10])
        file1.write("\n")    
        file2.write("\n")      
        for staid,sta in referencesurvey.output_dic.iteritems():
            file1.write("%4d "%sta[0])
            file2.write("%4d "%sta[0])
            for surv in survids:
                try:
                    currdata=DD.survey_dic[surv].output_dic[staid]          
                    if currdata:
                        grav=currdata[1]*1000
                        std=currdata[2]*1000
                        file1.write(" %11.1f"%grav)
                        file2.write(" %11.1f"%std)   
                    else:
                        file1.write("         NaN")
                        file2.write("         NaN")   
                except KeyError:
                    # means station is not available
                    file1.write("         NaN")
                    file2.write("         NaN") 
            file1.write("\n")    
            file2.write("\n")                      
        file1.close()
        file2.close()        
        
    def writeStationLocationFile(self,filename):
        """
        write a standard station location file for mcgravi, with (0,0)
        coordinates for all stations
        """
        file=open(filename,'w')
        stations=[]

        for keysurvey,survey in self.survey_dic.iteritems():            
            for keyloop,loop in self.survey_dic[keysurvey].loop_dic.iteritems():  
                for keysta,station in self.survey_dic[keysurvey].loop_dic[keyloop].station_dic.iteritems():                  
                    stations.append(station.station[0])
        
        stations=set(stations)
        #set([sta for sta in self.station])
        for sta in stations:
            file.write("%d 0 0\n"%(sta))
        file.close()      


    def saveProcessedData(self,output_root_dir):
        """
        only save modified c files
        
        Can be modified, for other format definitions, or for instance saving
        all data with an extra column with 1 or 0 depending on keepdata status
        
        A generic file containing all the data set hierarchy (survey, loops) is
        written, with the survey folder names and loops modified 'c' file names
        This file is further used to re-load the whole data set, which keeps
        the keepdata information and allow to retrieve the data selection 
        process left at some point.
        """
        Mainfile=open(output_root_dir+os.sep+'gravity_data_hierarchy.txt','w')
        
        datafordriftadj=deepcopy(self)        
        if not os.path.exists(output_root_dir):
            os.makedirs(output_root_dir)
            
        Mainfile.write("Directory %s\n"%output_root_dir)  
        for survid,surv in datafordriftadj.survey_dic.iteritems():
            if surv.keepitem==1:
                survdir=output_root_dir+os.sep+surv.name
                if not os.path.exists(survdir):
                    os.makedirs(survdir)
                #get number of loops
                nloops=0
                for loopid,loop in surv.loop_dic.iteritems():
                    if loop.keepitem==1:
                        nloops=nloops+1
                # save survey characteristics        
                Mainfile.write("Survey: %s nloops: %d directory: %s\n"%(surv.name,nloops,surv.name))     
                j=0
                file_list=[]
                for loopid,loop in surv.loop_dic.iteritems():
                    if loop.keepitem==1:
                        j=j+1
                        #get year and julian day:
                        stakeys=[key for key,sta in loop.station_dic.iteritems() if sta.keepitem==1]
                        datefirststa=loop.station_dic[stakeys[1]].t[0].timetuple()
                        yr=str(datefirststa.tm_year)
                        jday="%03d"%datefirststa.tm_yday
                        #filenameC=survdir+os.sep+"fn"+str(j)+"11c"+yr[2:4]+"."+jday+".txt"
                        filenameC="fn"+str(j)+"11c"+yr[2:4]+"."+jday+".txt"
                        file_list.append("fn"+str(j)+"11c"+yr[2:4]+"."+jday)
                        # getbase station: temporary formulation, should add a checkbox 
                        # in the tree object to select it manually for each loop
                        basestationkeys=[stakeys[i] for i in range(len(stakeys)) if stakeys[i][1]==2]
                        loop.writeModifCFile(survdir+os.sep+filenameC,basestationkeys[0])
                        Mainfile.write("Loop: %s filename: %s\n"%(loopid,filenameC))     
        Mainfile.close()
        
    def corrHeightBouguer(self,fname):
        """
        function for correcting height on all simple difference data
        """
        file=open(fname,'r')        
        lines=file.readlines()
                
        for survid,surv in self.survey_dic.iteritems():
            surv.height_bouguer_corr='ok'
            for sta_id,sta in surv.output_dic.iteritems():
                #initiate to NaNs
                sta=sta+(NaN,NaN,)
                for line in lines:
                    line = line.strip() 
                    vals=line.split()
                    if int(vals[0]==sta_id):
                        gcorr_height=0.3086*float(vals[1])+sta[1]
                        gcorr_Boug=-0.0419*2.65*float(vals[1])+gcorr_height                        
                        sta[3]=gcorr_height
                        sta[4]=gcorr_Boug
                        
        self.survey_dic[survid]=surv
        
    def tideCorrectionAgnew(self,lat,lon,alt):
        """
        
        code section still to be optimized
        """
        
        if any(self.t):
            #get tides and round to µgal level
            tides=np.round(np.array([earth_tide(lat,lon,t) for t in self.t])*1000)/1000 
            self.corr_g=[self.grav[i]-self.etc[i] +tides[i] for i in range(len(self.t))]
            self.grav=self.corr_g       
            self.etc=[tides[i] for i in range(len(self.t))]
        PBAR1 = ProgressBar(total=len(self.survey_dic.keys()),textmess='surveys')   
        PBAR1.show()
        i=1
        for keysurv,surv in self.survey_dic.iteritems():
            PBAR1.progressbar.setValue(i)
            i=i+1   
            if any(self.survey_dic[keysurv].t):
                #get tides and round to µgal level
                tides=np.round(np.array([earth_tide(lat,lon,t) for t in self.survey_dic[keysurv].t])*1000)/1000
                self.survey_dic[keysurv].corr_g=[self.survey_dic[keysurv].grav[i]-self.survey_dic[keysurv].etc[i]+tides[i] for i in range(len(self.survey_dic[keysurv].t))]           
                self.survey_dic[keysurv].grav=self.survey_dic[keysurv].corr_g
                self.survey_dic[keysurv].etc=[tides[i] for i in range(len(self.survey_dic[keysurv].t))]           
            PBAR2 = ProgressBar(total=len(self.survey_dic[keysurv].loop_dic.keys()),textmess='loops')   
            PBAR2.show()
            j=1
            for keyloop,loop in self.survey_dic[keysurv].loop_dic.iteritems():
                PBAR2.progressbar.setValue(j)
                j=j+1
                if any(self.survey_dic[keysurv].loop_dic[keyloop].t):
                    #get tides and round to µgal level
                    tides=np.round(np.array([earth_tide(lat,lon,t) for t in self.survey_dic[keysurv].loop_dic[keyloop].t])*1000)/1000
                    self.survey_dic[keysurv].loop_dic[keyloop].corr_g=[self.survey_dic[keysurv].loop_dic[keyloop].grav[i]-self.survey_dic[keysurv].loop_dic[keyloop].etc[i]+tides[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].t))] 
                    self.survey_dic[keysurv].loop_dic[keyloop].grav=self.survey_dic[keysurv].loop_dic[keyloop].corr_g
                    self.survey_dic[keysurv].loop_dic[keyloop].etc=[tides[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].t))] 
                PBAR3 = ProgressBar(total=len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic.keys()),textmess='stations')   
                PBAR3.show()
                k=1
                for keysta,sta in self.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
                    PBAR3.progressbar.setValue(k)
                    k=k+1
                    #get tides and round to µgal level
                    tides=np.round(np.array([earth_tide(lat,lon,t) for t in self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t])*1000)/1000
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g=[self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav[i]-self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].etc[i]+tides[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t))]
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav=self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g                                           
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].etc=[tides[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t))]                    
                    
    def tideCorrectionPredict(self,Tides):    
        """
        Apply tide corrections: populate/update the gcorr field of all gravity
        objects in self.campaigndata
        update all grav fields
        update all etc fields
        IMPORTANT: sign of the correction is important:
        - native etc from CG5 should be removed from .grav, and then predict etc
        should be removed too! Therefore the new etc should be '-' as well
        This is because native ETC CG5 value is not a synthetic tide but a 
        correction, hence '- synthetic tide'.
        Therefore, when removing etc to grav (grav-etc), this means that
        the synthetic tide is added back and the total signal is obtained.
        For further removing of another synthetic tide, one should therefore do
        new_corr_g = grav-old_etc-etc2
        new_etc=-etc2
        """
        #if data is loaded from already processed data, self.campaigndata time
        #series are not populated, and correction should not be applied at such
        #levels
        
        # CG5 data is acquired with a µgal resolution. therefore, apply 
        #corrections at the µgal level:
        Tides.d=np.round(np.array(Tides.d)*1000)/1000
        
        if any(self.t):
            Tides.interpolateOnGivenTimes(self.t)   
            self.corr_g=[self.grav[i]-self.etc[i] -Tides.d[i] for i in range(len(self.t))]
            self.grav=self.corr_g       
            self.etc=[-Tides.d[i] for i in range(len(self.t))]
        PBAR1 = ProgressBar(total=len(self.survey_dic.keys()),textmess='surveys')   
        PBAR1.show()
        i=1
        for keysurv,surv in self.survey_dic.iteritems():
            PBAR1.progressbar.setValue(i)
            i=i+1   
            tidetemp=deepcopy(Tides)
            if any(self.survey_dic[keysurv].t):
                tidetemp.interpolateOnGivenTimes(self.survey_dic[keysurv].t)   
                self.survey_dic[keysurv].corr_g=[self.survey_dic[keysurv].grav[i]-self.survey_dic[keysurv].etc[i]-tidetemp.d[i] for i in range(len(self.survey_dic[keysurv].t))]           
                self.survey_dic[keysurv].grav=self.survey_dic[keysurv].corr_g
                self.survey_dic[keysurv].etc=[-tidetemp.d[i] for i in range(len(self.survey_dic[keysurv].t))]           
            PBAR2 = ProgressBar(total=len(self.survey_dic[keysurv].loop_dic.keys()),textmess='loops')   
            PBAR2.show()
            j=1
            for keyloop,loop in self.survey_dic[keysurv].loop_dic.iteritems():
                PBAR2.progressbar.setValue(j)
                j=j+1
                tidetemp2=deepcopy(tidetemp)
                if any(self.survey_dic[keysurv].loop_dic[keyloop].t):
                    tidetemp2.interpolateOnGivenTimes(self.survey_dic[keysurv].loop_dic[keyloop].t)   
                    self.survey_dic[keysurv].loop_dic[keyloop].corr_g=[self.survey_dic[keysurv].loop_dic[keyloop].grav[i]-self.survey_dic[keysurv].loop_dic[keyloop].etc[i]-tidetemp2.d[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].t))] 
                    self.survey_dic[keysurv].loop_dic[keyloop].grav=self.survey_dic[keysurv].loop_dic[keyloop].corr_g
                    self.survey_dic[keysurv].loop_dic[keyloop].etc=[-tidetemp2.d[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].t))] 
                PBAR3 = ProgressBar(total=len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic.keys()),textmess='stations')   
                PBAR3.show()
                k=1
                for keysta,sta in self.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
                    PBAR3.progressbar.setValue(k)
                    k=k+1
                    tidetemp3=deepcopy(tidetemp2)
                    tidetemp3.interpolateOnGivenTimes(self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t)   
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g=[self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav[i]-self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].etc[i]-tidetemp3.d[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t))]
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav=self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g                                           
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].etc=[-tidetemp3.d[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t))]                    


    def oceanCorrectionAgnew(self,amp,phases,lon):
        """
        compute and apply ocean loading correction to the dataset:  
        
        update the gravity field but not the etc field for avoiding possible later conflict with earth tide correction        

        Arguments:
        amplitudes & phases: lists with amplitudes and phases. List order is:
        M2,S2,K1,O1,N2,P1,K2,Q1,Mf,Mm,Ssa      
        (amp,phases)
        lon:     site longitude                    
        """                         
             
        if any(self.t):
            #get tides and round to µgal level
            tides=np.round(np.array([oceanLoading(t,amp,phases,lon) for t in self.t])*1000)/1000 
            self.corr_g=[self.grav[i] +tides[i] for i in range(len(self.t))]
            self.grav=self.corr_g       
        PBAR1 = ProgressBar(total=len(self.survey_dic.keys()),textmess='surveys')   
        PBAR1.show()
        i=1
        for keysurv,surv in self.survey_dic.iteritems():
            PBAR1.progressbar.setValue(i)
            i=i+1   
            if any(self.survey_dic[keysurv].t):
                #get tides and round to µgal level
                tides=np.round(np.array([oceanLoading(t,amp,phases,lon) for t in self.survey_dic[keysurv].t])*1000)/1000
                self.survey_dic[keysurv].corr_g=[self.survey_dic[keysurv].grav[i]+tides[i] for i in range(len(self.survey_dic[keysurv].t))]           
                self.survey_dic[keysurv].grav=self.survey_dic[keysurv].corr_g
            PBAR2 = ProgressBar(total=len(self.survey_dic[keysurv].loop_dic.keys()),textmess='loops')   
            PBAR2.show()
            j=1
            for keyloop,loop in self.survey_dic[keysurv].loop_dic.iteritems():
                PBAR2.progressbar.setValue(j)
                j=j+1
                if any(self.survey_dic[keysurv].loop_dic[keyloop].t):
                    #get tides and round to µgal level
                    tides=np.round(np.array([oceanLoading(t,amp,phases,lon) for t in self.survey_dic[keysurv].loop_dic[keyloop].t])*1000)/1000
                    self.survey_dic[keysurv].loop_dic[keyloop].corr_g=[self.survey_dic[keysurv].loop_dic[keyloop].grav[i]+tides[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].t))] 
                    self.survey_dic[keysurv].loop_dic[keyloop].grav=self.survey_dic[keysurv].loop_dic[keyloop].corr_g
                PBAR3 = ProgressBar(total=len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic.keys()),textmess='stations')   
                PBAR3.show()
                k=1
                for keysta,sta in self.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
                    PBAR3.progressbar.setValue(k)
                    k=k+1
                    #get tides and round to µgal level
                    tides=np.round(np.array([oceanLoading(t,amp,phases,lon) for t in self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t])*1000)/1000
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g=[self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav[i]+tides[i] for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t))]
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav=self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g                                                                            
             
             
    def driftInversion(self,sd_factor=1,sd_add=0.005,drift_t=1,\
                        drift_temp=0,alpha=0.05,write_out_files='n',\
                        output_root_dir='./'):
         """
         calculate drift and find station gravity differences
         based on gravnet (Hwang, 2002) and to a lesser extent on mcgravi(Beilin, 2006)
         
         - update input SD based on sigma_factor and sigma_add
         - no ground reduction is applied, but this could be added easily
         - weighted mean using inverse of variance for each element (in MCGravi, weights are surprisingly SD values...)
         - also calculate SD of the sample. MCRAVI takes as input SE (standard 
            errors, i.e. SD/sqrt(nsample) for each sample and then fails in the
            addition of variance. SE is the standard error of taking the mean
            of the samples as the mean of the overall population. The more 
            samples are taken, the closer the sample mean is to the actual 
            distribution mean (but with no information on the distribution SD)
            see variance of weighted mean.
        - also calculate mean temperature by weighting with 1/SD² values
        - time is also a weighted mean...
            
        Parameters:
        sd_factor:          SD multiplication factor
        sd_add:             SD systematic error (mgal)
        drift_t:            polynomial degree of time drift
        drift_temp:         polynomial degree of temperature drift
        alpha:              Significance level for global model test
        write_out_files:    (y/n): write output files for drift adjustment
                            (similar to MCGravi output files)
        output_root_dir     directory for writting output files
         """
         
         #instanciate a Campaign object to store all observations prior to 
         #inversion
         LS_obs=Campaign()
         for keysurv,surv in self.survey_dic.iteritems():
             if surv.keepitem==1:
                 #empty survey
                 survtmp=Survey(ChannelList(),surv.name)             
                 LS_obs.survey_dic[keysurv]=survtmp
                 for keyloop,loop in self.survey_dic[keysurv].loop_dic.iteritems():
                     if loop.keepitem==1:
                         looptmp=Loop(ChannelList(),loop.name)
                         LS_obs.survey_dic[keysurv].loop_dic[keyloop]=looptmp
                         for keysta,sta in self.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():                     
                             if sta.keepitem==1:
                                 statmp=Station(ChannelList(),sta.stationName,sta.name)
                                 #update SD:
                                 # -get SE values:
                                 sdtmp=[sta.sd[i]/sqrt(sta.dur[i]) for i in range(len(sta.t))]
                                 #print np.mean(sdtmp)
                                 # -update with SD factor and additional terms:
                                 sdtmp=[sdtmp[i]*sd_factor for i in range(len(sta.t))] 
                                 sdtmp=[sqrt(sdtmp[i]*sdtmp[i]+sd_add*sd_add) for i in range(len(sta.t))] 
                                 
                                 gtmp=np.array([sta.grav[i] for i in range(len(sta.t)) if sta.keepdata[i]==1] )
            
                                 ttmp=np.array([date2num(sta.t[i]) for i in range(len(sta.t)) if sta.keepdata[i]==1])
                                 temptmp=np.array([sta.temp[i] for i in range(len(sta.t)) if sta.keepdata[i]==1])
            
                                 #get weighted mean g and variance and temp
                                 weights=np.array([1/(sdtmp[i]*sdtmp[i]) for i in range(len(sta.t)) if sta.keepdata[i]==1])
                                 sdmean=sqrt(1/sum(weights))
                                 #print sdmean
                                 
                                 #mcgravi version (seems wrong)
                                 #weights=np.array([sdtmp[i] for i in range(len(sta.t)) if sta.keepdata[i]==1])                     
                                 #mcgravi is wrong in this calculation, the result is pretty much this:                     
                                 #sdmean=sdtmp[-1]
                                 #sdmean=0.005 
                                 #sdmean=0
                                 #for i in range(len(sta.t)):
                                 #    if sta.keepdata[i]==1:
                                 #        sdmean=sdmean*sdmean+sdtmp[i]*sdtmp[i]
                                 #sdmean=sqrt(sdmean)
                                 
                                 gmean=sum(gtmp*weights)/sum(weights)
                                 tmean=sum(ttmp*weights)/sum(weights)
                                 tempmean=sum(temptmp*weights)/sum(weights)
                                 statmp.grav=gmean
                                 statmp.sd=sdmean   
                                 statmp.temp=tempmean
                                 statmp.t=tmean
                                 
                                 LS_obs.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta]=statmp
         
         #go to the survey level for drift estimates
         for keysurv,surv in self.survey_dic.iteritems(): 
             if surv.keepitem==1:
                 self.survey_dic[keysurv].calculateSimpleDiff(LS_obs.survey_dic[keysurv],\
                 drift_t,drift_temp,alpha,sd_factor,sd_add,write_out_files,output_root_dir)
                     
###############################################################################
class Survey(ChannelList):
    """Survey-type object
  
    Derived (inheritance) from a ChannelList-type object 
    all properties and functions from ChannelList apply
    
    Properties:
    loop_dic:             dictionary of Loop objects
    output_dic:           dictionary of output gravity values:
                          at keys = station numbers, 
                          a tuple (station_name,relative gravity value,SD)
                          when corrections are applied, the tuple if longer:
                          (station_name,rel. g.,SD, gcprr_height,gcorr_Bouguer)
    height_bouguer_corr:  'ok' or 'no' if height and Bouguer corrections have 
                          been applied on simple difference                          

    Functions:
    - populateLoopDic
    - read_from_mcgravi_output_file
    - writeOutputValues
    """

    loop_dic = {}
    output_dic={}
    height_bouguer_corr=None
    
    def __init__(self,k,name):
        
        """
        Initializes to an empty dictionary.
        Should be initialized with another ChannelList object (make sure the
        latter is deecopied first). Only way found to apply the extractsubset
        function to another ChannelList object, and having a Survey object 
        as an output...         
       
        """
        super(Survey,self).__init__()# call properties from the baseclass
        self.__dict__ = k.__dict__        
        self.loop_dic = {}
        self.output_dic={}   
        self.height_bouguer_corr='no'
        self.name=name
        
    def populateLoopDic(self,option,base_station,loop_number=1):
        """populate loop dictionary following user-defined rules
        

        Arguments:
        option:                 1: provide base station number, and number of 
                                   loops
                                2: automatic loop identification: identify base
                                   station based on gravity value (case where 
                                   station numbers have not been entered during 
                                   the survey). the base station is the one 
                                   which is repeated more often (loop number+1)
                                   
        base station:           base station number. If option = 1 

        loop number:            number of loops. If option = 2
        """            
        if option==1:
            #initialize "previous station"
            prev_sta=999999
            start_loop=0 
            first_base='n'
            nloop=1
            PBAR = ProgressBar(total=len(self.station),textmess='populate loop dictionary')   
            PBAR.show()
                
            #loop over samples:
            for i in range(len(self.station)) :
                PBAR.progressbar.setValue(i)
                curr_sta=self.station[i]
                if (curr_sta-prev_sta)!=0 :
                    if curr_sta==base_station and start_loop==1 and first_base=='y':
                        first_base='n'
                        ind_start_next=i
                    elif curr_sta==base_station and start_loop==0:
                        start_loop=1
                        first_base='y'
                        ind_start=i                        
                    elif curr_sta!=base_station and start_loop==1 and first_base=='n':
                        ind_end=i-1
                        temp_surv=deepcopy(self)
                        temp_loop=temp_surv.extract_subset(self.t[ind_start],self.t[ind_end])
                        # create a new Loop-type object:                                
                        temp_loop2=Loop(temp_loop,str(nloop))
                        temp_loop2.populateStationDic()
                        self.loop_dic[nloop]=temp_loop2
                        nloop=nloop+1
                        ind_start=ind_start_next
                        first_base='y'
                elif i == len(self.station)-1:
                    if curr_sta==base_station and start_loop==1 and first_base=='n':
                        ind_end=i
                        temp_surv=deepcopy(self)
                        temp_loop=temp_surv.extract_subset(self.t[ind_start],self.t[ind_end])
                        # create a new Loop-type object:                                
                        temp_loop2=Loop(temp_loop,str(nloop))
                        temp_loop2.populateStationDic()
                        self.loop_dic[nloop]=temp_loop2
                prev_sta=curr_sta            
        print 'number of loops: %d'%(nloop)        
        if option==2:
            raise UserWarning('option not implemented yet')        
          
    def read_from_mcgravi_output_file(self,mcgravi_filename):
        """
        read mcgravi output file, find the data at the end of the 
        file and record gravity values and SD relative to each
        station by filling a dictionnary
        """
    
        self.output_dic.clear()
        file=open(mcgravi_filename,'r')
        lines=file.readlines()
        
        # code for reading output in .gra file
        for line in lines:
          vals=line.split()
          if line.split() and vals[0].isdigit():
            self.output_dic[vals[0]]=(int(vals[0]),float(vals[1]),float(vals[2]))                
        
        # code for reading output in .lst file
        #test=0
        #for line in lines:
        #  vals=line.split()
        #  if test==1:
        #    test=2
        #  if line[0:23]==" Nom, Gravite compensee":
        #    test=1    
        #  if test==2:
        #    vals=line.split()
        #    if line.split() and vals[0].isdigit():
        #      self.output_dic[vals[0]]=(int(vals[0]),float(vals[1]),float(vals[2]))
        file.close()
           
        
    def writeOutputValues(self,filename):
        """
        write Output values, after processing and drift adjustment
        """
        file=open(filename,'w')
        file.write("station   g      sd\n")
        for station_id,sta in self.output_dic.iteritems():            
            file.write("%d %2.3f %2.3f\n"%(sta[0], sta[1], sta[2]))
        file.close()     
        
    def readSimpleDiff(self,filename):
        """
        Populate the output_dic dictionary by reading the given filename
        """
        file=open(filename,'r')
        lines=file.readlines()
        #remove header:
        lines.pop(0)
        self.output_dic.clear()        
        for line in lines:
            line = line.strip() 
            vals=line.split()
            self.output_dic[vals[0]]=(int(vals[0]),float(vals[1]),float(vals[2]))
        file.close()
        
        
    def calculateSimpleDiff(self,LS_obs,drift_t,drift_temp,alpha,\
                                sd_factor,sd_add,write_out_files,\
                                output_root_dir):
        """
        Calculate simple difference by LS inversion of dataset
        LS_obs is a SURVEY-type object which contain all the data for the inversion
        
        Parameters
        LS_obs:             a SURVEY-type object which contain all the data for
                            the inversion (single weighted mean values for each
                            stations)
        drift_t:            polynomial degree for time
        drift_k:            polynomial degree for temperature
        alpha:              Significance level for global model test
        sd_factor:          SD multiplication factor (for writting only)
        sd_add:             SD systematic error (mgal) (for writting only)
        write_out_files:    (y/n): write output files for drift adjustment
                            (similar to MCGravi output files)
        output_root_dir     directory for writting output files        
        """
        sta_dic_LS=LS_obs.getStationsIndices()
        #get number of relative observations:
        n_rel_Obs=0
        #number of absolute obs.
        n_obs_abs=1
        base_station_number=1
        for keyloop,loop in LS_obs.loop_dic.iteritems():
            n_rel_Obs=n_rel_Obs+len([1 for keysta,sta in LS_obs.loop_dic[keyloop].station_dic.iteritems()])-1
     
        nloops=sum([1 for keyloop,loop in LS_obs.loop_dic.iteritems()])
        
        #number of unknowns:
        nb_X=len(sta_dic_LS)+drift_t*nloops+drift_temp*nloops
        #model matrix:
        A=np.zeros((n_rel_Obs+n_obs_abs,nb_X))
        #weight matrix:
        P=np.zeros((n_rel_Obs+n_obs_abs,n_rel_Obs+n_obs_abs))#pas sur
        #observation matrix:
        Obs=np.zeros((n_rel_Obs+n_obs_abs,1))
        #datum-free constraint vector:
        S=np.zeros((nb_X,1))
        
        tot_stat=0
        tot_loop=0
        for keyloop,loop in LS_obs.loop_dic.iteritems():
            #check that same order is kept!!!!=> should be ok with OrderedDict
            g=np.array([sta.grav for keysta,sta in loop.station_dic.iteritems()])
            sd=np.array([sta.sd for keysta,sta in loop.station_dic.iteritems()])            
            temp=np.array([sta.temp for keysta,sta in loop.station_dic.iteritems()])            
            stanumber=np.array([keysta[0] for keysta,sta in loop.station_dic.iteritems()])   
            t=np.array([sta.t for keysta,sta in loop.station_dic.iteritems()])
            for i in range(len(g)-1):
                Obs[tot_stat+i]=g[i+1]-g[i]
                P[tot_stat+i,tot_stat+i]=1/(sd[i+1]**2+sd[i]**2)
                
                A[tot_stat+i,sta_dic_LS[stanumber[i+1]]]=1
                A[tot_stat+i,sta_dic_LS[stanumber[i]]]=-1
                
                S[sta_dic_LS[stanumber[i]]]=1
                S[sta_dic_LS[stanumber[i+1]]]=1                
                for j in range(drift_t):    
                    A[tot_stat+i,len(sta_dic_LS)+tot_loop*drift_t+j]=(t[i+1]-t[i])**(j+1)
                for j in range(drift_temp):                    
                    A[tot_stat+i,len(sta_dic_LS)+drift_t*nloops+tot_loop*drift_temp+j]=(temp[i+1]-temp[i])**(j+1)
                    
            tot_stat=tot_stat+len(g)-1
            tot_loop=tot_loop+1

            
        #add absolute observation (base station with known value)
        #should be looping over abs. values
        A[n_rel_Obs+n_obs_abs-1,sta_dic_LS[base_station_number]]=1
        P[n_rel_Obs+n_obs_abs-1,n_rel_Obs+n_obs_abs-1]=1/(0.001**2)
        
        Obs[n_rel_Obs+n_obs_abs-1]=0    
        
        # Do the inversion
        lsd=LSdata()
        lsd.A=A
        lsd.P=P
        lsd.Obs=Obs
        lsd.S=S
        lsd.dof=n_rel_Obs+n_obs_abs-nb_X
        lsd.lsInversion()
        for i in range(len(sta_dic_LS)):
            for key,val in sta_dic_LS.iteritems():
                if val==i:
                    keytmp=key
            self.output_dic[keytmp]=(keytmp,float(lsd.X[i]),sqrt(lsd.var[i]))
        

        #display results
        print "For survey: %s"%(self.name)
        print "Station , g (mgal), SD (mgal)"
        for tupid,tup in self.output_dic.iteritems():
            print "%d, %7.4f, %7.4f"%(tup)             
            
        #calculate and display statistics:            
        lsd.lsStatistics(alpha)
        
        #write output files
        if write_out_files=='y':            
            survdir=output_root_dir+os.sep+self.name
            if not os.path.exists(survdir):
                os.makedirs(survdir)
                
            #Write simple results (equivalent to simple diff)
            tday=datetime.now()  
            file=open(survdir+os.sep+\
            "LSresults_%4d%02d%02d_%02d%02d.dat"%(tday.year,\
            tday.month,tday.day,tday.hour,tday.minute),'w')
            file.write("station   g      sd\n")
            for station_id,sta in self.output_dic.iteritems():            
                file.write("%d %7.4f %7.4f\n"%(sta[0], sta[1], sta[2]))                      
            file.close()

            #Write complete results (similar to .lst MCGRAVI files)
            tday=datetime.now()  
            file=open(survdir+os.sep+\
            "LSresults_tot_%4d%02d%02d_%02d%02d.dat"%(tday.year,\
            tday.month,tday.day,tday.hour,tday.minute),'w')
            file.write(" ==============================================================================\n")
            file.write(" Modification of input SD values:\n")
            file.write(" Factor:        %5.2f\n"%(sd_factor))
            file.write(" Additive term: %6.4f\n"%(sd_add))            
            file.write(" ==============================================================================\n\n")
            file.write("Observations (weighted means):\n\n")
            file.write("station, gravity (mgal), SD(mgal), Temp, t\n\n")
            for keyloop,loop in LS_obs.loop_dic.iteritems():
                file.write("Loop: %s\n"%(keyloop))
                for keysta,sta in LS_obs.loop_dic[keyloop].station_dic.iteritems():
                    file.write("%04d    %9.4f %6.4f %5.2f %7.2f\n"%(int(keysta[0]),\
                    sta.grav,sta.sd,sta.temp,sta.t))
                file.write("\n")
            file.write(" ==============================================================================\n")
            file.write("Relative observations:\n\n")
            file.write("station1, station2, gravity (mgal), SD(mgal)\n\n")
            for keyloop,loop in LS_obs.loop_dic.iteritems():
                file.write("Loop: %s\n"%(keyloop))
                g=np.array([sta.grav for keysta,sta in loop.station_dic.iteritems()])
                sd=np.array([sta.sd for keysta,sta in loop.station_dic.iteritems()])            
                stanumber=np.array([keysta[0] for keysta,sta in loop.station_dic.iteritems()])                   
                for i in range(len(g)-1):
                    file.write("%04d  %04d    %9.4f %6.4f\n"%(stanumber[i],\
                    stanumber[i+1],g[i+1]-g[i],np.sqrt(sd[i+1]**2+sd[i]**2)))
                file.write("\n")    
            file.write(" ==============================================================================\n")                
            file.write("statistics\n\n")
            file.write("Number of stations:                 %d\n"%len(sta_dic_LS))                    
            file.write("Number of loops:                    %d\n"%tot_loop)                    
            file.write("Polynomial degree for time:         %d\n"%drift_t)                    
            file.write("Polynomial degree for temperature:  %d\n"%drift_temp)   
            file.write("\n")               
            file.write("Number of unknowns:                 %d\n"%nb_X)   
            file.write("Number of relative observations:    %d\n"%n_rel_Obs)   
            file.write("Number of absolute observations:    %d\n"%n_obs_abs)   
            file.write("Degree of freedom (nobs-nunknowns): %d\n"%lsd.dof)   
            file.write("\n")               
            file.write("SD a posteriori:                    %f\n"%float(np.sqrt(lsd.VtPV/lsd.dof)))   
            file.write("chi square value:                   %6.2f\n"%float(lsd.chi2))
            file.write("critical chi square value:          %6.2f\n"%float(lsd.chi2c))        
            if lsd.chi2<lsd.chi2c:
                print "Chi-test accepted\n"
            else:
                print "Chi-test rejected\n"    
                
            file.write("\n")    
            file.write(" ==============================================================================\n")                
            file.write("Final results\n\n")

            file.write("station   g      sd\n")
            for station_id,sta in self.output_dic.iteritems():            
                file.write("%d    %7.4f %7.4f\n"%(sta[0], sta[1], sta[2]))      
            file.write("\n")    
            file.write("Drift values:\n")
            for i in range(tot_loop):
                file.write("Loop %d:    gravi     SD\n"%(i+1))
                file.write("Time\n")
                for j in range(drift_t):       
                    ind=len(sta_dic_LS)+i*drift_t+j
                    file.write("Degree %d: %9.4f %9.4f\n"%(j+1,lsd.X[ind],\
                    lsd.var[ind]))
                file.write("Temperature\n")
                for j in range(drift_temp):       
                    ind=len(sta_dic_LS)+tot_loop*drift_t+i*drift_temp+j
                    file.write("Degree %d: %9.4f %9.4f\n"%(j+1,lsd.X[ind],\
                    lsd.var[ind]))    
                file.write("\n")    
                           
            file.close()
    def getStationsIndices(self):
        """
        return a dictionnary with unique (key=stationname,value=integer) pairs
        """
        
        station_list=[]
        sta_dic_LS={}
        for keyloop,loop in self.loop_dic.iteritems():
            #sta[1] est l'occurence...
            station_list=np.append(station_list,[keysta[0] for keysta,sta in self.loop_dic[keyloop].station_dic.iteritems()])
        station_list=np.sort(np.unique(station_list))
        for i in range(len(station_list)):
            sta_dic_LS[station_list[i]]=i
        return sta_dic_LS
###############################################################################
class Loop(ChannelList):
    """Loop-type object
  
    Derived (inheritance) from a ChannelList-type object 
    all properties and functions from ChannelList apply  
  
    Properties:
  
    station_dic		dictionary of Station objects. Keys are tuples: 
                       (station number, number of reoccupations)


    Functions:
    - getStationReoccupation
    - populateStationDic
    - writeCgxCFile
    - writeModifCFile
    - writeAssociatedCgxSFile
    
    """

    #station_dic = {}
    station_dic = OrderedDict()#Dictionary that remembers insertion order
    
    def __init__(self,k,name):
        """
        Initializes data to an empty dictionary.
        Should be initialized with another ChannelList object (make sure the
        latter is deecopied first). Only way found to apply the extractsubset
        function to another ChannelList object, and having a Loop object 
        as an output...             
        
        """
        super(Loop,self).__init__()        # call properties from the baseclass       
        self.__dict__ = k.__dict__   
        #self.station_dic = {}
        self.station_dic = OrderedDict()
        self.name=name
        
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print 'length of field t: %d'%(len(self.t))       
        print 'first station: %s'%(self.station[0])   
        print 'last station: %s'%(self.station[len(self.station)-1]) 
        print 'station availables: %s'%(set(self.station))
        return "you're welcome"
        
        

    def getStationReoccupation(self,station):
        """get station reoccupation
        mostly to be used iteratively within th populateStationDic function
        """
        station_reoc=[sta[0] for sta,stakey in self.station_dic.iteritems() if sta[0]==station]
        return(len(station_reoc))
        
    def populateStationDic(self):
        """populate station dictionary
        
        """
        #initialize "previous station"
        prev_sta=self.station[0]      
        ind_start=0
        #loop over samples:        
        for i in range(len(self.station)) :
            curr_sta=self.station[i]
            if (curr_sta-prev_sta)!=0 :
                #enter in a new station                    
                ind_end=i-1
                temp_loop=deepcopy(self)
                # create a new Loop-type object:  
                temp_sta=temp_loop.extract_subset(self.t[ind_start],self.t[ind_end])
                keysta=(int(prev_sta),self.getStationReoccupation(prev_sta)+1)
                temp_sta2=Station(temp_sta,prev_sta,'station '+str(keysta[0])+' ('+str(keysta[1])+')')
                temp_sta2.keepitem=1
                self.station_dic[keysta]=temp_sta2
                ind_start=i
            if i==len(self.station)-1:
                # enter last data line
                ind_end=i
                temp_loop=deepcopy(self)
                # create a new Loop-type object:  
                temp_sta=temp_loop.extract_subset(self.t[ind_start],self.t[ind_end])
                keysta=(int(curr_sta),self.getStationReoccupation(curr_sta)+1)
                temp_sta2=Station(temp_sta,curr_sta,'station '+str(keysta[0])+' ('+str(keysta[1])+')')
                temp_sta2.keepitem=1                
                self.station_dic[keysta]=temp_sta2
            prev_sta=curr_sta
            
                
    def writeCgxCFile(self,filename,base_station):
        """
        write cgx type c file (without header)
        need to check which is corrected from tides
        re-write the sd as former SD/sqrt(duration)
        last column is the first gravity value of the 'o' file
        hence the first of the dictionnary, corrected for the tides
        
        should add a checkbox in the tree object to select the base station
        """
    
        offset=self.station_dic[base_station].grav[0]
        print "filename: %s"%filename
        file=open(filename,'w')
        
        for station_id,sta in sorted(self.station_dic.iteritems(), key=lambda x: x[1].t[1]):                  
        
        #for station_id,sta in self.station_dic.iteritems():
            if sta.keepitem==1:
                for i in range(len(sta.t)):
                    if sta.keepdata[i]==1:                
            	        file.write("%6d %4.3f %1.3f  %2d   %2d	%2.1f	%2.1f	%2.2f	%1.3f	%3d	%3.3f %02d%02d%02d %02d%02d%02d   0  0.000     %3.4f\n"%(sta.station[i],\
                            sta.grav[i]-sta.etc[i],sta.sd[i]/sqrt(sta.dur[i]),
                            sta.dur[i], sta.rej[i], sta.tiltx[i], sta.tilty[i],
                            sta.temp[i],sta.etc[i],sta.t[i].timetuple()[7],
                            sta.t[i].hour*60+sta.t[i].minute+sta.t[i].second/60,
                            sta.t[i].day,sta.t[i].month,sta.t[i].year-2000,
                            sta.t[i].hour,sta.t[i].minute,sta.t[i].second,
            	              sta.grav[i]-offset))
    
        file.close()

    def writeModifCFile(self,filename,base_station):
        """
        write modified cgx type c file (without header)
        need to check which is corrected from tides
        re-write the sd as former SD/sqrt(duration)
        last column is the first gravity value of the 'o' file
        hence the first of the dictionnary, corrected for the tides
        add an extra column with the keepdata value (0 or 1) for further 
        loading        
        
        should add a checkbox in the tree object to select the base station
        IMPORTANT: grav is corrected from tides (as in the classical raw data), 
        which is different than for classical c files
        and SD is SD, not SD/sqrt(dur) as for classical c files
        OR an other option could be to really write 'c' files and modify the 
        reading function (read_c_file)
        """
    
        offset=self.station_dic[base_station].grav[0]
        print "filename: %s"%filename
        file=open(filename,'w')
        
        for station_id,sta in sorted(self.station_dic.iteritems(), key=lambda x: x[1].t[1]):                  
        
        #for station_id,sta in self.station_dic.iteritems():
            if sta.keepitem==1:
                for i in range(len(sta.t)):               
        	         file.write("%6d %4.3f %1.3f  %2d   %2d	%2.1f	%2.1f	%2.2f	%1.3f	%3d	%3.3f %02d%02d%02d %02d%02d%02d   0  0.000     %3.4f    %1d\n"%(sta.station[i],\
                        sta.grav[i],sta.sd[i],
                        sta.dur[i], sta.rej[i], sta.tiltx[i], sta.tilty[i],
                        sta.temp[i],sta.etc[i],sta.t[i].timetuple()[7],
                        sta.t[i].hour*60+sta.t[i].minute+sta.t[i].second/60,
                        sta.t[i].day,sta.t[i].month,sta.t[i].year-2000,
                        sta.t[i].hour,sta.t[i].minute,sta.t[i].second,
        	              sta.grav[i]-offset, sta.keepdata[i]))
    
        file.close()
  
    def writeAssociatedCgxSFile(self,filename):
        """
        write cgx station 's' file associated to the dictionnary of recordings
       
        """
        file=open(filename,'w')
        for station_id,sta in sorted(self.station_dic.iteritems(), key=lambda x: x[1].t[1]):                          
#        for station_id,sta in self.station_dic.iteritems():
            if sta.keepitem==1:            
                for i in range(len(sta.t)):
                    if sta.keepdata[i]==1:  
                        file.write("      %2d %02d:%02d 0.000 9999.000 9999.000 9999.000\n"%(sta.station[i], sta.t[i].hour, sta.t[i].minute))
        file.close()
                
###############################################################################

class Station(ChannelList):
    """Station-type object
  
    Derived (inheritance) from a ChannelList-type object 
    all properties and functions from ChannelList apply 
    
    Properties:
  
    reading_dic		dictionary of Reading objects
  
  
    """


    reading_dic={}	
    station_name=None

    def __init__(self,k,stationName,name):
        """
        Initializes data to an empty dictionary.
        
        """
        super(Station,self).__init__()     # call properties from the baseclass   
        self.__dict__ = k.__dict__        
        self.reading_dic={}
        self.stationName=stationName
        self.name=name


    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        """
        print 'length of field t: %d'%(len(self.t))       
        print 'first station: %s'%(self.station[0])   
        print 'last station: %s'%(self.station[len(self.station)-1]) 
        print 'station availables: %s'%(set(self.station))
        return "you're welcome"


###############################################################################
#other functions and classes definitions

class timeSeries(object):
    """Time Series object
    used to store a simple time series
    could have been a ChannelList object
    Used for storing synthetic tides, or atmospheric time series for instance
    
    Properties:
    t:                      time vector: datetime format
    d:                      data
    
    Functions:
    - populateFromTsoftFile
    - populateFromEternaFile
    - interpolateOnGivenTimes
    """    
    t=[]
    d=[]
    name=None
    
  
    def __init__(self):
        """
        """
        self.t=[]
        self.d=[]

        
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print name    
        return "you're welcome"

    def populateFromTsoftFile(self,filename,channel):
        """
        load a .tsf file and populate the object's properties
        load data from channel
        Assume data is in nm/s² and convert it to mgal!
        """
        try:
            #essaye d'ouvrir le fichier
            fh = open(filename, 'r')
            i=0

            PBAR = ProgressBar(total=len([1 for line in  open(filename, 'r')]),textmess='Load data')   
            
            print "number of lines: %d"%len([1 for line in  open(filename, 'r')])
            PBAR.show()
            test=0
            for line in fh:    
                PBAR.progressbar.setValue(i)
                i+=1
                # Clean line
                line = line.strip()
                # Skip blank and comment lines
                if (not line) or (line[0] == '/') or (line[0] == 'L'): continue
        	     #parse string line first with respect to '/' caracters (used in the date format), 
        	     #then with ':' (used for the time display), eventually with the classic ' '
                vals=line.split()
                if test==1:
                    self.t.append(datetime(int(vals[0]),int(vals[1]),
                                           int(vals[2]),int(vals[3]),
                                           int(vals[4]),int(vals[5])))
                    #read nm/s² data and convert it to mgal:                        
                    self.d.append(float(vals[5+channel])/10000)                                           
                if vals[0]=="[DATA]":
                    test=1                                                                                      
        except IOError:
            #si ça ne marche pas, affiche ce message et continue le prog
            print 'No file : %s' %(filename)            
        except ValueError:
            print 'pb at line %d : Is it really .tsf (or .TSF) format? '%(i)
        except IndexError:
            print 'pb at line %d : check raw data file: possibly last line?'%(i)  
        fh.close()        
        
    def populateFromEternaFile(self,filename,channel):
        """
        load an eterna file and populate the object's properties
        load data from channel
        Assume data is in nm/s² and convert it to mgal!
        """
        try:
            #essaye d'ouvrir le fichier
            fh = open(filename, 'r')
            i=0

            PBAR = ProgressBar(total=len([1 for line in  open(filename, 'r')]),textmess='Load data')   
            
            print "number of lines: %d"%len([1 for line in  open(filename, 'r')])
            PBAR.show()
            test=0
            for line in fh:    
                PBAR.progressbar.setValue(i)
                i+=1
                # Clean line
                line = line.strip()
                # Skip blank and comment lines
                if (not line) or (line[0] == '/') or (line[0] == 'L'): continue
        	     #parse string line first with respect to '/' caracters (used in the date format), 
        	     #then with ':' (used for the time display), eventually with the classic ' '
                vals=line.split()
                if vals[0]=="99999999":
                    test=0     
           
                if test==1:
                    #padd with 0
                    ttemp=vals[1].zfill(6)
                    self.t.append(datetime(int(vals[0][0:4]),int(vals[0][4:6]),
                                           int(vals[0][6:8]),int(ttemp[0:2]),
                                           int(ttemp[2:4]),int(ttemp[4:6])))
                    #read nm/s² data and convert it to mgal:                       
                    self.d.append(float(vals[1+channel])/10000)                                           
                if vals[0]=="77777777":
                    test=1                                                                                      
        except IOError:
            #si ça ne marche pas, affiche ce message et continue le prog
            print 'No file : %s' %(filename)            
        except ValueError:
            print 'pb at line %d : Is it really eterna format? '%(i)
        except IndexError:
            print 'pb at line %d : check raw data file: possibly last line?'%(i)  
        fh.close() 
        
    def interpolateOnGivenTimes(self,t):
        """
        interpolate the time series on the user input time vector
        overlain the previous t and d fields
        """
        tord=[date2num(tmp) for tmp in self.t]
#        tord=date2num(self.t) # create np.array
        f=interp1d(tord, self.d, kind='linear',bounds_error=False)
        self.d=f([date2num(tmp) for tmp in t])    
#        self.d=f(date2num(t)) # create np.array...    
        self.t=t
        

def read_start_end_dates(filename):
    """ read user input start and end dates of each campaign
    current format is
    yyyy/mm/dd hh:mn:ss yyyy/mm/dd hh:mn:ss 
    yyyy/mm/dd hh:mn:ss yyyy/mm/dd hh:mn:ss 
    ...
    but can be changed anytime here
    """    
   
    try:
        #essaye d'ouvrir le fichier
        fh = open(filename, 'r')
        i=0
        start_end_dates=[]
        for line in fh:    
            i+=1
            # Clean line
            line = line.strip()
    	     #parse string line first with respect to '/' caracters (used in the date format), 
    	     #then with ':' (used for the time display), eventually with the classic ' '
            vals=line.replace('/',' ').replace(':',' ').split()
            start_end_dates.append((datetime(int(vals[0]),int(vals[1]),\
            int(vals[2]),int(vals[3]),int(vals[4]),int(vals[5])),\
            datetime(int(vals[6]),int(vals[7]),int(vals[8]),int(vals[9]),\
            int(vals[10]),int(vals[11]))))
            
        return start_end_dates             
    except IOError:
        #si ça ne marche pas, affiche ce message et continue le prog
        print 'No file : %s' %(filename)            
    except ValueError:
        print 'pb at line %d : check data file'%(i)
        
                       
        
class ProgressBar(QtGui.QWidget):
    """
    define progress bar
    """
    
    def __init__(self, parent=None, total=20, textmess='Progress'):
        super(ProgressBar, self).__init__(parent)
        self.progressbar = QtGui.QProgressBar()
        self.progressbar.setMinimum(1)
        self.progressbar.setMaximum(total)
        main_layout = QtGui.QGridLayout()
        main_layout.addWidget(self.progressbar, 0, 1)
        self.setLayout(main_layout)
        self.setWindowTitle(textmess)
        
        
def write_mcgravi_conf_file(filename,file_list,mode,drift_t,drift_k,write_list_fic,write_list_station,\
    write_obs,write_mat,write_resid,write_tau,write_grav,write_drift,sigma_factor,sigma_add,calf,create_r,outf,fcor,region):
    """
    write mcgravi configuration file, containing all the parameters for the program
   
    """
    file=open(filename,'w')
    file.write("MODE %1d\n"%mode)
    file.write("DRIFT_T %1d\n"%drift_t)
    file.write("DRIFT_k %1d\n"%drift_k)
    file.write("WRITE_LIST_FIC %s\n"%write_list_fic)
    file.write("WRITE_LIST_STATION %s\n"%write_list_station)
    file.write("WRITE_OBS %s\n"%write_obs)
    file.write("WRITE_MAT %s\n"%write_mat)
    file.write("WRITE_RESID %s\n"%write_resid)
    file.write("WRITE_TAU %s\n"%write_tau)
    file.write("WRITE_GRAV %s\n"%write_grav)
    file.write("WRITE_DRIFT %s\n"%write_drift)
    file.write("SIGMA_FACTOR %1.1f\n"%sigma_factor)
    file.write("SIGMA_ADD %s\n"%sigma_add)
    file.write("CALF %s\n"%calf)

    for ofilename in file_list:
        file.write("RELF %s 1.000 0.000\n"%ofilename)

    file.write("CREATE_R %s\n"%create_r)
    file.write("ABSF tunnel.abs 1.000 0.000\n")
    file.write("OUTF %s\n"%outf)
    file.write("FCOR %s\n"%fcor)
    file.write("REGION %s\n"%region)
    file.close()        
    
    
    
#########################################
    ############ PARTIE STATIQUE
    
class StaticDataSet(object):
    """
                
    """   

    Station_dic = {}
    def __init__(self):
        """
        Initializes to an empty dictionary.
    
        """
        self.Station_dic = {}

    def readSimpleDiff_static(self,filename):
        """
        Populate the output_dic dictionary by reading the given filename
        """
        file=open(filename,'r')
        lines=file.readlines()
        #remove header:
        lines.pop(0)
        for line in lines:
            line = line.strip() 
            vals=line.split()
            station_temp=Station_static()         
            station_temp.station=int(vals[0])
            station_temp.grav=float(vals[1])
            station_temp.sd=float(vals[2])
            print station_temp
            self.Station_dic[int(vals[0])]=station_temp
        file.close()        

    def readStation_Information(self,filename):
        """
        """
        file=open(filename,'r')
        lines=file.readlines()
        #remove header:
        lines.pop(0)
        self.output_dic.clear()        
        for line in lines:
            line = line.strip() 
            vals=line.split()
            self.station=int(vals[0])
            self.grav=float(vals[1])
            self.sd=float(vals[2])
        file.close() 
        
    def writePGraviFor3DInputFile(self,filename):
        """
        Write input files for PGraviFor3D
        """        
        file=open(filename,'w')
        file.write("Station_number  X(m)    Y(m)    Z(m)    lat(°)  lon(°)  g(mgal)\n")
        for station_id,sta in self.Station_dic.iteritems():                         
            file.write("%d  %f  %f  %f  %f  %f  %f\n"%(sta.station,\
                            sta.x,sta.y,sta.alt,sta.lat,sta.lon,sta.grav))
    
        file.close()        
        
        
class Station_static(object):     
    """Station_static-type object
    
    Properties:
    
    list objects (time series):
    station          station name
    alt              altitude (m)
    grav             gravity value (mgal)
    sd               Standard Deviation (mgal)
    lat
    lon
    x
    y

    functions:

                
    """
    station=int
    alt=float
    grav=float
    sd=float
    lat=float
    lon=float
    x=float
    y=float
  
    def __init__(self):
        """
        """
        self.station=0
        self.alt=0
        self.grav=0
        self.sd=0
        self.lat=48
        self.lon=7
        self.x=1000000
        self.y=2000000

    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print 'station name: %d'%(self.station)       
        print 'altitude (m): %f'%(self.alt)   
        print 'gravity (mgal): %f'%(self.grav) 
        print 'SD (mgal): %f'%(self.sd)
        return "you're welcome"




class LSdata(object):
    """LSdata-type object (Least-Square data)
  
    keep data for least square inversion    
    
    Properties:
    A                   model matrix (Nobs rel + NobsAbs)*Nunknowns
    P                   Weight matrix (Nobs rel + NobsAbs)*(Nobs rel + NobsAbs)
    Obs                 observation vector (Nobs rel + NobsAbs)
    S                   StX=0
    X                   Unknowns
    r                   residuals (V)
    var                 Diagonal of the a posteriori covariance matrix
    VtPV                
    dof                 degree of freedom (Nobs rel + NobsAbs - Nunknowns)            

    Functions:

    """
    
    A=np.array([])
    P=np.array([])
    Obs=np.array([])
    S=np.array([])
    X=np.array([])
    r=np.array([])
    var=np.array([])
    VtPV=np.array([])  
    SDaposteriori=0
    dof=0
    chi2=0
    chi2c=0
    
    def __init__(self):
        
        """
        Initializes        
       
        """
        self.A=np.array([])
        self.P=np.array([])
        self.Obs=np.array([])
        self.S=np.array([])
        self.X=np.array([])
        self.r=np.array([])
        self.var=np.array([])        
        self.VtPV=np.array([])  
        self.SDaposteriori=0
        self.dof=0
        self.chi2=0
        self.chi2C=0

    def lsInversion(self):
        """
        LS Inversion from Hwang et al (2002)
        """
        
        At=np.transpose(self.A)
        St=np.transpose(self.S)
        N=At.dot(self.P).dot(self.A)
        #solution:
        self.X=np.linalg.inv(N+self.S.dot(St)).dot(At).dot(self.P).dot(self.Obs)
        
        self.r=self.A.dot(self.X)-self.Obs
        rt=np.transpose(self.r)
        self.VtPV=rt.dot(self.P).dot(self.r)
        var_post_norm=self.VtPV/self.dof
        self.SDaposteriori=np.sqrt(var_post_norm)
        
        cov_post=np.linalg.inv(N)*var_post_norm
        self.var=np.diagonal(cov_post)        
        
    def lsStatistics(self,alpha):
        """
        a priori variance of unit weight = 1
        """
        self.chi2=self.VtPV
        
        t=np.sqrt(2*np.log(1/alpha))
        chi_1_alpha=t-(2.515517+0.802853*t+0.010328*t**2)/(1+1.432788*t+0.189269*t**2+0.001308*t**3)
        dof=float(self.dof)
        self.chi2c=dof*(chi_1_alpha*sqrt(2/(9*dof))+1-2/(9*dof))**3
        
        print "SD a posteriori: %6.2f"%float(self.SDaposteriori)
        print "chi square value: %6.2f"%float(self.chi2)
        print "critical chi square value: %6.2f"%float(self.chi2c)
        
        if self.chi2<self.chi2c:
            print "Chi-test accepted"
        else:
            print "Chi-test rejected"
        
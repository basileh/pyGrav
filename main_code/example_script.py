# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:26:25 2015

@author: Basile HECTOR


##########################
        pyGrav example
##########################
        
This is an example script for working with pyGrav functions without the
GUI, using only the data_objects.py file and not the pyGrav_main.py file.        
"""

from data_objects import *
import sys,os,shutil,subprocess,glob
from PyQt4 import QtGui,QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from copy import deepcopy
from datetime import *
import numpy as np
from matplotlib.dates import date2num,num2date

# still building a Qt application for displaying progress bars defined in the code.
qApp = QtGui.QApplication(sys.argv)

#options:
data_dir="../test_case/input_data/raw/"
output_dir="../test_case/output_data/"

#instanciate a Campaign object:
data=Campaign()

# 0) read raw data:
#filename=data_dir+"Atest_Raw_data.txt"
#data.readRawDataFile(filename)
#
#
## I.a) fill in the Campaign object :
## Options:
#opt=1                       # provide a list of starting and ending dates
#                            # for survey identification
#base_station=1              # base station number
#start_end_dates=read_start_end_dates("../test_case/input_data/raw/Atest_start_end_dates.txt")
##Fill:
#data.populateSurveyDic(opt,base_station,start_end_dates)

# I.b) alternatively do not read raw data but read pre-processed data
# in this case, comment parts 0), I.a) and III).
filename="../test_case/input_data/preprocessed/gravity_data_hierarchy.txt"
#populate a Campaign object: this populates stations and loops in the data structure
data.readModifCDataFilesAndPopulateSurveyDic(filename)

#once stations and loops are populated, populate the upstream 
#structures (i.e. surveys objects and the campaign object):
# iterate through the dictionary of surveys in the campaign object
for keysurvey,survey in data.survey_dic.iteritems():   
    # for each survey, populate its channels using the loop dictionary:
    # this means that now, each survey will also possess time series (as any
    # ChannelList object: stations, loops, surveys and campaign objects are all
    # derived from ChannelList objects). These times series are simply 
    # concatenation of loops time series.
    survey.populateFromSubDictionnaries(survey.loop_dic)
    #re-put the modified survey in its place in the dic.:
    data.survey_dic[keysurvey]=survey
# then do the same with the campaign object (data): populate from the survey 
# objects.
data.populateFromSubDictionnaries(data.survey_dic)
# now the whole dataset is stored several times: as time series in the
# campaign, survey, loops and stations objects.

# IIa) solid-earth tides corrections
lat=9.742
lon=1.606
alt=450
data.tideCorrectionAgnew(lat,lon,alt)

# IIb) ocean loading corrections
#read amplitudes and phases:
fh = open("../test_case/input_data/oceantidal.txt", 'r')
test=0
i=0
amp=[]
phases=[]
# read each line in the file
for line in fh:    
    i+=1
    # Clean line
    line = line.strip()
    # Skip blank and comment lines
    if (not line) or (line == '$$ END TABLE'): continue
    # separate fields in line:
    vals=line.split()
    if len(vals) >= 3:
        if vals[2]=='GRAV':                    
            test=1
            lon=float(vals[5])
            i=0
    if test==1 and i==1:
        for j in range(11):
            amp.append(float(vals[j]))
    if test==1 and i==4:
        for j in range(11):
            phases.append(float(vals[j]))
        break
                          
data.oceanCorrectionAgnew(amp,phases,lon)


#
## III) automatic data selection:
#g_threshold=4                   # keep |g|< g_threshold where |g| is the 
#                                # absolute mean of last 3 values of time series
#                                # (µgal)
#sd_threshold=20                 # keep values where standard deviation<sd_threshold
#                                # (µgal)
#tilt_threshold=5                # keep values with |tilts|<tilt_threshold
#                                # (arcsec)
#dur_threshold=60                # keep values with duration = dur_threshold
#                                # (sec)
#
#for keysurv,surv in data.survey_dic.iteritems():
#    for keyloop,loop in data.survey_dic[keysurv].loop_dic.iteritems():
#        for keysta,sta in data.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
#            # convert to array: needed for operations
#            g=np.array(data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav)*1000
#            g=g-g[len(g)-1]
#            #mean of the last three values:
#            stabilized_values=mean(g[len(g)-3:len(g)])
#            
#            #extract time series:
#            sd=data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].sd
#            tiltsx=data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tiltx
#            tiltsy=data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tilty
#            dur=data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].dur                
#            keepdata=data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata
#
#            for i in range(len(g)):
#                if abs(g[i])>g_threshold + stabilized_values:
#                    keepdata[i]=0
#            for i in range(len(sd)):
#                if 1000*sd[i]>sd_threshold:
#                    keepdata[i]=0        
#            for i in range(len(tiltsx)):
#                if abs(tiltsx[i])>tilt_threshold or abs(tiltsy[i])>tilt_threshold :
#                        keepdata[i]=0    
#            for i in range(len(dur)):
#                if dur[i]!=dur_threshold:
#                    keepdata[i]=0      
#            #update the keepdata time series:
#            data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=keepdata
#
## finish the data selection step, if a station has no remaining time series, it
## should be discarded (i.e., its keepitem property should be set to 0)
#for keysurv,surv in data.survey_dic.iteritems():
#    for keyloop,loop in data.survey_dic[keysurv].loop_dic.iteritems():
#        for keysta,sta in data.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
#            if np.sum(sta.keepdata)==0:
#                #all measurements of a single station have been discarded
#                data.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepitem=0


# IV) Drift adjustment:
# Options:
sd_factor=1              # constant factor
sd_add=0.005            # constant value to add (mgal)
drift_t=1               # polynomial degree for temporal drift
drift_temp=0            #^polynomial degree for temperature drift
alpha=0.05              # significance level for global model test
write_out_files='n'     # writing output files
data.driftInversion(sd_factor,sd_add,drift_t,drift_temp,alpha,write_out_files,\
                        output_dir)
 

# V) Save simple difference files:
# save the data hierarchy structure:
tday=datetime.now()
Mainfile=open(output_dir+os.sep+\
'simple_diff_data_hierarchy_%4d%02d%02d_%02d%02d%02d.txt'%(tday.year,\
tday.month,tday.day,tday.hour,tday.minute,tday.second),'w')    
# Then save all simple difference files:    
# iterate through the dictionary of surveys in the campaign object
for survid,surv in data.survey_dic.iteritems():
    if surv.keepitem==1:
        survdir=output_dir+os.sep+surv.name  
        if not os.path.exists(survdir):
            #create a folder by survey
            os.makedirs(survdir)                
        filename=survdir+os.sep+\
        "SimpleDifferences_%4d%02d%02d_%02d%02d%02d.dat"%(tday.year,\
        tday.month,tday.day,tday.hour,tday.minute,tday.second)
        #write the output file:
        Mainfile.write("%s %s %s\n"%(surv.name,survid,filename))                
        surv.writeOutputValues(filename)
Mainfile.close()
               
# VI) Calculate double differences
# select a reference survey:               
referenceDate=datetime(2013,9,15).toordinal()
# Survey keys (survid) of the survey dictionnary entries are ordinal dates
survids=[survid for survid,surv in data.survey_dic.iteritems() if surv.keepitem==1]
print survids
for surv in survids:
    if surv == referenceDate:
        print "reference date"
        print surv
        referenceSurvey=data.survey_dic[int(surv)]        

# (double difference options: 1=classic 2=network mean)
DD=data.calculateDoubleDifferences(referenceSurvey,1)

#Save double differences 
filenameGrav=output_dir+os.sep+"doubledifferencesGrav"
filenameSTD=output_dir+os.sep+"doubledifferencesSTD"        
data.saveProcessedDoubledifferencesData(filenameGrav,filenameSTD,DD,referenceSurvey)
    

# Close the Qt application (only for progress bars...)               
QCoreApplication.exit()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 08 15:56:40 2014

main Program, main window for GUI

@author: Basile HECTOR

##########################
        pyGrav GUI
##########################

Version 1.0
21 jan 2016

Please cite as:
Hector, B. and Hinderer, J.: pyGrav, a Python-based program for handling and 
processing relative gravity data, Computers & Geosciences, doi:10.1016/j.cageo.2016.03.010, 2016.


To do:
- compute a function for a statistical analysis of drifts: find all station 
repetitions in a survey, calculate drift (not from adjustment), and plot vs
time, temp, or distance to base station or whatever.

- change the gravity list handled by the program in the channelList object, 
should be the tide- and pressure- correcter one (gcorr?)

- update the channel assignation when loading eterna/tsoft files

- subclass the module for updating corr_g channel when applying tides or atmospheric corrections

- currently the code does not allow several surveys in one day

"""
import sys,os,shutil,subprocess,glob
from PyQt4 import QtGui,QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from data_objects import *
from synthetic_tides import *
from copy import deepcopy
from datetime import *
import numpy as np
from matplotlib.dates import date2num,num2date
import matplotlib.dates as md

#from model_Classes import *
from model_Classes_tree_and_table import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class mainProg(QtGui.QMainWindow):
    """
    Main window definition
    
    Properties:
    - campaigndata:                     whole dataset
    - data_path:                        input data directory
    - output_root_dir:                  output data directory: where all 
                                        output files are written
    - survey_selection_option:          Code for survey selection options
                                        (1=automatic survey selection)
                                        (2=load start/end dates files)
                                        (3=load start/end dates of single survey)
    - t_threshold:                      input parameter for the automatic 
                                        survey selection method
    - base_cycling_station:             Code (1/0) for the automatic 
                                        survey selection method: 1 means base
                                        station and cycling stations have the
                                        same name (number), 0 they have 
                                        different names
    - referencesurvey:                  reference survey for double diff.         
    - DD:                               Double difference result (Campaign object)                                        

    Options properties
    - startProjectAction
    - exitAction
    - openFile
    - openModifCFile
    - saveUnprocessedAction
    - saveProcessedAction
    - saveAction1 (simple diff)
    - saveAction2 (double diff)
    - tideCorrectionAction
    - oceanLoadingCorrectionAction
    - atmosphericCorrectionAction
    - dataSelectionAction
    - driftAdjustmentAction
    - computeDoubleDifferencesAction
          
    """    
    
    campaigndata=Campaign()
    static_dataset=StaticDataSet()
    data_path='./'
    output_root_dir='./'
    survey_selection_option=0
    t_threshold=0
    referencesurvey=Survey(ChannelList(),'reference survey')
    DD=Campaign()
    
    def __init__(self):
        """
        program instanciation
        """
        super(mainProg, self).__init__()
        self.initUI()
        data_path='./'
        output_root_dir='./'        
        
    def initUI(self):      
        """
        initialize user interface
        """
        self.create_menu()
        #initialize an empty window for further display of data structure
        #needed for the closing function
        self.create_status_bar("Please load a data file")
        self.setGeometry(50, 50, 350, 300)
        self.setWindowTitle('CG5 data processing')
        self.show()

    def startProject(self):
        """
        Ask the user for input and output directories
        """
        # note: comment and uncomment for faster prog. execution during 
        #development phases
        #self.data_path='C:/Users/Basile/Documents/These/mesures_gravi/donnees_nalohou/corrections_pygrav_20013_14_15/data'
        #self.output_root_dir='C:/Users/Basile/Documents/These/mesures_gravi/donnees_nalohou/corrections_pygrav_20013_14_15/output'
        #self.data_path=str(QFileDialog.getExistingDirectory(self, "Select Input Directory",'C:/Users/Basile/Documents/These/microgravi/code/tests'))            
        #self.output_root_dir=str(QFileDialog.getExistingDirectory(self, "Select Output Directory",'C:/Users/Basile/Documents/These/microgravi/code/tests'))        
        #self.output_root_dir='C:/Users/Basile/Documents/These/microgravi/code/tests/output_temp/'
           
        self.data_path=str(QFileDialog.getExistingDirectory(self, "Select Input Directory",'./'))            
        self.output_root_dir=str(QFileDialog.getExistingDirectory(self, "Select Output Directory",'./'))        
        
        #other options from the process menu are now available:          
        self.openFile.setEnabled(True)
        self.openModifCFile.setEnabled(True)
        self.openSimpleDiffFile.setEnabled(True)
        #self.openSimpleDiffFile_static.setEnabled(True)
        
    def openRawdata(self):
        """
        - Display a file opening window
        - Populate a Campaign object: read all raw data
        - Set a new window for survey selection options
        - link selection options to apropriate functions:
          3 options are currently available: an automatic selection, a
          a selection with a user input file containing start-end dates for
          each survey, and a single srvey selection with a single start-end
          date 
          Each option calls the appropriate function:
              automaticSurveySelection()
              load_start_end_dates()
              askUserSingleSurvey()              
        which then calls the generic self.baseStationSelection() function
        which eventually calls the data populating methods within the 
        data_objects module with the appropriate options.
        """
        # open file
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.data_path)                
        campaigndata=Campaign()
        #populate a Campaign object
        campaigndata.readRawDataFile(fname)
        self.campaigndata=campaigndata
        if fname:
            #create new window and set as central Widget:
            surveySelectionWin=QtGui.QWidget()
            self.setCentralWidget(surveySelectionWin)        
            self.statusBar().showMessage("Please choose survey selection method")
            # create buttons and actions
            surveySelectionWin.btn1 = QtGui.QPushButton('automatic survey selection', self)
            surveySelectionWin.btn1.clicked.connect(self.automaticSurveySelection)
            surveySelectionWin.btn2 = QtGui.QPushButton('Load survey dates file', self)
            surveySelectionWin.btn2.clicked.connect(self.load_start_end_dates)
            surveySelectionWin.btn3 = QtGui.QPushButton('Single survey selection', self)
            surveySelectionWin.btn3.clicked.connect(self.askUserSingleSurvey)                        
            #locations                
            grid = QtGui.QGridLayout()
            grid.addWidget(surveySelectionWin.btn1,0,0,1,1)
            grid.addWidget(surveySelectionWin.btn2,1,0,1,1)  
            grid.addWidget(surveySelectionWin.btn3,2,0,1,1)          
            surveySelectionWin.setLayout(grid)   
            surveySelectionWin.setWindowTitle('Survey selections')    
            surveySelectionWin.show()      
                      
           
    def openModifCData(self):
        """
        Option for loading already processed data. 
        Data should be arranged as when saved with the appropriate save 
        processed data option, that is to say: a hierarchical file containing 
        survey names and loop names, with the following format:
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
        loops is the number of loops included in each survey 

        Then, all loops should be stored in subfolders which names are given in
        the 'survey directory' of each Survey.
        
        Loop files are classical CGxTOOL 'c' files with one more column with 1
        and 0 values if line data are kept or not
        base station should always be the same and is requested to the user
        """
        # ask for hierarchical file
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.data_path)                        
        # ask for base station number
        text, ok = QtGui.QInputDialog.getText(self, 'Input parameters', 
            'Base station?')        
        if ok:
            self.base_station=int(text)    
        campaigndata=Campaign()
        #populate a Campaign object
        campaigndata.readModifCDataFilesAndPopulateSurveyDic(fname)
        self.campaigndata=campaigndata        
        
        #once stations and loops are populated, populate the upstream 
        #structures:
        for keysurvey,survey in self.campaigndata.survey_dic.iteritems():            
            survey.populateFromSubDictionnaries(survey.loop_dic)
            self.campaigndata.survey_dic[keysurvey]=survey
        self.campaigndata.populateFromSubDictionnaries(self.campaigndata.survey_dic)
                                   
        #other options from the process menu are now available:          
        self.tideCorrectionAction.setEnabled(True)
        self.oceanLoadingCorrectionAction.setEnabled(True)
        self.atmosphericCorrectionAction.setEnabled(True)
        self.dataSelectionAction.setEnabled(True)
        self.driftAdjustmentAction.setEnabled(True)        
        #self.saveUnprocessedAction.setEnabled(True)             
        self.saveProcessedAction.setEnabled(True)        
        self.correctRecordedTimeAction.setEnabled(True)
        #emptywin=QtGui.QWidget()
        #self.setCentralWidget(emptywin)
        self.dataSelection()   
        
    def openSimpleDiff(self):
        """
        Load simple differences file(s?)
        """
        
        #create new window and set as central Widget:
        simpleDiffSelectionWin=QtGui.QWidget()
        self.setCentralWidget(simpleDiffSelectionWin)        
        self.statusBar().showMessage("Please choose simple difference loading method")
        # create buttons and actions
        simpleDiffSelectionWin.btn1 = QtGui.QPushButton('Load single simple differences file', self)
        simpleDiffSelectionWin.btn1.clicked.connect(self.openSingleSimpleDiff)
        simpleDiffSelectionWin.btn2 = QtGui.QPushButton('Load several simple differences files', self)
        simpleDiffSelectionWin.btn2.clicked.connect(self.openSeveralSimpleDiff)                       
        #locations                
        grid = QtGui.QGridLayout()
        grid.addWidget(simpleDiffSelectionWin.btn1,0,0,1,1)
        grid.addWidget(simpleDiffSelectionWin.btn2,1,0,1,1)  
        simpleDiffSelectionWin.setLayout(grid)   
        simpleDiffSelectionWin.setWindowTitle('Simple differences Loading')    
        simpleDiffSelectionWin.show() 
        #self.displaySimpleDiffAction.setEnabled(True)           
        #self.heightBouguerCorrectionAction.setEnabled(False)        
        #self.terrainCorrectionAction.setEnabled(False)            
        #self.CorrectionAction.setEnabled(True)            
        
    def openSimpleDiff_static(self):
        """
        Load simple differences file(s?) for static objects
        """
        
        #create new window and set as central Widget:
        simpleDiffSelectionWin=QtGui.QWidget()
        self.setCentralWidget(simpleDiffSelectionWin)        
        self.statusBar().showMessage("Please choose simple difference loading method")
        # create buttons and actions
        simpleDiffSelectionWin.btn1 = QtGui.QPushButton('Load single simple differences file', self)
        simpleDiffSelectionWin.btn1.clicked.connect(self.openSingleSimpleDiff_static)                     
        #locations                
        grid = QtGui.QGridLayout()
        grid.addWidget(simpleDiffSelectionWin.btn1,0,0,1,1)
        simpleDiffSelectionWin.setLayout(grid)   
        simpleDiffSelectionWin.setWindowTitle('Simple differences Loading')    
        simpleDiffSelectionWin.show() 
        self.displaySimpleDiffAction.setEnabled(True)           
        self.heightBouguerCorrectionAction.setEnabled(False)        
        self.terrainCorrectionAction.setEnabled(False)            
        self.CorrectionAction.setEnabled(True)            
        
        
    def openSingleSimpleDiff(self):
        """
        In this case, the data can be stored in any object: the name is not 
        relevant
        """

        # ask for data file
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.output_root_dir)  
        survname=datetime(1990,02,11)
        survey=Survey(ChannelList(),str(survname.date()))
        survey.readSimpleDiff(fname)
        self.campaigndata.survey_dic[survname.toordinal()]=survey
        self.computeDoubleDifferencesAction.setEnabled(True)
        self.saveAction1.setEnabled(True)
        self.statusBar().showMessage("Save options and double difference computation enabled")        
        
        
    def openSingleSimpleDiff_static(self):
        """
        In this case, the data can be stored in any object: the name is not 
        relevant
        """

        # ask for data file
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.output_root_dir)  
        survname=datetime(1990,02,11)
        staticdataset=StaticDataSet()
        staticdataset.readSimpleDiff_static(fname)
        self.static_dataset=staticdataset


    def openSeveralSimpleDiff(self):
        """
        A file containing the list of all simple difference files to read is
        required, and used to populate the output_dic of each survey present
        in the survey dictionary of the dataset.
        """
     
        # ask for hierarchical file
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.output_root_dir)            
        # open hierarchical file
        f=open(fname,'r')  
        for line in f:
            # Clean line
            line = line.strip() 
            vals=line.split()        
            #instanciate a survey object
            survtemp=Survey(ChannelList(),vals[0])
            survtemp.keepitem=1
            keysurv=int(vals[1])
            survtemp.readSimpleDiff(vals[2])
            self.campaigndata.survey_dic[keysurv]=survtemp
            self.computeDoubleDifferencesAction.setEnabled(True)
            self.saveAction1.setEnabled(True)
            self.statusBar().showMessage("Save options and double difference computation enabled")     
                       
    def automaticSurveySelection(self):
        """
        Option for automatic selection of surveys among the raw data set.
        
        Algo:
        A time threshold is asked to the user, and used as a criteria to 
        separate different survey. Time intervals between measurement dates for
        which station number changes are compared to the threshold. 
        If it is higher, a new survey is considered.
        
        The base station is asked to the user
        
        Remarks:
        Only works when the base station (asked to the user) is allways the 
        same. Whether the base station number is the same as the cycling
        station number is asked to the user.
        when complicated loop geometries are used, or for specific survey 
        designs, this option is likely to fail.
        """
        #set the survey selection option required by the survey populating 
        #function:
        self.survey_selection_option=1
        text, ok = QtGui.QInputDialog.getText(self, 'Input parameters', 
            'time threshold (hr)')        
        if ok:
            self.t_threshold=int(text)
            text, ok2 = QtGui.QInputDialog.getText(self, 'Input parameters', 
            'base station=cycling station? (1=y/0=n)')        
            if ok2:
                self.base_cycling_station=int(text)
                #call the next generic step of the survey selection process
                self.baseStationSelection()

            
    def load_start_end_dates(self):
        """
        function for loading start and end dates of each survey we want to 
        process. File format is
        yyy/mm/dd hh:mn:ss yyy/mm/dd hh:mn:ss 
        yyy/mm/dd hh:mn:ss yyy/mm/dd hh:mn:ss 
        ...
        
        update the Campaign object stored in the program by populating its
        survey dictionary

        """
        #set the survey selection option required by the survey populating 
        #function:        
        self.survey_selection_option=2
        print 'hold on a sec'
        self.statusBar().showMessage("Hold on a sec")        
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.data_path)          
        self.surveydates=read_start_end_dates(fname)        
        #call the next generic step of the survey selection process        
        self.baseStationSelection()     


    def askUserSingleSurvey(self):
        """
        ask the user for start and end dates of a single survey
        """
        self.survey_selection_option=3
        chooseSingleSurvey=QtGui.QWidget()
        self.setCentralWidget(chooseSingleSurvey)      
        self.statusBar().showMessage("Please enter start/end dates of a survey")
        
        
        chooseSingleSurvey.yr=QtGui.QLabel('year')
        chooseSingleSurvey.month=QtGui.QLabel('month')
        chooseSingleSurvey.day=QtGui.QLabel('day')
        chooseSingleSurvey.hr=QtGui.QLabel('hr')
        chooseSingleSurvey.mn=QtGui.QLabel('mn')
        chooseSingleSurvey.ss=QtGui.QLabel('ss')
        chooseSingleSurvey.yrEdit=QtGui.QLineEdit()
        chooseSingleSurvey.monthEdit=QtGui.QLineEdit()
        chooseSingleSurvey.dayEdit=QtGui.QLineEdit()
        chooseSingleSurvey.hrEdit=QtGui.QLineEdit()
        chooseSingleSurvey.mnEdit=QtGui.QLineEdit()
        chooseSingleSurvey.ssEdit=QtGui.QLineEdit()
        chooseSingleSurvey.yrEditend=QtGui.QLineEdit()
        chooseSingleSurvey.monthEditend=QtGui.QLineEdit()
        chooseSingleSurvey.dayEditend=QtGui.QLineEdit()
        chooseSingleSurvey.hrEditend=QtGui.QLineEdit()
        chooseSingleSurvey.mnEditend=QtGui.QLineEdit()
        chooseSingleSurvey.ssEditend=QtGui.QLineEdit()        
        # create buttons and actions
        chooseSingleSurvey.btn1 = QtGui.QPushButton('ok', self)
        chooseSingleSurvey.btn1.clicked.connect(lambda : self.setSingleSurvey(chooseSingleSurvey))
                          
        #locations                
        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Start date'),1,1)
        grid.addWidget(QtGui.QLabel('End date'),1,2)        
        grid.addWidget(chooseSingleSurvey.yr,2,0)
        grid.addWidget(chooseSingleSurvey.yrEdit,2,1) 
        grid.addWidget(chooseSingleSurvey.yrEditend,2,2)         
        grid.addWidget(chooseSingleSurvey.month,3,0)
        grid.addWidget(chooseSingleSurvey.monthEdit,3,1)
        grid.addWidget(chooseSingleSurvey.monthEditend,3,2)        
        grid.addWidget(chooseSingleSurvey.day,4,0)
        grid.addWidget(chooseSingleSurvey.dayEdit,4,1)
        grid.addWidget(chooseSingleSurvey.dayEditend,4,2) 
        grid.addWidget(chooseSingleSurvey.hr,5,0)
        grid.addWidget(chooseSingleSurvey.hrEdit,5,1)
        grid.addWidget(chooseSingleSurvey.hrEditend,5,2)        
        grid.addWidget(chooseSingleSurvey.mn,6,0)
        grid.addWidget(chooseSingleSurvey.mnEdit,6,1)
        grid.addWidget(chooseSingleSurvey.mnEditend,6,2)        
        grid.addWidget(chooseSingleSurvey.ss,7,0)
        grid.addWidget(chooseSingleSurvey.ssEdit,7,1)   
        grid.addWidget(chooseSingleSurvey.ssEditend,7,2)   
        grid.addWidget(chooseSingleSurvey.btn1,8,0)
        chooseSingleSurvey.setLayout(grid)   
        chooseSingleSurvey.setWindowTitle('Survey selection')    
        chooseSingleSurvey.show()             
 
        
    def setSingleSurvey(self,swin):
        """
        When user chooses to enter manually start and end dates of a single
        survey. Fill in the caimpagndata property.
        """
        self.surveydates=[(datetime(int(swin.yrEdit.text()),int(swin.monthEdit.text()),
            int(swin.dayEdit.text()),int(swin.hrEdit.text()),
            int(swin.mnEdit.text()),int(swin.ssEdit.text())),
            datetime(int(swin.yrEditend.text()),int(swin.monthEditend.text()),
            int(swin.dayEditend.text()),int(swin.hrEditend.text()),
            int(swin.mnEditend.text()),int(swin.ssEditend.text())))]
        #call the next generic step of the survey selection process               
        self.baseStationSelection()     
     
    def baseStationSelection(self):
        """
        Ask for base station number
        """
        baseStationSelectionWin=QtGui.QWidget()
        
        self.setCentralWidget(baseStationSelectionWin)        
        self.statusBar().showMessage("Please enter base station number")
        
        text, ok = QtGui.QInputDialog.getText(self, 'Input parameters', 
            'Enter base station number')        
        if ok:
            self.populateSurveyDictionnary(int(text))       
        
        ############### other option for displaying a selection screen:
#        # create buttons and actions
#        baseStationSelectionWin.baseStationEdit=QtGui.QLineEdit()
#        baseStationSelectionWin.btn1 = QtGui.QPushButton('ok', self)
#        baseStationSelectionWin.btn1.clicked.connect(lambda : self.populateSurveyDictionnary(int(baseStationSelectionWin.baseStationEdit.text())))                       
#        #locations                
#        hbox = QtGui.QHBoxLayout()
#        hbox.addWidget(baseStationSelectionWin.baseStationEdit)
#        hbox.addWidget(baseStationSelectionWin.btn1)  
#        hbox.addStretch(1)
#        vbox = QtGui.QVBoxLayout()
#        vbox.addWidget(QtGui.QLabel('Enter base station number'))
#        vbox.addLayout(hbox)
#        vbox.addStretch(1)
#        baseStationSelectionWin.setLayout(vbox)   
#        baseStationSelectionWin.setWindowTitle('Base station selection')    
#        baseStationSelectionWin.show()     

        
    def populateSurveyDictionnary(self,base_station):
        """
        call the appropriate function in the data_objects processing module, 
        with the appropriate options based on the chosen survey selection 
        method and associated parameters
        """
        self.base_station=base_station
        # if selection is based on user-defined start/end dates:
        if self.survey_selection_option==2 or self.survey_selection_option==3:
            self.campaigndata.populateSurveyDic(1,base_station,self.surveydates)
        # if selection is automatic:
        elif self.survey_selection_option==1:
            self.campaigndata.populateSurveyDic(2,base_station,None,time_threshold=self.t_threshold,base_cycling_station=self.base_cycling_station)
    
        #emptywin=QtGui.QWidget()
        #self.setCentralWidget(emptywin)
        self.dataSelection()        
        #self.statusBar().showMessage("Process some data?")     
        
        #other options from the process menu are now available:          
        self.tideCorrectionAction.setEnabled(True)
        self.oceanLoadingCorrectionAction.setEnabled(True)
        self.atmosphericCorrectionAction.setEnabled(True)
        self.dataSelectionAction.setEnabled(True)
        self.driftAdjustmentAction.setEnabled(True)        
        #self.saveUnprocessedAction.setEnabled(True)        
        self.saveProcessedAction.setEnabled(True)
        self.correctRecordedTimeAction.setEnabled(True)
        
    def tideCorrection(self):
        """
        - set the Window for tide Correction options
        - link selection options to appropriate tide correction functions
        Three correction functions are computed:
            - useCG5TideCorr: uses the CG5 internal tide correction. only work
            if appropriate survey coordinates have been input in the CG5 screen
            prior to the survey
            - usePredict: calculate a synthetic tide based on tidal parameters
            and survey coordinates. it is currently not possible to apply 
            different tide corrections to different stations or loops.
            - loadTideTimeSeries: load a synthetic tide and use it as a 
            correction
        """
        tidecorrwin=QtGui.QWidget()
        self.statusBar().showMessage("Please choose tidal correction method")

        tidecorrwin.btn1 = QtGui.QPushButton('Use CG5 tide correction', self)
        tidecorrwin.btn1.clicked.connect(self.useCG5TideCorr)        
        tidecorrwin.btn2 = QtGui.QPushButton('Use synthetic tides from predict', self)        
        tidecorrwin.btn2.clicked.connect(self.usePredict)        
        tidecorrwin.btn3 = QtGui.QPushButton('Use synthetic tides from Agnew', self)        
        tidecorrwin.btn3.clicked.connect(self.useAgnew)             
        tidecorrwin.btn4 = QtGui.QPushButton('Load time series', self)
        tidecorrwin.btn4.clicked.connect(self.loadTideTimeSeries)
                
        #button locations                
        grid = QtGui.QGridLayout()
#        grid.setSpacing(10)
        grid.addWidget(tidecorrwin.btn1,0,0,1,1)
        grid.addWidget(tidecorrwin.btn2,1,0,1,1)  
        grid.addWidget(tidecorrwin.btn3,2,0,1,1)  
        grid.addWidget(tidecorrwin.btn4,3,0,1,1)  
        tidecorrwin.setLayout(grid)   
        
        tidecorrwin.setWindowTitle('Tide correction method')    
        tidecorrwin.setGeometry(50, 50, 350, 300)        
        #tidecorrwin.show()         
        self.popup=tidecorrwin
        self.popup.show()
        
    def useCG5TideCorr(self):       
        """
        Function for using internal CG5 tide correction option.
        update the Campaign().corr_g list
        """
        PBAR1 = ProgressBar(total=len(self.campaigndata.survey_dic.keys()),textmess='surveys')   
        PBAR1.show()

        self.campaigndata.corr_g=self.campaigndata.grav
        Tides=timeSeries()
        Tides.d=deepcopy(self.campaigndata.etc)
        #should be negative to match the further application of 
        #applyTideCorrection(Tides): native CG5 etc value is not the synthetic
        # tide that is subtracted but already the correction (hence - the synt
        #tide....)
        Tides.d=[-d for d in Tides.d]
        Tides.t=deepcopy(self.campaigndata.t)        
        self.applyTideCorrection(Tides)        
        
    def loadTideTimeSeries(self):
        """
        Load time series, and apply the correction to all the data.
        time series should be either a .TSF (Tsoft) or an eterna formatted file
        It is assumed that the synthetic tide is stored in the first data
        column.
        """
        self.statusBar().showMessage("Loading synthetic tide data")        
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', 
                'home')
                
        Tides=timeSeries()
        if fname[len(fname)-4:len(fname)]==".tsf" or fname[len(fname)-4:len(fname)]==".TSF":
            #tsoftfile
            Tides.populateFromTsoftFile(fname,1)
        else:
            #assume it's eterna: JPBoy loading calculations are usually on channel 3
            Tides.populateFromEternaFile(fname,3)

        self.applyTideCorrection(Tides)
        
    def usePredict(self):
        """
        generate synthetic tides from a list of gravimetric factors for tide
        corrections
        """
        
        enterCoordinates=QtGui.QWidget()
        self.statusBar().showMessage("Please enter survey coordinates")        
        
        enterCoordinates.lat=QtGui.QLabel('latitude')
        enterCoordinates.lon=QtGui.QLabel('longitude')
        enterCoordinates.alt=QtGui.QLabel('altitude')
        enterCoordinates.latEdit=QtGui.QLineEdit()
        enterCoordinates.lonEdit=QtGui.QLineEdit()
        enterCoordinates.altEdit=QtGui.QLineEdit()    
        # create buttons and actions
        enterCoordinates.btn1 = QtGui.QPushButton('Use default tidal parameters', self)
        enterCoordinates.btn1.clicked.connect(lambda : self.launchPredict(enterCoordinates,option=1))
        enterCoordinates.btn2 = QtGui.QPushButton('Load tidal parameters file', self)
        enterCoordinates.btn2.clicked.connect(lambda : self.launchPredict(enterCoordinates,option=2))        
                  
        #locations                
        grid = QtGui.QGridLayout()     
        grid.addWidget(enterCoordinates.lat,1,0)
        grid.addWidget(enterCoordinates.latEdit,1,1) 
        grid.addWidget(enterCoordinates.lon,2,0)
        grid.addWidget(enterCoordinates.lonEdit,2,1)
        grid.addWidget(enterCoordinates.alt,3,0)
        grid.addWidget(enterCoordinates.altEdit,3,1)
        grid.addWidget(enterCoordinates.btn1,4,0)
        grid.addWidget(enterCoordinates.btn2,4,1)        
        enterCoordinates.setLayout(grid)   
        enterCoordinates.setWindowTitle('Survey coordinates')  
        self.popup=enterCoordinates
        enterCoordinates.show()          
                
    def launchPredict(self,cwin,option):
        """        
        Launch predict program: write a standard project file
        IMPORTANT the format of the tide groups definition in the .ini is
        very important; specifications can be found in the Eterna doc if
        problems are encountered here.
        """
        self.lat=float(cwin.latEdit.text())
        self.lon=float(cwin.lonEdit.text())
        self.alt=float(cwin.altEdit.text())
                
        t=self.campaigndata.t
        dur=t[len(t)-1]-t[0] #dur is a timedelta object
        dur=int(dur.days*24+dur.seconds/60/60+24    )#
        projname="PR000000.ini"
        file=open(self.output_root_dir+os.sep+projname,'w')
        file.write("TEXTHEADER= MAREES THEORIQUES\n")
        file.write("TEXTHEADER= %d (1 minute)\n"%t[0].year)
        file.write("TEXTHEADER= FILTRE GGP2 (DLAG = 17.18 sec)\n\n")
        file.write("SENSORNAME=   CG5         #earth tide sensor name\n")
        file.write("SAMPLERATE=     60         #sampling interval in seconds \n")
        file.write("STATLATITU=   %f       #stations latitude  in degree\n"%self.lat)
        file.write("STATLONITU=    %f       #stations longitude in degree\n"%self.lon)
        file.write("STATELEVAT=  %f         #stations elevation in meter\n"%self.alt)
        file.write("STATGRAVIT=    0.0\n")
        file.write("TIDALCOMPO=      0         #tidal component, see manual\n")
        file.write("TIDALPOTEN=      7         #Earth Tide Potential\n")
        file.write("AMTRUNCATE=1.D-10\n")
        file.write("INITIALEPO=%d  %02d  %02d\n"%(t[0].year,t[0].month,t[0].day))
        file.write("PREDICSPAN=%d          #time span in hours\n "%dur)
        file.write("PRINTDEVEL=      1        #ANALYZE print param. for tidal development (1=yes)\n")
        file.write("SEARDATLIM=     -1        #ANALYZE search for data error threshold \n")
        file.write("NUMHIGPASS=      0        #ANALYZE highpass filtering = 1 \n")
        file.write("PRINTOBSER=      0        #ANALYZE print parameter for observations (1=yes)\n")
        file.write("RIGIDEARTH=      0        #ANALYZE parameter for rigid earth model (1=yes)\n")
        file.write("HANNWINDOW=      0        #ANALYZE parameter for Hann-window (1=yes)\n")
        file.write("QUICKLOOKA=      0        #ANALYZE parameter for quick look analysis (1=yes)\n")
        file.write("POLETIDCOR=      0        #ANALYZE parameter for pole corrections (1=yes)\n")
        file.write("LODTIDECOR=      0        #ANALYZE parameter for LOD corrections (1=yes)\n")
        file.write("STORENEQSY=      1\n\n")

        if option==1:
            #use default tidal parameters: load file
            tides=open("./eterna_files/200D.INI",'r')
            for line in tides:
                file.write(line)
                
            tides.close()
        if option==2:
            #Ask the user to provide tide files
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Load tidal parameters file', 
                self.data_path)        
            tides=open(fname,'r')
            for line in tides:
                # Clean line
                line = line.strip()      
                vals=line.split()
                file.write("TIDALPARAM=")
                file.write("%10.6f%10.6f%10.6f%10.6f %4s "%(float(vals[0]),
                float(vals[1]),float(vals[2]),float(vals[3]),vals[4]))
                file.write(" #ANALYZE wave group\n")                
            tides.close()                
        file.close()    
        file2=open(self.output_root_dir+os.sep+"project",'w')      
        file2.write(projname)
        file2.close()

        ### cp predict.exe to output directory
        shutil.copyfile("./eterna_files/predict.exe",self.output_root_dir+os.sep+"predict.exe")
               
        # run predict
        curr_dir=os.getcwd()
        os.chdir(self.output_root_dir)
        subprocess.call(["predict.exe"])
        os.chdir(curr_dir)      
                
        #load results and apply tide corrections
        Tides=timeSeries()
        Tides.populateFromEternaFile(self.output_root_dir+os.sep+"PR000000.prd",1)        
    
        self.campaigndata.tideCorrectionPredict(Tides)        
        
        self.popup.close()
        self.statusBar().showMessage("Process some data?")
                

    def useAgnew(self):
        """
        generate synthetic tides from the Agnew approach
        see 
        - Agnew, D.C., 2007, 3.06 - Earth Tides, in Schubert, G. ed., 
        Treatise on Geophysics, Amsterdam, Elsevier, p. 163–195.
        - Agnew, D.C., 2012, SPOTL: Some Programs for Ocean-Tide Loading
        """
        
        enterCoordinates=QtGui.QWidget()
        
#        self.setCentralWidget(enterCoordinates)      
        self.statusBar().showMessage("Please enter survey coordinates")        
        
        enterCoordinates.lat=QtGui.QLabel('latitude')
        enterCoordinates.lon=QtGui.QLabel('longitude')
        enterCoordinates.alt=QtGui.QLabel('altitude')
        enterCoordinates.latEdit=QtGui.QLineEdit()
        enterCoordinates.lonEdit=QtGui.QLineEdit()
        enterCoordinates.altEdit=QtGui.QLineEdit()    
        # create buttons and actions
        enterCoordinates.btn1 = QtGui.QPushButton('OK', self)
        enterCoordinates.btn1.clicked.connect(lambda : self.launchAgnew(enterCoordinates,option=1))
        #enterCoordinates.btn2 = QtGui.QPushButton('Load tidal parameters file', self)
        #enterCoordinates.btn2.clicked.connect(lambda : self.launchPredict(enterCoordinates,option=2))        
                  
        #locations                
        grid = QtGui.QGridLayout()     
        grid.addWidget(enterCoordinates.lat,1,0)
        grid.addWidget(enterCoordinates.latEdit,1,1) 
        grid.addWidget(enterCoordinates.lon,2,0)
        grid.addWidget(enterCoordinates.lonEdit,2,1)
        grid.addWidget(enterCoordinates.alt,3,0)
        grid.addWidget(enterCoordinates.altEdit,3,1)
        grid.addWidget(enterCoordinates.btn1,4,0)
        #grid.addWidget(enterCoordinates.btn2,4,1)        
        enterCoordinates.setLayout(grid)   
        enterCoordinates.setWindowTitle('Survey coordinates')    
        enterCoordinates.show()
        self.popup=enterCoordinates
        self.popup.show()
        
        
    def launchAgnew(self,cwin,option):
        """        
        Launch the Agnew tide correction
        """
        self.lat=float(cwin.latEdit.text())
        self.lon=float(cwin.lonEdit.text())
        self.alt=float(cwin.altEdit.text())
        
        #apply tide correction
        self.campaigndata.tideCorrectionAgnew(self.lat,self.lon,self.alt)
                        
        self.popup.close()
        self.statusBar().showMessage("Process some data?")


    def oceanLoadingCorrection(self):
        """
        get the ocean Loading tidal parameters, compute the correction and 
        apply to the dataset
        
        read amplitude and phases of ocean loading from standard BLQ file as 
        obtained from http://holt.oso.chalmers.se/loading/ for a single station.         
        
        amplitudes & phases: lists with amplitudes and phases. List order is:
        M2,S2,K1,O1,N2,P1,K2,Q1,Mf,Mm,Ssa      
        (amp,phases)
        lon:     site longitude              
        """

        # open file
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open tidal parameters file (BLQ)',self.data_path) 
        
        #read amplitudes and phases:
        fh = open(fname, 'r')
        test=0
        i=0
        amp=[]
        phases=[]
        for line in fh:    
            i+=1
            # Clean line
            line = line.strip()
            # Skip blank and comment lines
            if (not line) or (line == '$$ END TABLE'): continue
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
                                  
        self.campaigndata.oceanCorrectionAgnew(amp,phases,lon)
        
        
       
        
    def atmosphericCorrection(self):
        """
        load a time series and apply the correction
        time series should be either a .TSF (Tsoft) or an eterna formatted file
        It is assumed that the synthetic tide is stored in the first data
        column.
        populate/update the corrg field of all gravity object within 
        self.campaigndata and same for grav field
        IMPORTANT: atmospheric correction is currently not stored in the data
        structure, nor saved in an output file, this could be updated
        """
        self.statusBar().showMessage("Loading Atmospheric loading time series")        
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', 
                'home')
                        
        Atm=timeSeries()
        if fname[len(fname)-4:len(fname)]==".tsf" or fname[len(fname)-4:len(fname)]==".TSF":
            #tsoftfile
            Atm.populateFromTsoftFile(fname,1)
        else:
            #assume it's eterna
            Atm.populateFromEternaFile(fname,1)

        #if data is loaded from already processed data, self.campaigndata time
        #series are not populated, and correction should not be applied at such
        #levels
        if self.campaigndata.t:
            Atm.interpolateOnGivenTimes(self.campaigndata.t)           
            print len(self.campaigndata.corr_g)
            print len(Atm.d)
               
        self.campaigndata.corr_g=self.campaigndata.grav-Atm.d
        self.campaigndata.grav=self.campaigndata.corr_g
        
        i=1
        PBAR1 = ProgressBar(total=len(self.campaigndata.survey_dic.keys()),textmess='surveys')   
        PBAR1.show()        
        for keysurv,surv in self.campaigndata.survey_dic.iteritems():
            PBAR1.progressbar.setValue(i)
            i=i+1
            Atmtemp=deepcopy(Atm)
            if self.campaigndata.survey_dic[keysurv].t:
                Atmtemp.interpolateOnGivenTimes(self.campaigndata.survey_dic[keysurv].t)   
                self.campaigndata.survey_dic[keysurv].corr_g=self.campaigndata.survey_dic[keysurv].grav-Atmtemp.d 
                self.campaigndata.survey_dic[keysurv].grav=self.campaigndata.survey_dic[keysurv].corr_g
            PBAR2 = ProgressBar(total=len(self.campaigndata.survey_dic[keysurv].loop_dic.keys()),textmess='loops')   
            PBAR2.show()
            j=1
            for keyloop,loop in self.campaigndata.survey_dic[keysurv].loop_dic.iteritems():
                PBAR2.progressbar.setValue(j)
                j=j+1
                Atmtemp2=deepcopy(Atmtemp)
                if self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].t:
                    Atmtemp2.interpolateOnGivenTimes(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].t)   
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].corr_g=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].grav-Atmtemp2.d   
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].grav=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].corr_g
                PBAR3 = ProgressBar(total=len(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic.keys()),textmess='stations')   
                PBAR3.show()
                k=1
                for keysta,sta in self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
                    PBAR3.progressbar.setValue(k)
                    k=k+1
                    Atmtemp3=deepcopy(Atmtemp2)
                    Atmtemp3.interpolateOnGivenTimes(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t)   
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav-Atmtemp3.d   
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g   

        #emptywin=QtGui.QWidget()
        #self.setCentralWidget(emptywin)
        self.statusBar().showMessage("Process some data?")                        
    
    def dataSelection(self):
        """
        station data selection
        This is the main window of the program:

        on the left panel is a treeview of the stored data.
        TreeWidgetItems are defined in another module package (model_Classes_tree_and_table.py)
        Each Item has keysurv, keyloop (set to None for survey objects), and
        keysta (set to None for survey and loop objects) properties, which are
        their successive keys in the 
        Campaign.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta]-like
        hierarchy. This allows to easily access the original data when items on
        the tree view are checked/unchecked. Checking a checkbox changes the
        keepitem property (1 or 0) of the concerned object. This further 
        controls the data to be saved or used for drift adjustment.
        
        On the middle panel is a tableview of the selected station (within 
        a survey/loop, clicked on the left panel). Table items are also defined
        in another module package (model_Classes_tree_and_table.py)
        Each Item has keysurv, keyloop (set to None for survey objects), and
        keysta (set to None for survey and loop objects) properties, which are
        their successive keys in the 
        Campaign.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta]-like
        hierarchy. This allows to easily access the original data when items on
        the tree view are checked/unchecked. Checking a checkbox changes the
        keepdata property (1 or 0) of the given row. This further 
        control the data to be used for drift adjustment  
        
        On the right panel is a matplotlib-type figure of gravity values,std,
        tiltx and tilty channels. The button update plots has to be clicked to
        trigger the update plot function.
        
        Some options for automatic selection of data are also available, either
        to select data among the whole data set, or among the current displayed
        table.        
        
        To fix: 
        - clicking the table updates the plot
        
        """
        #### define the data selection window Widget:
        self.dataselectionwin=QtGui.QWidget()        
        mess="check surveys, loops and stations to keep for processing/saving (left panel) and select data on each station (middle panel). click the ""update plots"" button to refresh."
        self.statusBar().showMessage(mess)         
        # and the final layout:
        layout_final = QtGui.QGridLayout() 


        #### now define display behaviour for all panels       
        
        # Left panel: tree with data hierarchy (surveys, loops, stations)     
        #initialization: set the display to the first station
        keysurv=self.campaigndata.survey_dic.keys()
        keyloop=self.campaigndata.survey_dic[keysurv[0]].loop_dic.keys()
        keysta=self.campaigndata.survey_dic[keysurv[0]].loop_dic[keyloop[0]].station_dic.keys()
        #instanciate a MyTree model object:    
        self.tree = MyTree(self.campaigndata,parent=self)
        self.tree.connect(self.tree, QtCore.SIGNAL('itemClicked(QTreeWidgetItem*, int)'), self.onClick)
        self.tree.connect(self.tree, QtCore.SIGNAL('itemChanged(QTreeWidgetItem*, int)'), self.onChange)

        # center panel: table (station values)
        #instanciate a stationDataTableModel model object
        self.tablemodel = stationDataTableModel(self.campaigndata.survey_dic[keysurv[0]].loop_dic[keyloop[0]].station_dic[keysta[0]], self.dataselectionwin,keysurv[0],keyloop[0],keysta[0],)
        self.tableview = QtGui.QTableView()
        self.tableview.setModel(self.tablemodel)  
        self.tableview.resizeColumnsToContents()
        
        # right panel: plot and some options
        #figure
        self.main_frame = QWidget()        
        self.dpi = 100
        self.fig = Figure((3.0, 2.0), dpi=self.dpi)
        # make some room for the axes labels (set more place between subplots)
        self.fig.subplots_adjust(right=0.95,wspace=0.3,hspace=0.35)
        
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)        
        self.axes_grav = self.fig.add_subplot(221)
        self.axes_tiltx = self.fig.add_subplot(222)
        self.axes_tilty = self.fig.add_subplot(224)
        self.axes_temp = self.fig.add_subplot(223)
        
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        # define some buttons with actions (signal/slot events) and line edits
        updateplot_button = QtGui.QPushButton("&update plots",self)
        updateplot_button.clicked.connect(self.updatePlotTemp)
        checkselected_button = QtGui.QPushButton("&check selected",self)
        checkselected_button.clicked.connect(self.checkSelected)  
        uncheckselected_button = QtGui.QPushButton("&uncheck selected",self)
        uncheckselected_button.clicked.connect(self.uncheckSelected)        
        checkall_button = QtGui.QPushButton("&check all",self)
        checkall_button.clicked.connect(self.checkall)
        uncheckall_button = QtGui.QPushButton("&uncheck all",self)
        uncheckall_button.clicked.connect(self.uncheckall)  
        #apply_tide_corr_button = QtGui.QPushButton("&tide correction",self)
        #apply_tide_corr_button.clicked.connect(self.applyTideCorrSingleStation)          
        autoselec_tilts=QtGui.QWidget()        
        autoselec_sd=QtGui.QWidget()         
        autoselec_grav=QtGui.QWidget()   
        autoselec_dur=QtGui.QWidget()              
        autoselec_all=QtGui.QWidget()                    
        autoselec_tilts.button = QtGui.QPushButton("&auto uncheck tilts >",self)
        autoselec_tilts.button.clicked.connect(lambda : self.autoselectFct_tilt(autoselec_tilts))           
        autoselec_sd.button = QtGui.QPushButton("&auto uncheck SD >",self)
        autoselec_sd.button.clicked.connect(lambda : self.autoselectFct_sd(autoselec_sd))        
        autoselec_grav.button = QtGui.QPushButton("&auto uncheck g >",self)
        autoselec_grav.button.clicked.connect(lambda : self.autoselectFct_grav(autoselec_grav))   
        autoselec_dur.button = QtGui.QPushButton("&auto uncheck dur <>",self)
        autoselec_dur.button.clicked.connect(lambda : self.autoselectFct_dur(autoselec_dur))          
        autoselec_all.button = QtGui.QPushButton("&apply to all data",self)
        autoselec_all.button.clicked.connect(lambda : self.autoselectFct_all(autoselec_tilts.val,autoselec_sd.val,autoselec_grav.val,autoselec_dur.val))            
        autoselec_tilts.val=QtGui.QLineEdit()
        autoselec_sd.val=QtGui.QLineEdit()
        autoselec_grav.val=QtGui.QLineEdit()
        autoselec_dur.val=QtGui.QLineEdit()
        
        validate_button=QtGui.QPushButton("\n OK \n",self)
        validate_button.clicked.connect(self.finishDataSelection)
        
        #### add widgets to the layout
       
        # add left panel (tree)
        layout_final.addWidget(self.tree,0,0,1,1)
        
        # add center panel (table)
        layout_final.addWidget(self.tableview,0,1,1,1) 
        
        # for the right panel (options & display):        
        # create sublayouts (allow the line edit to be of finite extent)
        layout_options = QtGui.QGridLayout()         
        grid = QtGui.QGridLayout() 
        
        # fill subplayouts:
        grid.addWidget(self.canvas,3,0,19,4)
        grid.addWidget(self.mpl_toolbar,2,0,1,4)
        layout_options.addWidget(checkselected_button,0,0)
        layout_options.addWidget(uncheckselected_button,0,1)        
        layout_options.addWidget(autoselec_tilts.button,0,2)
        layout_options.addWidget(autoselec_tilts.val,0,3,1,1)    
        layout_options.addWidget(autoselec_grav.button,0,4)
        layout_options.addWidget(autoselec_grav.val,0,5,1,1)     
        layout_options.addWidget(autoselec_dur.button,1,4)
        layout_options.addWidget(autoselec_dur.val,1,5,1,1)             
        layout_options.addWidget(checkall_button,1,0)
        layout_options.addWidget(uncheckall_button,1,1)
        layout_options.addWidget(autoselec_sd.button,1,2)               
        layout_options.addWidget(autoselec_sd.val,1,3,1,1) 
        layout_options.addWidget(updateplot_button,0,6)        
        layout_options.addWidget(autoselec_all.button,1,6)
        layout_options.addWidget(validate_button,0,7,2,1)
        #layout_options.addWidget(apply_tide_corr_button,0,6)
        # add subplayouts to the main layout
        grid.addLayout(layout_options,0,0,1,2)
        layout_final.addLayout(grid,0,2,1,1)  
        #grid.setColumnMinimumWidth(1,750)       
        layout_final.setColumnMinimumWidth(0,200)
        layout_final.setColumnMinimumWidth(1,400)
        layout_final.setColumnMinimumWidth(2,700)

        #### set window geometry and set layout as the main layout      
        
        self.setGeometry(50, 50, 1300, 700)
        self.setWindowTitle('Data selection')    
        self.dataselectionwin.setLayout(layout_final)               
        self.setCentralWidget(self.dataselectionwin)
        self.dataselectionwin.show()       
        
        
    def autoselectFct_tilt(self,autoselec):
        """
        function for automatic selection of data based on simple thresholds:
        Tilts: absolute value higher than threshold are set to keepdata=0
        """
        tilt_threshold=float(autoselec.val.text())

        keysurv=self.tablemodel.keysurv
        keyloop=self.tablemodel.keyloop
        keysta=self.tablemodel.keysta
        
        tiltsx=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tiltx
        tiltsy=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tilty
        
         #keepdata=[1 for i in tiltsx]
        keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata
        for i in range(len(tiltsx)):
            if abs(tiltsx[i])>tilt_threshold or abs(tiltsy[i])>tilt_threshold :
                keepdata[i]=0
        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=keepdata
        self.updatePlotTemp()

    def autoselectFct_sd(self,autoselec):
        """
        function for automatic selection of data based on simple thresholds
        sd: sd values higher than threshold are set to keepdata=0
        """
        sd_threshold=float(autoselec.val.text())

        keysurv=self.tablemodel.keysurv
        keyloop=self.tablemodel.keyloop
        keysta=self.tablemodel.keysta
        
        sd=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].sd
        
        #keepdata=[1 for i in sd]
        keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata
        
        for i in range(len(sd)):
            if 1000*sd[i]>sd_threshold :
                keepdata[i]=0
        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=keepdata
        self.updatePlotTemp()
                
    def autoselectFct_dur(self,autoselec):
        """
        function for automatic selection of data based on simple thresholds
        dur: duration different than given value are set to keepdata=0
        """
        dur_threshold=float(autoselec.val.text())

        keysurv=self.tablemodel.keysurv
        keyloop=self.tablemodel.keyloop
        keysta=self.tablemodel.keysta
        
        dur=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].dur
        keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata
        
        for i in range(len(dur)):
            if dur[i]!=dur_threshold :
                keepdata[i]=0
        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=keepdata
        self.updatePlotTemp()                
                
    def autoselectFct_grav(self,autoselec):
        """
        function for automatic selection of data based on simple thresholds
        grav: absolute values higher than threshold offset with respect to the
        mean value from 3 last points are set to keepdata=0
        """
        g_threshold=float(autoselec.val.text())

        keysurv=self.tablemodel.keysurv
        keyloop=self.tablemodel.keyloop
        keysta=self.tablemodel.keysta
        
        # convert to array: needed for operations
        g=np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav)*1000
        g=g-g[len(g)-1]
        keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata
        #mean of the last three values:
        stabilized_values=np.mean(g[len(g)-3:len(g)])
        for i in range(len(g)):
            if abs(g[i])>g_threshold + stabilized_values:
                keepdata[i]=0
        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=keepdata
        self.updatePlotTemp()
               
    def autoselectFct_all(self,tilts_thrshld,sd_thrshld,grav_thrshld,dur_thrshld):
        """
        function for automatic selection of data based on simple thresholds   
        apply all selection critera which have been input by the user to all 
        the data set.
        """
        selec_grav='no';g_threshold=0
        selec_sd='no';sd_threshold=0
        selec_tilts='no';tilt_threshold=0
        selec_dur='no';dur_threshold=0
        if grav_thrshld.text():
            g_threshold=float(grav_thrshld.text());selec_grav='ok'                        
        if sd_thrshld.text():
            sd_threshold=float(sd_thrshld.text());selec_sd='ok'
        if tilts_thrshld.text():            
            tilt_threshold=float(tilts_thrshld.text());selec_tilts='ok'
        if dur_thrshld.text():            
            dur_threshold=float(dur_thrshld.text());selec_dur='ok'         
            
        for keysurv,surv in self.campaigndata.survey_dic.iteritems():
            for keyloop,loop in self.campaigndata.survey_dic[keysurv].loop_dic.iteritems():
                for keysta,sta in self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
                    # convert to array: needed for operations
                    g=np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav)*1000
                    g=g-g[len(g)-1]
                    #mean of the last three values:
                    stabilized_values=np.mean(g[len(g)-3:len(g)])
                    sd=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].sd
                    tiltsx=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tiltx
                    tiltsy=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tilty
                    dur=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].dur                
                    keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata

                    for i in range(len(g)):
                        if abs(g[i])>g_threshold + stabilized_values and selec_grav=='ok':
                            keepdata[i]=0
                    for i in range(len(sd)):
                        if 1000*sd[i]>sd_threshold and selec_sd=='ok':
                            keepdata[i]=0        
                    for i in range(len(tiltsx)):
                        if selec_tilts=='ok':
                            if abs(tiltsx[i])>tilt_threshold or abs(tiltsy[i])>tilt_threshold :
                                keepdata[i]=0    
                    for i in range(len(dur)):
                        if dur[i]!=dur_threshold and selec_dur=='ok':
                            keepdata[i]=0      
                        
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=keepdata
        self.updatePlotTemp()               
        
    def uncheckall(self):
        """
        uncheck all items in the displayed table when button is clicked
        """
        keysurv=self.tablemodel.keysurv
        keyloop=self.tablemodel.keyloop
        keysta=self.tablemodel.keysta
        keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata
        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=[0 for i in keepdata]
        # essai:        
        self.tablemodel.ChannelList_obj.keepdata=[0 for i in keepdata]
        
        self.updatePlotTemp()
        
    def checkall(self):        
        """
        check all items in the displayed table when button is clicked
        """        
        keysurv=self.tablemodel.keysurv
        keyloop=self.tablemodel.keyloop
        keysta=self.tablemodel.keysta
        keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata
        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata=[1 for i in keepdata]
        self.updatePlotTemp()        
    
    def checkSelected(self):
        """
        check all selected items (select the station column)
        """
        qmodelindex = self.tableview.selectedIndexes()
        for indx in list(qmodelindex):
            # because the setData function currently only works for column 0,
            # then pass the index of the element with same row as selected (if
            # the user does not select element from first column) and column=0
            indx2=indx.sibling(indx.row(),0)
            self.tablemodel.setData(indx2,QtCore.Qt.Checked,QtCore.Qt.CheckStateRole)
        self.updatePlotTemp()
            
    def uncheckSelected(self):
        """
        check all selected items (select the station column)
        """
        qmodelindex = self.tableview.selectedIndexes()
        for indx in list(qmodelindex):
            # because the setData function currently only works for column 0,
            # then pass the index of the element with same row as selected (if
            # the user does not select element from first column) and column=0
            indx2=indx.sibling(indx.row(),0)
            self.tablemodel.setData(indx2,QtCore.Qt.Unchecked,QtCore.Qt.CheckStateRole)
        self.updatePlotTemp()
        
    def onChange(self,item,column):
        """
        definition of the function triggered by a change in the tree:
        update the keepitem status of the considered object (survey, loop or
        station)
        """
        if item.keysurv is not None and item.keyloop is not None and item.keysta is not None:  
            # it's a station
            if not item.checkState(0): 
                self.campaigndata.survey_dic[item.keysurv].loop_dic[item.keyloop].station_dic[item.keysta].keepitem=0            
                print "station unchecked"
            else:   
                self.campaigndata.survey_dic[item.keysurv].loop_dic[item.keyloop].station_dic[item.keysta].keepitem=1            
                print "station checked"
        elif item.keysurv is not None and item.keyloop is not None and item.keysta is None:  
            # it's a loop
            if not item.checkState(0): 
                self.campaigndata.survey_dic[item.keysurv].loop_dic[item.keyloop].keepitem=0
                print "loop unchecked"            
            else:   
                self.campaigndata.survey_dic[item.keysurv].loop_dic[item.keyloop].keepitem=1
                print "loop checked"
        elif item.keysurv is not None and item.keyloop is None and item.keysta is None:  
            # it's a survey
            if not item.checkState(0): 
                self.campaigndata.survey_dic[item.keysurv].keepitem=0
                print "survey unchecked"            
            else:   
                self.campaigndata.survey_dic[item.keysurv].keepitem=1
                print "survey checked"
                                          
    def onClick(self, item, column):
        """
        definition of the function triggered by clicking a station on the tree:
        update the table on the middle panel, and the plots to show the chosen
        station)
        """        
        if item.keysurv is not None and item.keyloop is not None and item.keysta is not None:            
            self.updatePlot(item.keysurv,item.keyloop,item.keysta)

    def updatePlotTemp(self):
        """
        temp function because slot functions from .connect() does not take any
        argument
        
        to be fixed: use a lambda formulation
        """
        self.updatePlot(self.tablemodel.keysurv,self.tablemodel.keyloop,self.tablemodel.keysta)    
        
    def updatePlot(self,keysurv,keyloop,keysta): 
        """
        actual function for plot updates
        plot series for keepdata indices = 1
        """
        self.tablemodel.keysurv=keysurv
        self.tablemodel.keyloop=keyloop
        self.tablemodel.keysta=keysta
        self.tablemodel.layoutAboutToBeChanged.emit()
        self.tablemodel.createArrayData(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta])        
        self.tablemodel.layoutChanged.emit()
        
        t=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t        
        keepdata=self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepdata        
        t_selec=[t[i] for i in range(len(t)) if keepdata[i]==1]
        # gravity channel (convert to microgals for display)
        series = np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav)*1000        
        series=series-series[len(series)-1]
        series_selec=[series[i] for i in range(len(series)) if keepdata[i]==1]
        self.setPlot(self.axes_grav,t,series,t_selec,series_selec,'gravity','$\mu$gal')
        #tiltx channel
        series = np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tiltx)        
        series_selec=[series[i] for i in range(len(series)) if keepdata[i]==1]            
        self.setPlot(self.axes_tiltx,t,series,t_selec,series_selec,'tilt X','arcsec')
        #tilty channel
        series = np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].tilty)        
        series_selec=[series[i] for i in range(len(series)) if keepdata[i]==1]
        self.setPlot(self.axes_tilty,t,series,t_selec,series_selec,'tilt Y','arcsec')
        #SD channel
        series = np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].sd)*1000 
        series_selec=[series[i] for i in range(len(series)) if keepdata[i]==1]
        self.setPlot(self.axes_temp,t,series,t_selec,series_selec,'standard deviation','$\mu$gal')            
        return True     

    def setPlot(self,axe,seriex,seriey,seriex_selec,seriey_selec,serie_type,serie_unit):
        """
        plot a single station
        """
        axe.clear()
        axe.grid(True)
        if serie_type=='gravity' and seriey_selec:
            mean_g=np.mean(seriey_selec)
            axe.plot([seriex[0],seriex[len(seriex)-1]],[mean_g,mean_g],'o-',color='b',label=serie_type)        
            
        axe.plot(seriex,seriey,'o-',color='k',label=serie_type)
        axe.plot(seriex_selec,seriey_selec,'o-',color='b',label=serie_type)            
        axe.set_ylabel(serie_unit, size='x-small')
        axe.set_title(serie_type, size='x-small')
        labels = axe.get_xticklabels() + axe.get_yticklabels()
        for label in labels:
            label.set_size('x-small') 
        xfmt = md.DateFormatter('%H:%M')
        axe.xaxis.set_major_formatter(xfmt)            
        plt.setp(axe.get_xticklabels(), rotation=30, horizontalalignment='right')              
        self.canvas.draw()
        
    def finishDataSelection(self):
        """
        finish the data selection step, integrate all changes not dynamically 
        integrated (in the interaction of model/view vs data set), and get back
        to main window
        """
        
        for keysurv,surv in self.campaigndata.survey_dic.iteritems():
            for keyloop,loop in self.campaigndata.survey_dic[keysurv].loop_dic.iteritems():
                for keysta,sta in self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
                    if np.sum(sta.keepdata)==0:
                        #all measurements of a single station have been discarded
                        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].keepitem=0
        
        
        #tmp=QtGui.QWidget()
        #tmp.setGeometry(50, 50, 350, 300)
        #self.setCentralWidget(tmp)      
        self.statusBar().showMessage("Data selected")                
        
    def driftAdjustment(self):
        """
        Drift adjustment function
        """          
        driftAdjWin=QtGui.QWidget()
        driftAdjWin.setGeometry(50, 50, 350, 300)

        self.statusBar().showMessage("Please choose drift adjustment method")

        driftAdjWin.btn1 = QtGui.QPushButton('Use MCGravi', self)
        driftAdjWin.btn1.clicked.connect(self.useMCGraviDriftAdj)        
        driftAdjWin.btn2 = QtGui.QPushButton('Use datum-free least-square inversion', self)
        driftAdjWin.btn2.clicked.connect(self.useLSDriftAdj)               
                
        #button locations                
        grid = QtGui.QGridLayout()
#        grid.setSpacing(10)
        grid.addWidget(driftAdjWin.btn1,0,0,1,1)
        grid.addWidget(driftAdjWin.btn2,1,0,1,1)  
        driftAdjWin.setLayout(grid)   
        
        driftAdjWin.setWindowTitle('Drift Adjustment method')    
        #driftAdjWin.show()       
        self.popup=driftAdjWin
        self.popup.show()
      
    def useMCGraviDriftAdj(self):
        """
        Use MCGravi for drift adjustment (Beilin 2006):
        includes the datum-free and weighted constraint LS adjustment from
        Hwang et al (2002)
        """
        datafordriftadj=deepcopy(self.campaigndata)
        self.selectMCGraviOptions(datafordriftadj,self.output_root_dir)

        
    def useLSDriftAdj(self):
        """
        Use Datum-free Least Square inversion from Hwang et al (2002)
        """
        LSOptions=QtGui.QWidget()
        
        LSOptions.setGeometry(50, 50, 350, 300)
        self.statusBar().showMessage("enter LS options")
        
        LSOptions.drift_t=QtGui.QLabel('temporal drift polynomial?')
        LSOptions.drift_k=QtGui.QLabel('temperature drift polynomial?')        
        LSOptions.sigma_factor=QtGui.QLabel('SD factor to data')                
        LSOptions.sigma_add=QtGui.QLabel('SD add to data (mgal)')        
        LSOptions.alpha=QtGui.QLabel('Significance level for global model test')        
        LSOptions.woutfiles=QtGui.QLabel('Write Output Files (y/n)?')        
       
        LSOptions.drift_tEdit=QtGui.QLineEdit('1')
        LSOptions.drift_kEdit=QtGui.QLineEdit('0')
        LSOptions.sigma_factorEdit=QtGui.QLineEdit('1')
        LSOptions.sigma_addEdit=QtGui.QLineEdit('0.005')
        LSOptions.alphaEdit=QtGui.QLineEdit('0.05')
        LSOptions.woutfilesEdit=QtGui.QLineEdit('n')
             
        # create buttons and actions
        LSOptions.btn1 = QtGui.QPushButton('ok', self)
        LSOptions.btn1.clicked.connect(lambda : self.lsDriftAdj(LSOptions))
                          
        #locations                
        grid = QtGui.QGridLayout()       
        grid.addWidget(LSOptions.drift_t,1,0)
        grid.addWidget(LSOptions.drift_tEdit,1,1)
        grid.addWidget(LSOptions.drift_k,2,0)
        grid.addWidget(LSOptions.drift_kEdit,2,1)
        grid.addWidget(LSOptions.sigma_factor,3,0)
        grid.addWidget(LSOptions.sigma_factorEdit,3,1)       
        grid.addWidget(LSOptions.sigma_add,4,0)
        grid.addWidget(LSOptions.sigma_addEdit,4,1)
        grid.addWidget(LSOptions.alpha,5,0)
        grid.addWidget(LSOptions.alphaEdit,5,1) 
        grid.addWidget(LSOptions.woutfiles,6,0)
        grid.addWidget(LSOptions.woutfilesEdit,6,1)         
        grid.addWidget(LSOptions.btn1,7,0)  
        
        LSOptions.setLayout(grid)   
        LSOptions.setWindowTitle('LS options')    
        self.popup=LSOptions
        self.popup.show()
                             
    def lsDriftAdj(self,LSWin):
        """
        Do the LS drift adjustment    
        """
        self.campaigndata.driftInversion(float(LSWin.sigma_factorEdit.text()),\
        float(LSWin.sigma_addEdit.text()),int(LSWin.drift_tEdit.text()),\
        int(LSWin.drift_kEdit.text()),float(LSWin.alphaEdit.text()),\
        LSWin.woutfilesEdit.text(),self.output_root_dir)
 
        
#        self.displaySimpleDiffAction.setEnabled(True)
#        self.heightBouguerCorrectionAction.setEnabled(False)        
#        self.terrainCorrectionAction.setEnabled(False)           
        self.computeDoubleDifferencesAction.setEnabled(True)
        self.saveAction1.setEnabled(True)
        self.popup.close()      
        self.statusBar().showMessage("Drift adjusted. Save options and double difference computation enabled")        
        
    def saveUnprocessedData(self):
        """
        only save c files
        Copied from driftAdjustment hierarchy
        
        Can be modified, for other format definitions, or for instance saving
        all data with an extra column with 1 or 0 depending on keepdata status
        """
        datafordriftadj=deepcopy(self.campaigndata)       
        if not os.path.exists(self.output_root_dir):
            os.makedirs(self.output_root_dir)
            
        for survid,surv in datafordriftadj.survey_dic.iteritems():
            if surv.keepitem==1:
                survdir=self.output_root_dir+os.sep+surv.name
                if not os.path.exists(survdir):
                    os.makedirs(survdir)
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
                        filenameC=survdir+os.sep+"fn"+str(j)+"11c"+yr[2:4]+"."+jday
                        file_list.append("fn"+str(j)+"11c"+yr[2:4]+"."+jday)
                        # getbase station: temporary formulation, should add a checkbox 
                        # in the tree object to select it manually for each loop
                        basestationkeys=[stakeys[i] for i in range(len(stakeys)) if stakeys[i][1]==2]
                        loop.writeCgxCFile(filenameC,basestationkeys[0])


    def saveTemporaryProcessedData(self) :
        """
        """
        self.campaigndata.saveProcessedData(self.output_root_dir)
        
    def selectMCGraviOptions(self,datafordriftadj,output_root_dir):
        """

        - Other problem: extra files are needed (should maybe provide an empty
        template if the user does not want mcgravi gmt output plots. Or ask
        the user to provide one if needed)
        """
        mcGraviOptions=QtGui.QWidget()
        self.statusBar().showMessage("enter MCGravi options")
                
        mcGraviOptions.datafordriftadj=datafordriftadj
        mcGraviOptions.output_root_dir=output_root_dir
        
        mcGraviOptions.mode=QtGui.QLabel('Mode? (1=datum-free 2=weighted constraint')
        mcGraviOptions.drift_t=QtGui.QLabel('temporal drift polynomial?')
        mcGraviOptions.drift_k=QtGui.QLabel('temperature drift polynomial?')        
        mcGraviOptions.sigma_add=QtGui.QLabel('SD add to data (mgal)')        
        
        mcGraviOptions.modeEdit=QtGui.QLineEdit('1')
        mcGraviOptions.drift_tEdit=QtGui.QLineEdit('1')
        mcGraviOptions.drift_kEdit=QtGui.QLineEdit('0')
        mcGraviOptions.sigma_addEdit=QtGui.QLineEdit('0.005')
             
        # create buttons and actions
        mcGraviOptions.btn1 = QtGui.QPushButton('ok', self)
        mcGraviOptions.btn1.clicked.connect(lambda : self.writeMCGraviInputfiles(mcGraviOptions))
                          
        #locations                
        grid = QtGui.QGridLayout()       
        grid.addWidget(mcGraviOptions.mode,1,0)
        grid.addWidget(mcGraviOptions.modeEdit,1,1) 
        grid.addWidget(mcGraviOptions.drift_t,2,0)
        grid.addWidget(mcGraviOptions.drift_tEdit,2,1)
        grid.addWidget(mcGraviOptions.drift_k,3,0)
        grid.addWidget(mcGraviOptions.drift_kEdit,3,1)
        grid.addWidget(mcGraviOptions.sigma_add,4,0)
        grid.addWidget(mcGraviOptions.sigma_addEdit,4,1)
        grid.addWidget(mcGraviOptions.btn1,5,0)
        mcGraviOptions.setLayout(grid)   
        mcGraviOptions.setWindowTitle('MCGravi options')    
        self.popup=mcGraviOptions
        self.popup.show()
       
    def writeMCGraviInputfiles(self,MCwin):
        """
        - Other problem: extra files are needed (should maybe provide an empty
        template if the user does not want mcgravi gmt output plots. Or ask
        the user to provide one if needed)
        """       
        output_root_dir=MCwin.output_root_dir
        datafordriftadj=MCwin.datafordriftadj        
#        if not os.path.exists(output_root_dir):
#            os.makedirs(output_root_dir)            
        for survid,surv in datafordriftadj.survey_dic.iteritems():
            if surv.keepitem==1:
                survdir=output_root_dir+os.sep+surv.name
                if not os.path.exists(survdir):
                    os.makedirs(survdir)
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
                        filenameC=survdir+os.sep+"fn"+str(j)+"11c"+yr[2:4]+"."+jday
                        file_list.append("fn"+str(j)+"11c"+yr[2:4]+"."+jday)
                        filenameS=survdir+os.sep+"fn"+str(j)+"11s"+yr[2:4]+"."+jday
                        # getbase station: temporary formulation, should add a checkbox 
                        # in the tree object to select it manually for each loop
                        #basestationkeys=[stakeys[i] for i in range(len(stakeys)) if stakeys[i][1]==2]
                        # test:                        
                        basestationkeys=[stakeys[i] for i in range(len(stakeys)) if stakeys[i][0]==self.base_station]
                        
                        loop.writeCgxCFile(filenameC,basestationkeys[0])
                        loop.writeAssociatedCgxSFile(filenameS)
                
                        
                ############### MCGRAVI conf file ####################
                
                write_list_fic='N'
                write_list_station='N'
                write_obs='Y'
                write_mat='N'
                write_resid='N'
                write_tau='N'
                write_grav='Y'
                write_drift="N"
                sigma_factor=1.0
                calf="calib.cal"
                create_r='N'
                outf="mix_%s"%(surv.name)
                fcor="stations_coordinates.txt"
                region="1.602 9.7410 1.607 9.7460"
               
                write_mcgravi_conf_file("%s/conf_%s.conf"%(survdir,surv.name),\
                file_list,int(MCwin.modeEdit.text()),\
                int(MCwin.drift_tEdit.text()),int(MCwin.drift_kEdit.text()),\
                write_list_fic,write_list_station,write_obs,write_mat,\
                write_resid,write_tau,write_grav,write_drift,sigma_factor,\
                float(MCwin.sigma_addEdit.text()),calf,create_r,outf,fcor,region)
                  
                ### cp other mcgravi files to output directory    
                #shutil.copyfile("%s/calib.cal"%(self.data_path),"%s/calib.cal"%(survdir))
                #shutil.copyfile("%s/tunnel.abs"%(self.data_path),"%s/tunnel.abs"%(survdir))
                shutil.copyfile("./MCGravi_files/calib.cal","%s/calib.cal"%(survdir))
                shutil.copyfile("./MCGravi_files/tunnel.abs","%s/tunnel.abs"%(survdir))                
                shutil.copyfile("C:/mcgravi/mcgravi.exe","%s/mcgravi.exe"%(survdir))
                
                self.campaigndata.writeStationLocationFile("%s/stations_coordinates.txt"%(survdir))              
                
                #write tunnel file (absolute value of the base station)
                file=open("%s/tunnel.abs"%(survdir),'w')
                file.write("%d 0 0.001\n"%(self.base_station))
                file.close()                   
                
                # run mcgravi
                curr_dir=os.getcwd()
                os.chdir("%s"%(survdir))
                subprocess.call(["mcgravi","conf_%s.conf"%(surv.name)])
                os.chdir(curr_dir)
                # read mcgravi ouputs
#                self.campaigndata.survey_dic[survid].read_from_mcgravi_output_file(mcgravi_filename)
                #other options from the process menu are now available:          
        self.readMCGraviOutputfiles(datafordriftadj,self.output_root_dir)        
#        self.displaySimpleDiffAction.setEnabled(True)
#        self.heightBouguerCorrectionAction.setEnabled(False)        
#        self.terrainCorrectionAction.setEnabled(False)    
        self.computeDoubleDifferencesAction.setEnabled(True)
        self.saveAction1.setEnabled(True)
        self.popup.close()
        self.statusBar().showMessage("Drift adjusted. Save options and double difference computation enabled")  

    def readMCGraviOutputfiles(self,datafordriftadj,output_root_dir,pattern="mix*",output_file_pattern="*.gra"):
        """
        Read output files *.lst files from the mix_ folders outputs from mcgravi
        
        populate output_dic dictionnary from each survey.
        At station numbers keys of this dictionary, one finds the following tuple:
        (station, gravity, std)
        """
        for survid,surv in datafordriftadj.survey_dic.iteritems():
            if surv.keepitem==1:
                survdir=output_root_dir+os.sep+surv.name     
                #identify every folders matching the given pattern
                folders=glob.glob(survdir+os.sep+pattern)
                #sort the folder list (mcgravi outputs folder names with prog execution date...)
                folders.sort()
                #get the last one
                folder=folders.pop()
                
                mcgravi_filename=glob.glob(folder+os.sep+output_file_pattern)

                #fill in the survey object:
                self.campaigndata.survey_dic[survid].read_from_mcgravi_output_file(mcgravi_filename[0])
                print "For survey: %s"%(num2date(survid))
                print "Station , g (mgal), SD (mgal)"                
                for tupid,tup in self.campaigndata.survey_dic[survid].output_dic.iteritems():
                    print "%d, %7.4f, %7.4f"%(tup)
               
                                                     
    def saveSimpleDiff(self):
        """
        Save simple differences:
        in each selected survey-folder, write a SimpleDifference.dat file
        containing a 1 line header, a station number column, and grav and std
        """

        tday=datetime.now()
        Mainfile=open(self.output_root_dir+os.sep+\
        'simple_diff_data_hierarchy_%4d%02d%02d_%02d%02d%02d.txt'%(tday.year,\
        tday.month,tday.day,tday.hour,tday.minute,tday.second),'w')        
        for survid,surv in self.campaigndata.survey_dic.iteritems():
            print survid
            if surv.keepitem==1:
                survdir=self.output_root_dir+os.sep+surv.name  
                if not os.path.exists(survdir):
                    os.makedirs(survdir)                
                filename=survdir+os.sep+\
                "SimpleDifferences_%4d%02d%02d_%02d%02d%02d.dat"%(tday.year,\
                tday.month,tday.day,tday.hour,tday.minute,tday.second)
                Mainfile.write("%s %s %s\n"%(surv.name,survid,filename))                
                surv.writeOutputValues(filename)
        Mainfile.close()
        

    def computeDoubleDifferences(self):
        """
        subfunction for actually chosing the reference survey
        stored in self.referencesurvey. The computing step is actually 
        executed when saving data...
        """
        #axis=0 to sort on first axis
        survnames=np.sort([(survid,surv.name) for survid,surv in self.campaigndata.survey_dic.iteritems() if surv.keepitem==1],axis=0)
        surveySelectList=QtGui.QListWidget()
        self.statusBar().showMessage("Please choose reference survey (double-click)")
        
        for surv in survnames:
            QListWidgetItem(str(surv[1]),surveySelectList)
                
        surveySelectList.itemActivated.connect(self.referenceSurveySelected)
        self.setWindowTitle('survey selection')  
        self.popup=surveySelectList         
        self.popup.show()
#        self.dataselectionwin.show()   
    
    def referenceSurveySelected(self,item):
        """
        once the reference survey is selected, choose either to
        - compute classic double differences
        - compute double differences from network mean values        
        """             
        survnames=np.sort([(survid,surv.name) for survid,surv in self.campaigndata.survey_dic.iteritems() if surv.keepitem==1],axis=0)
        for surv in survnames:
            if surv[1] == item.text():
                #print self.campaigndata.survey_dic[surv[0]].name
                print surv[1]
                self.referencesurvey=self.campaigndata.survey_dic[int(surv[0])]        
        
        #create new window and set as central Widget:
        ddSelectionWin=QtGui.QWidget()
        self.statusBar().showMessage("Please choose double differences method")
        
        # create buttons and actions
        ddSelectionWin.btn1 = QtGui.QPushButton('Classic double differences', self)
        ddSelectionWin.btn1.clicked.connect(self.ddClassic)
        ddSelectionWin.btn2 = QtGui.QPushButton('Double differences from network mean', self)
        ddSelectionWin.btn2.clicked.connect(self.ddFromNetworkMean)                        
        #locations                
        grid = QtGui.QGridLayout()
        grid.addWidget(ddSelectionWin.btn1,0,0,1,1)
        grid.addWidget(ddSelectionWin.btn2,1,0,1,1)  
        ddSelectionWin.setLayout(grid)   
        ddSelectionWin.setWindowTitle('Double differences') 
        ddSelectionWin.setGeometry(50, 50, 350, 300)
        self.popup=ddSelectionWin
        self.popup.show()
        
    def ddClassic(self):
        """       
        the save option is now enabled
        """
        #double difference option 1=classic 2=network mean
        DD=self.campaigndata.calculateDoubleDifferences(self.referencesurvey,1)
        self.DD=DD
        print "Done: save option for double difference is now activated"
        self.statusBar().showMessage("Done: save option for double difference is now activated")                
        self.saveAction2.setEnabled(True) 
        self.popup.close()
        self.setWindowTitle('CG5 data processing')

    def ddFromNetworkMean(self):        
        """       
        the save option is now enabled
        """
        #double difference option 1=classic 2=network mean
        DD=self.campaigndata.calculateDoubleDifferences(self.referencesurvey,2)
        self.DD=DD
        print "Done: save option for double difference is now activated"
        self.statusBar().showMessage("Done: save option for double difference is now activated")                
        self.saveAction2.setEnabled(True)     
        self.popup.close()
        self.setWindowTitle('CG5 data processing')
        self.show()        
        
    def saveProcessedDoubledifferences(self):
        """
        Save double differences:
        in the output folder, write a doubledifferencesGrav.dat file and a
        doubledifferencesSTD.dat file containing a 1 line header, 
        a survey date column, and grav or std values for each station
        """
        filenameGrav=self.output_root_dir+os.sep+"doubledifferencesGrav"
        filenameSTD=self.output_root_dir+os.sep+"doubledifferencesSTD"        
        self.campaigndata.saveProcessedDoubledifferencesData(filenameGrav,
            filenameSTD,self.DD,self.referencesurvey)
            

        
    def correctRecordedTime(self):
        """
        Correct all times from an offset: when GMT time entered in CG5 is bad.
        """
        #ask for time difference to apply
        text, ok = QtGui.QInputDialog.getText(self, 'Input parameters', 
            'time offset to apply (hr)?')        
        if ok:
            self.t_offset=int(text)                           
            for keysurv,surv in self.campaigndata.survey_dic.iteritems():
                t=np.array(self.campaigndata.survey_dic[keysurv].t)+timedelta(self.t_offset/24,0,0)                              
                self.campaigndata.survey_dic[keysurv].t=t                 
                for keyloop,loop in self.campaigndata.survey_dic[keysurv].loop_dic.iteritems():
                    t=np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].t)+timedelta(self.t_offset/24,0,0)                              
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].t=t                    
                    for keysta,sta in self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic.iteritems():
                        # convert to array: needed for operations
                        t=np.array(self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t)+timedelta(self.t_offset/24,0,0)                              
                        self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].t=t    
            
        #once stations and loops are populated, populate the upstream 
        #structures:
        #for keysurvey,survey in self.campaigndata.survey_dic.iteritems():            
        #    survey.populateFromSubDictionnaries(survey.loop_dic)
        #    self.campaigndata.survey_dic[keysurvey]=survey
        #self.campaigndata.populateFromSubDictionnaries(self.campaigndata.survey_dic)        
    def displaySimpleDiff(self):
        """
        """
        #### define the data selection window Widget:
        gdispwin=QtGui.QWidget()        
        mess="select a survey to see simple differences plots"
        self.statusBar().showMessage(mess)         
        # and the final layout:
        layout_final = QtGui.QGridLayout() 


        #### now define display behaviour for all panels       
        
        # Left panel: 
        survnames=[(survid,surv.name) for survid,surv in self.campaigndata.survey_dic.iteritems() if surv.keepitem==1]
        self.surveySelectList=QtGui.QListWidget()
        self.statusBar().showMessage("Please choose a survey (double-click)")        
        for surv in survnames:
            QListWidgetItem(str(surv[1]),self.surveySelectList)                
        self.surveySelectList.itemActivated.connect(self.surveySelected)
     
        # right panel: plot and some options
        #figure
        self.main_frame = QWidget()        
        self.dpi = 100
        self.fig = Figure((3.0, 2.0), dpi=self.dpi)
        # make some room for the axes labels (set more place between subplots)
        self.fig.subplots_adjust(right=0.95,wspace=0.3)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)        
        self.axes_grav = self.fig.add_subplot(221)
        self.axes_tiltx = self.fig.add_subplot(222)
        self.axes_tilty = self.fig.add_subplot(224)
        self.axes_temp = self.fig.add_subplot(223)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)      
       
        #### add widgets to the layout
       
        # add left panel (tree)
        layout_final.addWidget(self.surveySelectList,0,0,1,1)
        
        # for the right panel (options & display):        
        # create sublayouts (allow the line edit to be of finite extent)
        layout_options = QtGui.QGridLayout()         
        grid = QtGui.QGridLayout() 
        
        # fill subplayouts:
        grid.addWidget(self.canvas,3,0,19,4)
        grid.addWidget(self.mpl_toolbar,2,0,1,4)   
        
        # add subplayouts to the main layout
        grid.addLayout(layout_options,0,0,1,2)
        layout_final.addLayout(grid,0,1,1,1)  
        #grid.setColumnMinimumWidth(1,750)       
        layout_final.setColumnMinimumWidth(0,200)
        layout_final.setColumnMinimumWidth(1,700)

        #### set window geometry and set layout as the main layout      
        
        self.setGeometry(50, 50, 1300, 700)
        self.setWindowTitle('Data selection')    
        gdispwin.setLayout(layout_final)               
        self.setCentralWidget(gdispwin)
        gdispwin.show()       
     
    def surveySelected(self,item):
        """
        once the survey is selected, it can be plotted
        """             
        survnames=[(survid,surv.name) for survid,surv in self.campaigndata.survey_dic.iteritems() if surv.keepitem==1]
        for surv in survnames:
            if surv[1] == item.text():
                print self.campaigndata.survey_dic[surv[0]].name
                simpleDiffDic=self.campaigndata.survey_dic[surv[0]].output_dic
                stanames=[station[0] for keysta,station in sorted(simpleDiffDic.iteritems(), key=lambda x: int(x[0]))] 
                gval=[station[1] for keysta,station in sorted(simpleDiffDic.iteritems(), key=lambda x: int(x[0]))] 
                STD=[station[2] for keysta,station in sorted(simpleDiffDic.iteritems(), key=lambda x: int(x[0]))] 
                if surv.height_bouguer_corr=='ok':
                    gcorr_height=[station[3] for keysta,station in sorted(simpleDiffDic.iteritems(), key=lambda x: int(x[0]))] 
                    gcorr_Bouguer=[station[4] for keysta,station in sorted(simpleDiffDic.iteritems(), key=lambda x: int(x[0]))] 
                    
                    self.setSimpleDiffPlot(self.axes_grav,stanames,gval,STD,gcorr_height,gcorr_Bouguer)
                else:
                    self.setSimpleDiffPlot(self.axes_grav,stanames,gval,STD)
                    

    def setSimpleDiffPlot(self,axe,seriex,seriey,std,serie_corr1=None,serie_corr2=None,serie_corr3=None):
        """
        plot simple differences
        """
        axe.clear()
        axe.grid(True)        
        axe.plot(seriex,seriey,linestyle="dashed", marker="o",color='k',label='g (mgal)')        
        axe.errorbar(seriex,seriey,yerr=std,linestyle="None", marker="None",color='k',label='g (mgal)')
        #axe.plot(seriex_selec,seriey_selec,'o-',color='b',label=serie_type)            
        axe.set_ylabel('g (mgal)', size='x-small')
        #axe.set_title(serie_type, size='x-small')
        labels = axe.get_xticklabels() + axe.get_yticklabels()
        for label in labels:
            label.set_size('x-small')            
            
        self.canvas.draw()                
        
    def heightBouguerCorrection(self):
        """
        """
         # ask for data file
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.output_root_dir)  
        self.campaigndata.corrHeightBouguer(fname)        
        
        
    def terrainCorrection(self):
        """
        """
        raise "not implemented"
              
    def Correction(self):
        """
        """
        filename=self.output_root_dir+os.sep+'data_input.dat'
        self.static_dataset.writePGraviFor3DInputFile(filename)
        
                
                #write_mcgravi_conf_file("%s/conf_%s.conf"%(survdir,surv.name),file_list,mode,drift_t,drift_k,write_list_fic,\
                #write_list_station,write_obs,write_mat,write_resid,write_tau,write_grav,write_drift,sigma_factor,sigma_add,calf,create_r,outf,fcor,region)
                  
                ### cp other mcgravi files to output directory    
                #shutil.copyfile("./calib.cal","%s/calib.cal"%(survdir))
                #shutil.copyfile("./tunnel.abs","%s/tunnel.abs"%(survdir))                
                #shutil.copyfile("C:/mcgravi/mcgravi.exe","%s/mcgravi.exe"%(survdir))
                
                #self.campaigndata.writeStationLocationFile("%s/stations_coordinates.txt"%(survdir))              
                
                #write tunnel file (absolute value of the base station)
                #file=open("%s/tunnel.abs"%(survdir),'w')
                #file.write("%d 0 0.001\n"%(self.base_station))
                #file.close()                   
                
                # run mcgravi
                #curr_dir=os.getcwd()
                #os.chdir("%s"%(survdir))
                #subprocess.call(["mcgravi","conf_%s.conf"%(surv.name)])
                #os.chdir(curr_dir)                
        
    def close_windows(self):
        """ function for closing all windows
        """
#        self.displaydatawin.close()
        self.close()
        
    def create_menu(self):     
        """
        Menu creation
        """         
        ########################## Create a File Menu:
        self.file_menu = self.menuBar().addMenu("&File")             
        ############ define actions: 
        # for project starting
        startProjectAction = self.create_action("&Start project",
            shortcut="Ctrl+D", slot=self.startProject, tip="Start new project") 
        # for file opening
        self.openFile = self.create_action("&Load raw data",
            shortcut="Ctrl+O", slot=self.openRawdata, tip="Open new File",
            enabled=False)       
        # for processed file opening
        self.openModifCFile = self.create_action("&Load processed data",
            shortcut="Ctrl+P",slot=self.openModifCData, 
            tip="Open modified c-type data",enabled=False)           
        # for simple differences file opening
        self.openSimpleDiffFile = self.create_action("&Load simple differences",
            slot=self.openSimpleDiff, 
            tip="Open simple difference file",enabled=False)    
#        # for simple differences file opening
#        self.openSimpleDiffFile_static = self.create_action("&Load simple differences (static)",
#            slot=self.openSimpleDiff_static, 
#            tip="Open simple difference file (static)",enabled=False)               
        # for menu closing
        exitAction=self.create_action("&Exit",
            shortcut="Ctrl+Q", slot=self.close_windows, tip="Exit App") 
            #shortcut="Ctrl+Q", slot=sys.exit(), tip="Exit App")         
#        # for simple save
#        self.saveUnprocessedAction=self.create_action("&Save 'c' files",
#            slot=self.saveUnprocessedData, 
#            tip="Save unprocessed data",enabled=False)    
        # for simple save
        self.saveProcessedAction=self.create_action("&Save processed data (modified 'c' files)",
            slot=self.saveTemporaryProcessedData,shortcut="Ctrl+S", 
            tip="Save processed data",enabled=False)              
        # for simple save of processed data
        self.saveAction1=self.create_action("&Save simple differences",
            slot=self.saveSimpleDiff, 
            tip="Save processed simple differences",enabled=False)       
        # for double diff save
        self.saveAction2=self.create_action("&Save double differences",
            slot=self.saveProcessedDoubledifferences, 
            tip="Save processed double differences surveys",enabled=False)              
        ############ add actions to menu
#        self.add_actions(self.file_menu,(startProjectAction,self.openFile,
#                        self.openModifCFile,self.openSimpleDiffFile, 
#                        self.openSimpleDiffFile_static, 
#                        self.saveUnprocessedAction,self.saveProcessedAction,
#                        self.saveAction1,self.saveAction2,exitAction))
        self.add_actions(self.file_menu,(startProjectAction,self.openFile,
                        self.openModifCFile,self.openSimpleDiffFile,  
                        self.saveProcessedAction,
                        self.saveAction1,self.saveAction2,exitAction))        
        
        ########################## Create a Process Menu:
        self.process_menu = self.menuBar().addMenu("&Process")             
        ############ define actions:
        # for tide correction
        self.tideCorrectionAction = self.create_action("&Tide correction",
            shortcut="Ctrl+T", slot=self.tideCorrection, 
            tip="Choose tide correction method",enabled=False)    
        # for ocean loading correction
        self.oceanLoadingCorrectionAction = self.create_action("&Ocean loading correction",
            shortcut="Ctrl+L", slot=self.oceanLoadingCorrection, 
            tip="Ocean loading corrections",enabled=False)              
        # for atmospheric correction
        self.atmosphericCorrectionAction = self.create_action("&Atmospheric correction",
            shortcut="Ctrl+P", slot=self.atmosphericCorrection, 
            tip="Apply atmospheric correction by loading a data file",enabled=False)               
        # for data selection
        self.dataSelectionAction = self.create_action("&Data selection",
            shortcut="Ctrl+d", slot=self.dataSelection, tip="Select data",
            enabled=False)  
        # for drift retrieval
        self.driftAdjustmentAction = self.create_action("&Drift adjustment",
            shortcut="Ctrl+A", slot=self.driftAdjustment,
            tip="Least-square drift adjustment",enabled=False)
        # for computing double differences
        self.computeDoubleDifferencesAction = self.create_action("&Compute double differences",
            slot=self.computeDoubleDifferences,
            tip="Compute double differences",enabled=False)          
        # for correcting time 
        self.correctRecordedTimeAction = self.create_action("&Correct recorded time",
            slot=self.correctRecordedTime,
            tip="Correct recorded time",enabled=False)                 
        ############ add actions to menu
        self.add_actions(self.process_menu,
            (self.tideCorrectionAction,self.oceanLoadingCorrectionAction, self.atmosphericCorrectionAction,
             self.dataSelectionAction,self.driftAdjustmentAction,
             self.computeDoubleDifferencesAction,self.correctRecordedTimeAction))
        
#        ########################## Create a Static corrections Menu:
#        self.static_menu = self.menuBar().addMenu("&Static corrections")             
#        ############ define actions:
#        # for diplaying simple differences
#        self.displaySimpleDiffAction = self.create_action("&Display simple diff.",
#            slot=self.displaySimpleDiff,tip=" ",enabled=False)             
#        # for height corrections
#        self.heightBouguerCorrectionAction = self.create_action("&Height-Bouguer correction",
#            slot=self.heightBouguerCorrection, 
#            tip=" ",enabled=False)                   
#        # for terrain selection
#        self.terrainCorrectionAction = self.create_action("&terrain correction",
#            slot=self.terrainCorrection, tip=" ",
#            enabled=False)          
#        # for CORRECTIONS
#        self.CorrectionAction = self.create_action("&corrections totale",
#            slot=self.Correction, tip=" ",
#            enabled=False)                 
#        ############ add actions to menu
#        self.add_actions(self.static_menu,
#            ( self.displaySimpleDiffAction,self.heightBouguerCorrectionAction,
#             self.terrainCorrectionAction,self.CorrectionAction))        
             
             
    def create_status_bar(self,status_string):
        self.status_text = QtGui.QLabel(status_string)
        self.statusBar().addWidget(self.status_text, 1)
                
    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()",enabled=True):
        """
        Simplify action creation
        """                            
                            
        action = QtGui.QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, QtCore.SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        if enabled is not True:
            action.setEnabled(False)
        return action                            
        
    def add_actions(self, target, actions):
        """
        Add actions to a target
        """  
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action) 

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
        
        
def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = mainProg()
    app.exec_()
#    sys.exit(app.exec_())        
if __name__ == '__main__':
    main()
   

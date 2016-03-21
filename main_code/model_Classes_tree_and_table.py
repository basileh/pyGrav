# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 14:38:09 2014

@author: Basile


Package containing classes for models definitions for model/view programming in
the GUI:

Models contain the generic display information of an item which is not stored 
in the model

Models defined here:
- myTree: Tree-type model for representing the data hierarchical 
    structures
- stationDataTableModel: Table-type model for representing channels of a single
    station (and eventually use it for data selection)   

Inspired from http://www.apileofgrains.nl/qtreeview-qmodel-example-pyqt/
Tutorial for model/view programing: http://www.yasinuludag.com/blog/?p=98
Other example: 
http://rowinggolfer.blogspot.fr/2010/05/qtreeview-and-qabractitemmodel-example.html
"""


from PyQt4 import QtCore,QtGui
#from data_objects import *
#from item_classes_campaign import *
import numpy as np
from matplotlib.dates import date2num,num2date


class MyTreeItem(QtGui.QTreeWidgetItem):
    """
    Tree item
    each item stores its coordinates within the grav_obj hierarchy, to easily
    access the object instance from the main program (and modify it!)
    
    the QTreeWidgetItem is a standardized tree widget item, which implicitely
    associate a standard tree model. It therefore defines the view associated 
    to a standard model (of tree type)
    
    properties:
    s:                    a string: name to display on the tree (object name)
    keysurv:              the key of a survey object within the grav_obj
                          hierarchy
    keyloop:              the key of a loop object within the grav_obj
                          hierarchy
    keysta:               the key of a station object within the grav_obj
                          hierarchy
                          
                          
    """
    def __init__(self, s, keysurv=None, keyloop=None, keysta=None,parent = None):

        super(MyTreeItem, self).__init__(parent, [s])
        self.keysurv=keysurv
        self.keyloop=keyloop
        self.keysta=keysta
        self.setFlags(self.flags() | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsUserCheckable) 
        self.setCheckState(0, QtCore.Qt.Checked) # 0 is the column number
        self.setCheckState(1, QtCore.Qt.Checked) # 1 is the column number
        self.setCheckState(2, QtCore.Qt.Checked) # 2 is the column number


class MyTree(QtGui.QTreeWidget):
    """
    a standard tree widget, which provides a standard tree view, implicitely
    associated with a predefined tree model

    The tree widget items are populated based on the data hierarchy, provided
    by a Campaign-type object.
    
    properties:
    campaigndata:         a Campaign object, containing all the data
    
    """

    def __init__(self,campaigndata,parent = None):

        super(MyTree, self).__init__(parent)

        # create some data members, these will be set from the outside and trigger a model change
        self.campaigndata = campaigndata          
        for keysurvey,survey in self.campaigndata.survey_dic.iteritems():            
            # Create a survey tree item
            survey_item = MyTreeItem(survey.name,keysurv=keysurvey,parent=self)
            # for all the loops attached to the survey create a loop tree item
            for keyloop,loop in self.campaigndata.survey_dic[keysurvey].loop_dic.iteritems():  
                # get loop name: date of the first measure of the first station
                #after base station (the latter being potentially acquired the 
                # previous day)
                station_keys=[keysta for keysta,station in sorted(self.campaigndata.survey_dic[keysurvey].loop_dic[keyloop].station_dic.iteritems(), key=lambda x: x[1].t[1])]                  
                loopname=loop.station_dic[station_keys[1]].t[1].date()
                loop_item = MyTreeItem("loop" + loop.name+' ('+ str(loopname)+')',keysurv=keysurvey,keyloop=keyloop,parent=survey_item)                
                
                # the following sorts dictionary items following the first value of the 't'-list attribute
                for keysta,station in sorted(self.campaigndata.survey_dic[keysurvey].loop_dic[keyloop].station_dic.iteritems(), key=lambda x: x[1].t[1]):                  
                
                #for keysta,station in self.campaigndata.survey_dic[keysurvey].loop_dic[keyloop].station_dic.iteritems():                    
                    # create the station item
#                    name_to_disp="sta" +
                    
                    station_item = MyTreeItem(station.name,keysurv=keysurvey,keyloop=keyloop,keysta=keysta, parent=loop_item)
                    

# création du modèle de table pour la sélection des données de station
class stationDataTableModel(QtCore.QAbstractTableModel):
    """
    a table model object. Derived from an abstract model, so all main functions
    have to be reimplemented (rowCount,columnCount,data)
    see http://www.yasinuludag.com/blog/?p=98 for rich tutorials on model/view
    programing
    
    The station position in the data hierarchy are stored, so that if a 
    modification is triggered, the original data can be accessed and changed
    accordingly (keysurv,keyloop,keysta)
    
    If other columns are to be chosed, the following have to be modified here:
    __headers, arraydata construction protocol
    
    by default, all table entries are checked (this can be modified to allow
    pre-check based on user criteria (tiltx,tilty,...)). Then, if one is 
    unchecked, the keepdata property of the ChannelList object at the table row
    position is set to 0
    
    properties:
    __headers:            table headers        
    unchecked:            a dictionnary of unchecked items. Keys are item 
                          indexes, entries are states
    ChannelList_obj:      an object of ChannelList-type: used to store the 
                          table data as the structured data
    arraydata:            an array representation of the data from the
                          ChannelList_obj
    keysurv:              the key of a survey object within the grav_obj
                          hierarchy
    keyloop:              the key of a loop object within the grav_obj
                          hierarchy
    keysta:               the key of a station object within the grav_obj
                          hierarchy                                                          
    """
    def __init__(self,ChannelList_obj,headers = [], keysurv=None, keyloop=None, keysta=None,  parent = None):
        self.__headers=["station",u"g (\u00b5gal)",u"sd (\u00b5gal)","tiltx","tilty","temp (K)","dur (s)","rej","t (mn)","date/time"]
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.unchecked = {}
        self.createArrayData(ChannelList_obj)
        self.keysurv=keysurv
        self.keyloop=keyloop
        self.keysta=keysta
      
    def createArrayData(self,ChannelList_obj)  :
        """
        Create the np array data for table display, and update the 
        ChannelList_obj. This function can be called from outside to update the
        table display
        """
        self.ChannelList_obj = ChannelList_obj
        self.arraydata=np.concatenate((ChannelList_obj.station,
            np.array(ChannelList_obj.grav)*1000,np.array(ChannelList_obj.sd)*1000,ChannelList_obj.tiltx,
            ChannelList_obj.tilty,ChannelList_obj.temp,ChannelList_obj.dur,
            ChannelList_obj.rej,
            (date2num(ChannelList_obj.t)-date2num(ChannelList_obj.t[0]))*24*60,
            np.array(ChannelList_obj.t))).reshape(len(ChannelList_obj.t),10,order='F')
            
    def rowCount(self, parent):
        return len(self.ChannelList_obj.t)

    def columnCount(self, parent):
        return len(self.__headers)
        
    def flags(self, index):
        """
        allow the checkboxes (Qt.ItemIsUserCheckable)
        """
        return QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable 
        
    def data(self, index, role):
        """
        display data definition
        """
        if not index.isValid():
            return None
        if role == QtCore.Qt.DisplayRole:
            # view definition
            row=index.row()
            column=index.column()
            strdata=str(self.arraydata[row][column])
            if self.__headers[column]=="t (mn)":
                strdata="%2.1f"%self.arraydata[row][column]
            if self.__headers[column]=="rej":
                strdata="%2.0f"%self.arraydata[row][column]    
            if self.__headers[column]=="station" or self.__headers[column]=="dur (s)":
                strdata="%3.0f"%self.arraydata[row][column]  
            if self.__headers[column]==u"g (\u00b5gal)" or self.__headers[column]==u"sd (\u00b5gal)":
                strdata="%8.0f"%self.arraydata[row][column]                  
            return strdata
            
        if role == QtCore.Qt.CheckStateRole: 
            # check status definition
            if index.column() == 0:
                return self.checkState(index)
                 
            
    def checkState(self, index):
        """
        by default, everything is checked. If keepdata property from the 
        ChannelList object is 0, it is unchecked
        """
        if self.ChannelList_obj.keepdata[index.row()]==0:
            self.unchecked[index]= QtCore.Qt.Unchecked           
            return self.unchecked[index]
        else: 
            return QtCore.Qt.Checked
        
         
            
    def setData(self, index, value, role):
        """
        if a row is unchecked, update the keepdata value to 0
        setData launched when role is acting
        value is QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if (role == QtCore.Qt.CheckStateRole and index.column() == 0):
            if value == QtCore.Qt.Checked:
                self.ChannelList_obj.keepdata[index.row()]=1
            elif value == QtCore.Qt.Unchecked:    
                self.unchecked[index] = value
                self.ChannelList_obj.keepdata[index.row()]=0
            return True 
        return QtCore.QAbstractTableModel.setData(self, index, value, role)
            
    def headerData(self, section, orientation, role):
        """
        display header
        """
        if role == QtCore.Qt.DisplayRole:            
            if orientation == QtCore.Qt.Horizontal:                
                if section < len(self.__headers):
                    return self.__headers[section]
                else:
                    return "not implemented"        
              
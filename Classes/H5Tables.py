# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 12:17:47 2014

@author: szordan
"""
import numpy as np, h5py
from tables import *
import time
import collections



class TimeStampData(IsDescription):
    timeStamp = Int32Col(pos=1)   # 16-character String
    idUnitStim = Int32Col(pos=2)     # Signed 64-bit integer
    idShape = Int32Col(pos=3)
    
class NeuronMatrixData(IsDescription):
    id  = Int32Col(pos=1)  
    neuronName  = StringCol(256,pos=2)   
    posX  = Int32Col(pos=3)     
    posY  = Int32Col(pos=4)
    unitNum = Int32Col(pos=5)
    marked = Int32Col(pos=6)
    
class StimulusMatrixData(IsDescription):
    id  = Int32Col(pos=1)
    idSequence  = Int32Col(pos=2)   
    idFrame  = Int32Col(pos=3)     
    idImage  = Int32Col(pos=4)

class StimuliParameterData(IsDescription):
    id  = Int32Col(pos=1)
    pathSequence  = StringCol(512,pos=2) 
    nRefreshStimulus = Int32Col(pos=3) 
    nStimuliToBeDisplayed =Int32Col(pos=4) 
    nOfStimuliDisplayed =Int32Col(pos=5) 
    duration =Int32Col(pos=6) 
    modality=Int32Col(pos=7) 
    nOfImage = Int32Col(pos=8) 

class RecordingData(IsDescription):
    spikeTrainMat  = StringCol(512,pos=1)    
    startFrame  = Int32Col(pos=2)   
    stopFrame  = Int32Col(pos=3)     
    samplingFreq  = Float64Col(pos=4)
    
class ImageData(IsDescription):
    id  = Int32Col(pos=1)
    pathImage = StringCol(512,pos=2)
    
class ExpPathData(IsDescription):
    
    report = StringCol(512,pos=1)
    trigger =StringCol(512,pos=2)
    rootStimuli = StringCol(512,pos=3)
    pathReplace = StringCol(512,pos=4) 
    maskingFile=StringCol(512,pos=5) 

class H5():
    
    def __init__(self,pathfile,mode):
        
        self.path=pathfile
        self.mode=mode
        self.filehdf5=open_file(pathfile,mode=self.mode)
        
        
    '''
    ***************************************************************************
    Method for write structure in a hdf5
    ***************************************************************************
    '''
    def addStructures(self,Dstruct,Position):
        
        '''
        Add the matrixes structures contain in a dictionary in root if posiotion= r
        or in a group named as position value
        
        '''
        if Position=='r':
            '''
            info=open('..\\InfoFile.txt','r')
            listLine=info.readlines()
            info.close()
            dset=self.filehdf5.create_array("/","InfoFile",obj=listLine)
            '''
            for matrix in Dstruct:
                if matrix=="TimeStampMatrix":
                    
                    
                    dset=self.filehdf5.create_table("/",matrix,TimeStampData,'Matrix of time stamp event')
                    particle =dset.row
                    
                    for i in range(0,len(Dstruct[matrix][0])):
                        particle['timeStamp']=Dstruct[matrix][0][i]
                        particle['idUnitStim']=Dstruct[matrix][1][i]
                        particle['idShape']=Dstruct[matrix][2][i]
                        particle.append()
                    #dset.cols.timeStamp.create_index()
                    #dset.cols.idUnitStim.create_index()
                elif matrix=="StimuliMatrix":
                    
                    
                    dset=self.filehdf5.create_table("/",matrix,StimulusMatrixData,'Matrix of stimulus event')
                    particle =dset.row
                    
                    for i in range(0,len(Dstruct[matrix][0])):
                        particle['id']=i+1
                        particle['idSequence']=Dstruct[matrix][1][i]
                        particle['idFrame']=Dstruct[matrix][2][i]
                        particle['idImage']=Dstruct[matrix][3][i]
                        particle.append()
                        
                elif matrix=="NeuronMatrix":
                    
                    
                    dset=self.filehdf5.create_table("/",matrix,NeuronMatrixData,'Neuron Information')
                    particle =dset.row
                    
                    for i in range(0,len(Dstruct[matrix][0])):
                        particle['id']=Dstruct[matrix][0][i]
                        particle['neuronName']=Dstruct[matrix][1][i]
                        particle['posX']=Dstruct[matrix][2][i]
                        particle['posY']=Dstruct[matrix][3][i]
                        particle['unitNum']=Dstruct[matrix][4][i]
                        particle['marked']=Dstruct[matrix][5][i]
                        particle.append()
                
                
                elif matrix=="StimuliSequenceMatrix":
                    dset=self.filehdf5.create_table("/",matrix,StimuliParameterData,'Stimuli sequence data')
                    particle =dset.row
                
                    for i in range(0,len(Dstruct[matrix][0])):
                        particle['id']=Dstruct[matrix][0][i]
                        particle['pathSequence']=str(Dstruct[matrix][1][i])
                        particle['nRefreshStimulus']=Dstruct[matrix][2][i]
                        particle['nStimuliToBeDisplayed']=Dstruct[matrix][3][i]
                        particle['nOfStimuliDisplayed']=Dstruct[matrix][4][i]
                        particle['duration']=Dstruct[matrix][5][i]
                        particle['modality']=Dstruct[matrix][6][i]
                        particle['nOfImage']=Dstruct[matrix][7][i]
                        particle.append()
                        
                elif matrix=="StimuliImageArray":
                    
                    dset=self.filehdf5.create_table("/",matrix,ImageData,'Image path data')
                    particle =dset.row
                    for i in range(0,len(Dstruct[matrix][0])):
                        particle['id']=i+1
                        particle['pathImage']=str(Dstruct[matrix][1][i])
                        particle.append()
                        
                    
                elif matrix=="NeuronNameArray":
                    pass
                elif matrix=="StimuliSequenceArray":
                    pass
                
                    
                else:
                    print matrix
                    if isinstance(Dstruct[matrix],dict):
                    
                        dset=self.filehdf5.create_array("/",matrix,obj=Dstruct[matrix]['data'])
                        
                        for index,attribute in enumerate(Dstruct[matrix]['Attributes']):
                            dset.attrs[attribute]=index
                        
                    else:
                         dset=self.filehdf5.create_array("/",matrix,obj=Dstruct[matrix])
                    
                    
                    
        else:
            
            if Position in self.filehdf5.root.__members__:
                group="/"+Position
                
            else:
                group = self.filehdf5.create_group("/", Position,'Global Parameter')
            
            for matrix in Dstruct:
                
                if matrix == 'RecordingParameter':
                    
                    dset=self.filehdf5.create_table(group,matrix,RecordingData,'Recording parameter data')
                    particle =dset.row
                    
                    for i in range(0,len(Dstruct[matrix][0])):
                        
                        particle['spikeTrainMat']=str(Dstruct[matrix][0][i])
                        particle['startFrame']=Dstruct[matrix][1][i]
                        particle['stopFrame']=Dstruct[matrix][2][i]
                        particle['samplingFreq']=Dstruct[matrix][3][i]
                        particle.append()
                    
                elif matrix == 'ExpPath':
                    
                    dset=self.filehdf5.create_table(group,matrix,ExpPathData,'Path sorce file')
                    particle =dset.row
                    for i in range(0,len(Dstruct[matrix][0])):
                        #particle['spikeTrainMat']=Dstruct[matrix][0]
                        particle['trigger']=Dstruct[matrix][0][i]
                        particle['report']=Dstruct[matrix][1][i]
                        particle['rootStimuli']=Dstruct[matrix][2][i]
                        particle['pathReplace']=Dstruct[matrix][3][i]
                        particle['maskingFile']=Dstruct[matrix][4][i]
                        particle.append()
                    
                
                        
                else:
                    if isinstance(Dstruct[matrix],dict):
                    
                        dset=self.filehdf5.create_array(group,matrix,obj=Dstruct[matrix]['data'])
                        
                        for index,attribute in enumerate(Dstruct[matrix]['Attributes']):
                            dset.attrs[attribute]=index
                        
                    else:
                       # if matrix in self.filehdf5.root.Roi._
                         dset=self.filehdf5.create_array(group,matrix,obj=Dstruct[matrix])
                         
                         
                         
    
    '''
    ***************************************************************************
    Method for read structure in a hdf5
    ***************************************************************************
    '''           
    
    def getSpike2xn(self):
        spike=[row['timeStamp'] for row in self.filehdf5.root.TimeStampMatrix.where('idUnitStim>0')]
        channel=[row['idUnitStim'] for row in self.filehdf5.root.TimeStampMatrix.where('idUnitStim>0')]
        
        samplingFreq=self.filehdf5.root.GlobalParameter.RecordingParameter.col('samplingFreq')[0]
        
        dt=1000 / samplingFreq
        
        TotalTimeRecording= self.filehdf5.root.GlobalParameter.RecordingParameter.col('stopFrame')[0]*dt
        spike2xn=np.row_stack((np.asarray(channel),np.asarray(spike)*dt))
        return spike2xn
        
    def getStructureList(self):
        return self.filehdf5.keys()
    
    def getStructureDict(self):
        
        dictS={}
        for key in self.filehdf5.root.__members__:
            print key
        
            if isinstance(self.filehdf5.root._v_children[key],Group):
                dictS[key]={}
                for element in self.filehdf5.root._v_children[key].__members__:
                    
                    if isinstance(self.filehdf5.root._v_children[key]._v_children[element],Table):
                        matrixList=[]
                        for colname in self.filehdf5.root._v_children[key]._v_children[element].colnames:
                            matrixList.append(self.filehdf5.root._v_children[key]._v_children[element].col(colname))
                        dictS[key][element]=matrixList
                    else:
                        dictS[key][element]=np.asarray(self.filehdf5.root._v_children[key]._v_children[element])
                    
                    
            else:
                
                if isinstance(self.filehdf5.root._v_children[key],Table):
                    matrixList=[]
                    for colname in self.filehdf5.root._v_children[key].colnames:
                        matrixList.append(self.filehdf5.root._v_children[key].col(colname))
                    dictS[key]=matrixList
                else:
                    dictS[key]=np.asarray(self.filehdf5.root._v_children[key])
            
        return dictS
        
    def getMfrMatrix(self,thr):
        
        matrixMfr=np.zeros((64,64))
        matrixActiveChannel=np.zeros((64,64))-1
        timeStampMatrix=self.filehdf5.root.TimeStampMatrix
        NeuronMatrix=self.filehdf5.root.NeuronMatrix.cols
        print len(NeuronMatrix.neuronName)
        
        mtr=[]
        for i in range(0,64):
            line=[]
            for j in range(0,64):
              line.append([])
            mtr.append(line)
       
        
        
        startTime= time.time()
        index=timeStampMatrix.read(field='idUnitStim')   
        counter=collections.Counter(index)
        #print counter
        for i in range(0,len(NeuronMatrix.neuronName)):
        #for i in range(0,100) :
            rowX=int(NeuronMatrix.posX[i])-1
            
            
            colY=int(NeuronMatrix.posY[i])-1
            
            name=NeuronMatrix.neuronName[i]
            
            #whereStr='(idevent=='+str(i)+')'
            #spikeTrain= [row['frame'] for row in timeStampMatrix.where(whereStr)]
            #spikeTrain=timeStampMatrix.get_where_list(whereStr)            
            #mfr=int(len(spikeTrain))
            
            #if mfr > matrixMfr[rowX,colY]:
            #print i
            matrixMfr[rowX,colY] += int(counter[i])
            if counter[i] > thr:
                matrixActiveChannel[rowX,colY]= i
            mtr[rowX][colY].append(counter[i])
            
        print time.time()-startTime    
        return matrixMfr,matrixActiveChannel,mtr
     
    def checkStructure(self,name):
        if name in self.filehdf5.keys():
            return True
        else:
            return False
    def getStructure(self,name):
        return np.asarray(self.filehdf5[name])
        
    def getMaskingPath(self):
        maskingFile=str(self.filehdf5.root.GlobalParameter.ExpPath.cols.maskingFile[0])
        return maskingFile
    def getInfoFileString(self):
        
        report=str(self.filehdf5.root.GlobalParameter.ExpPath.cols.report[0])
        spikeTrain=str(self.filehdf5.root.GlobalParameter.RecordingParameter.cols.spikeTrainMat[0])
        trigger=str(self.filehdf5.root.GlobalParameter.ExpPath.cols.trigger[0])
        rootStimuli=str(self.filehdf5.root.GlobalParameter.ExpPath.cols.rootStimuli[0])
        
        
        InfoStr="SOURCE FILE \n"+ "Report: " + report + "\n" + "Trigger: "+trigger+"\n"+ "Stimuli: " + rootStimuli +"\n" + "SpikeTrain"+ spikeTrain + "\n"
        return InfoStr
                
                
    def getRoiList_Matrix(self):
        Droi={}
        DroiMatrix={}
        NeuronMatrix=self.filehdf5.root.NeuronMatrix.cols
        
        for roiName in self.filehdf5.root.Roi.__members__ :
            roiMatrix=np.zeros((64,64))            
            Droi[roiName]=np.asarray(self.filehdf5.root.Roi._v_children[roiName])
            for neuronIndex in Droi[roiName]:
                row=NeuronMatrix.posX[int(neuronIndex)]-1
                col=NeuronMatrix.posY[int(neuronIndex)]-1
                roiMatrix[row,col]=1
            DroiMatrix[roiName]=roiMatrix
        return Droi,DroiMatrix
                    
    def close(self):
        self.filehdf5.close()
           
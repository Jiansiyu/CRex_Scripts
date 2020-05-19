'''
Read in the beam information and parser
Get the average etc
'''
import os
import sys
import matplotlib.pyplot as plt
import statistics
import numpy as np
import pandas as pd
runList={21626:'C_12 0%',21632:'C_12 1%',21641:'c_12 -1%',21642:'C_12 -2%'}
class beamInfor(object):
    def __init__(self):
        self.DataDict={}
        pass

    def ReadRawFile(self,Path='',nameTemplate="RHRS_{}_BeamE.txt"):
        for runID in runList:
            fullPath=os.path.join(Path,nameTemplate.format(runID))
            if os.path.isfile(fullPath):
                with open (fullPath) as flileio:
                    lines=flileio.readlines()
                    content=[x.strip() for x in lines]
                    beamEArray=[float(line.split()[-1]) for line in content]
                    self.DataDict[runList[runID]]=beamEArray
    def BeamEPlot(self):
        self.dataframe=pd.DataFrame()
        for item in self.DataDict:
            df=pd.DataFrame(self.DataDict[item],index=[x for x in range(0,len(self.DataDict[item]))],columns=[str(item)])
            mean=statistics.mean(self.DataDict[item])
            meandf=pd.DataFrame([mean for x in range(0,len(self.DataDict[item]))],index=[x for x in range(0,len(self.DataDict[item]))],columns=['mean'])
            df=df.join(meandf)
            
            print(df)
            # df.plot()
            plt.plot(df)
            plt.text(10,mean,"Mean:{}".format(mean),size=18)
            plt.title('EPICS BeamE {}'.format(item))
            plt.show()
            
        

    def test(self):
        self.ReadRawFile(Path='/home/newdriver/Storage/Server/JLabTempStorage/EPICS_BeamE')
        print(self.DataDict)
        self.BeamEPlot()
            

def BeamE(fname=''):
    if not os.path.isfile(fname):
        print("Can Not Find File {}".format(fname))
    
    with open(fname) as fileio:
        lines=fileio.readlines()
        content=[x.strip() for x in lines]
        beamEArray=[float(line.split()[-1]) for line in content]
        #print(beamEArray)
        print(statistics.mean(beamEArray))           

if __name__ == "__main__":
    test=beamInfor()
    test.test()
    # for i in sys.argv[1:]:
    #     #print("process File:{}".format(i))
    #     if i.endswith(".txt"):
    #         print("process File:{}".format(i))
    #         BeamE(i)
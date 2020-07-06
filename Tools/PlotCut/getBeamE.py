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

#runList={21626:'C_12 0%',21632:'C_12 1%',21641:'c_12 -1%',21642:'C_12 -2%',21739:'H2O_0%',21762:'H2O_1%',21789:'H2O_-1%'}
# runList={ 
#     21641:'C_12 -1%_0',
#     21640:'C_12 -1%_1',
#     21639:'C_12 -1%_2',
#     21638:'C_12 -1%_3',
#     21637:'C_12 -1%_4'
#     }

# runList={ 
#     21784:'C_12 -1%_0',
#     21783:'C_12 -1%_1',
#     21782:'C_12 -1%_2',
#     21781:'C_12 -1%_3',
#     21785:'C_12 -1%_4'
#     }

# runList={ 
#     21739:'H2O_ 0%_0',
#     21740:'H2O_ 0%_1',
#     21789:'H2O -1%_0',
#     21790:'H20 -1%_1'
#     }


runList={
    2566:'C_12 -2%',
    2565:'C_12 -1%',
    2556:'C_12 +1%',
    2550:'C_12  0%',
    2674:'H2O_0%',
    2697:'H2O_1%',
    2726:'H2O_-1%'
    }

# runList={
#     20804:'H2O -1%',
#     20805:'H2O +1%',
# #    20808:'H2O 0%',
#     20825:'C12 1%',
#     20826:'C12 -1%',
#     20827:'C12 0%'
#     }

class beamInfor(object):
    def __init__(self):
        self.DataDict={}
        pass

    def ReadRawFile(self,Path='',nameTemplate="LHRS_{}_BeamE.txt"):
        for runID in runList:
            fullPath=os.path.join(Path,nameTemplate.format(runID))
            if os.path.isfile(fullPath):
                with open (fullPath) as flileio:
                    lines=flileio.readlines()
                    content=[x.strip() for x in lines]
                    #beamEArray=[(line.split()[-1]) for line in content]
                    beamEArray=[float(line.split()[-1]) for line in content]
                    # print(beamEArray)
                    beamEArrayFilter=[]
                    for x in beamEArray:
                        if x > 940:
                            beamEArrayFilter.append(x)
                    self.DataDict[runID]=beamEArrayFilter#[0:len(beamEArray)//4]
    def BeamEPlot(self):
        self.dataframe=pd.DataFrame()
        for item in self.DataDict:
            df=pd.DataFrame(self.DataDict[item],index=[x for x in range(0,len(self.DataDict[item]))],columns=[str(item)])
            mean=statistics.mean(self.DataDict[item])
            stdv=statistics.stdev(self.DataDict[item])
            meandf=pd.DataFrame([mean for x in range(0,len(self.DataDict[item]))],index=[x for x in range(0,len(self.DataDict[item]))],columns=['mean'])
            df=df.join(meandf)
            print(df)
            # df.plot()
            plt.plot(df)
            plt.text(10,mean,"Mean:{}, stdv{}".format(mean,stdv),size=18)
            plt.title('EPICS BeamE {}'.format(item))
            plt.savefig("/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/BeamE/BeamE{}.jpg".format(item))
            plt.show()
            

    def test(self):
        self.ReadRawFile(Path='/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/BeamE')
        #self.ReadRawFile(Path='/home/newdriver/Storage/Server/JLabFarm/PRex/BeamE/')
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
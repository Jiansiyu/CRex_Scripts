import os
import sys
from fpdf import FPDF
import csv
from dateutil import tz
from dateutil.parser import parse
import datetime
import PIL
def run(runID=0):
    pass
def decodeRunFile(RunListFileName):
    '''
    read the file list from the file 
    add support for multi format support
    formate can be 
    run nmae (or run number) start run ID (optional) end run ID(optional)
    prexRHRS_1111.dat,     0 , 1
    or
    prexRHRS_1111.dat
    or 
    1111     0 , 1
    '''
    runListArray=[]
    FileListArray=[]
    runIDListArray=[]
    if ".txt" in RunListFileName:
        with open(RunListFileName) as fileio:
            for line in fileio:
                if '#' not in line:
                    runListArray.append([item for item in line.strip().split(',')])
                else:
                    runListArray.append([item for item in line.strip().split('#')[0].strip().split(',')])
    for line in runListArray:
        if len(line) == 2 and len(line[0]) > 2:
            print('Currently unsupport')
        elif len(line) == 3 and len(line[0]) > 2:
            print('Currently unsupport')
        elif len(line) ==1 and len(line[0]) > 2:
            runIDListArray.append(line[0])
    return runIDListArray
def GetRunStartTime(runID=0,runLogFile="./runList.csv"):

    if os.path.isfile(runLogFile):
        with open(runLogFile) as f:
            reader=csv.DictReader(f)
            logList=list(reader)
    else:
        print('CAN NOT FIND FILE {}'.format(runLogFile))

    ET = tz.gettz('US/Eastern')
    CT = tz.gettz('US/Central')
    MT = tz.gettz('US/Mountain')
    PT = tz.gettz('US/Pacific')

    us_tzinfos = {'CST': CT, 'CDT': CT,
                    'EST': ET, 'EDT': ET,
                    'MST': MT, 'MDT': MT,
                    'PST': PT, 'PDT': PT}


    for item in logList:
        if runID == int(item["runID"]):
            startTime=item["StartTimestamp"]
            dt_pst = parse(startTime)
            return dt_pst.strftime("%Y-%m-%d %H-%M-%S")
    return None
    
    
    # print(dt_pst.strftime("%Y-%m-%d %H-%M-%S"))
def GetRunTarget(runID=0,runLogFile="./runList.csv"):
    if os.path.isfile(runLogFile):
        with open(runLogFile) as f:
            reader=csv.DictReader(f)
            logList=list(reader)
    else:
        print('CAN NOT FIND FILE {}'.format(runLogFile))
    
    for item in logList:
        if runID == int(item["runID"]):
            return item["production_target_type"]
    return None

def getImageSize(imageFile):
    cover=PIL.Image.open(imageFile)
    width, height=cover.size
    return width,height


if __name__ == "__main__":
    runList=[]
    if len(sys.argv) == 2:
        if os.path.isfile(sys.argv[1]):
            for runID in decodeRunFile(RunListFileName=sys.argv[1]):
                runList.append(runID)
    
    runList.sort()
    for runID in runList:
        runCommand='analyzer -b -q \'TargetCheck.C (\'{runID}\') \''.format(runID=runID)
        os.system(runCommand)
    
    # generate the report
    pdf=FPDF(orientation = 'L', unit = 'mm', format='A4')
    for runID in runList:
        imageFileName='./carbonCheck/target_{}.jpg'.format(runID)
        if not os.path.isfile(imageFileName):
            print("[WORNING]:: Can NOT find {}".format(imageFileName))
            continue
        width, height =getImageSize(imageFileName)
        width, height = float(width * 0.264583), float(height * 0.264583)
        pdf_size = {'P': {'w': 210, 'h': 297}, 'L': {'w': 297, 'h': 210}}
        orientation = 'P' if width < height else 'L'
        width = width if width < pdf_size[orientation]['w'] else pdf_size[orientation]['w']
        height = height if height < pdf_size[orientation]['h'] else pdf_size[orientation]['h']
        pdf.add_page(orientation=orientation)
        pdf.image(imageFileName,0, 0, width, height)

        # ADD THE TARGET information
        target=GetRunTarget(runID=int(runID))
        print(target)
        pdf.set_font('Arial', 'B', 20)
        pdf.cell(width,0 , "Target:{} Time: {}".format(GetRunTarget(runID=int(runID)),GetRunStartTime(runID=int(runID))))
    
    print("Creating the PDF file")

    #Generate FileName
    HRS='LHRS'
    targetUsed=''
    if int(runList[0]) > 20000:
        HRS='RHRS'
    for item in runList:
        if GetRunTarget(runID=int(item)):
            targetUsed=GetRunTarget(runID=int(item))
            break

    pdf.output("./TargetReport_{}_{}.pdf".format(HRS,targetUsed.replace(' ','')),"F")



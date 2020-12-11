# Scripts instruction

Most of the code are root or Hall A analyzer micros

1. [Optics automatic cut tool](#cutpro-c)
2. [Pointing Measurement Script](#pointing-measurement-script)


# Usage

#### Optics Dp check Scripts
Used for check the Carbon with the Water Target Pointing Measurement

#### cutPro.C
Used for cut the Sieve Holes. It also can search for the ground state peak, and first excited state peak, and apply cut on it. It also have the scripts to generate the average focal plane variables.

```c
analyzer 
.L cutPro.C
cuPro()

```

#### Get the HRS Pointing Measurement Result 
Optional requrements:
1. Get the average Beam E information [file](https://github.com/Jiansiyu/GeneralScripts/blob/master/halog/beamE.txt). 

```shell script
OpticsGraphicCutPro()
```

#### Plot the focal Plane Variables with Cut 

```
OpticsFocalVariableCheck.C
```

* Get the Focal Plane x y theta and phi plot

```

```



#### plot the beam E information (Moved to [scripts](https://github.com/Jiansiyu/GeneralScripts/tree/master/halog))
![beamE](https://github.com/Jiansiyu/GeneralScripts/blob/master/halog/result/BeamE1696.jpg)
* Before 2020, the beam E information is not write to the data file. To extract the beam E information
on apar@aonl:

```
myget -b "2019-12-16 04:11:14" -e "2019-12-16 04:14:12" -c HALLA:p > RHRS_21739_BeamE.txt

```
--- 

* To plot the information(need to change the path maybe):

```
python3 getBeamE.py
```

### Pointing Measurement Script
#### Regular mode 

Measure the Central Sieve Pointing Angle (HRS pointing angle)
```c++
OpticsGraphicCutH2O.C +
OpticsGraphicCutProH20(UInt_t runID)

```

#### Calculate All the Sieve holes for single Pointing Run 

```c++
OpticsGraphicCutH2O.C +

getAllSievePointing(UInt_t runID)
```


## PRexCRex_cut branch 

modify the cutPro code, to handle both PRex and CRex data

Other change 
* plot the origional cut if exist
* after the user click the central, select a more accurate center before pass to the cunter selection 

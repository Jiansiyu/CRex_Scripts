# Target Check Scripts 
Used for Check PRex/CRex Target thickness.

## Get the Target Average thickness 
### runList.txt 

```c
id, runNumber # ID is the run label, run Number
              # id =0 will be take as the reference run 
```
Get the average thickness:
```
$ analyzer
analyzer [0]: TargetThicknessCal("runlis.txt")
```

## Get the run Charge plot
### Get the target accumulate charge
Those command need to run on 'apar@aonl1'. The timestamp are used for get the charge for each run. The output structure CAN NOT change.
```c

```

# Contact me
Siyu Jian @ UVa [<jiansiyu@gmail.com>](mailto:sj9va@virginia.edu)

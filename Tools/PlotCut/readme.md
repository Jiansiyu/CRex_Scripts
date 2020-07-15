# Scripts instruction

Most of the code are root or Hall A analyzer micros

## Plot the focal Plane Variables with Cut 
```
OpticsFocalVariableCheck.C
```

* Get the Focal Plane x y theta and phi plot

```

```



## plot the beam E information 
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



##PRexCRex_cut branch 

modify the cutPro code, to handle both PRex and CRex data

Other change 
* plot the origional cut if exist
* after the user click the central, select a more accurate center before pass to the cunter selection 








# Other reference 

### Target Thickness Check 
* Logbook[https://prex.jlab.org/wiki/index.php/Target_density_variation_check]
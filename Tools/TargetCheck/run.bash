yourfilenames=`ls ./runList/*.txt`
for eachfile in $yourfilenames
do
   python3 run.py $eachfile
done
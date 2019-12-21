#
# bashScripts used for analysis the GEM and VDC data
#
RAWFILENAME=$1



pkill -f tmux
rm ./replay/finish.txt
./scripts/tmuxthreadkiller.sh &
# start tmux session
tmux has-session -t  mlabviet -d

if [ $? != 0 ]; then
	
	tmux new-session -s mlabviet -n os -d
	tmux split-window -h -t mlabviet
	tmux split-window -v -t mlabviet:os.0
	tmux send-keys -t mlabviet:os.0 'cd /home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/ClumAlign/replay && source  loaddb.sh && analyzer' C-m
	sleep 3
	tmux send-keys -t mlabviet:os.0 ".x replay.C (${RAWFILENAME})" C-m 
	tmux send-keys -t mlabviet:os.2 "watch -n 1 tail log.log" C-m
fi
sleep 3
tmux send-keys -t mlabviet:os.1 "nano temp/temp.log" C-m
tmux a
#if [ test -f "finish.txt" ]; then 
	tmux send-keys -t mlabviet:os.1 " End of run" C-m
#fi
#pkill tmux

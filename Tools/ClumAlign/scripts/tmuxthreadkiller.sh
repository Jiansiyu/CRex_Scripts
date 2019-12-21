#!/bin/bash

cd ./replay

while ! test -f "finish.txt" ;
 do 
 	
	tmux send-keys -t mlabviet:os.1 " waiting for the analyze process ...." C-m
	sleep 5
done

tmux send-keys -t mlabviet:os.1 " Analyze process done, kill tmux ....." C-m
pkill tmux

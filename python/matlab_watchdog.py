# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 18:25:48 2020

@author: observer
"""

# reference: https://stackoverflow.com/questions/89228/calling-an-external-command-from-python#2251026
# reference: https://datatofish.com/python-script-windows-scheduler/
#reference: https://stackoverflow.com/questions/4249542/run-a-task-every-x-minutes-with-windows-task-scheduler

import subprocess
import psutil

pid = None
process = None

list = [p for p in psutil.process_iter()]

if "MATLAB.exe" in (p.name() for p in list): # check if matlab is in the list of running processes
    # print("MATLAB is running...")
    process = next((p for p in list if p.name()=='MATLAB.exe'), None)
    pid = process.pid
else:
    subprocess.Popen('matlab -useStartupFolderPref') # the optional argument tells matlab to start at the preferred startup folder defined in the preferences menu
    

# add a check of the logfiles and kill matlab if it is not updated often enough
if pid is not None:
    print("MATLAB process id is " + str(pid))


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
import os
from datetime import datetime, timedelta
import math
import time

pid = None
process = None

list = [p for p in psutil.process_iter()]

if "MATLAB.exe" in (p.name() for p in list): # check if matlab is in the list of running processes
    # print("MATLAB is running...")
    process = next((p for p in list if p.name()=='MATLAB.exe'), None)
    pid = process.pid
else:
    print('MATLAB is not working! Starting a new instance... ')
    subprocess.Popen('matlab -useStartupFolderPref') # the optional argument tells matlab to start at the preferred startup folder defined in the preferences menu
    

if pid is not None:
    print("MATLAB process id is " + str(pid))

    log_file = os.environ.get('LOGFILE', None)
    
    if log_file is not None:
        t_folder = datetime.utcnow()
        if t_folder.hour<12:
            t_folder -= timedelta(days=1)  # until noon UTC we count it as last night 
    
        folder_date_str = str(t_folder.date())
    
        log_file = os.path.join(os.environ["DATA"], 
                                'WFAST/logfiles', 
                                folder_date_str, 
                                folder_date_str+'_'+log_file+'.txt')
    
    if os.path.exists(log_file):
        
        with open(log_file) as f:
            lines = f.readlines()
        
        lines = [l.strip() for l in lines if len(l)>0]
        lines.reverse()
        
        t0 = None
        today_str = str(datetime.today().date())
        for l in lines:
            try:
                t0 = datetime.strptime(today_str + ' ' + l[0:12] + 'UTC', "%Y-%m-%d %H:%M:%S.%f%Z")
                break
            except ValueError:
                pass
        
        if t0 is not None: 
            dt = (datetime.utcnow()-t0).seconds
            if dt>900:  # more than 15 minutes 
                print(f'MATLAB has not updated in more than {math.floor(dt/60)} minutes!')
                
                process.terminate()  # kill the MATLAB! 
                
                time.sleep(10)
                
                # start a new MATLAB instance... 
                print('Starting up a new instance of MATLAB.')
                subprocess.Popen('matlab -useStartupFolderPref') # the optional argument tells matlab to start at the preferred startup folder defined in the preferences menu
    
                
            else: 
                print(f'MATLAB has updated "{os.environ["LOGFILE"]}" logfile {math.floor(dt/60)} minutes ago...')
        
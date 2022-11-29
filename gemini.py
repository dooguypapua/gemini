'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : gemini.py                       |
description      : in silico python companion      |
author           : dooguypapua                     |
lastmodification : 20210628                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import sys
import os
import shutil
import time
import geminiset
from datetime import datetime
from geminini import printcolor, init_gemini_dico, listfct, infofct, build_gemini_dico, title
from geminini import fct_params_from_args, gemini_header, gemini_fct_header, get_sys_info, printDico
from geminiannot import *
from geminicluster import *
from geminidl import *
from geminiparse import *
from geminiphage import *
from geminiphylo import *
from geminiplot import *


if __name__ == "__main__":
    dicoGemini = init_gemini_dico()
    if len(sys.argv) == 1 or sys.argv[1] in ["-h", "-help", "--help"]:
        listfct()
    else:
        startTime = time.time()
        fctName = sys.argv[1]
        # function dispatcher
        if fctName not in dicoGemini:
            listfct()
            printcolor("[ERROR: invalid task] \""+fctName+"\" does not exist\n", 1, "212;64;89", "None", True)
            exit("\n")
        # Display help
        if len(sys.argv) == 2 or len(set(["-h", "-help", "--help"]) & set(sys.argv)) > 0:
            infofct(fctName)
            sys.exit(0)
        # Check and send params to function
        build_gemini_dico(sys.argv)
        lstArgValue = fct_params_from_args(fctName)
        gemini_header()
        gemini_fct_header(fctName)
        # Check CPU
        slurmBool, cpu, memMax, memMin = get_sys_info()
        if "--cpu" in sys.argv:
            cpu = sys.argv[sys.argv.index("--cpu")+1]
        # Specific environment import
        if slurmBool is True:
            sys.path.append('/home/umr8227/gv/dgoudenege/script/python_lib/')
            sys.path.append('/home/umr8227/gv/dgoudenege/.local/lib/python3.9/site-packages/')
            sys.path.append('/shared/software/miniconda/envs/python-pytorch-tensorflow-3.9-1.11.0-2.6.2/lib/python3.9/site-packages/')
        # Display execution infos
        printcolor("⌛ "+str(datetime.fromtimestamp(startTime).strftime("%H:%M:%S [%d-%m-%y]"))+"\n")
        printcolor("➰ SLURM: "+str(slurmBool)+"\n")
        printcolor("➰ PID  : "+str(os.getpid())+"\n")
        printcolor("➰ CPU  : "+str(cpu)+"\n")
        printcolor("➰ MEM  : "+str(memMin)+"-"+str(memMax)+"Go\n")
        if dicoGemini[fctName]['pathTMP'] != "":
            printcolor("➰ TMP  : "+dicoGemini[fctName]['pathTMP']+"\n")
        # Set gemini variables
        geminiset.fctName = fctName
        geminiset.title = dicoGemini[fctName]['title']
        geminiset.slurmBool = slurmBool
        geminiset.cpu = cpu
        geminiset.memMax = memMax
        geminiset.memMin = memMin
        geminiset.pathTMP = dicoGemini[fctName]['pathTMP']
        # Launch function
        returnObj = globals()[fctName](*lstArgValue)
        if returnObj is not None and "-j" not in sys.argv:
            printcolor(returnObj)
        printcolor("✅ Completed"+"\n")
        # Delete temp folder
        if dicoGemini[fctName]['pathTMP'] != "":
            shutil.rmtree(dicoGemini[fctName]['pathTMP'])
        # Execution time
        elapseDay = int((time.time()-startTime)/86400)
        execTime = str(elapseDay)+"days"+time.strftime("%Hh%Mm%Ss", time.gmtime((time.time()-startTime)))
        printcolor("⌛ "+str(execTime)+"\n")
        title(text="✅", pbar=None)
    sys.exit(0)

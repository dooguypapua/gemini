'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ phylo
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminiphylo.py                  |
description      : gemini phylogenetic functions   |
author           : dooguypapua                     |
lastmodification : 20210830                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import os
import sys
import re
import shutil
import geminiset
import xlsxwriter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from yaspin import yaspin
from typing import Tuple
from yaspin.spinners import Spinners
from ete3 import Tree, TreeStyle, NodeStyle
from geminini import printcolor, fct_checker, get_input_files, path_converter, get_gemini_path, get_sys_info
from geminini import dump_json, load_json, launch_threads, linear_gradient, reverse_complement, title, exit_gemini
from geminicluster import mmseqs_rbh, make_rbhcluster_dict
from geminiparse import make_fasta_dict, gbk_to_fna


@fct_checker
def mash_matrix(pathIN: str, pathOUT: str, sketchSize: int = 10000, ext: str = ".fna") -> Tuple[str, str, int, str]:
    '''
     ------------------------------------------------------------
    |                    MASH DISTANCE MATRIX                    |
    |------------------------------------------------------------|
    |                    Mash distance matrix                    |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN    : path of input fna files or folder (required)|
    |    pathOUT   : path of output sketch folder (required)     |
    |    sketchSize: fastaANI fragments length (default=10000)   |
    |    ext       : extension of input files (default=.fna)     |
     ------------------------------------------------------------
    |  Sketch size = number of min-hashes that are kept.         |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "mash_matrix", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: mash_matrix]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    if pathIN[-1] == "/":
        pathINreplace = pathIN
    else:
        pathINreplace = pathIN+"/"
    pathOUTSKETCH = pathOUT+"/sketch_"+str(sketchSize)
    os.makedirs(pathOUTSKETCH, exist_ok=True)
    pathOUTDIST = pathOUT+"/dist_"+str(sketchSize)
    os.makedirs(pathOUTDIST, exist_ok=True)
    PATHJSON = pathOUT+"/mash_matrix_"+str(sketchSize)+".json"
    PATHPNG = pathOUT+"/mash_matrix_"+str(sketchSize)+".png"
    PATHSVG = pathOUT+"/mash_matrix_"+str(sketchSize)+".svg"
    PATHXLSX = pathOUT+"/mash_matrix_"+str(sketchSize)+".xlsx"
    lstFiles.sort()
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    if "." not in ext:
        ext = "."+ext
    # ***** INIT matrix ***** #
    printcolor("♊ Init matrix"+"\n")
    if os.path.isfile(PATHJSON):
        dicoMATRIX = load_json(PATHJSON)
    else:
        dicoMATRIX = {}
    lstDistOut = os.listdir(pathOUTDIST)
    pbar = tqdm(total=len(lstDistOut), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for distOut in lstDistOut:
        orgName = distOut.replace(".tsv", "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        pathDIST = pathOUTDIST+"/"+distOut
        TMP = open(pathDIST, 'r')
        lstLines = TMP.read().split("\n")[: -1]
        TMP.close()
        for line in lstLines:
            splitLine = line.split("\t")
            try:
                dicoMATRIX[orgName][splitLine[0]] = float(splitLine[1])
            except KeyError:
                dicoMATRIX[orgName] = {splitLine[0]: float(splitLine[1])}
        pbar.update(1)
        title("Init matrix", pbar)
    pbar.close()
    # Save matrix
    dump_json(dicoMATRIX, PATHJSON)
    # ***** SKETCH genomes ***** #
    printcolor("♊ Sketch genomes"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFNA in lstFiles:
        orgName = os.path.basename(pathFNA).replace(ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        pathSKETCH = pathOUTSKETCH+"/"+orgName+".msh"
        if not os.path.isfile(pathSKETCH) or os.path.getsize(pathSKETCH) == 0:
            cmdSKETCH = dicoGeminiPath['TOOLS']['mash']+" sketch -o "+pathSKETCH.replace(".msh", "")+" -p "+str(cpu)+" -s "+str(sketchSize)+" "+pathFNA+" > /dev/null 2>&1"
            os.system(cmdSKETCH)
        pbar.update(1)
        title("Sketch genomes", pbar)
    pbar.close()
    # ***** INIT distance ***** #
    printcolor("♊ Init distance"+"\n")
    dicoThread = {}
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFNA1 in lstFiles:
        orgName1 = os.path.basename(pathFNA1).replace(ext, "")
        pbar.set_description_str(orgName1+" ".rjust(maxpathSize-len(orgName1)))
        pathSKETCH1 = pathOUTSKETCH+"/"+orgName1+".msh"
        pathDIST = pathOUTDIST+"/"+orgName1+".tsv"
        pathTMPLST = pathTMP+"/"+orgName1+".lst"
        strpathSKETCH2 = ""
        for pathFNA2 in lstFiles:
            orgName2 = os.path.basename(pathFNA2).replace(ext, "")
            pathSKETCH2 = pathOUTSKETCH+"/"+orgName2+".msh"
            if orgName1 not in dicoMATRIX or orgName2 not in dicoMATRIX[orgName1]:
                strpathSKETCH2 += pathSKETCH2+"\n"
        if strpathSKETCH2 != "":
            TMPLIST = open(pathTMPLST, 'w')
            TMPLIST.write(strpathSKETCH2)
            TMPLIST.close()
            cmdDIST = dicoGeminiPath['TOOLS']['mash']+" dist -s "+str(sketchSize)+" -l "+pathSKETCH1+" "+pathTMPLST + \
                " | cut -f 2, 3,4, 5 | awk \'{gsub(\""+pathINreplace+"\", \"\");gsub(\""+ext+"\", \"\")}1\' >> "+pathDIST
            dicoThread[orgName1] = {"cmd": cmdDIST, "returnstatut": None, "returnlines": []}
        pbar.update(1)
        title("Init dist", pbar)
    pbar.close()
    # ***** COMPUTE distance ***** #
    if len(dicoThread) != 0:
        spinner = yaspin(Spinners.aesthetic, text="♊ Compute distance", side="right")
        spinner.start()
        title("Compute dist", None)
        launch_threads(dicoThread, "mash_matrix", cpu, pathTMP)
        spinner.stop()
        printcolor("♊ Compute distance"+"\n")
        # ***** CREATE matrix ***** #
        printcolor("♊ Update Matrix"+"\n")
        lstDistOut = os.listdir(pathOUTDIST)
        pbar = tqdm(total=len(lstDistOut), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        for distOut in lstDistOut:
            orgName = distOut.replace(".tsv", "")
            pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
            pathDIST = pathOUTDIST+"/"+distOut
            TMP = open(pathDIST, 'r')
            lstLines = TMP.read().split("\n")[: -1]
            TMP.close()
            for line in lstLines:
                splitLine = line.split("\t")
                try:
                    dicoMATRIX[orgName][splitLine[0]] = float(splitLine[1])
                except KeyError:
                    dicoMATRIX[orgName] = {splitLine[0]: float(splitLine[1])}
            pbar.update(1)
            title("Update matrix", pbar)
        pbar.close()
        # Save matrix
        dump_json(dicoMATRIX, PATHJSON)
    # ***** CLUSTERING matrix ***** #
    printcolor("♊ Clustering"+"\n")
    # Clustering
    df = pd.DataFrame(dicoMATRIX)
    figsize = len(df) / np.sqrt(len(df))
    cg = sns.clustermap(df, cmap='Spectral', figsize=(figsize, figsize), tree_kws={'linewidths': 1.5}, dendrogram_ratio=0.05, xticklabels=True, yticklabels=True, linewidths=0.0, rasterized=True)
    # Plot clustered heatmap
    cg.ax_cbar.remove()
    plt.tick_params(axis='both', which='major', labelsize=1.5, width=0.2)
    fig = plt.gcf()  # or by other means, like plt.subplots
    figsize = fig.get_size_inches()
    fig.set_size_inches(figsize * 1.5)  # scale current size by 1.5
    plt.savefig(PATHPNG, dpi=300)
    plt.savefig(PATHSVG)
    # Retrieve ordered ticks label
    newColums = df.columns[cg.dendrogram_col.reordered_ind]
    newIndexs = df.index[cg.dendrogram_row.reordered_ind]
    newData = df.loc[newIndexs, newColums]
    orderedOrg = list(newData.keys())
    # Write output excel
    workbook = xlsxwriter.Workbook(PATHXLSX)
    worksheet = workbook.add_worksheet()
    headerFormat1 = workbook.add_format({'align': 'center', 'valign': 'bottom', 'border': 1, 'font_size': 11})
    headerFormat2 = workbook.add_format({'align': 'center', 'valign': 'bottom', 'border': 1, 'font_size': 11})
    headerFormat1.set_rotation(90)
    defaultRowFormat = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 8, 'text_wrap': True, 'bg_color': "#FFFFFF"})
    HEX_list = linear_gradient("#FFFFFF", "#D40000", 101)
    lstRowFormat = []
    for i in range(101):
        lstRowFormat.append(defaultRowFormat)
    for i in range(len(HEX_list)):
        rowFormat = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 8, 'text_wrap': True, 'bg_color': HEX_list[i].upper()})
        lstRowFormat[i] = rowFormat
    dicoWidth = {0: 0}
    dicoHeight = {0: 0}
    worksheet.write(0, 0, "", headerFormat2)
    col = 1
    for orgName1 in orderedOrg:
        worksheet.write(0, col, orgName1, headerFormat1)
        dicoHeight[0] = max(dicoHeight[0], len(orgName1))
        col += 1
    row = 1
    for orgName1 in orderedOrg:
        worksheet.write(row, 0, orgName1, headerFormat2)
        dicoWidth[0] = max(dicoWidth[0], len(orgName1))
        col = 1
        for orgName2 in orderedOrg:
            value = dicoMATRIX[orgName1][orgName2]
            worksheet.write(row, col, value, lstRowFormat[int(value)])
            try:
                dicoWidth[col] = max(dicoWidth[col], len(str(value)))
            except KeyError:
                dicoWidth[col] = len(str(value))
            col += 1
        row += 1
    # Adjust row height and column width
    for row in dicoHeight:
        worksheet.set_row(row, dicoHeight[row]*2.5)
    for col in dicoWidth:
        worksheet.set_column(col, col, dicoWidth[col])
    workbook.close()


@fct_checker
def fastani_db(pathIN: str, pathIN2: str, pathJSON: str, fragLen: int = 3000, ext: str = ".fna") -> Tuple[str, str, str, int, str]:
    '''
     ------------------------------------------------------------
    |                         FASTANI DB                         |
    |------------------------------------------------------------|
    |                Create fastANI JSON database                |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN   : path of input fna files or folder (required) |
    |    pathIN2  : path of in/output fastANI folder (required)  |
    |    pathJSON : path of output files (required)              |
    |    fragLen  : fastaANI fragments length (default=3000)     |
    |    ext      : extension of input files (default=.fna)      |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "fastani_db", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: fastani_db]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathIN2 == "":
        printcolor("[ERROR: fastani_db]\nMissing '-i2'pathIN2\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathIN2 = path_converter(pathIN2)
    os.makedirs(pathIN2, exist_ok=True)
    lstFiles2, maxpathSize2 = get_input_files(pathIN2, "fastani_db", [".txt"])
    lstFiles.sort()
    lstFiles2.sort()
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    if pathJSON == "":
        printcolor("[ERROR: fastani_db]\nMissing '-j'pathJSON\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathJSON = path_converter(pathJSON)
    # ***** PARSE fastANI files ***** # (orgHit<tab>ani)
    printcolor("♊ Parse fastANI files"+"\n")
    dicoANI = {}
    pbar = tqdm(total=len(lstFiles2), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFASTANI in lstFiles2:
        orgName1 = os.path.basename(pathFASTANI).replace(".txt", "")
        pbar.set_description_str(orgName1+" ".rjust(maxpathSize-len(orgName1)))
        dicoANI[orgName1] = {}
        IN = open(pathFASTANI, 'r')
        lstLines = IN.read().split("\n")
        IN.close()
        for line in lstLines:
            if line != "":
                splitLine = line.split("\t")
                orgName2 = splitLine[0]
                ani = float(splitLine[1])
                dicoANI[orgName1][orgName2] = ani
        pbar.update(1)
        title("Parse prev", pbar)
    pbar.close()
    # ***** CHECK missing fastANi couple ***** #
    printcolor("♊ Check missing results"+"\n")
    # Retrieve organism name list
    setOrgName = set()
    for fna in lstFiles:
        setOrgName.add(os.path.basename(fna).replace("."+ext, "").replace(ext, ""))
    # Retrieve missing ani results
    dicoThread = {}
    for orgName1 in setOrgName:
        for orgName2 in setOrgName:
            if orgName1 not in dicoANI or orgName2 not in dicoANI[orgName1]:
                pathTMPOUT = pathTMP+"/"+orgName1+"#####"+orgName2+".txt"
                cmdFASTANI = dicoGeminiPath['TOOLS']['fastani']+" -q "+pathIN+"/"+orgName1+".fna -r "+pathIN+"/"+orgName2+".fna -o "+pathTMPOUT+" --threads 1 --fragLen "+str(fragLen)+" > /dev/null 2>&1"
                dicoThread[orgName1+"_"+orgName2] = {"cmd": cmdFASTANI, "returnstatut": None, "returnlines": []}
    printcolor("⏩ "+str(len(dicoThread))+" missing couples"+"\n")
    # ***** LAUNCH missing fastANI ***** #
    if len(dicoThread) != 0:
        spinner = yaspin(Spinners.aesthetic, text="♊ Launch missing fastANI", side="right")
        spinner.start()
        title("fastANI", None)
        launch_threads(dicoThread, "fastani_db", cpu, pathTMP)
        spinner.stop()
        printcolor("♊ Launch missing fastANI"+"\n")
        # ***** Parse new results ***** #
        printcolor("♊ Read new results"+"\n")
        pbar = tqdm(total=len(dicoThread), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
        for file in os.listdir(pathTMP):
            if "#####" in file:
                pathFILE = pathTMP+"/"+file
                orgName1 = file.split("#####")[0]
                orgName2 = file.split("#####")[1].replace(".txt", "")
                # Parse temp results files
                IN = open(pathFILE, 'r')
                try:
                    ani = float(IN.read().split("\n")[0].split("\t")[2])
                except IndexError:
                    ani = -1
                IN.close()
                # Add to dictionnary
                dicoANI[orgName1][orgName2] = ani
                # Cat to previous file
                pathORGRES = pathIN2+"/"+orgName1+".txt"
                if not os.path.isfile(pathORGRES):
                    os.system("echo \""+orgName2+"\t"+str(ani)+"\""+" > "+pathORGRES)
                else:
                    os.system("echo \""+orgName2+"\t"+str(ani)+"\""+" >> "+pathORGRES)
                pbar.update(1)
                title("Parse new", pbar)
        pbar.close()
    # ***** SAVE fastANI JSON file ***** #
    printcolor("♊ Dump fastANI JSON"+"\n")
    dump_json(dicoANI, pathJSON)


@fct_checker
def best_gene_tree_topology(pathIN1: str, pathIN2: str, pathIN3: str, pathOUT: str, outgroup: str, idThrClust: int = 80, covThrClust: int = 80, ext1: str = ".ffn", ext2: str = ".faa") -> Tuple[str, str, str, str, str, int, int, str, str]:
    '''
     ------------------------------------------------------------
    |                   BEST GENE TREE TOPOLOGY                  |
    |------------------------------------------------------------|
    |        Construct tree and search best group topology       |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN1    : path of FFN files or folder (required)     |
    |    pathIN2    : path of FAA files or folder (required)     |
    |    pathIN3    : path of TSV org\tgrp\tcolor (required)     |
    |    pathOUT    : path of output files (required)            |
    |    outgroup   : outgroup organism (required)               |
    |    idThrClust : %id clustering threshold (default=80)      |
    |    covThrClust: %cov clustering threshold (default=80)     |
    |    ext1        : extension of input files (default=.ffn)   |
    |    ext2        : extension of input files (default=.faa)   |
     ------------------------------------------------------------
    '''
    lstFilesFFN, maxpathSize1 = get_input_files(pathIN1, "best_gene_tree_topology", [ext1])
    if len(lstFilesFFN) == 0:
        printcolor("[ERROR: best_gene_tree_topology]\nAny input FFN files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    lstFilesFAA, maxpathSize2 = get_input_files(pathIN2, "best_gene_tree_topology", [ext2])
    if len(lstFilesFAA) == 0:
        printcolor("[ERROR: best_gene_tree_topology]\nAny input FAA files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathIN3 = path_converter(pathIN3)
    if not os.path.isfile(pathIN3) or os.path.getsize(pathIN3) == 0:
        printcolor("[ERROR: best_gene_tree_topology]\nAny input TSV file found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Check same both FFN and FAA
    if "." not in ext1:
        ext1 = "."+ext1
    if "." not in ext2:
        ext2 = "."+ext2
    lstMissing = []
    for fileFFN in lstFilesFFN:
        orgName = os.path.basename(fileFFN).replace(ext1, "")
        if not os.path.isfile(pathIN2+"/"+orgName+ext2):
            lstMissing.append("FAA for \""+orgName+"\"")
    for fileFAA in lstFilesFAA:
        orgName = os.path.basename(fileFAA).replace(ext2, "")
        if not os.path.isfile(pathIN1+"/"+orgName+ext1):
            lstMissing.append("FFN for \""+orgName+"\"")
    if len(lstMissing) > 0:
        printcolor("[ERROR: best_gene_tree_topology]\nMissing input files found\n", 1, "212;64;89", "None", True)
        for miss in lstMissing:
            printcolor("  "+miss+"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Initialization
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    if pathOUT == "":
        printcolor("[ERROR: best_tree_topology]\nMissing '-o'pathOUT\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if outgroup == "":
        printcolor("[ERROR: best_tree_topology]\nMissing '-outgrp'outgroup\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if not os.path.isfile(pathIN1+"/"+outgroup+"."+ext1.replace(".", "")):
        printcolor("[ERROR: best_tree_topology]\nAny outgroup input file found \""+pathIN1+"/"+outgroup+"."+ext1.replace(".", "")+"\"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    pathOUTfasta = pathOUT+"/core_gene"
    pathOUTalign = pathOUT+"/core_align"
    pathOUTtree = pathOUT+"/core_tree"
    pathOUTpng = pathOUT+"/core_besttree_png"
    os.makedirs(pathOUT, exist_ok=True)
    os.makedirs(pathOUTfasta, exist_ok=True)
    os.makedirs(pathOUTalign, exist_ok=True)
    os.makedirs(pathOUTtree, exist_ok=True)
    # ***** Get organism groups ***** #
    printcolor("♊ Read groups"+"\n")
    TSV = open(pathIN3, 'r')
    lstLines = TSV.read().split("\n")
    TSV.close()
    dicoOrgToGroup = {}
    dicoOrgToColor = {}
    dicoGroupToOrg = {}
    for line in lstLines:
        if line != "" and line[0] != "#":
            try:
                orgName = line.split("\t")[0]
                group = line.split("\t")[1]
                color = line.split("\t")[2]
            except IndexError:
                printcolor("[ERROR: best_tree_topology]\nParsing TSV file\n", 1, "212;64;89", "None", True)
                exit_gemini()
            dicoOrgToGroup[orgName] = group
            dicoOrgToColor[orgName] = color
            if group not in dicoGroupToOrg:
                dicoGroupToOrg[group] = set()
            dicoGroupToOrg[group].add(orgName)
    # ***** Get gene sequence and Keep organism name per locusTag ***** #
    printcolor("♊ Read sequences"+"\n")
    dicoLTtoOrg = {}
    dicoLTtoGeneSeq = {}
    pbar = tqdm(total=len(lstFilesFFN), ncols=50+maxpathSize1, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for fileFFN in lstFilesFFN:
        orgName = os.path.basename(fileFFN).replace(ext1, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize1-len(orgName)))
        dicoFFN = make_fasta_dict(fileFFN)
        for key in dicoFFN:
            lt = key.split("|")[0]
            dicoLTtoOrg[lt] = orgName
            dicoLTtoGeneSeq[lt] = dicoFFN[key]
        pbar.update(1)
        title("Read seq", pbar)
    pbar.close()
    # ***** Protein clustering ***** #
    pathRBHJSON = pathOUT+"/rbh.json"
    if os.path.isfile(pathRBHJSON) and os.path.getsize(pathRBHJSON) > 0:
        dicoCLUSTER = load_json(pathRBHJSON)
    else:
        # Launch rbh
        mmseqs_rbh(pathIN=pathIN2, pathOUT=pathOUT, ref=outgroup, idThrClust=idThrClust, covThrClust=covThrClust, ext=ext2)
        # Parse rbh results
        printcolor("♊ Read rbh clusters"+"\n")
        dicoCLUSTER = {}
        MMSEQS = open(pathOUT+"/"+outgroup+".rbh", 'r')
        lstLines = MMSEQS.read().split("\n")[: -1]
        MMSEQS.close()
        pbar = tqdm(total=len(lstLines), ncols=50+maxpathSize2, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
        # Parse and assign cluster
        for line in lstLines:
            splitLine = line.split("\t")
            ltQuery = splitLine[0].split("|")[0]
            orgQuery = dicoLTtoOrg[ltQuery]
            query = ltQuery+"_____"+orgQuery
            ltTarget = splitLine[1].split("|")[0]
            orgTarget = dicoLTtoOrg[ltTarget]
            target = ltTarget+"_____"+orgTarget
            findCluster = False
            for clusterNum in dicoCLUSTER.keys():
                if query in dicoCLUSTER[clusterNum]:
                    dicoCLUSTER[clusterNum].add(target)
                    findCluster = True
                    break
                if target in dicoCLUSTER[clusterNum]:
                    dicoCLUSTER[clusterNum].add(query)
                    findCluster = True
                    break
            if findCluster is False:
                dicoCLUSTER["cluster"+str(len(dicoCLUSTER)+1).zfill(8)] = set([query, target])
            pbar.update(1)
            title("Parse clusters", pbar)
        pbar.close()
        # Transform set to list
        for clusterNum in dicoCLUSTER:
            dicoCLUSTER[clusterNum] = list(dicoCLUSTER[clusterNum])
        dump_json(dicoCLUSTER, pathRBHJSON)
    # ***** Define core ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Define core genes", side="right")
    spinner.start()
    title("Core genes", None)
    pathCOREJSON = pathOUT+"/core.json"
    if os.path.isfile(pathCOREJSON) and os.path.getsize(pathCOREJSON) > 0:
        dicoCORE = load_json(pathCOREJSON)
    else:
        dicoCORE = {}
        for clusterNum in dicoCLUSTER:
            setOrg = set()
            outLT = ""
            for prot in dicoCLUSTER[clusterNum]:
                lt = prot.split("_____")[0]
                orgName = prot.split("_____")[1]
                setOrg.add(orgName)
                if orgName == outgroup:
                    outLT = lt
            # Filter core cluster
            if len(setOrg) == len(lstFilesFFN) and outLT != "":
                dicoCORE[outLT] = {}
                for prot in dicoCLUSTER[clusterNum]:
                    lt = prot.split("_____")[0]
                    orgName = prot.split("_____")[1]
                    dicoCORE[outLT][orgName] = lt
        dump_json(dicoCORE, pathCOREJSON)
    spinner.stop()
    printcolor("♊ Define core genes"+"\n")
    printcolor("⏩ Found "+str(len(dicoCORE))+" core genes"+"\n")
    # ***** Align core genes ***** #
    printcolor("♊ Align core genes"+"\n")
    pbar = tqdm(total=len(dicoCORE), ncols=50+maxpathSize1, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for outLT in dicoCORE:
        pbar.set_description_str(outLT+" ".rjust(maxpathSize1-len(outLT)))
        # Make gene FASTA
        pathCoreGeneFFN = pathOUTfasta+"/"+outLT+".ffn"
        if not os.path.isfile(pathCoreGeneFFN) or os.path.getsize(pathCoreGeneFFN) == 0:
            OUTcoreGene = open(pathCoreGeneFFN, 'w')
            for orgName in dicoCORE[outLT]:
                lt = dicoCORE[outLT][orgName]
                OUTcoreGene.write(">"+orgName+"\n"+dicoLTtoGeneSeq[lt]+"\n")
            OUTcoreGene.close()
        # Align gene FASTA
        pathCoreGeneALIGN = pathOUTalign+"/"+outLT+"_align.ffn"
        if not os.path.isfile(pathCoreGeneALIGN) or os.path.getsize(pathCoreGeneALIGN) == 0:
            os.system(dicoGeminiPath['TOOLS']['famsa']+" -t "+str(cpu)+" "+pathCoreGeneFFN+" "+pathCoreGeneALIGN+" > /dev/null 2>&1")
        pbar.update(1)
        title("Align core", pbar)
    pbar.close()
    # ***** Make individual core gene tree ***** #
    printcolor("♊ Make genes trees"+"\n")
    pbar = tqdm(total=len(dicoCORE), ncols=50+maxpathSize1, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for outLT in dicoCORE:
        pbar.set_description_str(outLT+" ".rjust(maxpathSize1-len(outLT)))
        pathCoreGeneALIGN = pathOUTalign+"/"+outLT+"_align.ffn"
        pathCoreGeneTREE = pathOUTtree+"/"+outLT+"_GTR_fasttree.nwk"
        if not os.path.isfile(pathCoreGeneTREE) or os.path.getsize(pathCoreGeneTREE) == 0:
            os.system(dicoGeminiPath['TOOLS']['fasttree']+" -quiet -nt -gtr "+pathCoreGeneALIGN+" > "+pathCoreGeneTREE+" 2>/dev/null")
        pbar.update(1)
        title("Make trees", pbar)
    pbar.close()
    # ***** ANALYSE TREE TOPOLOGY ***** #
    printcolor("♊ Analyse topologies"+"\n")
    circularStyle = TreeStyle()
    circularStyle.mode = "c"
    circularStyle.scale = 20
    cptBestTree = 0
    pbar = tqdm(total=len(dicoCORE), ncols=50+maxpathSize1, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for outLT in dicoCORE:
        pbar.set_description_str(outLT+" ".rjust(maxpathSize1-len(outLT)))
        pathCoreGeneTREE = pathOUTtree+"/"+outLT+"_GTR_fasttree.nwk"
        tree = Tree(pathCoreGeneTREE)
        tree.set_outgroup(outgroup)
        lstBadGroup = []
        for group in dicoGroupToOrg:
            groupAncestor = tree.get_common_ancestor(dicoGroupToOrg[group])
            for node in groupAncestor.iter_descendants("postorder"):
                if node.name in dicoOrgToGroup and dicoOrgToGroup[node.name] != group:
                    lstBadGroup.append(group)
        # Draw best tree
        if len(lstBadGroup) == 0:
            pbar.set_description_str(outLT+" (plot)".rjust(maxpathSize1-len(outLT+" (plot)")))
            cptBestTree += 1
            os.makedirs(pathOUTpng, exist_ok=True)
            pathTreePNG = pathOUTpng+"/"+outLT+".png"
            for node in tree.traverse("postorder"):
                if node.name in dicoOrgToColor:
                    nst = NodeStyle()
                    nst["bgcolor"] = dicoOrgToColor[node.name]
                    node.set_style(nst)
            tree.render(pathTreePNG, w=183, units="mm", tree_style=circularStyle)
        pbar.update(1)
        title("Topology", pbar)
    pbar.close()
    if cptBestTree == 0:
        printcolor("⏩ Any gene groups tree topology found"+"\n")
    else:
        printcolor("⏩ Found "+str(cptBestTree)+" gene groups tree topology"+"\n")


@fct_checker
def specific_kmers(pathIN: str, pathIN2: str, pathOUT: str, kmerLen: int = 25, ext: str = ".fna") -> Tuple[str, str, str, int, str]:
    '''
     ------------------------------------------------------------
    |                   SEARCH SPECIFIC KMERS                    |
    |------------------------------------------------------------|
    |         Search kmers specific to a organism group          |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN     : path of input files or folder (required)   |
    |    pathIN2    : path of file with organism list (required) |
    |    pathOUT    : path of output files (required)            |
    |    kmerLen    : kmer length (default=25)                   |
    |    ext        : extension of input files (default=.fna)    |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "specific_kmers", [ext])
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    if "." not in ext:
        ext = "."+ext
    if len(lstFiles) == 0:
        printcolor("[ERROR: specific_kmers]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathIN2 == "":
        printcolor("[ERROR: specific_kmers]\nMissing '-i2'pathIN2\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathIN2 = path_converter(pathIN2)
    # Check organism specific list
    IN2 = open(pathIN2, 'r')
    lstLines = IN2.read().split("\n")
    IN2.close()
    setSelectedOrg = set()
    lstNotFound = []
    for line in lstLines:
        if line != "" and not line[0] == "#":
            pathOrgFile = pathIN+"/"+line+ext
            if not os.path.isfile(pathOrgFile):
                lstNotFound.append(pathOrgFile)
            else:
                setSelectedOrg.add(line)
    if len(lstNotFound) > 0:
        printcolor("[ERROR: specific_kmers]\nSome specific organism not found\n"+"\n".join(lstNotFound)+"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Make output folder
    if pathOUT == "":
        printcolor("[ERROR: specific_kmers]\nMissing '-o'pathOUT\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    pathOUTprimers = pathOUT+"primers.txt"
    # ***** MAKE Kmers ***** #
    printcolor("♊ Make Kmers"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        pathKMERS = pathOUT+"/"+orgName+".k"+str(kmerLen)+".kmers"
        if not os.path.isfile(pathKMERS) or os.path.getsize(pathKMERS) == 0:
            cmdWORDCOUNT = dicoGeminiPath['TOOLS']['wordcount']+" -sequence "+pathFile+" -wordsize "+str(kmerLen)+" -outfile "+pathKMERS+" > /dev/null 2>&1"
            os.system(cmdWORDCOUNT)
        pbar.update(1)
        title("Make kmers", pbar)
    pbar.close()
    # ***** SELECTED ORGANISM KMERS ***** #
    printcolor("♊ Selected organism Kmers"+"\n")
    dicoKmers = {}
    pbar = tqdm(total=len(setSelectedOrg), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        if orgName in setSelectedOrg:
            pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
            pathKMERS = pathOUT+"/"+orgName+".k"+str(kmerLen)+".kmers"
            dicoKmers[orgName] = set()
            KMERS = open(pathKMERS, 'r')
            lstLines = KMERS.read().split("\n")
            KMERS.close()
            for line in lstLines:
                if line != "":
                    kmer = line.split("\t")[0]
                    dicoKmers[orgName].add(kmer)
                    dicoKmers[orgName].add(reverse_complement(kmer))
            pbar.update(1)
            title("Select kmers", pbar)
    pbar.close()
    # ***** SELECTED ORGANISM COMMON KMERS ***** #
    printcolor("♊ Common Kmers"+"\n")
    pbar = tqdm(total=len(setSelectedOrg), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    pivotOrg = list(setSelectedOrg)[0]
    setIntersectKmers = dicoKmers[pivotOrg]
    for orgName in setSelectedOrg:
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        setIntersectKmers = setIntersectKmers.intersection(dicoKmers[orgName])
        pbar.update(1)
        title("Common kmers", pbar)
    pbar.close()
    printcolor("⏩ Found "+str(len(setIntersectKmers))+" common Kmers"+"\n")
    # ***** ANALYSE Kmers ***** #
    printcolor("♊ Specific Kmers"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFile in lstFiles:
        if "Vibrio" not in pathFile:
            continue
        if len(setIntersectKmers) == 0:
            break
        orgName = os.path.basename(pathFile).replace(ext, "")
        if orgName not in setSelectedOrg:
            pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName))+"/"+str(len(setIntersectKmers))+"/")
            pathKMERS = pathOUT+"/"+orgName+".k"+str(kmerLen)+".kmers"
            KMERS = open(pathKMERS, 'r')
            lstLines = KMERS.read().split("\n")
            KMERS.close()
            for line in lstLines:
                if line != "":
                    kmer = line.split("\t")[0]
                    rckmer = reverse_complement(kmer)
                    if kmer in setIntersectKmers:
                        setIntersectKmers.remove(kmer)
                    if rckmer in setIntersectKmers:
                        setIntersectKmers.remove(rckmer)
        pbar.update(1)
        title("Specific kmers", pbar)
    pbar.close()
    printcolor("⏩ Found "+str(len(setIntersectKmers))+" common & specific Kmers"+"\n")
    # ***** SEARCH primers ***** #
    printcolor("♊ Get primers"+"\n")
    pbar = tqdm(total=len(setIntersectKmers), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    OUT = open(pathOUTprimers, 'w')
    for kmer1 in setIntersectKmers:
        for kmer2 in setIntersectKmers:
            for pathFNA in lstFiles:
                orgName = os.path.basename(pathFNA).replace(ext, "")
                if orgName in setSelectedOrg:
                    seq = ""
                    dicoFNA = make_fasta_dict(pathFNA)
                    for key in dicoFNA:
                        seq += dicoFNA[key]
                    search = re.search(kmer1+"(.+)"+kmer2, seq)
                    if search:
                        OUT.write(kmer1+";"+kmer2+";"+str(len(search.group(1)))+"n")
        pbar.update(1)
        title("Get primers", pbar)
    pbar.close()
    OUT.close()


@fct_checker
def core_prot_tree(pathIN: str, pathOUT: str, idThr: int = 30, covThr: int = 80, ext: str = ".faa") -> Tuple[str, str, int, int, str]:
    '''
     ------------------------------------------------------------
    |               CORE PROTEINS PHYLOGENETIC TREE              |
    |------------------------------------------------------------|
    |                  Make core proteins tree                   |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathOUT : path of output files (required)               |
    |    idThr   : %identity threshold (default=30)              |
    |    covThr  : %coverage threshold (default=80)              |
    |    ext     : extension of input files (default=.faa)       |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "core_prot_tree", [ext])
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    if "." not in ext:
        ext = "."+ext
    if len(lstFiles) == 0:
        printcolor("[ERROR: core_prot_tree]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if os.path.isdir(pathOUT+"/core_align"):
        shutil.rmtree(pathOUT+"/core_align")
    os.makedirs(pathOUT+"/core_align", exist_ok=True)
    # ***** MMSEQS RBH ***** #
    mmseqs_rbh(pathIN=pathIN, pathOUT=pathOUT+"/rbh", ref="None", idThrClust=idThr, covThrClust=covThr, boolNucl=True, ext=".faa")
    # ***** RETRIEVE PROTEIN SEQUENCE ***** #
    printcolor("♊ Get proteins"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoLTtoOrg = {}
    dicoLTtoSeq = {}
    lstOrgName = []
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        lstOrgName.append(orgName)
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        dicoFAA = make_fasta_dict(pathFile)
        for key in dicoFAA:
            if "[" in key:
                lt = key.split(" [")[0]
            else:
                lt = key
            dicoLTtoOrg[lt] = orgName
            dicoLTtoSeq[lt] = dicoFAA[key]
        pbar.update(1)
        title("Get proteins", pbar)
    pbar.close()
    # ***** MAKE RBH CLUSTERS ***** #
    printcolor("♊ RBH clustering"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoCluster = {}
    for rbhFile in os.listdir(pathOUT+"/rbh"):
        pbar.set_description_str(rbhFile.replace(".rbh", "")+" ".rjust(maxpathSize-len(rbhFile.replace(".rbh", ""))))
        pathRBH = pathOUT+"/rbh/"+rbhFile
        RBH = open(pathRBH, 'r')
        lstLines = RBH.read().split("\n")[: -1]
        RBH.close()
        for line in lstLines:
            splitLine = line.split("\t")
            query = splitLine[0]
            target = splitLine[1]
            findCluster = False
            for key in dicoCluster:
                if query in dicoCluster[key] or target in dicoCluster[key]:
                    findCluster = True
                    dicoCluster[key].update([query, target])
            if findCluster is False:
                dicoCluster[len(dicoCluster)] = set([query, target])
        pbar.update(1)
        title("RBH clustering", pbar)
    pbar.close()
    # Final core dictionnary
    dicoCore = {}
    for clusterNum in dicoCluster:
        if len(dicoCluster[clusterNum]) == len(lstFiles):
            dicoCore[len(dicoCore)] = dicoCluster[clusterNum]
    printcolor("⏩ Found "+str(len(dicoCore))+" core genes"+"\n")
    # ***** ALIGN CORE PROTEINS ***** #
    printcolor("♊ Align core proteins"+"\n")
    pbar = tqdm(total=len(dicoCore), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    lstAlignFiles = []
    for coreProtNum in dicoCore:
        orgName = os.path.basename(pathFile).replace(ext, "")
        pbar.set_description_str("Core protein "+str(coreProtNum)+" ".rjust(maxpathSize-len("Core protein "+str(coreProtNum))))
        # Create protein fasta
        pathTMPFAA = pathTMP+"/coreprot_"+str(coreProtNum)+".FAA"
        TMPFAA = open(pathTMPFAA, 'w')
        for lt in dicoCore[coreProtNum]:
            TMPFAA.write(">"+lt+" ["+dicoLTtoOrg[lt]+"]\n"+dicoLTtoSeq[lt]+"\n")
        TMPFAA.close()
        # Launch muscle
        pathALIGN = pathOUT+"/core_align/coreprot_"+str(coreProtNum)+"_align.fasta"
        lstAlignFiles.append(pathALIGN)
        cmdMUSCLE = dicoGeminiPath['TOOLS']['muscle']+" -in "+pathTMPFAA+" -out "+pathALIGN+" -quiet"
        os.system(cmdMUSCLE)
        pbar.update(1)
        title("Align core", pbar)
    pbar.close()
    # Cat all genes
    dicoAlign = {}
    for alignFile in os.listdir(pathOUT+"/core_align"):
        dicoFASTA = make_fasta_dict(pathOUT+"/core_align/"+alignFile)
        for key in dicoFASTA:
            orgName = key.split("[")[1].replace("]", "")
            try:
                dicoAlign[orgName] += dicoFASTA[key]
            except KeyError:
                dicoAlign[orgName] = dicoFASTA[key]
    ALL = open(pathOUT+"/core_align/all_core_align.fasta", 'w')
    for orgName in dicoAlign:
        ALL.write(">"+orgName+"\n"+dicoAlign[orgName]+"\n")
    ALL.close()
    # ***** MAKE CORE TREE ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Make core tree", side="right")
    spinner.start()
    title("Core tree", None)
    cmdTree = dicoGeminiPath['TOOLS']['iqtree2']+" -s "+pathOUT+"/core_align/all_core_align.fasta -T "+str(cpu)+" --mem "+str(memMax)+"GB --quiet -m LG -B 1000 --seqtype AA --prefix "+pathOUT+"/iqtree2"
    os.system(cmdTree)
    shutil.copyfile(pathOUT+"/iqtree2.treefile", pathOUT+"/all_core.nwk")
    os.system("rm -f "+pathOUT+"/iqtree2*")
    spinner.stop()
    printcolor("♊ Make core tree"+"\n")


@fct_checker
def individual_core_tree(pathIN: str, pathOUT: str, idThr: int = 30, covThr: int = 80, boolNucl: bool = False, ext: str = ".faa") -> Tuple[str, str, int, int, bool, str]:
    '''
     ------------------------------------------------------------
    |      INDIVIDUAL CORE GENES/PROTEINS PHYLOGENETIC TREE      |
    |------------------------------------------------------------|
    |        Make individual and core genes/proteins tree        |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathOUT : path of output files (required)               |
    |    idThr   : %identity threshold (default=30)              |
    |    covThr  : %coverage threshold (default=80)              |
    |    ext     : extension of input files (default=.ffn)       |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "individual_core_tree", [ext])
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    if "." not in ext:
        ext = "."+ext
    if len(lstFiles) == 0:
        printcolor("[ERROR: individual_core_tree]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if os.path.isdir(pathOUT+"/core_align"):
        shutil.rmtree(pathOUT+"/core_align")
    os.makedirs(pathOUT+"/core_align", exist_ok=True)
    if os.path.isdir(pathOUT+"/core_tree"):
        shutil.rmtree(pathOUT+"/core_tree")
    os.makedirs(pathOUT+"/core_tree", exist_ok=True)
    # ***** MMSEQS RBH ***** #
    mmseqs_rbh(pathIN=pathIN, pathOUT=pathOUT+"/rbh", ref="None", idThrClust=idThr, covThrClust=covThr, boolNucl=True, ext=ext)
    # ***** RETRIEVE GENE SEQUENCE ***** #
    if boolNucl is True:
        printcolor("♊ Get genes"+"\n")
    else:
        printcolor("♊ Get proteins"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoLTtoOrg = {}
    dicoLTtoSeq = {}
    lstOrgName = []
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        lstOrgName.append(orgName)
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        dicoFFN = make_fasta_dict(pathFile)
        for key in dicoFFN:
            if "[" in key:
                lt = key.split(" [")[0]
            else:
                lt = key
            dicoLTtoOrg[lt] = orgName
            dicoLTtoSeq[lt] = dicoFFN[key]
        pbar.update(1)
        title("Get genes", pbar)
    pbar.close()
    # ***** MAKE RBH CLUSTERS ***** #
    printcolor("♊ RBH clustering"+"\n")
    dicoCluster = make_rbhcluster_dict(pathIN=pathOUT+"/rbh", pathIN2=pathIN, pathJSON="None", idThrClust=idThr, covThrClust=covThr, ext=".rbh", ext2=ext)
    # Final core dictionnary
    dicoCore = {}
    for clusterNum in dicoCluster:
        setOrg = set()
        for gene in dicoCluster[clusterNum]:
            orgName = dicoLTtoOrg[gene.split(" ")[0]]
            setOrg.add(orgName)
        if len(setOrg) == len(lstFiles) == len(dicoCluster[clusterNum]):  # to avoid paralogous or split gene
            dicoCore[len(dicoCore)] = dicoCluster[clusterNum]
    if boolNucl is True:
        printcolor("⏩ Found "+str(len(dicoCore))+" core genes"+"\n")
        printcolor("♊ Align core genes"+"\n")
    else:
        printcolor("⏩ Found "+str(len(dicoCore))+" core proteins"+"\n")
        printcolor("♊ Align core proteins"+"\n")
    # ***** ALIGN CORE GENES/PROTEINS ***** #
    pbar = tqdm(total=len(dicoCore), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    lstAlignFiles = []
    for coreNum in dicoCore:
        if boolNucl is True:
            pbar.set_description_str("Core gene "+str(coreNum)+" ".rjust(maxpathSize-len("Core gene "+str(coreNum))))
            pathTMPFASTA = pathTMP+"/coregene_"+str(coreNum)+".fasta"
            pathALIGN = pathOUT+"/core_align/coregene_"+str(coreNum)+"_align.fasta"        
        else:            
            pbar.set_description_str("Core protein "+str(coreNum)+" ".rjust(maxpathSize-len("Core protein "+str(coreNum))))
            pathTMPFASTA = pathTMP+"/coreprot_"+str(coreNum)+".fasta"
            pathALIGN = pathOUT+"/core_align/coreprot_"+str(coreNum)+"_align.fasta"        
        # Create gene fasta
        TMPFASTA = open(pathTMPFASTA, 'w')
        for header in dicoCore[coreNum]:
            lt = header.split(" ")[0]
            TMPFASTA.write(">"+dicoLTtoOrg[lt]+"\n"+dicoLTtoSeq[lt]+"\n")
        TMPFASTA.close()
        # Launch muscle    
        lstAlignFiles.append(pathALIGN)
        cmdMUSCLE = dicoGeminiPath['TOOLS']['muscle']+" -in "+pathTMPFASTA+" -out "+pathALIGN+" -quiet"
        os.system(cmdMUSCLE)
        pbar.update(1)
        title("Align core", pbar)
    pbar.close()
    # Cat all genes
    dicoAlign = {}
    for alignFile in os.listdir(pathOUT+"/core_align"):
        dicoFASTA = make_fasta_dict(pathOUT+"/core_align/"+alignFile)
        for orgName in dicoFASTA:
            try:
                dicoAlign[orgName] += dicoFASTA[orgName]
            except KeyError:
                dicoAlign[orgName] = dicoFASTA[orgName]
    ALL = open(pathOUT+"/core_align/all_core_align.fasta", 'w')
    for orgName in dicoAlign:
        ALL.write(">"+orgName+"\n"+dicoAlign[orgName]+"\n")
    ALL.close()
    # ***** MAKE TREES ***** #
    # indiviual trees
    printcolor("♊ Make trees"+"\n")
    pbar = tqdm(total=len(dicoCore), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    lstAlignFiles = []
    for coreNum in dicoCore:
        # Launch
        if boolNucl is True:
            pbar.set_description_str("Core gene "+str(coreNum)+" ".rjust(maxpathSize-len("Core gene "+str(coreNum))))
            pathALIGN = pathOUT+"/core_align/coregene_"+str(coreNum)+"_align.fasta"
            pathTREE = pathOUT+"/core_tree/coregene_"+str(coreNum)+"_iqtree2"
            cmdTree = dicoGeminiPath['TOOLS']['iqtree2']+" -s "+pathALIGN+" -T "+str(cpu)+" --mem "+str(memMax)+"GB --quiet -m GTR -B 1000 --seqtype DNA --prefix "+pathTREE+" --keep-ident > /dev/null 2>&1"
            os.system(cmdTree)
            os.system("rm -f "+pathTREE+".splits.nex "+pathTREE+".log "+pathTREE+".mldist "+pathTREE+".bionj "+pathTREE+".ckp.gz "+pathTREE+".contree "+pathTREE+".iqtree")
        else:
            pbar.set_description_str("Core protein "+str(coreNum)+" ".rjust(maxpathSize-len("Core protein "+str(coreNum))))            
            pathALIGN = pathOUT+"/core_align/coreprot_"+str(coreNum)+"_align.fasta"
            pathTREE = pathOUT+"/core_tree/coreprot_"+str(coreNum)+"_iqtree2"
            cmdTree = dicoGeminiPath['TOOLS']['iqtree2']+" -s "+pathALIGN+" -T "+str(cpu)+" --mem "+str(memMax)+"GB --quiet -m LG -B 1000 --seqtype AA --prefix "+pathTREE+" --keep-ident > /dev/null 2>&1"
            os.system(cmdTree)
            os.system("rm -f "+pathTREE+".splits.nex "+pathTREE+".log "+pathTREE+".mldist "+pathTREE+".bionj "+pathTREE+".ckp.gz "+pathTREE+".contree "+pathTREE+".iqtree")
        pbar.update(1)
        title("Make trees", pbar)
    pbar.close()
    # core tree
    spinner = yaspin(Spinners.aesthetic, text="♊ Make core tree", side="right")
    spinner.start()
    title("Core tree", None)
    pathALIGN = pathOUT+"/core_align/all_core_align.fasta"
    pathTREE = pathOUT+"/core_tree/all_core_align_iqtree2"
    cmdTree = dicoGeminiPath['TOOLS']['iqtree2']+" -s "+pathALIGN+" -T "+str(cpu)+" --mem "+str(memMax)+"GB --quiet -m LG -B 1000 --seqtype AA --prefix "+pathTREE+" --keep-ident > /dev/null 2>&1"
    os.system(cmdTree)
    os.system("rm -f "+pathTREE+".splits.nex "+pathTREE+".log "+pathTREE+".mldist "+pathTREE+".bionj "+pathTREE+".ckp.gz "+pathTREE+".contree "+pathTREE+".iqtree")
    spinner.stop()
    printcolor("♊ Make core tree"+"\n")


@fct_checker
def protein_similarity_matrix(pathIN: str, pathOUT: str, locusTag: str, idThr: int = 30, covThr: int = 80, ext: str = ".faa") -> Tuple[str, str, str, int, int, str]:
    '''
     ------------------------------------------------------------
    |                 Protein similarity matrix                  |
    |------------------------------------------------------------|
    |          Create similarity matrix for one protein          |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN   : path of input folder (required)              |
    |    pathOUT  : path of output xslx file (required)          |
    |    locusTag : protein locus tag (required)                 |
    |    idThr    : min % identity (default=30)                  |
    |    covThr   : min % coverage (default=80)                  |
    |    ext      : extension of input files (default=.faa)      |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "protein_similarity_matrix", [ext])
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    pathQueryFASTA = pathTMP+"/query.faa"
    pathSubjectFASTA = pathTMP+"/subject.faa"
    pathSubjectDMND = pathTMP+"/subject.dmnd"
    pathDiamondOUT = pathTMP+"/diamond.out"
    findQuery = False
    maxTargetSeq = len(lstFiles)*10
    if "." not in ext:
        ext = "."+ext
    if len(lstFiles) == 0:
        printcolor("[ERROR: protein_similarity_matrix]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if locusTag == "":
        printcolor("[ERROR: protein_similarity_matrix]\nAny input locusTag\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # ***** GET SEQUENCES *****""
    printcolor("♊ Get sequences"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    QUERY = open(pathQueryFASTA, 'w')
    SUBJECT = open(pathSubjectFASTA, 'w')
    dicoSeq = {}
    setOrg = set()
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        setOrg.add(orgName)
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        dicoFAA = make_fasta_dict(pathFile)
        for key in dicoFAA:
            currentLT = key.split("|")[0]
            dicoSeq[currentLT] = dicoFAA[key]
            SUBJECT.write(">"+currentLT+"_____"+orgName+"\n"+dicoFAA[key]+"\n")
            if locusTag.lower() == currentLT.lower():
                QUERY.write(">"+locusTag+"_____"+orgName+"\n"+dicoFAA[key]+"\n")
                findQuery = True
        pbar.update(1)
        title("Get seq", pbar)
    pbar.close()
    QUERY.close()
    SUBJECT.close()
    if findQuery is False:
        printcolor("[ERROR: protein_similarity_matrix]\nLocusTag \""+locusTag+"\" not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # ***** MAKE DIAMOND DATABASE ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Make database", side="right")
    spinner.start()
    title("makeblastdb", None)
    cmdMakeDB = dicoGeminiPath['TOOLS']['diamond']+" makedb --in "+pathSubjectFASTA+" -d "+pathSubjectDMND+" --quiet --threads "+str(cpu)
    os.system(cmdMakeDB)
    spinner.stop()
    printcolor("♊ Make database"+"\n")
    # ***** LAUNCH BLASTP ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Launch blastp", side="right")
    spinner.start()
    title("blastP", None)
    cmdBlastP = dicoGeminiPath['TOOLS']['diamond']+" blastp -d "+pathSubjectDMND+" -q "+pathQueryFASTA+" -o "+pathDiamondOUT + \
        " --id "+str(idThr)+" --query-cover "+str(covThr)+" --subject-cover "+str(covThr)+" --max-target-seqs "+str(maxTargetSeq)+" --quiet --threads "+str(cpu)
    os.system(cmdBlastP)
    spinner.stop()
    printcolor("♊ Launch blastp"+"\n")
    # ***** GET ORTHOLOGOUS ***** #
    printcolor("♊ Get orthologous"+"\n")
    dicoBlastP = {}
    IN = open(pathDiamondOUT, 'r')
    lstLines = IN.read().split("\n")[: -1]
    IN.close()
    for line in lstLines:
        splitLine = line.split("\t")
        qseqid = splitLine[0]
        sseqid = splitLine[1]
        ltHit = sseqid.split("_____")[0]
        orgHit = sseqid.split("_____")[1]
        pident = float(splitLine[2])
        evalue = float(splitLine[10])
        if orgHit not in dicoBlastP or dicoBlastP[orgHit]['evalue'] > evalue:
            dicoBlastP[orgHit] = {'lt': ltHit, 'pident': pident, 'evalue': evalue}
    # ***** COMPUTE pairwise orthologous similarity ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Orthologous similarity", side="right")
    spinner.start()
    title("Orthosim", None)
    # Write FASTA
    FASTA = open(pathQueryFASTA, 'w')
    for orgHit in dicoBlastP:
        lt = dicoBlastP[orgHit]['lt']
        FASTA.write(">"+lt+"_____"+orgHit+"\n"+dicoSeq[lt]+"\n")
    FASTA.close()
    shutil.copyfile(pathQueryFASTA, pathSubjectFASTA)
    # Make database
    cmdMakeDB = dicoGeminiPath['TOOLS']['diamond']+" makedb --in "+pathSubjectFASTA+" -d "+pathSubjectDMND+" --quiet --threads "+str(cpu)
    os.system(cmdMakeDB)
    # Launch Blastp
    cmdBlastP = dicoGeminiPath['TOOLS']['diamond']+" blastp -d "+pathSubjectDMND+" -q "+pathQueryFASTA+" -o "+pathDiamondOUT + \
        " --id "+str(idThr)+" --query-cover "+str(covThr)+" --subject-cover "+str(covThr)+" --max-target-seqs "+str(maxTargetSeq)+" --quiet --threads "+str(cpu)
    os.system(cmdBlastP)
    # Parse output
    dicoPairwise = {}
    IN = open(pathDiamondOUT, 'r')
    lstLines = IN.read().split("\n")[: -1]
    IN.close()
    for line in lstLines:
        splitLine = line.split("\t")
        qseqid = splitLine[0]
        # ltSubject = qseqid.split("_____")[0]
        orgSubject = qseqid.split("_____")[1]
        sseqid = splitLine[1]
        ltHit = sseqid.split("_____")[0]
        orgHit = sseqid.split("_____")[1]
        pident = float(splitLine[2])
        evalue = float(splitLine[10])
        if orgSubject not in dicoPairwise:
            dicoPairwise[orgSubject] = {}
        if orgHit not in dicoPairwise[orgSubject] or dicoPairwise[orgSubject][orgHit]['evalue'] > evalue:
            dicoPairwise[orgSubject][orgHit] = {'lt': ltHit, 'pident': pident, 'evalue': evalue}
    spinner.stop()
    printcolor("♊ Orthologous similarity"+"\n")
    # ***** CONSTRUCT MATRIX ***** #
    printcolor("♊ Create clustering matrix"+"\n")
    # Create matrix dictionnaries
    dicoMatrixID = {}
    dicoMatrixLT = {}
    minId = 100
    for orgName1 in setOrg:
        dicoMatrixID[orgName1] = {}
        dicoMatrixLT[orgName1] = {}
        for orgName2 in setOrg:
            if orgName1 not in dicoPairwise or orgName2 not in dicoPairwise[orgName1]:
                dicoMatrixID[orgName1][orgName2] = 0.0
                dicoMatrixLT[orgName1][orgName2] = "None"
            else:
                dicoMatrixID[orgName1][orgName2] = dicoPairwise[orgName1][orgName2]['pident']
                dicoMatrixLT[orgName1][orgName2] = dicoPairwise[orgName1][orgName2]['lt']
                minId = min(minId, dicoPairwise[orgName1][orgName2]['pident'])
    # Clustering
    df = pd.DataFrame(dicoMatrixID)
    cg = sns.clustermap(df, cmap='crest', figsize=(50, 50), tree_kws={'linewidths': 2.5}, dendrogram_ratio=0.01, annot_kws={"size": 35 / np.sqrt(len(df))}, cbar_kws={'label': 'similarity %'})
    # Retrieve ordered ticks label
    newColums = df.columns[cg.dendrogram_col.reordered_ind]
    newIndexs = df.index[cg.dendrogram_row.reordered_ind]
    newData = df.loc[newIndexs, newColums]
    orderedOrg = list(newData.keys())
    # Write output excel
    workbook = xlsxwriter.Workbook(pathOUT)
    worksheet = workbook.add_worksheet()
    headerFormat1 = workbook.add_format({'align': 'center', 'valign': 'bottom', 'border': 1, 'font_size': 11})
    headerFormat2 = workbook.add_format({'align': 'center', 'valign': 'bottom', 'border': 1, 'font_size': 11})
    headerFormat1.set_rotation(90)
    defaultRowFormat = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 8, 'text_wrap': True, 'bg_color': "#FFFFFF"})
    HEX_list = linear_gradient("#FFFFFF", "#D40000", 101-int(minId))
    lstRowFormat = []
    for i in range(101):
        lstRowFormat.append(defaultRowFormat)
    for i in range(len(HEX_list)):
        value = int(minId)+i
        rowFormat = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 8, 'text_wrap': True, 'bg_color': HEX_list[i].upper()})
        lstRowFormat[value] = rowFormat
    dicoWidth = {0: 0}
    dicoHeight = {0: 0}
    worksheet.write(0, 0, "", headerFormat2)
    col = 1
    for orgName1 in orderedOrg:
        worksheet.write(0, col, orgName1, headerFormat1)
        dicoHeight[0] = max(dicoHeight[0], len(orgName1))
        col += 1
    row = 1
    for orgName1 in orderedOrg:
        worksheet.write(row, 0, orgName1, headerFormat2)
        dicoWidth[0] = max(dicoWidth[0], len(orgName1))
        col = 1
        for orgName2 in orderedOrg:
            value = dicoMatrixID[orgName1][orgName2]
            worksheet.write(row, col, value, lstRowFormat[int(value)])
            try:
                dicoWidth[col] = max(dicoWidth[col], len(str(value)))
            except KeyError:
                dicoWidth[col] = len(str(value))
            col += 1
        row += 1
    # Adjust row height and column width
    for row in dicoHeight:
        worksheet.set_row(row, dicoHeight[row]*2.5)
    for col in dicoWidth:
        worksheet.set_column(col, col, dicoWidth[col])
    workbook.close()


@fct_checker
def panacota_flexible_tree(pathIN: str, pathOUT: str, filterOrg: str = "None") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                PanACoTA FLEXIBLE GENES TREE                |
    |------------------------------------------------------------|
    |       Make flexible genes tree form PanACoTA results       |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN    : path of input PanACoTA folder (required)    |
    |    pathOUT   : path of output folder (required)            |
    |    filterOrg : organism name filter (default=None)         |
     ------------------------------------------------------------
    '''
    pathPanacotaDir = path_converter(pathIN)
    if not os.path.isdir(pathPanacotaDir):
        printcolor("[ERROR: panacota_flexible_tree]\nPanACoTA folder not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathPanacotaGenesDir = pathPanacotaDir+"/Genes"
    if not os.path.isdir(pathPanacotaGenesDir):
        printcolor("[ERROR: panacota_flexible_tree]\nPanACoTA \"Genes\" folder not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathPanacotaPangenome = pathPanacotaDir+"/PanGenome.lst"
    if not os.path.isfile(pathPanacotaPangenome):
        printcolor("[ERROR: panacota_flexible_tree]\nPanACoTA \"PanGenome.lst\" file not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    pathJSONmean = pathOUT+"/flexible_mean_similarities.json"
    pathJSONgene = pathOUT+"/flexible_genes_similarities.json"
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    dicoGeneToOrg = {}
    dicoGenes = {}
    setOrg = set()
    # ***** FILTERING by organism name ***** #
    lstGenesFiles = []
    maxOrgSize = 0
    for genesFile in os.listdir(pathPanacotaGenesDir):
        orgName = genesFile.replace(".gen", "")
        if filterOrg == "None" or filterOrg in orgName:
            setOrg.add(orgName)
            lstGenesFiles.append(genesFile)
            maxOrgSize = max(maxOrgSize, len(orgName))
    # ***** BROWSE all genes FASTA ***** #
    printcolor("♊ Get genes"+"\n")
    pbar = tqdm(total=len(lstGenesFiles), ncols=50+maxOrgSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for genesFile in lstGenesFiles:
        orgName = genesFile.replace(".gen", "")
        pbar.set_description_str(orgName+" ".rjust(maxOrgSize-len(orgName)))
        dicoFASTA = make_fasta_dict(pathPanacotaGenesDir+"/"+genesFile)
        for key in dicoFASTA:
            lt = key.split(" ")[0]
            dicoGeneToOrg[lt] = orgName
            dicoGenes[lt] = dicoFASTA[key]
        pbar.update(1)
        title("Get genes", pbar)
    pbar.close()
    # ***** PARSE pangenome file ***** #
    printcolor("♊ Read pangenome file"+"\n")
    PAN = open(pathPanacotaPangenome, 'r')
    lstLines = PAN.read().split("\n")[: -1]
    PAN.close()
    dicoPan = {}
    for line in lstLines:
        splitLine = line.split(" ")
        numPan = int(splitLine[0])
        dicoPan[numPan] = {}
        for i in range(1, len(splitLine), 1):
            lt = splitLine[i]
            if lt in dicoGeneToOrg:
                orgName = dicoGeneToOrg[lt]
                try:
                    dicoPan[numPan][orgName].append(lt)
                except KeyError:
                    dicoPan[numPan][orgName] = [lt]
    # ***** WRITE FLEXIBLE genes file ***** #
    printcolor("♊ Get flexibles"+"\n")
    pbar = tqdm(total=len(dicoPan), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    lstFASTA = []
    lstPanNumCore = []
    for numPan in dicoPan:
        formatNumPan = str(numPan).zfill(len(str(len(dicoPan))))
        pbar.set_description_str(formatNumPan)
        if len(dicoPan[numPan]) > 1 and len(dicoPan[numPan]) != len(setOrg):
            pathFlexibleFASTA = pathTMP+"/pan_gene_"+str(numPan)+".fasta"
            lstFASTA.append(pathFlexibleFASTA)
            lstPanNumCore.append(numPan)
            FASTA = open(pathFlexibleFASTA, 'w')
            for orgName in dicoPan[numPan]:
                for lt in dicoPan[numPan][orgName]:
                    FASTA.write(">"+lt+"\n"+dicoGenes[lt]+"\n")
            FASTA.close()
        pbar.update(1)
        title("Get flexibles", pbar)
    pbar.close()
    # ***** ALIGN FLEXIBLE genes ***** #
    printcolor("♊ Align flexibles"+"\n")
    pbar = tqdm(total=len(lstFASTA), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    lstALIGN = []
    for pathFASTA in lstFASTA:
        pbar.set_description_str(os.path.basename(pathFASTA).replace(".fasta", ""))
        pathALIGN = pathFASTA.replace(".fasta", "_align.fasta")
        cmdFAMSA = dicoGeminiPath['TOOLS']['famsa']+" -t "+str(cpu)+" "+pathFASTA+" "+pathALIGN+" > /dev/null 2>&1"
        os.system(cmdFAMSA)
        lstALIGN.append(pathALIGN)
        pbar.update(1)
        title("Align flexibles", pbar)
    pbar.close()
    # ***** Init dict ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Init dict", side="right")
    spinner.start()
    title("Init", None)
    dicoMeanFlexibleIdentity = {}
    dicoFlexibleIdentityGene = {}
    for numPan in lstPanNumCore:
        dicoFlexibleIdentityGene[numPan] = {}
    for orgName1 in setOrg:
        dicoMeanFlexibleIdentity[orgName1] = {}
        for numPan in lstPanNumCore:
            dicoFlexibleIdentityGene[numPan][orgName1] = {}
        for orgName2 in setOrg:
            if orgName1 != orgName2:
                dicoMeanFlexibleIdentity[orgName1][orgName2] = []
    spinner.stop()
    printcolor("♊ Init dict"+"\n")
    # ***** COMPUTE flexible genes similarity ***** #
    printcolor("♊ Flexibles similarity"+"\n")
    pbar = tqdm(total=len(lstALIGN), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathALIGN in lstALIGN:
        numPan = int(os.path.basename(pathALIGN).replace("_align.fasta", "").replace("pan_gene_", ""))
        formatNumPan = str(numPan).zfill(len(str(len(dicoPan))))
        pbar.set_description_str(formatNumPan)
        dicoALIGN = make_fasta_dict(pathALIGN)
        for lt1 in dicoALIGN:
            orgName1 = dicoGeneToOrg[lt1]
            seq1 = dicoALIGN[lt1]
            for lt2 in dicoALIGN:
                orgName2 = dicoGeneToOrg[lt2]
                seq2 = dicoALIGN[lt2]
                if orgName1 != orgName2:
                    nbMatch = 0
                    for i in range(len(seq1)):
                        if seq1[i] == seq2[i]:
                            nbMatch += 1
                    pident = (nbMatch*100)/len(seq1)
                    if orgName2 not in dicoFlexibleIdentityGene[numPan][orgName1] or pident >= dicoFlexibleIdentityGene[numPan][orgName1][orgName2]:
                        dicoFlexibleIdentityGene[numPan][orgName1][orgName2] = pident
        # Add to mean dictionnary
        for orgName1 in dicoFlexibleIdentityGene[numPan]:
            for orgName2 in dicoFlexibleIdentityGene[numPan][orgName1]:
                dicoMeanFlexibleIdentity[orgName1][orgName2].append(dicoFlexibleIdentityGene[numPan][orgName1][orgName2])
        pbar.update(1)
        title("Flexible sim", pbar)
    pbar.close()
    # ***** COMPUTE MEAN flexible genes similarity ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Compute mean flexible similarity", side="right")
    spinner.start()
    title("Flexible meansim", None)
    for orgName1 in dicoMeanFlexibleIdentity:
        for orgName2 in dicoMeanFlexibleIdentity[orgName1]:
            if len(dicoMeanFlexibleIdentity[orgName1][orgName2]) > 0:
                dicoMeanFlexibleIdentity[orgName1][orgName2] = np.mean(dicoMeanFlexibleIdentity[orgName1][orgName2])
            else:
                dicoMeanFlexibleIdentity[orgName1][orgName2] = 0.0
    spinner.stop()
    printcolor("♊ Compute mean flexible similarity"+"\n")
    # ***** DUMP JSON ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Dump results JSON", side="right")
    spinner.start()
    title("Dump", None)
    dump_json(dicoMeanFlexibleIdentity, pathJSONmean)
    dump_json(dicoFlexibleIdentityGene, pathJSONgene)
    spinner.stop()
    printcolor("♊ Dump results JSON"+"\n")


@fct_checker
def snippy(pathIN: str, pathOUT: str, ext: str = ".gbk") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                   Snippy - SNPs calling                    |
    |------------------------------------------------------------|
    |              Call pairwise SNPs using Snippy               |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN   : path of input PanACoTA folder (required)     |
    |    pathOUT  : path of output folder (required)             |
    |    ext      : extension of input files (default=.faa)      |
     ------------------------------------------------------------
    # dicoSNIPPY[org1][org2] = {'count':{'total':0,'intergenic':0,'cds':0,'cds_notsyn':0,'gene':0,'gene_notsyn':0},
    #                           'dicoSNP': {lt:{pos:{'contig','hgvs_n','hgvs_p','effect','gene','product'}}}
    #                         }
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "snippy", [ext, ext+".gz"])
    maxpathSize = (maxpathSize*2)+4
    pathOUT = path_converter(pathOUT)
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    os.makedirs(pathOUT, exist_ok=True)
    pathJSON = pathOUT+"/snippy.json"
    pathLOG = pathOUT+"/log.out"
    if "." not in ext:
        ext = "."+ext
    if len(lstFiles) == 0:
        printcolor("[ERROR: snippy]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if len(lstFiles) == 1:
        printcolor("[ERROR: snippy]\nSnippy required a minimum of 2 input files\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # ***** CONVERT gbk to contigs ***** #
    printcolor("♊ GBK to FNA"+"\n")
    lstFNAfiles = []
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathGBK in lstFiles:
        orgName = os.path.basename(pathGBK).replace(ext, "").replace(".gz", "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        pathFNA = pathTMP+"/"+orgName+".fna"
        lstFNAfiles.append(pathFNA)
        gbk_to_fna(pathIN=pathGBK, pathOUT=pathFNA)
        pbar.update(1)
        title("gbk2fna", pbar)
    pbar.close()
    # ***** LAUNCH Snippy ***** #
    printcolor("♊ Launch"+"\n")
    lstSnippyOUT = []
    pbar = tqdm(total=len(lstFiles)*(len(lstFiles)-1), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathGBK in lstFiles:
        # decompress if required
        if pathGBK[-3:] == ".gz":
            os.system("gzip -c -d "+pathGBK+" > "+pathTMP+"/"+os.path.basename(pathGBK).replace(".gz", ""))
            pathGBK = pathTMP+"/"+os.path.basename(pathGBK).replace(".gz", "")
        orgName1 = os.path.basename(pathGBK).replace(ext, "")
        # Launch against each contigs FNA
        for pathFNA in lstFNAfiles:
            orgName2 = os.path.basename(pathFNA).replace(".fna", "")
            if orgName1 != orgName2:
                pbar.set_description_str(orgName1+" VS "+orgName2.rjust(maxpathSize-len(orgName1+" VS "+orgName2)))
                pathTMPOUT = pathTMP+"/snippy_"+orgName1+"_vs_"+orgName2
                pathSNPOUT = pathOUT+"/snippy_"+orgName1+"_vs_"+orgName2+".tab"
                lstSnippyOUT.append(pathSNPOUT)
                if not os.path.isfile(pathSNPOUT):
                    cmdSnippy = dicoGeminiPath['TOOLS']['snippy']+" --outdir "+pathTMPOUT+" --ref "+pathGBK+" --ctgs "+pathFNA + \
                        " --cpus "+str(cpu)+" --ram "+str(memMax)+" --tmpdir "+pathTMP+" >> "+pathLOG+" 2>&1"
                    os.system(cmdSnippy)
                    shutil.copyfile(pathTMPOUT+"/snps.tab", pathSNPOUT)
                    shutil.rmtree(pathTMPOUT)
                pbar.update(1)
                title("Launch", pbar)
    pbar.close()
    # ***** PARSE Snippy ***** #
    printcolor("♊ Parse"+"\n")
    dicoSNIPPY = {}
    pbar = tqdm(total=len(lstSnippyOUT), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathOUT in lstSnippyOUT:
        pbar.set_description_str(os.path.basename(pathOUT).replace(".tab", "").rjust(maxpathSize-len(os.path.basename(pathOUT).replace(".tab", ""))))
        orgName1 = os.path.basename(pathOUT).split("_vs_")[0].replace("snippy_", "")
        orgName2 = os.path.basename(pathOUT).split("_vs_")[1].replace(".tab", "").replace("snippy_", "")
        if orgName1 not in dicoSNIPPY:
            dicoSNIPPY[orgName1] = {}
        if orgName2 not in dicoSNIPPY[orgName1]:
            dicoSNIPPY[orgName1][orgName2] = {'count': {'total': 0, 'intergenic': 0, 'cds': 0, 'cds_notsyn': 0, 'gene': 0, 'gene_notsyn': 0}, 'dicoSNP': {}}
        OUT = open(pathOUT, 'r')
        lstLines = OUT.read().split("\n")[1: -1]
        OUT.close()
        for line in lstLines:
            splitLine = line.split("\t")
            contig = splitLine[0]
            pos = splitLine[1]
            # typeSNP = splitLine[2]
            ref = splitLine[3]
            alt = splitLine[4]
            # evidence = splitLine[5]
            ftype = splitLine[6]
            # strand = splitLine[7]
            # nt_pos = splitLine[8]
            # aa_pos = splitLine[9]
            effect = splitLine[10]
            lt = splitLine[11]
            gene = splitLine[12]
            product = splitLine[13]
            # HGVS
            hgvs_n = str(pos)+ref+">"+alt
            hgvs_p = effect.split(" ")[-1]
            if "p." not in hgvs_p:
                hgvs_p = ""
            # Normalize effect
            effect = effect.split(" ")[0]
            # Count entries
            dicoSNIPPY[orgName1][orgName2]['count']['total'] += 1
            if ftype == "CDS":
                dicoSNIPPY[orgName1][orgName2]['count']['cds'] += 1
            else:
                dicoSNIPPY[orgName1][orgName2]['count']['intergenic'] += 1
                lt = "intergenic"
            # Fill dicoSNP
            if lt not in dicoSNIPPY[orgName1][orgName2]['dicoSNP']:
                dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt] = {}
            dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt][pos] = {'contig': contig, 'hgvs_n': hgvs_n, 'effect': effect, 'gene': gene, 'product': product}
            if "synonymous_variant" not in effect:
                dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt][pos]['hgvs_p'] = hgvs_p
                dicoSNIPPY[orgName1][orgName2]['count']['cds_notsyn'] += 1
        # Count number of SNPs per gene
        for lt in dicoSNIPPY[orgName1][orgName2]['dicoSNP']:
            if lt != "intergenic":
                dicoSNIPPY[orgName1][orgName2]['count']['gene'] += 1
                notSyn = False
                for pos in dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt]:
                    if "synonymous_variant" not in dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt][pos]['effect']:
                        notSyn = True
                        break
                if notSyn is True:
                    dicoSNIPPY[orgName1][orgName2]['count']['gene_notsyn'] += 1
            else:
                for pos in dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt]:
                    del(dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt][pos]['effect'])
                    del(dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt][pos]['gene'])
                    del(dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt][pos]['product'])
                    del(dicoSNIPPY[orgName1][orgName2]['dicoSNP'][lt][pos]['hgvs_p'])
        pbar.update(1)
        title("Parse", pbar)
    pbar.close()
    # ***** SAVE SNPs JSON file ***** #
    printcolor("♊ Dump SNPs JSON"+"\n")
    dump_json(dicoSNIPPY, pathJSON)


@fct_checker
def wgrr_matrix(pathIN: str, pathOUT: str, ext: str = ".faa") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                    wGRR DISTANCE MATRIX                    |
    |------------------------------------------------------------|
    |                    wGRR distance matrix                    |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN    : path of input faa files or folder (required)|
    |    pathOUT   : path of output sketch folder (required)     |
    |    ext       : extension of input files (default=.faa)     |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "wgrr_matrix", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: wgrr_matrix]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    lstFiles.sort()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    if "." not in ext:
        ext = "."+ext
    # PATHS
    PATHJSON = pathOUT+"/wGRR.json"
    PATHPNG = pathOUT+"/wGRR.png"
    PATHSVG = pathOUT+"/wGRR.svg"
    PATHXLSX = pathOUT+"/wGRR.xlsx"
    # ***** REFORMAT FASTA ***** #
    printcolor("♊ Reformat FASTA"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        # Reformat FAA
        pathTMPFAA = geminiset.pathTMP+"/"+orgName+".faa"
        OUT = open(pathTMPFAA, 'w')
        dicoFAA = make_fasta_dict(pathFile)
        for key in dicoFAA:
            OUT.write(">"+key+"_____"+orgName+"\n"+dicoFAA[key]+"\n")
        OUT.close()
        pbar.update(1)
        title("Reformat FASTA", pbar)
    pbar.close()
    # ***** COMPUTE RBH ***** #
    mmseqs_rbh(pathIN=geminiset.pathTMP, pathOUT=pathOUT, ref="None", idThrClust=20, covThrClust=50, boolNucl=False, ext=".faa")
    # ***** GET Number of proteins per file ***** #
    printcolor("♊ Count proteins"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoNbProt = {}
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        dicoFAA = make_fasta_dict(geminiset.pathTMP+"/"+orgName+".faa")
        dicoNbProt[os.path.basename(pathFile).replace(".faa", "")] = len(dicoFAA)
        pbar.update(1)
        title("Count proteins", pbar)
    pbar.close()
    # ***** GET percentage identity between all best hits (evalue<10-5) ***** #
    printcolor("♊ Get RBH pident"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoRBH = {}
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        pathRBH = pathOUT + "/" + orgName + ".rbh"
        RBH = open(pathRBH, 'r')
        lstLines = RBH.read().split("\n")[:-1]
        RBH.close()
        dicoRBH[orgName] = {}
        for line in lstLines:
            splitLine = line.split("\t")
            target = splitLine[1]
            orgTarget = target.split("_____")[1]
            pident = float(splitLine[2])
            evalue = float(splitLine[3])
            if evalue < 0.00001:
                try:
                    dicoRBH[orgName][orgTarget] += pident
                except KeyError:
                    dicoRBH[orgName][orgTarget] = pident
        pbar.update(1)
        title("Get RBH pident", pbar)
    pbar.close()
    # ***** COMPUTE wGRR *****#
    printcolor("♊ Compute wGRR"+"\n")
    pbar = tqdm(total=len(dicoRBH), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoWGRR = {}
    for orgName1 in dicoRBH:
        pbar.set_description_str(orgName1+" ".rjust(maxpathSize-len(orgName1)))
        dicoWGRR[orgName1] = {}
        for orgName2 in dicoRBH:
            if orgName1 == orgName2:
                dicoWGRR[orgName1][orgName2] = 100.0
            else:
                dicoWGRR[orgName1][orgName2] = dicoRBH[orgName1][orgName2] / min(dicoNbProt[orgName1], dicoNbProt[orgName2])
        pbar.update(1)
        title("Compute wGRR", pbar)
    pbar.close()
    # Save matrix
    dump_json(dicoWGRR, PATHJSON)
    # ***** CREATE wGRR matrix *****#
    printcolor("♊ wGRR matrix"+"\n")
    # Make matrix
    df = pd.DataFrame(dicoWGRR)
    figsize = len(df) / np.sqrt(len(df))
    cg = sns.clustermap(df, cmap='Spectral', figsize=(figsize, figsize), tree_kws={'linewidths': 1.5}, dendrogram_ratio=0.05, xticklabels=True, yticklabels=True, linewidths=0.0, rasterized=True)
    cg.ax_cbar.remove()
    plt.tick_params(axis='both', which='major', labelsize=1.5, width=0.2)
    fig = plt.gcf()  # or by other means, like plt.subplots
    figsize = fig.get_size_inches()
    fig.set_size_inches(figsize * 1.5)  # scale current size by 1.5
    plt.savefig(PATHPNG, dpi=300)
    plt.savefig(PATHSVG)
    # Retrieve ordered ticks label
    newColums = df.columns[cg.dendrogram_col.reordered_ind]
    newIndexs = df.index[cg.dendrogram_row.reordered_ind]
    newData = df.loc[newIndexs, newColums]
    orderedOrg = list(newData.keys())
    # Write output excel
    workbook = xlsxwriter.Workbook(PATHXLSX)
    worksheet = workbook.add_worksheet()
    headerFormat1 = workbook.add_format({'align': 'center', 'valign': 'bottom', 'border': 1, 'font_size': 11})
    headerFormat2 = workbook.add_format({'align': 'center', 'valign': 'bottom', 'border': 1, 'font_size': 11})
    headerFormat1.set_rotation(90)
    defaultRowFormat = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 8, 'text_wrap': True, 'bg_color': "#FFFFFF"})
    HEX_list, RBG_list = linear_gradient("#FFFFFF", "#D40000", 101)
    lstRowFormat = []
    for i in range(101):
        lstRowFormat.append(defaultRowFormat)
    for i in range(len(HEX_list)):
        rowFormat = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 8, 'text_wrap': True, 'bg_color': HEX_list[i].upper()})
        lstRowFormat[i] = rowFormat
    dicoWidth = {0: 0}
    dicoHeight = {0: 0}
    worksheet.write(0, 0, "", headerFormat2)
    col = 1
    for orgName1 in orderedOrg:
        worksheet.write(0, col, orgName1, headerFormat1)
        dicoHeight[0] = max(dicoHeight[0], len(orgName1))
        col += 1
    row = 1
    for orgName1 in orderedOrg:
        worksheet.write(row, 0, orgName1, headerFormat2)
        dicoWidth[0] = max(dicoWidth[0], len(orgName1))
        col = 1
        for orgName2 in orderedOrg:
            value = dicoWGRR[orgName1][orgName2]
            worksheet.write(row, col, value, lstRowFormat[int(value)])
            try:
                dicoWidth[col] = max(dicoWidth[col], len(str(value)))
            except KeyError:
                dicoWidth[col] = len(str(value))
            col += 1
        row += 1
    # Adjust row height and column width
    for row in dicoHeight:
        worksheet.set_row(row, dicoHeight[row]*2.5)
    for col in dicoWidth:
        worksheet.set_column(col, col, dicoWidth[col])
    workbook.close()

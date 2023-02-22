'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ cluster
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminicluster.py                |
description      : gemini clustering functions     |
author           : dooguypapua                     |
lastmodification : 20210715                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''
import os
import geminiset
import sys
import shutil
import glob
from yaspin import yaspin
from yaspin.spinners import Spinners
from tqdm import tqdm
from typing import Tuple
from geminini import printcolor, path_converter, get_gemini_path, get_sys_info, dump_json, load_json, exit_gemini
from geminini import fct_checker, get_input_files, launch_threads, read_group, linear_gradient, title
from geminiparse import make_fasta_dict, gbk_to_fna, make_gbk_dict


@fct_checker
def mmseqs_easycluster(pathIN: str, pathOUT: str, idThr: int = 30, maxLRthr: int = 80, ext: str = ".faa") -> Tuple[str, str, int, int, str]:
    '''
     ------------------------------------------------------------
    |                     MMSEQS easycluster                     |
    |------------------------------------------------------------|
    |          MMSEQS easycluster: sensitive clustering          |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathOUT : path of output files (required)               |
    |    idThr   : % identity threshold (default=30)             |
    |    maxLRthr: % maxLrap threshold (default=80)              |
    |    ext     : extension of input files (default=.faa)       |
     ------------------------------------------------------------
    |TOOLS: mmseqs                                               |
     ------------------------------------------------------------
    '''
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "mmseqs_easycluster", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: mmseqs_easycluster]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if len(lstFiles) == 1:
        printcolor("[ERROR: mmseqs_easycluster]\nMMseqs require more than one FASTA file\n", 1, "212;64;89", "None", True)
        exit_gemini()
    os.makedirs(pathOUT, exist_ok=True)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'mmseqs' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['mmseqs'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    printcolor("♊ easycluster"+"\n")
    # Launch command
    pathMMSEQS = pathOUT+"/mmseqs"
    if not os.path.isfile(pathMMSEQS+"_cluster.tsv") or os.path.getsize(pathMMSEQS+"_cluster.tsv") == 0:
        cmdMmseqs = dicoGeminiPath['TOOLS']['mmseqs']+" easy-cluster "+pathIN+"/*"+ext+" "+pathMMSEQS+" "+geminiset.pathTMP+" -c "+str(idThr/100)+" --cov-mode 0 --min-seq-id "+str(maxLRthr/100)+" --threads "+str(cpu)+" -v 0"
        os.system(cmdMmseqs)
        os.remove(pathMMSEQS+"_rep_seq.fasta")
        os.remove(pathMMSEQS+"_all_seqs.fasta")


@fct_checker
def mmseqs_rbh(pathIN: str, pathOUT: str, ref: str = "None", idThrClust: int = 80, covThrClust: int = 80, boolNucl: bool = False, ext: str = ".faa") -> Tuple[str, str, str, int, int, bool, str]:
    '''
     ------------------------------------------------------------
    |                         MMSEQS rbh                         |
    |------------------------------------------------------------|
    |            MMSEQS rbh: find reciprocal best hit            |
    |              (createdb + rbh + convertalis)                |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN     : path of input files or folder (required)   |
    |    pathOUT    : path of output files (required)            |
    |    ref        : ref organism (default="None")              |
    |    idThrClust : %identity clustering threshold (default=80)|
    |    covThrClust: %coverage clustering threshold (default=80)|
    |    boolNucl   : bool for nucleotide input (default=False)  |
    |    ext        : extension of input files (default=.faa)    |
     ------------------------------------------------------------
    |TOOLS: mmseqs                                               |
     ------------------------------------------------------------
    '''
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "mmseqs_easyrbh", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: mmseqs_easyrbh]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if len(lstFiles) == 1:
        printcolor("[ERROR: mmseqs_easyrbh]\nMMseqs require more than one FASTA file\n", 1, "212;64;89", "None", True)
        exit_gemini()
    os.makedirs(pathOUT, exist_ok=True)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'mmseqs' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['mmseqs'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    # ***** Check if some results are missing ***** #
    boolMissing = False
    for pathFile in lstFiles:
        orgName = os.path.basename(pathFile).replace(ext, "").replace("."+ext, "")
        pathOUTrbh = pathOUT+"/"+orgName+".rbh"
        if (ref == "None" and not os.path.isfile(pathOUTrbh)) or (ref != "None" and orgName == ref and not os.path.isfile(pathOUTrbh)):
            boolMissing = True
    if boolMissing is False:
        printcolor("♊ rbh (already done)"+"\n")
    else:
        # Paths
        pathTMPdirDB = geminiset.pathTMP+"/mmseqsdb"
        os.makedirs(pathTMPdirDB, exist_ok=True)
        pathTMPres = geminiset.pathTMP+"/mmseqsres"
        os.makedirs(pathTMPres, exist_ok=True)
        # ***** Create all required mmseqs db files ***** #
        printcolor("♊ createdb"+"\n")
        pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        for pathFile in lstFiles:
            file = os.path.basename(pathFile)
            orgName = file.replace(ext, "").replace("."+ext, "")
            pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
            if boolNucl is False:
                cmdCreateDB = dicoGeminiPath['TOOLS']['mmseqs']+" createdb "+pathFile+" "+pathTMPdirDB+"/"+orgName+".db -v 0 --dbtype 1"
            else:
                cmdCreateDB = dicoGeminiPath['TOOLS']['mmseqs']+" createdb "+pathFile+" "+pathTMPdirDB+"/"+orgName+".db -v 0 --dbtype 2"
            os.system(cmdCreateDB)
            pbar.update(1)
            title("createdb", pbar)
        pbar.close()
        # ***** Launch mmseqs rbh and convertalis ***** # [multithread]
        if ref == "None":
            printcolor("♊ rbh"+"\n")
            pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        else:
            spinner = yaspin(Spinners.aesthetic, text="♊ rbh", side="right")
            spinner.start()
            title("RBH", None)
        for pathFAA1 in lstFiles:
            orgName1 = os.path.basename(pathFAA1).replace(ext, "").replace("."+ext, "")
            if ref == "None":
                pbar.set_description_str(orgName1+" ".rjust(maxpathSize-len(orgName1)))
            if ref == "None" or (ref != "" and orgName1 == ref):
                pathOUTrbh = pathOUT+"/"+orgName1+".rbh"
                pathTMPdb1 = pathTMPdirDB+"/"+orgName1+".db"
                if not os.path.isfile(pathOUTrbh):
                    pathTMPresOrg = pathTMPres+"/"+orgName1
                    os.makedirs(pathTMPresOrg, exist_ok=True)
                    # Construct threads
                    dicoThread = {}
                    for pathFAA2 in lstFiles:
                        if pathFAA1 != pathFAA2:
                            orgName2 = os.path.basename(pathFAA2).replace(ext, "").replace("."+ext, "")
                            pathTMPdb2 = pathTMPdirDB+"/"+orgName2+".db"
                            pathTMPresPair = pathTMPresOrg+"/"+orgName2
                            os.makedirs(pathTMPresPair, exist_ok=True)
                            # Launch mmseqs (use thread = 1 because multithread)
                            if boolNucl is False:
                                cmdMmseqs = dicoGeminiPath['TOOLS']['mmseqs']+" rbh "+pathTMPdb1+" "+pathTMPdb2+" "+pathTMPresPair+"/mmseqs_ali.db "+pathTMPresPair + \
                                            " -c "+str(covThrClust/100)+" --search-type 1 --cov-mode 0 --min-seq-id "+str(idThrClust/100)+" -v 2 ; " + \
                                            dicoGeminiPath['TOOLS']['mmseqs']+" convertalis "+pathTMPdb1+" "+pathTMPdb2+" "+pathTMPresPair+"/mmseqs_ali.db "+pathTMPresPair+"/mmseqs.tsv " + \
                                            "--format-output query,target,pident,evalue,alnlen,qstart,qend,qlen,tstart,tend,tlen -v 2 --threads 1"
                            else:
                                cmdMmseqs = dicoGeminiPath['TOOLS']['mmseqs']+" rbh "+pathTMPdb1+" "+pathTMPdb2+" "+pathTMPresPair+"/mmseqs_ali.db "+pathTMPresPair + \
                                            " -c "+str(covThrClust/100)+" --search-type 3 --cov-mode 0 --min-seq-id "+str(idThrClust/100)+" -v 2 --threads 1 ; " + \
                                            dicoGeminiPath['TOOLS']['mmseqs']+" convertalis "+pathTMPdb1+" "+pathTMPdb2+" "+pathTMPresPair+"/mmseqs_ali.db "+pathTMPresPair+"/mmseqs.tsv " +\
                                            "--format-output query,target,pident,evalue,alnlen,qstart,qend,qlen,tstart,tend,tlen -v 2 --threads 1"
                            dicoThread[orgName1+"_"+orgName2] = {"cmd": cmdMmseqs, "returnstatut": None, "returnlines": []}
                    # Launch threads
                    launch_threads(dicoThread, "mmseqs_rbh", cpu, pathTMPresOrg)
                    # Merge TSV pairwise files
                    cmdMergeTSV = "cat "+pathTMPresOrg+"/*/mmseqs.tsv > "+pathOUTrbh
                    os.system(cmdMergeTSV)
                    # Delete TSV pairwise files
                    shutil.rmtree(pathTMPresOrg)
            if ref == "None":
                pbar.update(1)
                title("rbh", pbar)
        if ref == "None":
            pbar.close()
        else:
            spinner.stop()
            printcolor("♊ rbh"+"\n")


@fct_checker
def make_rbhcluster_dict(pathIN: str, pathIN2: str, pathJSON: str = "None", idThrClust: int = 80, covThrClust: int = 80, ext: str = ".rbh", ext2: str = ".faa") -> Tuple[str, str, str, int, int, str, str]:
    '''
     ------------------------------------------------------------
    |             MAKE MMSEQS RBH CLUSTER DICTIONNARY            |
    |------------------------------------------------------------|
    |      Parse MMSEQS rbh  file and create a dictionnary       |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN     : path of input rbh files/folder (required)  |
    |    pathIN2    : path of input fasta files/folder (required)|
    |    pathJSON   : path of output JSON (required)             |
    |    idThrClust : %identity threshold (default=80)           |
    |    covThrClust: %coverage clustering threshold (default=80)|
    |    ext        : extension of input files (default=.rbh)    |
    |    ext2       : extension of input files (default=.faa)    |
    |RETURN                                                      |
    |    dicoCLUSTER: { clusterNum: [lt, ..., lt], ... }         |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathIN2 = path_converter(pathIN2)
    pathJSON = path_converter(pathJSON)
    dicoLTtoOrg = {}
    dicoCLUSTER = {}
    dicoLTtoCluster = {}
    setAllheader = set()
    # rbh files
    lstFiles, maxpathSize = get_input_files(pathIN, "make_rbhcluster_dict", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_rbhcluster_dict]\nAny rbh input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    maxpathSize += 15
    # faa files
    lstFiles2, maxpathSize2 = get_input_files(pathIN2, "make_rbhcluster_dict", [ext2])
    if len(lstFiles2) == 0:
        printcolor("[ERROR: make_rbhcluster_dict]\nAny fasta input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # link rbh and faa
    dicoRBHtoFAA = {}
    for pathRBH in lstFiles:
        orgName = os.path.basename(pathRBH).replace(ext, "").replace("."+ext, "")
        for pathFAA in lstFiles2:
            orgName2 = os.path.basename(pathFAA).replace(ext2, "").replace("."+ext2, "")
            if orgName == orgName2:
                dicoRBHtoFAA[pathRBH] = pathFAA
                break
    printcolor("♊ LT to Org"+"\n")
    pbar = tqdm(total=len(lstFiles2), ncols=75+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFAA in lstFiles2:
        file = os.path.basename(pathFAA)
        orgName = file.replace(ext2, "").replace("."+ext2, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        dicoFAA = make_fasta_dict(pathFAA)
        for header in dicoFAA:
            lt = header.split(" [")[0]
            dicoLTtoOrg[lt] = orgName
            setAllheader.add(header)
        pbar.update(1)
        title("LT to org", pbar)
    pbar.close()
    printcolor("♊ Clustering"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=75+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathRBH in lstFiles:
        file = os.path.basename(pathRBH)
        orgName = file.replace(ext, "").replace("."+ext, "")
        pbar.set_description_str(orgName+" ("+str(len(dicoCLUSTER))+" clusters)"+" ".rjust(maxpathSize-len(orgName+" ("+str(len(dicoCLUSTER))+" clusters)")))
        setAllOrgLT = set()
        dicoTMP = make_fasta_dict(dicoRBHtoFAA[pathRBH])
        for key in dicoTMP:
            setAllOrgLT.add(key.split(" [")[0]+" ["+orgName+"]")
        if os.path.getsize(pathRBH) != 0:
            # Read mmseqs RBH file
            MMSEQS = open(pathRBH, 'r')
            lstLines = MMSEQS.read().split("\n")[: -1]
            MMSEQS.close()
            # Parse and assign cluster
            for line in lstLines:
                splitLine = line.split("\t")
                query = splitLine[0]+" ["+dicoLTtoOrg[splitLine[0]]+"]"
                target = splitLine[1]+" ["+dicoLTtoOrg[splitLine[1]]+"]"
                pident = float(splitLine[2])
                # case of not percent identity
                if pident <= 1:
                    pident = pident*100
                align_len = int(splitLine[4])
                qlen = int(splitLine[7])
                tlen = int(splitLine[10])
                minLrap = (align_len/min(qlen, tlen))*100
                maxLrap = (align_len/max(qlen, tlen))*100
                if pident >= idThrClust and minLrap >= covThrClust and maxLrap >= covThrClust:
                    if query in dicoLTtoCluster and target in dicoLTtoCluster:
                        continue
                    elif query in dicoLTtoCluster:
                        clusterNum = dicoLTtoCluster[query]
                        dicoCLUSTER[clusterNum].add(target)
                        dicoLTtoCluster[target] = clusterNum
                    elif target in dicoLTtoCluster:
                        clusterNum = dicoLTtoCluster[target]
                        dicoCLUSTER[clusterNum].add(query)
                        dicoLTtoCluster[query] = clusterNum
                    else:
                        clusterNum = "cluster"+str(len(dicoCLUSTER)+1).zfill(4)
                        dicoCLUSTER[clusterNum] = set([query, target])
                        dicoLTtoCluster[query] = clusterNum
                        dicoLTtoCluster[target] = clusterNum
        # ADD specific protein
        for orgLT in setAllOrgLT:
            if orgLT not in dicoLTtoCluster:
                dicoCLUSTER["cluster"+str(len(dicoCLUSTER)+1).zfill(4)] = set([orgLT])
        pbar.update(1)
        title("Clustering", pbar)
    pbar.close()
    # Convert set to list before dump
    for clusterNum in dicoCLUSTER.keys():
        dicoCLUSTER[clusterNum] = list(dicoCLUSTER[clusterNum])
    # Dump and return dictionnary
    if pathJSON != "None":
        dump_json(dicoCLUSTER, pathJSON)
    else:
        return dicoCLUSTER


@fct_checker
def make_group_core_align(pathIN: str, pathJSON: str, pathGROUP: str, pathOUT: str, boolGene: bool, boolProt: bool, extN: str = ".ffn", extP: str = ".faa") -> Tuple[str, str, str, str, str, str]:
    '''
     ------------------------------------------------------------
    |              CONSTRUCT GROUP CORE ALIGNEMENT               |
    |------------------------------------------------------------|
    |      Construct group core alignment from rbh clusters      |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN   : path of input fasta folder (required)        |
    |    pathJSON : path of input rbh JSON (required)            |
    |    pathGROUP: path of tabulated input group file (required)|
    |    pathOUT  : path of output folder (required)             |
    |    boolGene : bool for gene alignment (required)           |
    |    boolProt : bool for protein alignment (required)        |
    |    extN     : ext of gene fasta files (default=.ffn)       |
    |    extP     : ext of prot fasta files (default=.faa)       |
     ------------------------------------------------------------
    |TOOLS: muscle                                               |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    pathGROUP = path_converter(pathGROUP)
    pathOUT = path_converter(pathOUT)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'muscle' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['muscle'])
    printcolor("♊ RBH to core align"+"\n")
    if pathJSON != "" and os.path.isfile(pathJSON):
        dicoCLUSTER = load_json(pathJSON)
    else:
        printcolor("[ERROR: make_group_core_align]\nRBH cluster JSON file not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathGROUP != "" and os.path.isfile(pathGROUP):
        dicoGroup, maxGroupStrSize = read_group(pathGROUP)
    else:
        printcolor("[ERROR: make_group_core_align]\nGROUP TSV file not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # ***** Retrieve all sequences (nucleic and proteic) ***** #
    dicoGene = {}
    dicoProt = {}
    maxpathSize = 0
    for file in os.listdir(pathIN):
        pathFile = pathIN+"/"+file
        if extN in file:
            dicoGene.update(make_fasta_dict(pathFile))
            maxpathSize = max(len(os.path.basename(file).replace(extN, "")), maxpathSize)
        if extP in file:
            dicoProt.update(make_fasta_dict(pathFile))
            maxpathSize = max(len(os.path.basename(file).replace(extP, "")), maxpathSize)
    # ***** Read Cluster ***** #
    pbar = tqdm(total=len(dicoGroup)*len(dicoCLUSTER), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for group in dicoGroup:
        dicoAlignCoreGene = {}
        dicoAlignCoreProt = {}
        pathOUTgroup = pathOUT+"/"+group
        os.makedirs(pathOUTgroup, exist_ok=True)
        pathOUTmergeGene = pathOUTgroup+"/all_gene.align"
        pathOUTmergeProt = pathOUTgroup+"/all_prot.align"
        pbar.set_description_str(group+" ".rjust(maxGroupStrSize-len(group)))
        for cluster in dicoCLUSTER:
            lstLT = []
            for header in dicoCLUSTER[cluster]:
                org = header.split("[")[1].replace("]", "")
                if org in dicoGroup[group]:
                    lstLT.append(header)
            # Search core cluster
            if len(lstLT) == len(dicoGroup[group]):
                # Construct paths
                pathOUTgroupFastaGene = pathOUTgroup+"/"+cluster+"_gene.ffn"
                pathOUTgroupFastaProt = pathOUTgroup+"/"+cluster+"_prot.faa"
                pathOUTgroupAlignGene = pathOUTgroup+"/"+cluster+"_gene.align"
                pathOUTgroupAlignProt = pathOUTgroup+"/"+cluster+"_prot.align"
                if boolGene is True:
                    OUTgene = open(pathOUTgroupFastaGene, 'w')
                if boolProt is True:
                    OUTprot = open(pathOUTgroupFastaProt, 'w')
                # Write FASTA
                for header in lstLT:
                    seqGene = dicoGene[header]
                    seqProt = dicoProt[header]
                    if boolGene is True:
                        OUTgene.write(">"+header+"\n"+seqGene+"\n")
                    if boolProt is True:
                        OUTprot.write(">"+header+"\n"+seqProt+"\n")
                if boolGene is True:
                    OUTgene.close()
                if boolProt is True:
                    OUTprot.close()
                # Align sequence
                if boolGene is True:
                    os.system(dicoGeminiPath['TOOLS']['muscle']+" -in "+pathOUTgroupFastaGene+" -out "+pathOUTgroupAlignGene+" > /dev/null 2>&1")
                if boolProt is True:
                    os.system(dicoGeminiPath['TOOLS']['muscle']+" -in "+pathOUTgroupFastaProt+" -out "+pathOUTgroupAlignProt+" > /dev/null 2>&1")
                # Parse align
                if boolGene is True:
                    dicoAlign = make_fasta_dict(pathOUTgroupAlignGene)
                    for key in dicoAlign:
                        org = key.split("[")[1].replace("]", "")
                        if org not in dicoAlignCoreGene:
                            dicoAlignCoreGene[org] = ""
                        dicoAlignCoreGene[org] += dicoAlign[key]
                if boolProt is True:
                    dicoAlign = make_fasta_dict(pathOUTgroupAlignProt)
                    for key in dicoAlign:
                        org = key.split("[")[1].replace("]", "")
                        if org not in dicoAlignCoreProt:
                            dicoAlignCoreProt[org] = ""
                        dicoAlignCoreProt[org] += dicoAlign[key]
            pbar.update(1)
            title("RBH2align", pbar)
        # Write merge core align
        if boolGene is True:
            MERGE = open(pathOUTmergeGene, 'w')
            for org in dicoAlignCoreGene:
                MERGE.write(">"+org+"\n"+dicoAlignCoreGene[org]+"\n")
            MERGE.close()
        if boolProt is True:
            MERGE = open(pathOUTmergeProt, 'w')
            for org in dicoAlignCoreProt:
                MERGE.write(">"+org+"\n"+dicoAlignCoreProt[org]+"\n")
            MERGE.close()
    pbar.close()


@fct_checker
def pan_core_group(pathJSON: str, pathDIST: str, pathGROUP: str, pathOUT: str = "stdout") -> Tuple[str, str, str, str]:
    '''
     ------------------------------------------------------------
    |               PAN-CORE GENOME GROUP ANALYSIS               |
    |------------------------------------------------------------|
    |    Compute pan, core, variable genome by organism group    |
    | PAN       = NBCLUSTER TOTAL + NB SPECIFIC                  |
    | CORE      = NBCLUSTER CORE                                 |
    | VARIABLE  = NBCLUSTER NOT CORE + NB SPECIFIC               |
    | ACCESSORY = NBCLUSTER NOT CORE                             |
    | SPECIFIC  = NB SPECIFIC                                    |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathJSON : path of input rbh JSON (required)            |
    |    pathDIST : path of input JSON distance file (required)  |
    |    pathGROUP: path of tabulated input group file (required)|
    |    pathOUT  : path of output results file (default=stdout) |
     ------------------------------------------------------------
    '''
    pathJSON = path_converter(pathJSON)
    pathDIST = path_converter(pathDIST)
    pathGROUP = path_converter(pathGROUP)
    printcolor("♊ Pan-Core group"+"\n")
    if pathJSON != "" and os.path.isfile(pathJSON):
        dicoCLUSTER = load_json(pathJSON)
    else:
        printcolor("[ERROR: pan_core_group]\nRBH cluster JSON file not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathDIST != "" and os.path.isfile(pathDIST):
        dicoDIST = load_json(pathDIST)
    else:
        printcolor("[ERROR: pan_core_group]\nDIST TSV file not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathGROUP == "" or not os.path.isfile(pathGROUP):
        printcolor("[ERROR: pan_core_group]\nGROUP TSV file not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # ***** PARSE input file with group definition ***** #
    dicoGroup, maxGroupStrSize = read_group(pathGROUP)
    # ***** BROWSE group ***** #
    toDisplay = "Group\tOrganism\tDistance\tCore\tVariable\tAccessory\tSpecific\n"
    for groupName in dicoGroup:
        dicoOrgCompGenom = {}
        # Browse clusters
        for cluster in dicoCLUSTER:
            setOrg = set()
            setHeader = set()
            for header in dicoCLUSTER[cluster]:
                orgName = header.split(" [")[1].replace("]", "")
                if orgName in dicoGroup[groupName]:
                    setOrg.add(orgName)
                    setHeader.add(header)
            # Apply cluster type
            if len(setOrg) == len(dicoGroup[groupName]):
                clusterType = "core"
            else:
                if len(setOrg) == 1:
                    clusterType = "specific"
                else:
                    clusterType = "accessory"
            # Apply cluster type to each LT
            for header in setHeader:
                orgName = header.split(" [")[1].replace("]", "")
                if orgName not in dicoOrgCompGenom:
                    dicoOrgCompGenom[orgName] = {"core": set(), "variable": set(), "accessory": set(), "specific": set(), "dist": ""}
                if clusterType in ["specific", "accessory"]:
                    dicoOrgCompGenom[orgName]['variable'].add(header)
                dicoOrgCompGenom[orgName][clusterType].add(header)
        # Count DIST value
        minDIST = 100
        maxDIST = 0
        for orgName1 in dicoGroup[groupName]:
            for orgName2 in dicoGroup[groupName]:
                if orgName1 != orgName2:
                    minDIST = min(minDIST, dicoDIST[orgName1][orgName2])
                    maxDIST = max(minDIST, dicoDIST[orgName1][orgName2])
        for orgName in dicoGroup[groupName]:
            dicoOrgCompGenom[orgName]['dist'] = str(round(minDIST, 2))+"-"+str(round(maxDIST, 2))
            if len(dicoGroup[groupName]) == 1:
                toDisplay += groupName+"\t"+orgName+"\tNA\tNA\tNA\tNA\tNA\n"
            else:
                toDisplay += groupName+"\t"+orgName+"\t"+str(dicoOrgCompGenom[orgName]['dist'])+"\t" + \
                             str(len(dicoOrgCompGenom[orgName]['core']))+"\t"+str(len(dicoOrgCompGenom[orgName]['variable']))+"\t" + \
                             str(len(dicoOrgCompGenom[orgName]['accessory']))+"\t"+str(len(dicoOrgCompGenom[orgName]['specific']))+"\n"
    # ***** OUTPUT ***** #
    if pathOUT == "stdout":
        print(toDisplay)
    else:
        OUT = open(pathOUT, 'w')
        OUT.write(toDisplay)
        OUT.close()


@fct_checker
def make_vibrio_core(pathIN: str, pathIN2: str, pathOUT: str, ref: str, idThrClust: int = 30, covThrClust: int = 80, maxMash: float = 0.3, cutN: int = 5, maxContig: int = 1000, maxL90: int = 100, persRatio: float = 0.9, minPersPart: float = 0.75, minsize: float = 1.0, ext: str = ".faa") -> Tuple[str, str, str, str, int, int, float, int, int, int, float, float, float, str]:
    '''
     ------------------------------------------------------------
    |                VIBRIO CORE PROTEIN FROM FAA                |
    |------------------------------------------------------------|
    |        Construct protein core from vibrio faa files        |
    |------------------------------------------------------------|
    | PARAMETERS                                                 |
    |    pathIN     : path of input faa folder (required)        |
    |    pathIN2    : path of input gzipped gbk folder (required)|
    |    pathOUT    : path of output results folder (required)   |
    |    ref        : reference organism (required)              |
    |    idThrClust : %identity clustering threshold (default=30)|
    |    covThrClust: %coverage clustering threshold (default=80)|
    |    maxMash    : max mash distance (default=0.3)            |
    |    cutN       : number of N to split (default=5)           |
    |    maxContig  : max contig number (default=1000)           |
    |    maxL90     : max L90 value (default=100)                |
    |    persRatio  : max L90 value (default=0.9)                |
    |    minPersPart: max L90 value (default=0.75)               |
    |    minsize    : min genome length (Mb) (default=1.0)       |
    |    ext        : ext of prot fasta files (default=.faa)     |
     ------------------------------------------------------------
    |TOOLS: mash, famsa                                          |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathIN2 = path_converter(pathIN2)
    pathOUT = path_converter(pathOUT)
    slurmBool, cpu, memMax, memMin = get_sys_info()
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'mash' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['mash'])
    if 'famsa' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['famsa'])
    lstFiles, maxpathSize = get_input_files(pathIN, "make_vibrio_core", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_vibrio_core]\nAny input fasta files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if len(lstFiles) == 1:
        printcolor("[ERROR: make_vibrio_core]\nMaking core require more than one FASTA file\n", 1, "212;64;89", "None", True)
        exit_gemini()
    lstFiles2, maxpathSize2 = get_input_files(pathIN2, "make_vibrio_core", [".gz"])
    if len(lstFiles2) == 0:
        printcolor("[ERROR: make_vibrio_core]\nAny input gbk files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if len(lstFiles2) == 1:
        printcolor("[ERROR: make_vibrio_core]\nMaking core require more than one GBK file\n", 1, "212;64;89", "None", True)
        exit_gemini()
    printcolor("♊ Vibrio supertree"+"\n")
    os.makedirs(pathOUT, exist_ok=True)
    pathRBH = pathOUT+"/"+ref+"_id"+str(idThrClust)+"_cov"+str(covThrClust)+".rbh"
    pathDirFASTAOUT = pathOUT+"/FASTA_id"+str(idThrClust)+"_cov"+str(covThrClust)+"_mash"+str(cutN)+"-"+str(maxMash)+"_ctg"+str(maxContig)
    os.makedirs(pathDirFASTAOUT, exist_ok=True)
    pathMashOUT = pathOUT+"/mash.out"
    pathFilterContigJSON = pathOUT+"/filter_contigs.json"
    persRatio = 0.9
    minPersPart = 0.75
    # ***** CREATE ORGANISM LIST ***** #
    setAllOrg = set()
    for pathGBK in lstFiles2:
        orgName2 = os.path.basename(pathGBK).replace(".gbff", "").replace(".gbk", "").replace(".gb", "").replace(".gz", "").replace("_genomic", "")
        for pathFAA in lstFiles:
            orgName1 = os.path.basename(pathFAA).replace(ext, "").replace("."+ext, "")
            if orgName1 == orgName2:
                setAllOrg.add(orgName1)
                break
    printcolor("⏩ Found "+str(len(setAllOrg))+" organisms"+"\n")
    # ***** LAUNCH mmseqs RBH ***** #
    relaunchRBH = True
    if os.path.isfile(pathRBH) and os.path.getsize(pathRBH) > 0:
        relaunchRBH = False
    else:
        for pathPrevRBH in glob.glob(pathOUT+"/"+ref+"_id*_cov*.rbh"):
            prevId = float(pathPrevRBH.split("_id")[1].split("_cov")[0])
            prevCov = float(pathPrevRBH.split("_id")[1].split("_cov")[1].split(".rbh")[0])
            if prevId <= idThrClust and prevCov <= covThrClust and os.path.getsize(pathPrevRBH) > 0:
                pathRBH = pathPrevRBH
                relaunchRBH = False
                break
    if relaunchRBH is True:
        mmseqs_rbh(pathIN=pathIN, pathOUT=pathOUT, pathTMP=geminiset.pathTMP, ref=ref, idThrClust=idThrClust, covThrClust=covThrClust, ext=ext)
        shutil.move(pathOUT+"/"+ref+".rbh", pathRBH)
    # ***** PARSE mmseqs RBH ***** #
    printcolor("♊ Parse rbh"+"\n")
    # Read mmseqs RBH file
    pathJSONdicoPairwiseOrg = pathOUT+"/"+ref+"_id"+str(idThrClust)+"_cov"+str(covThrClust)+"_dicoPairwiseOrg.json"
    pathJSONdicoPerProt = pathOUT+"/"+ref+"_id"+str(idThrClust)+"_cov"+str(covThrClust)+"_dicoPerProt.json"
    if not os.path.isfile(pathJSONdicoPairwiseOrg) or os.path.getsize(pathJSONdicoPairwiseOrg) == 0 or not os.path.isfile(pathJSONdicoPairwiseOrg) or os.path.getsize(pathJSONdicoPairwiseOrg) == 0:
        dicoPairwiseOrg = {}
        dicoPerProt = {}
        MMSEQS = open(pathRBH, 'r')
        lstLines = MMSEQS.read().split("\n")
        MMSEQS.close()
        # Parse and assign cluster
        for line in lstLines:
            if line != "":
                splitLine = line.split("\t")
                query = splitLine[0]
                orgQuery = query.split("|")[2]
                target = splitLine[1]
                orgTarget = target.split("|")[2]
                pident = float(splitLine[2])*100
                align_len = int(splitLine[4])
                qlen = int(splitLine[7])
                tlen = int(splitLine[10])
                minLrap = (align_len/min(qlen, tlen))*100
                maxLrap = (align_len/max(qlen, tlen))*100
                if orgQuery != orgTarget and orgQuery in setAllOrg and pident >= idThrClust and minLrap >= covThrClust and maxLrap >= covThrClust:
                    try:
                        dicoPerProt[query].append(target)
                    except KeyError:
                        dicoPerProt[query] = [target]
                    try:
                        dicoPairwiseOrg[orgTarget] += 1
                    except KeyError:
                        dicoPairwiseOrg[orgTarget] = 1
        dump_json(dicoPairwiseOrg, pathJSONdicoPairwiseOrg)
        dump_json(dicoPerProt, pathJSONdicoPerProt)
    else:
        dicoPairwiseOrg = load_json(pathJSONdicoPairwiseOrg)
        dicoPerProt = load_json(pathJSONdicoPerProt)
    # ***** MASH distance ***** #
    printcolor("♊ Mash"+"\n")
    if not os.path.isfile(pathMashOUT) or os.path.getsize(pathMashOUT) == 0:
        # Make FNA
        pbar = tqdm(total=len(lstFiles2), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        TMPLIST = open(geminiset.pathTMP+"/mash_genome.txt", 'w')
        for pathGBK in lstFiles2:
            orgName = os.path.basename(pathGBK).replace(".gbff", "").replace(".gbk", "").replace(".gb", "").replace(".gz", "").replace("_genomic", "")
            pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
            gbk_to_fna(pathIN=pathGBK, pathOUT=geminiset.pathTMP+"/"+orgName+".fna")
            TMPLIST.write(geminiset.pathTMP+"/"+orgName+".fna\n")
            pbar.update(1)
            title("Mash", pbar)
        pbar.close()
        TMPLIST.close()
        # Launch MASH
        cmdMASH = dicoGeminiPath['TOOLS']['mash']+" dist -l -p "+str(cpu)+" "+geminiset.pathTMP+"/"+ref+".fna "+geminiset.pathTMP+"/mash_genome.txt | sed 's, "+geminiset.pathTMP+"/, ,g' | sed 's, .fna, ,g' > "+pathMashOUT+" > /dev/null 2>&1"
        os.system(cmdMASH)
    # ***** Contigs filters ***** #
    printcolor("♊ Contigs stats"+"\n")
    if not os.path.isfile(pathFilterContigJSON) or os.path.getsize(pathFilterContigJSON) == 0:
        dicoFilterContig = {}
        pbar = tqdm(total=len(lstFiles2), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        for pathGBK in lstFiles2:
            orgName = os.path.basename(pathGBK).replace(".gbff", "").replace(".gbk", "").replace(".gb", "").replace(".gz", "").replace("_genomic", "")
            pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
            totalSize = 0
            lstContigLength = []
            dicoGBK = make_gbk_dict(pathGBK)[orgName]
            for contig in dicoGBK:
                splitNcontig = dicoGBK[contig]['seq'].split("N"*cutN)
                for splitSeq in splitNcontig:
                    totalSize += len(splitSeq)
                    lstContigLength.append(len(splitSeq))
            lstContigLength.sort(reverse=True)
            l90 = 0
            incrSize = 0
            for length in lstContigLength:
                incrSize += length
                l90 += 1
                if incrSize >= (totalSize*0.9):
                    break
            dicoFilterContig[orgName] = {'totalSize': totalSize, 'nbContigs': len(dicoGBK), 'nbSplitNContig': len(lstContigLength), 'L90': l90}
            pbar.update(1)
            title("Contigs stats", pbar)
        pbar.close()
        dump_json(dicoFilterContig, pathFilterContigJSON)
    else:
        dicoFilterContig = load_json(pathFilterContigJSON)
    # ***** FILTER organism ***** #
    printcolor("♊ Assembly filtering"+"\n")
    setExcludeOrg = set()
    # Parse mash distance
    MASH = open(pathMashOUT, 'r')
    lstLines = MASH.read().split("\n")
    MASH.close()
    dicoMASH = {}
    for line in lstLines:
        if line != "":
            splitLine = line.split("\t")
            orgName = splitLine[1]
            distance = float(splitLine[2])
            dicoMASH[orgName] = distance
            if orgName in setAllOrg and distance >= maxMash:
                setExcludeOrg.add(orgName)
    # Filtering based on contigs stats
    for orgName in dicoFilterContig:
        if orgName in setAllOrg and (dicoFilterContig[orgName]['totalSize'] < minsize*1000000 or dicoFilterContig[orgName]['nbSplitNContig'] >= maxContig or dicoFilterContig[orgName]['L90'] >= maxL90):
            setExcludeOrg.add(orgName)
    # Filtering JCM genome
    for orgName in setAllOrg:
        if orgName in setAllOrg and ("JCM" in orgName.upper() or orgName == "Vibrio_metoecus_08-2459"):
            setExcludeOrg.add(orgName)
    # ***** RETRIEVE LIST OF PERSISTENT CORE PROTEINS ***** #
    nbCoreOrg = int((len(setAllOrg)-len(setExcludeOrg))*persRatio)
    printcolor("♊ Persistent core"+"\n")
    lstCoreProt = []
    for prot in dicoPerProt:
        setOrg = set([ref])
        for target in dicoPerProt[prot]:
            orgTarget = target.split("|")[2]
            if orgTarget in setAllOrg and orgTarget not in setExcludeOrg:
                setOrg.add(orgTarget)
        if len(setOrg) >= nbCoreOrg:
            lstCoreProt.append(prot)
    # Filter organism without more of "minPersPart" persistent proteins
    printcolor("♊ Persistent filtering"+"\n")
    dicoPersistentPerOrg = {}
    setFinalOrgList = set()
    for prot in lstCoreProt:
        for target in dicoPerProt[prot]:
            orgTarget = target.split("|")[2]
            if orgTarget in setAllOrg and orgTarget not in setExcludeOrg:
                try:
                    dicoPersistentPerOrg[orgTarget] += 1
                except KeyError:
                    dicoPersistentPerOrg[orgTarget] = 1
    for orgName in dicoPersistentPerOrg:
        if dicoPersistentPerOrg[orgName] < len(lstCoreProt)*minPersPart:
            setExcludeOrg.add(orgName)
        else:
            setFinalOrgList.add(orgName)
    printcolor("⏩ Found "+str(len(lstCoreProt))+" core proteins"+"\n")
    printcolor("⏩ "+str(len(setFinalOrgList))+" retained organisms"+"\n")
    printcolor("⏩ "+str(len(setExcludeOrg))+" filtered organisms"+"\n")
    # ***** CREATE CORE FAA FILES ***** #
    # Retrieve all proteins sequence
    printcolor("♊ Read proteins"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoProteins = {}
    for pathFAA in lstFiles:
        file = os.path.basename(pathFAA)
        orgName = file.replace(ext, "").replace("."+ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        dicoProteins[orgName] = make_fasta_dict(pathFAA)
        pbar.update(1)
        title("Read prot", pbar)
    pbar.close()
    # Create individual core protein fasta
    printcolor("♊ Protein FASTAs"+"\n")
    pbar = tqdm(total=len(lstCoreProt), ncols=50, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for prot in lstCoreProt:
        ltRef = prot.split("|")[0]
        pathFAA = pathDirFASTAOUT+"/"+ltRef+".faa"
        pbar.set_description_str(ltRef+" ".rjust(maxpathSize-len(ltRef)))
        if not os.path.isfile(pathFAA) or os.path.getsize(pathFAA) == 0:
            lstOrg = []
            FAA = open(pathFAA, 'w')
            FAA.write(">"+ref+"\n"+dicoProteins[ref][prot]+"\n")
            for target in dicoPerProt[prot]:
                orgTarget = target.split("|")[2]
                if orgTarget in setAllOrg and orgTarget not in setExcludeOrg and orgTarget not in lstOrg:
                    FAA.write(">"+orgTarget+"\n"+dicoProteins[orgTarget][target]+"\n")
                    lstOrg.append(orgTarget)
            FAA.close()
        pbar.update(1)
        title("Prot FASTAs", pbar)
    pbar.close()
    # ***** ALIGN INDIVIDUAL CORE ***** #
    printcolor("♊ Protein Alignments"+"\n")
    pbar = tqdm(total=len(lstCoreProt), ncols=50, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for prot in lstCoreProt:
        ltRef = prot.split("|")[0]
        pathFAA = pathDirFASTAOUT+"/"+ltRef+".faa"
        pathALIGN = pathDirFASTAOUT+"/"+ltRef+".align"
        pbar.set_description_str(ltRef+" ".rjust(maxpathSize-len(ltRef)))
        if not os.path.isfile(pathALIGN) or os.path.getsize(pathALIGN) == 0:
            cmdFAMSA = dicoGeminiPath['TOOLS']['famsa']+" -t "+str(cpu)+" "+pathFAA+" "+pathALIGN+" > /dev/null 2>&1"
            os.system(cmdFAMSA)
        pbar.update(1)
        title("Prot align", pbar)
    pbar.close()
    # ***** CONCAT INDIVIDUAL CORE ***** #
    printcolor("♊ Cat core align"+"\n")
    pathCoreFAA = pathDirFASTAOUT+"/allCoreProt.align"
    if not os.path.isfile(pathCoreFAA) or os.path.getsize(pathCoreFAA) == 0:
        # Read all alignment
        dicoCoreAlign = {}
        pbar = tqdm(total=len(lstCoreProt), ncols=50, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
        for prot in lstCoreProt:
            ltRef = prot.split("|")[0]
            pathALIGN = pathDirFASTAOUT+"/"+ltRef+".align"
            pathTMP = geminiset.pathTMP+"/"+ltRef+".fasta"
            TMP = open(pathTMP, 'w')
            dicoAlign = make_fasta_dict(pathALIGN)
            pbar.set_description_str(ltRef+" ".rjust(maxpathSize-len(ltRef)))
            for key in dicoAlign:
                TMP.write(">"+key+"\n"+dicoAlign[key]+"\n")
                try:
                    dicoCoreAlign[key] += dicoAlign[key]
                except KeyError:
                    dicoCoreAlign[key] = dicoAlign[key]
                size = len(dicoAlign[key])
            # Add missing organism
            for orgName in setFinalOrgList:
                if orgName not in dicoAlign:
                    TMP.write(">"+orgName+"\n"+"-"*size+"\n")
                    try:
                        dicoCoreAlign[orgName] += "-"*size
                    except KeyError:
                        dicoCoreAlign[orgName] = "-"*size
            TMP.close()
            pbar.update(1)
            title("Cat core align", pbar)
        pbar.close()
        # Concatenate per organism
        CORE = open(pathCoreFAA, 'w')
        for org in dicoCoreAlign:
            CORE.write(">"+org+"\n"+dicoCoreAlign[org]+"\n")
        CORE.close()


@fct_checker
def ppanggolin(pathIN: str, pathIN2: str, pathOUT: str, maxRGP: int = -1, prefix: str = "None", ext: str = ".gbk.gz") -> Tuple[str, str, str, int, str, str]:
    '''
     ------------------------------------------------------------
    |                         PPanGGOLiN                         |
    |------------------------------------------------------------|
    |                  PPanGGOLiN RGP analysis                   |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input files or folder (required)       |
    |    pathIN2: path of ordered organism file (required)       |
    |    pathOUT: path of output files (required)                |
    |    maxRGP : maximum number of common RGP (default=-1)      |
    |    prefix : organism prefix (default=None)                 |
    |    ext    : ext of prot fasta files (default=.gbk.gz)      |
     ------------------------------------------------------------
    |TOOLS: ppanggolin, circos                                   |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "ppanggolin", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: ppanggolin]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'ppanggolin' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['ppanggolin'])
    if 'circos' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['circos'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathDIRTMP = geminiset.pathTMP
    # ***** Ordered organisms list ***** #
    pathIN2 = path_converter(pathIN2)
    IN2 = open(pathIN2, 'r')
    lstLines = IN2.read().split("\n")
    IN2.close()
    lstOrderedOrg = []
    for line in lstLines:
        if line != "" and line[0] != "#":
            lstOrderedOrg.append(line)
    # ***** PPanGGOLiN panRGP ***** #
    pathOUTppanggolin = pathOUT+"/ppanggolinRGP"
    pathLOGppanggolin = pathOUT+"/ppanggolinRGP.log"
    if os.path.isfile(pathOUTppanggolin+"/plastic_regions.tsv"):
        printcolor("♊ PPanGGOLiN panRGP"+"\n")
    else:
        spinner = yaspin(Spinners.aesthetic, text="♊ PPanGGOLiN panRGP", side="right")
        spinner.start()
        title("PPanGGOLiN", None)
        # Create genbank list file
        pathTMPlst = pathDIRTMP+"/lst_gbk.txt"
        LST = open(pathTMPlst, 'w')
        for pathFile in lstFiles:
            LST.write(os.path.basename(pathFile).replace(ext, "")+"\t"+pathFile+"\n")
        LST.close()
        # Launch PPanGGOLiN
        cmdPPANGGOLINRGP = dicoGeminiPath['TOOLS']['ppanggolin']+" panrgp --anno "+pathTMPlst+" --output "+pathOUTppanggolin+" --tmpdir "+pathDIRTMP + \
            " --cpu "+str(cpu)+" --disable_prog_bar --log "+pathLOGppanggolin+" > "+pathLOGppanggolin+" 2>&1"
        os.system(cmdPPANGGOLINRGP)
        os.system("mv "+pathLOGppanggolin+" "+pathOUTppanggolin)
        spinner.stop()
        printcolor("♊ PPanGGOLiN panRGP"+"\n")
    # ***** PPanGGOLiN parser ***** #
    printcolor("♊ Parse PPanGGOLiN RGPs"+"\n")
    # Parse plastic_regions.tsv
    RGP = open(pathOUTppanggolin+"/plastic_regions.tsv", 'r')
    lstLines = RGP.read().split("\n")[1: -1]
    RGP.close()
    dicoRGPtoOrg = {}
    for line in lstLines:
        splitLine = line.split("\t")
        rgpName = splitLine[0]
        orgName = splitLine[1]
        dicoRGPtoOrg[rgpName] = orgName
    # Parse spots.tsv
    SPOT = open(pathOUTppanggolin+"/spots.tsv", 'r')
    lstLines = SPOT.read().split("\n")[1: -1]
    SPOT.close()
    dicoRGP = {}
    for line in lstLines:
        splitLine = line.split("\t")
        spotId = splitLine[0]
        rgpName = splitLine[1]
        orgName = dicoRGPtoOrg[rgpName]
        try:
            dicoRGP[spotId].add(orgName)
        except KeyError:
            dicoRGP[spotId] = set([orgName])
    printcolor("♊ Make PPanGGOLiN RGP matrix"+"\n")
    # ***** Make RGP matrix ***** #
    # Construct symmetrix matrix dictionnary
    dicoRGPmatrix = {}
    for rgp in dicoRGP:
        for orgName1 in dicoRGP[rgp]:
            if orgName1 not in dicoRGPmatrix:
                dicoRGPmatrix[orgName1] = {}
            for orgName2 in dicoRGP[rgp]:
                if orgName2 not in dicoRGPmatrix[orgName1]:
                    dicoRGPmatrix[orgName1][orgName2] = 1
                else:
                    dicoRGPmatrix[orgName1][orgName2] += 1
    # Write to output file
    pathOUTmatrix = pathOUT+"/rgp_matrix.csv"
    MATRIX = open(pathOUTmatrix, 'w')
    MATRIX.write("Organism;"+";".join(lstOrderedOrg)+"\n")
    for orgName1 in lstOrderedOrg:
        line = orgName1
        for orgName2 in lstOrderedOrg:
            if orgName2 in dicoRGPmatrix[orgName1]:
                line += ";"+str(dicoRGPmatrix[orgName1][orgName2])
            else:
                line += ";"
        MATRIX.write(line+"\n")
    MATRIX.close()
    # ***** Make circos plot ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Make Circos RGP plot", side="right")
    spinner.start()
    title("Plotting", None)
    # List number of commons RGPs
    lstNbCommon = set()
    for orgName1 in dicoRGPmatrix:
        for orgName2 in dicoRGPmatrix[orgName1]:
            if orgName1 != orgName2:
                try:
                    lstNbCommon.add(dicoRGPmatrix[orgName1][orgName2])
                except KeyError:
                    pass
    if maxRGP == -1:
        maxRGP = max(lstNbCommon)
    # Circos conf header
    CIRCOSCONF = open(pathDIRTMP+"/circos.conf", 'w')
    CIRCOSCONF.write("# KARYOTYPE\n")
    CIRCOSCONF.write("karyotype = "+pathDIRTMP+"/karyotype.txt\n")
    CIRCOSCONF.write("# CHROMOSOMES\n")
    CIRCOSCONF.write("chromosomes_units = 1000\n\n")
    CIRCOSCONF.write("<links>\n")
    CIRCOSCONF.write("z = 0\n")
    CIRCOSCONF.write("radius = 0.99r\n")
    CIRCOSCONF.write("crest = 0.5\n")
    CIRCOSCONF.write("bezier_radius = 0.25r\n")
    CIRCOSCONF.write("bezier_radius_purity = 1\n\n")
    # Color gradient
    lstColor = linear_gradient("#FFFFFF", "#000080", maxRGP)['hex']
    # Circos conf links blocks
    for i in range(1, maxRGP+1, 1):
        if i in lstNbCommon:
            CIRCOSCONF.write("<link>\n")
            CIRCOSCONF.write("color = "+lstColor[i-1].replace("#", "")+"\n")
            CIRCOSCONF.write("thickness = 2\n")
            CIRCOSCONF.write("file="+pathDIRTMP+"/links_"+str(i)+"_rgps.txt\n")
            CIRCOSCONF.write("record_limit = "+str(len(lstFiles)*len(lstFiles))+"\n")
            CIRCOSCONF.write("</link>\n\n")
    CIRCOSCONF.write("</links>\n\n")
    # Create links files
    dicoLinks = {}
    cptId1 = 1
    for orgName1 in lstOrderedOrg:
        cptId2 = 1
        for orgName2 in lstOrderedOrg:
            if orgName1 != orgName2:
                try:
                    dicoLinks[str(dicoRGPmatrix[orgName1][orgName2])] += "hs"+str(cptId1)+" 100 900 hs"+str(cptId2)+" 100 900\n"
                except KeyError:
                    dicoLinks[str(dicoRGPmatrix[orgName1][orgName2])] = "hs"+str(cptId1)+" 100 900 hs"+str(cptId2)+" 100 900\n"
            cptId2 += 1
        cptId1 += 1
    for num in dicoLinks:
        LINKFILE = open(pathDIRTMP+"/links_"+num+"_rgps.txt", 'w')
        LINKFILE.write(dicoLinks[num])
        LINKFILE.close()
    # Circos conf footer
    CIRCOSCONF.write("<<include "+pathDIRTMP+"/ideogram.conf>>\n\n")
    CIRCOSCONF.write("<<include "+pathDIRTMP+"/ticks.conf>>\n")
    CIRCOSCONF.write("<image>\n")
    CIRCOSCONF.write("dir = "+pathOUT+"\n")
    CIRCOSCONF.write("file=rgp_circos.png\n")
    CIRCOSCONF.write("png = yes\n")
    CIRCOSCONF.write("svg = yes\n")
    CIRCOSCONF.write("radius = 1500p\n")
    CIRCOSCONF.write("angle_offset = -90\n")
    CIRCOSCONF.write("auto_alpha_colors = yes\n")
    CIRCOSCONF.write("auto_alpha_steps  = 5\n")
    CIRCOSCONF.write("</image>\n")
    CIRCOSCONF.write("<<include colors_fonts_patterns.conf>>\n")
    CIRCOSCONF.write("<<include housekeeping.conf>>\n")
    CIRCOSCONF.write("track_defaults* = undef\n")
    CIRCOSCONF.close()
    # Circos other conf files
    TICKS = open(pathDIRTMP+"/ticks.conf", 'w')
    TICKS.write("show_ticks = no\nshow_tick_labels = yes\nshow_grid = yes\n\n<ticks>\nradius = 1.15r\ncolor = black\nthickness = 2p\nmultiplier = 1\nformat = %d\ngrid_start = 0.7r\ngrid_end = 1.1r\ngrid_color = grey\ngrid_thickness = 2p\nsize = 20p\nspacing = 100u\nlabel_size = 30p\nshow_label = yes\nlabel_offset = 10p\ngrid = yes\n\n<tick>\nspacing = 1u\n</tick>\n\n<tick>\nshow_label = yes\nposition = start\nlabel = -\nformat = %s\n</tick>\n\n</ticks>")
    TICKS.close()
    IDEOGRAM = open(pathDIRTMP+"/ideogram.conf", 'w')
    IDEOGRAM.write("<ideogram>\n\n<spacing>\ndefault=0.005r\n</spacing>\n\nradius = 0.77r\nthickness = 40p\nfill = yes\n\nshow_label = yes\nlabel_font = bold\nlabel_radius = 1.02r\nlabel_center = no\nlabel_size = 25\nlabel_parallel = no\n\n<<include "+pathDIRTMP+"/bands.conf>>\n</ideogram>\n")
    IDEOGRAM.close()
    BANDS = open(pathDIRTMP+"/bands.conf", 'w')
    BANDS.write("show_bands = yes\nfill_bands = yes\nband_transparency = 0\nband_stroke_thickness = 2p\nband_stroke_color = dgrey\n")
    BANDS.close()
    KARYOTYPE = open(pathDIRTMP+"/karyotype.txt", 'w')
    cpt = 1
    for orgName1 in lstOrderedOrg:
        if prefix != 'None':
            KARYOTYPE.write("chr - hs"+str(cpt)+" "+orgName1.replace(prefix, "")+" 0 1000 black\n")
        else:
            KARYOTYPE.write("chr - hs"+str(cpt)+" "+orgName1[-10:]+" 0 1000 black\n")
        cpt += 1
    KARYOTYPE.close()
    # Launch circos
    cmdCIRCOS = dicoGeminiPath['TOOLS']['circos']+" -conf "+pathDIRTMP+"/circos.conf > /dev/null 2>&1"
    os.system(cmdCIRCOS)
    spinner.stop()
    printcolor("♊ Make Circos RGP plot"+"\n")

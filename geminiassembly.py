'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ assembly
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminiassembly.py               |
description      : gemini assembly functions       |
author           : dooguypapua                     |
lastmodification : 20220805                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import os
import sys
import shutil
import geminiset
from tabulate import tabulate
from typing import Tuple
from tqdm import tqdm
from yaspin import yaspin
from yaspin.spinners import Spinners
from geminini import path_converter, get_input_files, printcolor, get_gemini_path, get_sys_info, fct_checker
from geminini import title, exit_gemini, read_file, launch_threads
from geminiparse import unwrap_fasta, make_blast_dict, make_fasta_dict, check_circular


@fct_checker
def phage_assembly(pathIN: str, pathOUT: str) -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |                       PHAGE ASSEMBLY                       |
     ------------------------------------------------------------
    |    Trim FASTQ & de novo assembly using spades (-careful)   |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathOUT : path of output files (required)               |
     ------------------------------------------------------------
    |TOOLS: trimmomatic, spades                                  |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "phage_assembly", [".fastq.gz", ".fastq", ".fq.gz", ".fq"])
    if len(lstFiles) == 0:
        printcolor("[ERROR: phage_assembly]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'trimmomatic' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['trimmomatic'])
    if 'spades' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['spades'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    if pathOUT == "":
        printcolor("[ERROR: phage_assembly]\nMissing '-o'pathOUT\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    trimmomatic_param = "LEADING: 3 TRAILING: 3 SLIDINGWINDOW: 4:15 MINLEN: 36"
    spades_param = "--careful --cov-cutoff auto -k 21, 33, 55, 77 -m 10"
    pbar = tqdm(total=int(len(lstFiles)/2), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFASTQ in lstFiles:
        file = os.path.basename(pathFASTQ)
        if "_R1" in file:
            sampleName = os.path.basename(pathFASTQ).split("_R1")[0]
            pathFASTA = pathOUT+"/"+sampleName+".fasta"
            if not os.path.isfile(pathFASTA):
                # Construct requierd path
                pathFASTQr1 = pathFASTQ
                pathFASTQr2 = pathFASTQ.replace("_R1", "_R2")
                if not os.path.isfile(pathFASTQr2):
                    printcolor("[ERROR: phage_assembly]\nUnable to find R2 for \""+pathFASTQr1+"\"\n", 1, "212;64;89", "None", True)
                    exit_gemini()
                pathTRIMr1 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr1).replace("_R1", "_trimmed_R1")
                pathTRIMr2 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr2).replace("_R2", "_trimmed_R2")
                pathUNTRIMr1 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr1).replace("_R1", "_untrimmed_R1")
                pathUNTRIMr2 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr2).replace("_R2", "_untrimmed_R2")
                pathUNTRIM = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr1).replace("_R1", "_untrimmed")
                pathSPADES = geminiset.pathTMP+"/spades_"+sampleName
                pathSPADESlog = geminiset.pathTMP+"/spades_"+sampleName+".log"
                # Launch command
                cmdTrim = "java -jar "+dicoGeminiPath['TOOLS']['trimmomatic']+" PE -threads "+str(cpu)+" -quiet "+pathFASTQr1+" "+pathFASTQr2+" "+pathTRIMr1+" "+pathUNTRIMr1+" "+pathTRIMr2+" "+pathUNTRIMr2+" "+trimmomatic_param
                pbar.set_description_str(sampleName+": trimmomatic"+" ".rjust(maxpathSize-len(sampleName)))
                os.system(cmdTrim)
                cmdCat = "cat "+pathUNTRIMr1+" "+pathUNTRIMr2+" > "+pathUNTRIM
                os.system(cmdCat)
                cmdSpades = dicoGeminiPath['TOOLS']['spades']+" -1 "+pathTRIMr1+" -2 "+pathTRIMr2+" -s "+pathUNTRIM+" -o "+pathSPADES+" --threads "+str(cpu)+" "+spades_param+" > "+pathSPADESlog
                pbar.set_description_str(sampleName+": spades"+" ".rjust(maxpathSize-len(sampleName)+5))
                os.system(cmdSpades)
                # Copy output contigs FASTA file
                shutil.copyfile(pathSPADES+"/contigs.fasta", pathFASTA)
            # Unwrap FASTA
            unwrap_fasta(pathFASTA)
            pbar.update(1)
            title("Phage ass", pbar)
    pbar.close()


@fct_checker
def filter_phage_assembly(pathIN: str, pathOUT: str, minLen: int = 100, minCov: int = 2, ext: str = ".fasta") -> Tuple[str, str, int, int, str]:
    '''
     ------------------------------------------------------------
    |                     FILTER PHAGE ASSEMBLY                  |
     ------------------------------------------------------------
    |               Filter phage assembly contigs                |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    pathIN : path of input files or folder (required)       |
    |    pathOUT: path of output files (required)                |
    |    minLen : minimum contig length (default=100)            |
    |    minCov : minimum contig coverage (default=2)            |
    |    ext    : extension of input files (default=.fasta)      |
     ------------------------------------------------------------
    |TOOLS: blastn                                               |
     ------------------------------------------------------------
    '''
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "filter_phage_assembly", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: filter_phage_assembly]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    os.makedirs(pathOUT, exist_ok=True)
    tabulateTable = []
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'blastn' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['blastn'])
    printcolor("♊ Filter assembly"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathAssembly in lstFiles:
        file = os.path.basename(pathAssembly)
        orgName = file.replace(ext, "").replace("."+ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        # Blast contigs against contamination sequence
        cmdBlastn = dicoGeminiPath['TOOLS']['blastn']+" -task blastn -query "+pathAssembly+" -db "+dicoGeminiPath['DATABASES']['contaPhageDB']+" -out "+geminiset.pathTMP+"/blast.xml -outfmt 5 -perc_identity 80 -qcov_hsp_perc 30 -max_hsps 10 -max_target_seqs 5"
        os.system(cmdBlastn)
        dicoBlast = make_blast_dict(pathIN=geminiset.pathTMP, dbName="blast", pathJSON="None", idThr=30, minLRthr=0, maxLRthr=0, ext=".xml")
        lstContamContig = list(dicoBlast.keys())
        # Read contigs FASTA
        dicoAssembly = make_fasta_dict(pathAssembly)
        dicoFiltered = {}
        totalSize = 0
        for contig in dicoAssembly:
            contigCov = float(contig.split("_cov_")[1])
            if contigCov >= minCov and len(dicoAssembly[contig]) >= minLen and contig not in lstContamContig:
                oneNuclContig = False
                for nucl in ["A", "T", "G", "C"]:
                    if dicoAssembly[contig].count(nucl) == len(dicoAssembly[contig]):
                        oneNuclContig = True
                        break
                if not oneNuclContig:
                    overlap = check_circular(dicoAssembly[contig])
                    if overlap == "":
                        dicoFiltered["CONTIG_"+str(len(dicoFiltered)+1)+"|cov: "+str(contigCov)] = dicoAssembly[contig]
                    else:
                        dicoFiltered["CONTIG_"+str(len(dicoFiltered)+1)+"|cov: "+str(contigCov)+"|circular: "+overlap] = dicoAssembly[contig]
                    totalSize += len(dicoAssembly[contig])
        # Write new filtered FASTA
        pathOUTassembly = pathOUT+"/"+file
        OUT = open(pathOUTassembly, 'w')
        for key in dicoFiltered:
            OUT.write(">"+key+"\n"+dicoFiltered[key]+"\n")
        OUT.close()
        # Results table
        tabulateTable.append([orgName, str(len(dicoAssembly)), str(len(dicoFiltered)), overlap != "", str(totalSize)])
        pbar.update(1)
        title("Filter assembly", pbar)
    pbar.close()
    printcolor("♊ Results\n")
    printcolor(tabulate(tabulateTable, headers=["orgName", "Initial contigs", "Final contigs", "Circular", "Final size"], tablefmt="pretty")+"\n")


@fct_checker
def replicon_distribution(pathIN: str, pathREF: str, pathOUT: str, idThr: int = 80, extIN: str = ".fasta", extREF: str = ".fasta") -> Tuple[str, str, str, int, str, str]:
    '''
     ------------------------------------------------------------
    |                     REPLICON DISTRIBUTION                  |
     ------------------------------------------------------------
    |            Search replicons in assembly contigs            |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    pathIN : path of input files or folder (required)       |
    |    pathREF: path of reference files (required)             |
    |    pathOUT: path of output tsv (required)                  |
    |    idThr  : %identity threshold (default=80)               |
    |    extIN  : extension of input files (default=.fasta)      |
    |    extREF : extension of reference files (default=.fasta)  |
     ------------------------------------------------------------
    '''
    pathOUT = path_converter(pathOUT)
    lstQuery, maxpathSize1 = get_input_files(pathIN, "replicon_distribution", [extIN])
    if len(lstQuery) == 0:
        printcolor("[ERROR: replicon_distribution]\nAny query files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    lstRef, maxpathSize2 = get_input_files(pathREF, "replicon_distribution", [extREF])
    if len(lstQuery) == 0:
        printcolor("[ERROR: replicon_distribution]\nAny reference files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    # Variables
    maxpathSize = max(maxpathSize1, maxpathSize2)
    setRefRepliconPath = set()
    setRefRepliconName = set()
    dicoRepliconCov = {}
    blastoutfmt = "6 delim=; qseqid sseqid sstart send bitscore slen"
    # Split reference replicons in distinct FASTA files (required header ">repliconType [orgName]")
    printcolor("♊ Split replicons"+"\n")
    pbar = tqdm(total=len(lstRef), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFile in lstRef:
        fileName = os.path.basename(pathFile).replace(extREF, "").replace("."+extREF, "")
        pbar.set_description_str(fileName+" ".rjust(maxpathSize-len(fileName)))
        dicoFASTA = make_fasta_dict(pathFile)
        for key in dicoFASTA:
            repliconType = key.replace("|", " ").split(" ")[0]
            orgName = key.split("[")[1].replace("]", "").replace(" ", "_")
            pathREPLICON = geminiset.pathTMP+"/"+repliconType+".fasta"
            setRefRepliconPath.add(pathREPLICON)
            setRefRepliconName.add(repliconType)
            REPLICON = open(pathREPLICON, 'a')
            REPLICON.write(">"+orgName+"\n"+dicoFASTA[key]+"\n")
            REPLICON.close()
        pbar.update(1)
        title("Split replicons", pbar)
    pbar.close()
    # Blast each query on each replicon fasta
    spinner = yaspin(Spinners.aesthetic, text="♊ Blast replicons", side="right")
    spinner.start()
    title("Blast", None)
    dicoThread = {}
    for pathFile in lstQuery:
        orgName = os.path.basename(pathFile).replace(extIN, "").replace("."+extIN, "")
        # Foreach replicon FASTA
        for pathREPLICON in setRefRepliconPath:
            repliconName = os.path.basename(pathREPLICON).replace(".fasta", "")
            # Launch blastn
            pathREPLICONOUT = geminiset.pathTMP+"/"+orgName+"_____"+repliconName+".out"
            cmdBlastn = dicoGeminiPath["TOOLS"]["blastn"]+" -query "+pathFile+" -subject "+pathREPLICON+" -out "+pathREPLICONOUT+" -outfmt \""+blastoutfmt+"\" -perc_identity "+str(idThr)
            dicoThread[orgName+"_"+repliconName] = {"cmd": cmdBlastn, "returnstatut": None, "returnlines": []}
    # Launch threads
    if len(dicoThread) > 0:
        launch_threads(dicoThread, "blast", geminiset.cpu, geminiset.pathTMP, spinner)
    spinner.stop()
    printcolor("♊ Blast replicons"+"\n")
    # Parse blast results
    printcolor("♊ Search replicons"+"\n")
    pbar = tqdm(total=len(lstQuery), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFile in lstQuery:
        orgName = os.path.basename(pathFile).replace(extIN, "").replace("."+extIN, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        dicoQuery = {}
        dicoSubject = {}
        dicoRepliconCov[orgName] = {}
        for pathREPLICON in setRefRepliconPath:
            repliconName = os.path.basename(pathREPLICON).replace(".fasta", "")
            pathREPLICONOUT = geminiset.pathTMP+"/"+orgName+"_____"+repliconName+".out"
            if os.path.isfile(pathREPLICONOUT):
                lstLine = read_file(pathREPLICONOUT, yaspinBool=False)
                for line in lstLine:
                    splitLine = line.split(";")
                    query = splitLine[0]
                    bitscore = float(splitLine[4])
                    # Keep best for each contig query
                    if query not in dicoQuery:
                        dicoQuery[query] = {}
                    if repliconName not in dicoQuery[query]:
                        dicoQuery[query][repliconName] = {'bestScore': bitscore, 'line': []}
                    dicoQuery[query][repliconName]['bestScore'] = max(dicoQuery[query][repliconName]['bestScore'], bitscore)
                    dicoQuery[query][repliconName]['line'].append(line)
        # Compute coverage for each bestScore replicon
        for query in dicoQuery:
            bestScore = 0
            bestReplicon = ""
            for repliconName in dicoQuery[query]:
                if dicoQuery[query][repliconName]['bestScore'] > bestScore:
                    bestScore = dicoQuery[query][repliconName]['bestScore']
                    bestReplicon = repliconName
            if bestReplicon not in dicoSubject:
                dicoSubject[bestReplicon] = {}
            # Parse correspunding blast line
            for line in dicoQuery[query][bestReplicon]['line']:
                splitLine = line.split(";")
                query = splitLine[0]
                subject = splitLine[1]
                sstart = int(splitLine[2])
                send = int(splitLine[3])
                bitscore = float(splitLine[4])
                slen = int(splitLine[5])
                # Add covered position
                for i in range(min(sstart, send), max(sstart, send)+1, 1):
                    try:
                        dicoSubject[bestReplicon][subject]['cov_pos'][i] = 1
                    except KeyError:
                        dicoSubject[bestReplicon][subject] = {'cov_pos': {i: 1}, 'len': slen}
        # Look for coverage per replicon
        for repliconName in dicoSubject:
            bestCov = 0
            for subject in dicoSubject[repliconName]:
                percentCov = len(dicoSubject[repliconName][subject]['cov_pos'])*100/dicoSubject[repliconName][subject]['len']
                bestCov = max(bestCov, percentCov)
            dicoRepliconCov[orgName][repliconName] = bestCov
        pbar.update(1)
        title("Search replicons", pbar)
    pbar.close()
    # Write output
    printcolor("♊ Write output"+"\n")
    OUT = open(pathOUT, 'w')
    OUT.write("organism")
    # Sorted replicon name with chromosome as first
    lstSortedRepliconChrom = []
    lstSortedRepliconOthers = []
    for repliconName in setRefRepliconName:
        if "chr" in repliconName or "chrom" in repliconName or "chromosome" in repliconName:
            lstSortedRepliconChrom.append(repliconName)
        else:
            lstSortedRepliconOthers.append(repliconName)
    lstSortedRepliconChrom.sort()
    lstSortedRepliconOthers.sort()
    lstSortedReplicon = lstSortedRepliconChrom+lstSortedRepliconOthers
    # Write header
    for repliconName in lstSortedReplicon:
        OUT.write("\t"+repliconName)
    OUT.write("\n")
    # Write % coverage
    for orgName in dicoRepliconCov:
        OUT.write(orgName)
        for repliconName in lstSortedReplicon:
            if repliconName in dicoRepliconCov[orgName]:
                OUT.write("\t"+str(round(dicoRepliconCov[orgName][repliconName], 1)).replace(".", ","))
            else:
                OUT.write("\t0")
        OUT.write("\n")
    OUT.close()

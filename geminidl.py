'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ dl
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminidl.py                     |
description      : gemini download functions       |
author           : dooguypapua                     |
lastmodification : 20210713                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import os
import re
import sys
import shutil
import geminiset
import time
import hashlib
from yaspin import yaspin
from yaspin.spinners import Spinners
from tqdm import tqdm
from ftplib import FTP
from ete3 import NCBITaxa
from typing import Tuple
import requests
import subprocess
from geminini import printcolor, fct_checker, path_converter, get_gemini_path, get_sys_info, title, exit_gemini


'''
-------------------------------------------------------------------------------------------------------
                                          REQUESTS FUNCTIONS
-------------------------------------------------------------------------------------------------------
'''


# ***** GET DISTANT FILE ***** #
def get_distant_file(url, path, pbar):
    try:
        request = requests.get(url, timeout=30)
        with open(path, 'wb') as file:
            file.write(request.content)
    except Exception as e:
        printcolor("[ERROR: Download \""+url+"\"] "+str(e)+"\n", 1, "212;64;89", "None", True)
        exit_gemini()


# ***** GET DISTANT FILE WITH PROGRESSBAR ***** #
def get_distant_file_tqdm(url, path, filesize, chunksize):
    try:
        with requests.get(url, stream=True) as r, open(path, "wb") as f, tqdm(unit="B", unit_scale=True, leave=False, unit_divisor=1024, total=filesize, ncols=75, file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{rate_fmt}{postfix}]") as progress:
            for chunk in r.iter_content(chunk_size=chunksize):
                # download the file chunk by chunk
                datasize = f.write(chunk)
                # on each chunk update the progress bar.
                progress.update(datasize)
    except Exception as e:
        printcolor("[ERROR: Download \""+url+"\"] "+str(e)+"\n", 1, "212;64;89", "None", True)
        exit_gemini()


# ***** GET DISTANT FILE SIZE ***** #
def get_ftp_file_size(urlFTP, urlFTPSummary):
    try:
        ftp = FTP(urlFTP)
        ftp.login()
        ftp.voidcmd('TYPE I')
        filesize = ftp.size(urlFTPSummary)
        ftp.quit()
    except Exception as e:
        printcolor("[ERROR: FTP summary size request] "+str(e)+"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    return filesize


# ***** GET DISTANT FILE ETAG ***** #
def get_etag(pathEtag, url):
    # local file etag
    if url is None:
        if not os.path.isfile(pathEtag):
            return None
        else:
            with open(pathEtag, 'r') as f:
                return f.read()
    # distant file etag
    else:
        etag = ""
        exception = ""
        nbAttempt = 0
        while etag == "" and nbAttempt < 10:
            try:
                with requests.Session() as s:
                    etag = s.head(url).headers['ETag'].replace("\"", "")
            except Exception as e:
                exception = e
                pass
            nbAttempt += 1
            time.sleep(5)
        if etag == "":
            return "ERROR any etag "+str(exception)
        return etag


'''
-------------------------------------------------------------------------------------------------------
                                           OTHERS FUNCTIONS
-------------------------------------------------------------------------------------------------------
'''


# ***** CONTROL MD5 against md5sum file ***** #
def file_as_bytes(file):
    with file:
        return file.read()


def checkMD5(pathFile, pathMD5):
    if not os.path.isfile(pathFile) or os.path.getsize(pathFile) == 0:
        return False
    if not os.path.isfile(pathMD5) or os.path.getsize(pathMD5) == 0:
        return False
    # Make file md5
    fileMD5 = hashlib.md5(file_as_bytes(open(pathFile, 'rb'))).hexdigest()
    # Read md5sumfile
    MD5SUM = open(pathMD5, 'r')
    lstLines = MD5SUM.read().split("\n")
    MD5SUM.close()
    findMD5 = False
    for line in lstLines:
        if line != "":
            splitLine = re.sub(r"\s+", "\t", line).split("\t")
            md5 = splitLine[0]
            filename = splitLine[1]
            if os.path.basename(pathFile) in filename and md5 == fileMD5:
                findMD5 = True
    return findMD5


'''
-------------------------------------------------------------------------------------------------------
                                           GENBANK DOWNLOAD
-------------------------------------------------------------------------------------------------------
'''


@fct_checker
def dl_genbank_bacteria(section: str, taxonomyID: int, pathOUT: str, chunkSize: int = 100) -> Tuple[str, int, str, int]:
    '''
     ------------------------------------------------------------
    |                   DOWNLOAD BACTERIA GBFF                   |
     ------------------------------------------------------------
    |        Download all bacteria Genbank .gbff.gz files        |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    section   : refseq or genbank section (required)        |
    |    taxonomyID: taxonomy id (required)                      | bacteria => 2 , vibrionaceae => 641
    |    pathOUT   : path of output files (required)             |
    |    chunkSize : download chunk size (default=100)           |
     ------------------------------------------------------------
    '''
    pathOUT = path_converter(pathOUT)
    if section not in ["refseq", "genbank"]:
        printcolor("[ERROR: dl_genbank_bacteria]\nSection must be \"refseq\" or \"genbank\"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT += "/"+section
    os.makedirs(pathOUT, exist_ok=True)
    urlSUMMARY = "https: //ftp.ncbi.nlm.nih.gov/genomes/"+section+"/assembly_summary_"+section+".txt"
    pathSummary = pathOUT+"/assembly_summary_"+section+".txt"
    pathSummaryEtag = pathOUT+"/assembly_summary_"+section+".etag"
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathTMP = geminiset.pathTMP
    pathLOG = pathTMP+"/wget.log"
    # ***** RETRIEVE SUMMARY ***** #
    # Get current local/distant etag (or None if empty)
    spinner = yaspin(Spinners.aesthetic, text="♊ Retrieve summary ETAG", side="right")
    spinner.start()
    title("Etag", None)
    previousEtag = get_etag(pathSummaryEtag, None)
    newEtag = get_etag(None, urlSUMMARY)
    spinner.stop()
    printcolor("♊ Retrieve summary ETAG"+"\n")
    if newEtag == "ERROR":
        printcolor("[ERROR: Etag request for \""+urlSUMMARY+"\"]\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Download if new distant etag
    if previousEtag == newEtag:
        printcolor("♊ "+section.capitalize()+" summary [up-to-date]"+"\n")
    else:
        printcolor("♊ "+section.capitalize()+" summary [downloading]"+"\n")
        # Get filesize for progress bar
        filesize = get_ftp_file_size("ftp.ncbi.nlm.nih.gov", urlSUMMARY.replace("https: //ftp.ncbi.nlm.nih.gov/", ""))
        # Download file
        get_distant_file_tqdm(urlSUMMARY, pathSummary, filesize, 1024)
    # Write new etag
    ETAG = open(pathSummaryEtag, 'w')
    ETAG.write(newEtag)
    ETAG.close()
    # ***** GENBANK SUMMARY PARSER ***** #
    printcolor("♊ Get available"+"\n")
    # Init NCBITaxa
    ncbi = NCBITaxa()
    # Count summary line
    p = subprocess.Popen("wc -l "+pathSummary, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    nbLines = int(output.decode().split(" ")[0])
    # Browse summary
    dicoSummary = {}
    lstMissing = []
    with open(pathSummary) as infile:
        pbar = tqdm(total=nbLines, ncols=50, leave=False, desc="", file=sys.stdout, bar_format="   Parsing {percentage: 3.0f}%|{bar}|")
        for line in infile:
            if line[0] != "#":
                splitLine = line[: -1].split("\t")
                gcaAccession = splitLine[0]
                taxID = int(splitLine[5])
                specieID = int(splitLine[6])
                genomeStatut = splitLine[13]
                ftpPath = splitLine[19]
                # Check taxonomy
                try:
                    lineage = ncbi.get_lineage(taxID)
                except ValueError:
                    try:
                        lineage = ncbi.get_lineage(specieID)
                    except ValueError:
                        continue
                if taxonomyID in lineage and genomeStatut != "Partial" and ("GCA" in gcaAccession or "GCF" in gcaAccession):
                    dicoSummary[gcaAccession] = ftpPath
                    pathGBK = pathOUT+"/"+gcaAccession+".gbk.gz"
                    pathGBKmd5 = pathOUT+"/"+gcaAccession+".md5"
                    md5Check = checkMD5(pathGBK, pathGBKmd5)
                    if not os.path.isfile(pathGBK) or os.path.getsize(pathGBK) == 0 or md5Check is False:
                        lstMissing.append(gcaAccession)
            pbar.update(1)
            title("Get available", pbar)
        pbar.close()
    # Check deprecated genomes
    nbDeprecated = 0
    for file in os.listdir(pathOUT):
        if ".gbk.gz" in file and not file.replace(".gbk.gz", "") in dicoSummary:
            nbDeprecated += 1
            os.remove(pathOUT+"/"+file)
    # Display statistics
    printcolor("⏩ Found "+str(len(dicoSummary))+" accession files"+"\n")
    printcolor("⏩ with "+str(len(lstMissing))+" missing"+"\n")
    printcolor("⏩ with "+str(nbDeprecated)+" deprecated"+"\n")
    # ***** DOWNLOAD MISSING FILES ***** #
    printcolor("♊ Download missing genomes"+"\n")
    lstDLerror = []
    # Chunks missing genomes
    chunkedLstMissing = [lstMissing[i: i+chunkSize] for i in range(0, len(lstMissing), chunkSize)]
    pbar = tqdm(total=len(lstMissing), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for subLstMissing in chunkedLstMissing:
        TMP = open(pathTMP+"/chunks.txt", 'w')
        for gcaAccession in subLstMissing:
            urlGBKFTP = dicoSummary[gcaAccession]+"/"+os.path.basename(dicoSummary[gcaAccession])+"_genomic.gbff.gz"
            TMP.write(urlGBKFTP+"\n")
        TMP.close()
        cmdDL = dicoGeminiPath['TOOLS']['python']+" "+dicoGeminiPath['TOOLS']['multidl']+" -c "+pathTMP+"/chunks.txt -o "+pathTMP+" -n "+str(cpu)+" > /dev/null 2>&1"
        os.system(cmdDL)
        for gcaAccession in subLstMissing:
            pathTmpGBK = pathTMP+"/"+os.path.basename(dicoSummary[gcaAccession])+"_genomic.gbff.gz"
            if not os.path.isfile(pathTmpGBK):
                lstDLerror.append(gcaAccession)
            elif os.path.getsize(pathTmpGBK) == 0:
                os.remove(pathTmpGBK)
                lstDLerror.append(gcaAccession)
            else:
                pathGBK = pathOUT+"/"+gcaAccession+".gbk.gz"
                shutil.move(pathTmpGBK, pathGBK)
                pathMD5 = pathOUT+"/"+gcaAccession+".md5"
                urlMD5FTP = dicoSummary[gcaAccession]+"/md5checksums.txt"
                cmdCurlMD5 = "curl -s "+urlMD5FTP+" | grep \"_genomic.gbff.gz\" | sed -E s/\"\\..+_genomic.gbff.gz\"/\""+gcaAccession+".gbk.gz\"/ > "+pathMD5
                os.system(cmdCurlMD5)
        pbar.update(chunkSize)
    pbar.close()
    if len(lstDLerror):
        printcolor("⛔ Download error for "+",".join(lstDLerror)+"\n")
    pbar = tqdm(total=len(lstMissing[: 2]), ncols=50, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for gcaAccession in lstMissing[: 2]:
        pbar.set_description_str(gcaAccession)
        pathGBK = pathOUT+"/"+gcaAccession+".gbk.gz"
        urlGBKFTP = dicoSummary[gcaAccession]+"/"+os.path.basename(dicoSummary[gcaAccession])+"_genomic.gbff.gz"
        cmdWGET = "wget --header = \"accept-encoding: gzip\" -O "+pathGBK+" "+urlGBKFTP+" >> "+pathLOG+" 2>&1"
        pathMD5 = pathOUT+"/"+gcaAccession+".md5"
        urlMD5FTP = dicoSummary[gcaAccession]+"/md5checksums.txt"
        cmdCurlMD5 = "curl -s "+urlMD5FTP+" | grep \"_genomic.gbff.gz\" | sed -E s/\"\\..+_genomic.gbff.gz\"/\""+gcaAccession+".gbk.gz\"/ > "+pathMD5
        os.system(cmdWGET)
        os.system(cmdCurlMD5)
        pbar.update(1)
        title("Download", pbar)
    pbar.close()

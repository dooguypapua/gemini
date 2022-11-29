'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ plot
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminiphage.py                  |
description      : gemini phage functions          |
author           : dooguypapua                     |
lastmodification : 20210720                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import sys
import re
import os
import shutil
from tqdm import tqdm
from typing import Tuple
from yaspin import yaspin
from yaspin.spinners import Spinners
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from svgpathtools import svg2paths
from cairosvg import svg2png
from openpyxl import Workbook, load_workbook, styles
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Alignment
import geminiset
from geminini import fct_checker, get_input_files, printcolor, path_converter, load_json, get_gemini_path
from geminini import random_hex_color, linear_gradient, title, exit_gemini, cat_lstfiles, to_ranges
from geminicluster import mmseqs_rbh
from geminiparse import make_gff_dict, make_gbk_dict, gbk_to_faa, make_fasta_dict


@fct_checker
def gff_to_linear_geneplot(pathIN: str, pathOUT: str, pathLT: str = "None", length: int = -1, ext: str = ".gff") -> Tuple[str, str, str, int, str]:
    '''
     ------------------------------------------------------------
    |                  GFF3 TO LINEAR GENE PLOT                  |
    |------------------------------------------------------------|
    |           Create linear gene plot from GFF3 file           |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input files or folder (required)       |
    |    pathOUT: path of output files (required)                |
    |    pathLT : path of input LT list file (default=None)      |
    |    length : representation length (default=-1 use sequence)|
    |    ext    : extension of input files (default=.gff)        |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "gff_to_linear_geneplot", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: gff_to_linear_geneplot]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    # Retrieve LT list if required
    if pathLT != "None":
        setLT = set()
        pathLT = path_converter(pathLT)
        LTFILE = open(pathLT, 'r')
        lstLines = LTFILE.read().split("\n")
        LTFILE.close()
        for line in lstLines:
            if line != "" and not line[0] == "#":
                setLT.add(line)
    # Browse GFF3
    printcolor("♊ Plotting"+"\n")
    pbar = tqdm(total=int(len(lstFiles)), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathGFF in lstFiles:
        file = os.path.basename(pathGFF)
        orgName = file.replace(ext, "").replace("."+ext, "")
        pathPNG = pathOUT+"/"+orgName+".png"
        pathSVG = pathOUT+"/"+orgName+".svg"
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        features = []
        # ***** PARSE GFF ***** #
        dicoGFF = make_gff_dict(pathIN=pathGFF, ext=ext)
        startRegion = 0
        endRegion = 0
        # ***** BROWSE GENES ***** #
        for geneType in dicoGFF[orgName]:
            if geneType != 'length':
                for geneEntry in dicoGFF[orgName][geneType]:
                    if geneType == "CDS":
                        color = "#2a7fff"
                    elif geneType == "tRNA":
                        color = "#37c8ab"
                    else:
                        continue
                    if 'locus_tag' in geneEntry['attributes']:
                        geneFeature = GraphicFeature(start=geneEntry['start'], end=geneEntry['end'], strand=int(geneEntry['strand']+"1"), color=color, linewidth=0)
                    else:
                        geneFeature = GraphicFeature(start=geneEntry['start'], end=geneEntry['end'], strand=int(geneEntry['strand']+"1"), color=color, linewidth=0)
                    if pathLT == "None":
                        features.append(geneFeature)
                    elif geneEntry['attributes']['locus_tag'] in setLT:
                        if startRegion == 0:
                            startRegion = geneEntry['start']
                        endRegion = geneEntry['end']
                    features.append(geneFeature)
        # ***** PLOT GENES ***** #
        if length != -1:
            seqLen = length
        elif pathLT == "None":
            seqLen = dicoGFF[orgName]['length']
        else:
            seqLen = endRegion-startRegion
        record = GraphicRecord(sequence_length=seqLen+int(seqLen/10), features=features, first_index=startRegion-100)
        ax, _ = record.plot(figure_width=50)
        ax.figure.savefig(pathPNG, dpi=300)
        ax.figure.savefig(pathSVG)
        plt.close('all')
        pbar.update(1)
        title("Plotting", pbar)
    pbar.close()


@fct_checker
def gff_to_linear_group_geneplot(pathIN: str, pathCLUSTER: str, pathGROUP: str, pathOUT: str, ext: str = ".gff") -> Tuple[str, str, str, str, str]:
    '''
     ------------------------------------------------------------
    |             GFF3 TO LINEAR GENE PLOT PER GROUP             |
    |------------------------------------------------------------|
    | Create linear gene plot from GFF3 file and merge per group |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN     : path of input files or folder (required)   |
    |    pathOUT    : path of output files (required)            |
    |    pathCLUSTER: path of JSON rbh cluster (required)        |
    |    pathGROUP  : path of input tab group file (required)    |
    |    ext        : extension of input files (default=.gff)    |
     ------------------------------------------------------------
     input tab group: <name>\t<group>\t<subgroup color>
     singleton => black
     core => grey
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "gff_to_linear_group_geneplot", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: gff_to_linear_group_geneplot]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathCLUSTER = path_converter(pathCLUSTER)
    pathOUT = path_converter(pathOUT)
    pathGROUP = path_converter(pathGROUP)
    os.makedirs(pathOUT, exist_ok=True)
    # ***** PARSE tabular group file ***** #
    printcolor("♊ Read groups"+"\n")
    dicoOrgToGroup = {}
    dicoGroupToOrg = {}
    dicoOrgToSubGroupColor = {}
    GROUP = open(pathGROUP, 'r')
    lstLines = GROUP.read().split("\n")
    GROUP.close()
    for line in lstLines:
        if line != "":
            orgName = line.split("\t")[0]
            groupName = line.split("\t")[1]
            subGroupColor = line.split("\t")[2]
            dicoOrgToGroup[orgName] = groupName
            dicoOrgToSubGroupColor[orgName] = subGroupColor
            try:
                dicoGroupToOrg[groupName].append(orgName)
            except KeyError:
                dicoGroupToOrg[groupName] = [orgName]
    # ***** CLUSTERING genes ***** #
    printcolor("♊ Read clustering"+"\n")
    dicoCLUSTER = load_json(pathCLUSTER)
    printcolor("♊ Group clustering"+"\n")
    dicoGFF = make_gff_dict(pathIN=pathIN, ext=ext)
    dicoLTtoOrg = {}
    dicoGroupToOrgLT = {}
    dicoLTcolor = {}
    seqLen = 0
    for orgName in dicoOrgToGroup:
        for geneType in dicoGFF[orgName]:
            if geneType == 'length':
                seqLen = max(seqLen, dicoGFF[orgName][geneType])
            if geneType != 'length':
                for geneEntry in dicoGFF[orgName][geneType]:
                    if not geneType[0] == "#":
                        lt = geneEntry['attributes']['locus_tag']
                        dicoGroupToOrgLT[lt] = dicoOrgToGroup[orgName]
                        dicoLTtoOrg[lt] = orgName
    for group in dicoGroupToOrg:
        for cluster in dicoCLUSTER:
            setOrg = set()  # To avoid paralogous
            setOrgSubGroupColor = set()
            for header in dicoCLUSTER[cluster]:
                lt = header.split(" [")[0].split("|")[0]
                if lt in dicoGroupToOrgLT and dicoGroupToOrgLT[lt] == group:
                    setOrg.add(dicoLTtoOrg[lt])
                    setOrgSubGroupColor.add(dicoOrgToSubGroupColor[dicoLTtoOrg[lt]])
            # Singleton gene in black
            if len(setOrg) == 1:
                color = "#1a1a1a"
            # If group core gene => grey
            elif len(setOrg) == len(dicoGroupToOrg[group]):
                color = "#cccccc"
            # Look for core subgroup
            elif len(setOrgSubGroupColor) == 1 and list(setOrgSubGroupColor)[0].lower() != "none":
                color = list(setOrgSubGroupColor)[0]
            # Else random color
            else:
                color = random_hex_color(False)
            # Apply to all LT
            for header in dicoCLUSTER[cluster]:
                lt = header.split(" [")[0].split("|")[0]
                dicoLTcolor[lt] = color
    # ***** BROWSE input GFF files ***** #
    dicoGroupToOrgPlot = {}
    printcolor("♊ GraphicRecords"+"\n")
    pbar = tqdm(total=int(len(lstFiles)), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathGFF in lstFiles:
        file = os.path.basename(pathGFF)
        orgName = file.replace(ext, "").replace("."+ext, "")
        if orgName in dicoOrgToGroup:
            pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
            features = []
            # ***** PARSE GFF ***** #
            dicoGFF = make_gff_dict(pathIN=pathGFF, ext=ext)
            # ***** BROWSE GENES ***** #
            for geneType in dicoGFF[orgName]:
                if geneType != 'length':
                    for geneEntry in dicoGFF[orgName][geneType]:
                        if geneType in ["CDS", "tRNA"]:
                            color = dicoLTcolor[geneEntry['attributes']['locus_tag']]
                        else:
                            continue
                        geneFeature = GraphicFeature(start=geneEntry['start'], end=geneEntry['end'],
                                                     strand=int(geneEntry['strand']+"1"), color=color, linewidth=0)
                        features.append(geneFeature)
            # ***** PLOT GENES ***** #
            record = GraphicRecord(sequence_length=seqLen, features=features)
            try:
                dicoGroupToOrgPlot[dicoOrgToGroup[orgName]].append((orgName, record))
            except KeyError:
                dicoGroupToOrgPlot[dicoOrgToGroup[orgName]] = [(orgName, record)]
        pbar.update(1)
        title("GraphicRecords", pbar)
    pbar.close()
    # ***** Unify group plot ***** #
    printcolor("♊ Unify plots"+"\n")
    pbar = tqdm(total=int(len(dicoGroupToOrgPlot)), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    dicoGroupPNG = {}
    setSVGfiles = set()
    for group in dicoGroupToOrg:
        if len(dicoGroupToOrg[group]) >= 2:
            pathOUTgroup = pathOUT+"/"+group
            os.makedirs(pathOUTgroup, exist_ok=True)
            pbar.set_description_str("Group = "+group+" ".rjust(maxpathSize-len("Group = "+group)))
            dicoGroupPNG[group] = []
            # Get longest sequence for current group and enlarge to 10%
            maxLen = 0
            for (orgName, record) in dicoGroupToOrgPlot[group]:
                maxLen = max(record.sequence_length, maxLen)
            maxLen = maxLen+int(maxLen/10)
            # Get GraphicFeature
            for (orgName, record) in dicoGroupToOrgPlot[group]:
                pathPNG = pathOUTgroup+"/"+orgName+".png"
                pathSVG = pathOUTgroup+"/"+orgName+".svg"
                ax, _ = record.plot(figure_width=50)
                ax.figure.savefig(pathPNG, dpi=300)
                ax.figure.savefig(pathSVG, backend="SVG")
                setSVGfiles.add(pathSVG)
                plt.close('all')
        pbar.update(1)
        title("Unify", pbar)
    pbar.close()
    # ***** Reformat SVG output  ***** # (selectable color)
    printcolor("♊ Reformat SVG"+"\n")
    pbar = tqdm(total=len(setSVGfiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathSVG in setSVGfiles:
        pbar.set_description_str(os.path.basename(pathSVG)+" ".rjust(maxpathSize-len(os.path.basename(pathSVG))))
        # Read initial SVG
        SVG = open(pathSVG, 'r')
        dataSVG = SVG.read()
        splitG = dataSVG.split("<g")
        SVG.close()
        dicoPatchToColor = {}
        for partG in splitG:
            if "id=\"patch_" in partG:
                splitLine = partG.split("\n")
                patchID = ""
                color = ""
                for line in splitLine:
                    if "id=\"patch_" in line:
                        patchID = line.split("\"")[1]
                    if "style=\"fill:" in line:
                        color = line.split(":")[1].split(";")[0].replace("\"", "")
                dicoPatchToColor[patchID] = color
        # Apply modifications
        SVG = open(pathSVG, 'w')
        for patchID in dicoPatchToColor:
            dataSVG = dataSVG.replace("<g id=\""+patchID+"\">", "<g id=\""+patchID+"\"\nstyle=\"fill:"+dicoPatchToColor[patchID]+"\">")
            dataSVG = dataSVG.replace(dicoPatchToColor[patchID], dicoPatchToColor[patchID]+";fill-opacity:1")
        SVG.write(dataSVG)
        SVG.close()
        pbar.update(1)
        title("Reformat", pbar)
    pbar.close()


@fct_checker
def svg_dna_transform(pathIN: str, pathOUT: str) -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |             TRANSFORM dna_features_viewer SVG              |
    |------------------------------------------------------------|
    |         Transform SVG plot from dna_features_viewer        |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN     : path of input/output SVG file              |
    |    pathOUT    : path of output PNG file (required)         |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    # ***** DELETE background rectangle and middle line ***** #
    SVG = open(pathIN, 'r')
    lstLines = SVG.read().split("\n")
    SVG.close()
    NEWSVG = open(pathIN, 'w')
    findBg = False
    findHline = False
    for i in range(len(lstLines)):
        if "id = \"patch_1\"" in lstLines[i]:
            findBg = True
        if findBg and "</g>" in lstLines[i]:
            findBg = False
            continue
        if "id = \"line2d_1\"" in lstLines[i]:
            findHline = True
        if findHline and "</g>" in lstLines[i]:
            findHline = False
            continue
        if findBg is False and findHline is False:
            if "<g" in lstLines[i] and "<g" in lstLines[i-1]:
                continue
            NEWSVG.write(lstLines[i]+"\n")
    NEWSVG.close()
    # Parse SVG
    paths, attributes = svg2paths(pathIN)
    # Search max y pos
    maxAllY = 0
    for attr in attributes:
        if "style" in attr:
            searchallY = re.findall(r"([lL]) (\S+) (\S+)\s+Q", attr['d'])
            for searchY in searchallY:
                maxAllY = max(maxAllY, float(searchY[2]))
    # Check required y position transformation
    dicoReplace = {}
    for attr in attributes:
        if "style" in attr:
            maxY = 0
            searchallY = re.findall(r"([lL]) (\S+) (\S+)\s+Q", attr['d'])
            for searchY in searchallY:
                maxY = max(maxY, float(searchY[2]))
            searchDLine = re.search(r"([mM]) (\S+) (\S+)\s+Q", attr['d'])
            if maxY < maxAllY:
                dicoReplace[searchDLine.group(1)+" "+searchDLine.group(2)+" "+searchDLine.group(3)] = maxAllY-float(maxY)
            else:
                dicoReplace[searchDLine.group(1)+" "+searchDLine.group(2)+" "+searchDLine.group(3)] = None
    # Apply y position transformation
    SVG = open(pathIN, 'r')
    lstLines = SVG.read().split("\n")
    SVG.close()
    lstToWrite = []
    for i in range(len(lstLines)):
        for key in dicoReplace:
            if key in lstLines[i]:
                lstToWrite[-1] = lstToWrite[-1].replace(">", "")
                if dicoReplace[key] is None:
                    lstToWrite.append("     style = \"stroke: #000000;stroke-opacity: 1\">")
                else:
                    lstToWrite.append("     style = \"stroke: #000000;stroke-opacity: 1\"")
                    lstToWrite.append("     transform = \"translate(0, "+str(dicoReplace[key])+")\">")
                break
        lstToWrite.append(lstLines[i])
    # Convert to PNG
    svg2png(bytestring="\n".join(lstToWrite), write_to=pathOUT, dpi=300)
    os.system("convert "+pathOUT+" -gravity South -background white -splice 0x1 -background black -splice 0x1 -trim +repage -chop 0x1" +
              " -gravity North -background white -splice 0x1 -background black -splice 0x1 -trim +repage -chop 0x1" +
              " -gravity East -background white -splice 5x0 -background black -splice 5x0 -trim +repage -chop 5x " +
              pathOUT.replace(".png", ".bak"))
    shutil.move(pathOUT.replace(".png", ".bak"), pathOUT)


@fct_checker
def xlsx_to_heatmap(pathIN: str, pathOUT: str, colorStart: str = "FFFFFF", colorEnd: str = "FF0000", headerRow: int = -1, headerCol: int = -1) -> Tuple[str, str, str, str, int, int]:
    '''
     ------------------------------------------------------------
    |                      XLSX TO HEATMAP                       |
    |------------------------------------------------------------|
    |            Convert xlsx table to heatmap table             |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN    : path of input files or folder (required)    |
    |    pathOUT   : path of output files (required)             |
    |    colorStart: hex start color (default=FFFFFF)            |
    |    colorEnd  : hex end color (default=FF0000)              |
    |    headerRow : row header line (default=-1 no header)      |
    |    headerCol : column header line (default=-1 no header)   |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    printcolor("♊ XLSX to heatmap"+"\n")
    # Create color gradient list
    HEX_list, RBG_list = linear_gradient("#"+colorStart, "#"+colorEnd, 101)
    # Open input/output XLSX
    spinner = yaspin(Spinners.aesthetic, text="♊ Loading input XLSX", side="right")
    spinner.start()
    title("Loading", None)
    wbIN = load_workbook(pathIN)
    wsIN = wbIN.active
    spinner.stop()
    printcolor("♊ Load XLSX"+"\n")
    wbOUT = Workbook(write_only=True)
    wsOUT = wbOUT.create_sheet("heatmap")
    # Browse input XLSX
    printcolor("♊ Heatmap creation"+"\n")
    cptRow = 0
    pbar = tqdm(total=wsIN.max_row, ncols=50, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for rowIN in wsIN.rows:
        rowOUT = []
        cptCol = 0
        for cellIN in rowIN:
            # For header copy value
            if cptRow == headerRow or cptCol == headerCol:
                cellOUT = WriteOnlyCell(wsOUT, value=cellIN.value)
            # For data apply color without value
            else:
                cellOUT = WriteOnlyCell(wsOUT, value="")
                value = int(round(float(str(cellIN.value).replace(", ", ".")), 0))
                cellColor = styles.colors.Color(rgb=HEX_list[value].replace("#", ""))
                if cellColor != "#FFFFFF":
                    cellFill = styles.fills.PatternFill(patternType='solid', fgColor=cellColor)
                    cellOUT.fill = cellFill
            cellOUT.alignment = Alignment(horizontal='center')
            rowOUT.append(cellOUT)
            cptCol += 1
        wsOUT.append(rowOUT)
        cptRow += 1
        pbar.update(1)
        title("Heatmap", pbar)
    pbar.close()
    # Close the workbook after reading
    wbIN.close()
    wbOUT.save(pathOUT)


@fct_checker
def circos_plot(pathIN1: str, pathOUT: str, pathIN2: str = "None", pident: int = 30, cov: int = 80) -> Tuple[str, str, str, int, int]:
    '''
     ------------------------------------------------------------
    |                        CIRCOS PLOT                         |
    |------------------------------------------------------------|
    |           Circular circos plot from genbank file           |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN1 : path of reference genbank file (required)     |
    |    pathOUT : path of output file (required)                |
    |    pathIN2 : path of query FAA file/folder (default=None)  |
    |    pident  : min %id for orthologuous (default=30)         |
    |    cov     : min %cov for orthologuous (default=80)        |
     ------------------------------------------------------------
    '''
    pathIN1 = path_converter(pathIN1)
    pathOUT = path_converter(pathOUT)
    if pathIN2 != "None":
        lstQueryFAAFiles, maxpathSize = get_input_files(pathIN2, "circos_plot", [".faa"])

    # ***** READ reference GBK ***** #
    printcolor("♊ Parse GBK"+"\n")
    dicoGBK = list(make_gbk_dict(pathIN1).values())[0]

    # ***** BLASTP query FAA ***** #
    if pathIN2 != "None":
        # Init dicoRBH
        dicoRBH = {}
        for contig in dicoGBK:
            for lt in dicoGBK[contig]['dicoLT']:
                dicoRBH[lt] = {}
                for FAAFile in lstQueryFAAFiles:
                    dicoRBH[lt][os.path.basename(FAAFile).replace(".faa", "")] = 0
        # Make FAA and create symlinks
        gbk_to_faa(pathIN=pathIN1, pathOUT=geminiset.pathTMP+"/ref.faa")
        for FAAFile in lstQueryFAAFiles:
            orgQuery = os.path.basename(FAAFile).replace(".faa", "")
            os.makedirs(geminiset.pathTMP+"/mmseqs")
            os.symlink(FAAFile, geminiset.pathTMP+"/mmseqs/"+os.path.basename(FAAFile))
            os.symlink(geminiset.pathTMP+"/ref.faa", geminiset.pathTMP+"/mmseqs/ref.faa")
            # Launch mmseqs_rbh
            mmseqs_rbh(pathIN=geminiset.pathTMP+"/mmseqs", pathOUT=geminiset.pathTMP+"/mmseqs", ref="ref", idThrClust=pident, covThrClust=cov, boolNucl=False, ext=".faa")
            shutil.copyfile(geminiset.pathTMP+"/mmseqs/ref.rbh", geminiset.pathTMP+"/"+orgQuery+".rbh")
            os.system("rm -rf "+geminiset.pathTMP+"/mmseqs*")
            # Parse mmseqs_rbh
            TSV = open(geminiset.pathTMP+"/"+orgQuery+".rbh", 'r')
            lstLines = TSV.read().split("\n")[:-1]
            TSV.close()
            for line in lstLines:
                splitLine = line.split("\t")
                ltRef = splitLine[0].replace("#", "|").split("|")[0]
                dicoRBH[ltRef][orgQuery] = float(splitLine[2])

    # ***** KARYOTYPE ***** #
    KARYOTYPE = open(geminiset.pathTMP+"/karyotype.txt", 'w')
    lstColorBand = ["black", "white", "grey"]
    start = 0
    bandIndex = 0
    totalSize = 0
    LstTowrite = []
    for contig in dicoGBK:
        ctgLen = len(dicoGBK[contig]['seq'])
        totalSize += ctgLen
        LstTowrite.append("band chr1 1.1 1.1 "+str(start)+" "+str(start+ctgLen)+" "+lstColorBand[bandIndex])
        bandIndex += 1
        if bandIndex == 3:
            bandIndex = 0
        start = start+ctgLen
    # Add total size at first karyotype line
    LstTowrite.insert(0, "chr - chr1 1 0 "+str(totalSize)+" black")
    KARYOTYPE.write("\n".join(LstTowrite))
    KARYOTYPE.close()
    # Define ticks correspunding to total size
    if totalSize >= 4000000:
        ticksMultiplier = 0.000001
        majorTicksSpacing = "1000000u"
        minorTicksSpacing = "100000u"
        majorTicksSuffix = "Mb"
    elif totalSize >= 100000:
        ticksMultiplier = 0.00001
        majorTicksSpacing = "100000u"
        minorTicksSpacing = "10000u"
        majorTicksSuffix = "Kb"

    # ***** IDEOGRAM ***** #
    IDEOGRAM = open(geminiset.pathTMP+"/ideogram.conf", 'w')
    IDEOGRAM.write("<ideogram>\n")
    IDEOGRAM.write("<spacing>\n")
    IDEOGRAM.write("default = 0u\n")
    IDEOGRAM.write("break   = 0u\n")
    IDEOGRAM.write("</spacing>\n")
    IDEOGRAM.write("thickness             = 20p\n")  # thickness (px) of chromosome ideogram
    IDEOGRAM.write("stroke_thickness      = 2\n")
    IDEOGRAM.write("stroke_color          = black\n")  # ideogram border color
    IDEOGRAM.write("fill                  = yes\n")
    IDEOGRAM.write("fill_color            = black\n")  # the default chromosome color is set here and any value, defined in the karyotype file overrides it
    IDEOGRAM.write("radius                = 0.85r\n")  # fractional radius position of chromosome ideogram within image
    IDEOGRAM.write("show_label            = no\n")
    IDEOGRAM.write("label_font            = default\n")
    IDEOGRAM.write("label_radius          = dims(ideogram,radius) + 0.05r\n")
    IDEOGRAM.write("label_size            = 36\n")
    IDEOGRAM.write("label_parallel        = yes\n")
    IDEOGRAM.write("label_case            = upper\n")
    IDEOGRAM.write("band_stroke_thickness = 2\n")  # cytogenetic bands
    IDEOGRAM.write("show_bands            = yes\n")  # show_bands determines whether the outline of cytogenetic bands will be seen
    IDEOGRAM.write("fill_bands            = yes\n")  # in order to fill the bands with the color defined in the karyotype file you must set fill_bands
    IDEOGRAM.write("</ideogram>\n")
    IDEOGRAM.close()

    # ***** TICKS ***** #
    TICKS = open(geminiset.pathTMP+"/ticks.conf", 'w')
    TICKS.write("show_ticks          = yes\n")
    TICKS.write("show_tick_labels    = yes\n")
    TICKS.write("show_grid          = no\n")
    TICKS.write("grid_start         = dims(ideogram,radius_inner)-0.5r\n")
    TICKS.write("grid_end           = dims(ideogram,radius_inner)\n")
    TICKS.write("<ticks>\n")
    TICKS.write("skip_first_label     = yes\n")
    TICKS.write("skip_last_label      = no\n")
    TICKS.write("radius               = dims(ideogram,radius_outer)\n")
    TICKS.write("tick_separation      = 2p\n")
    TICKS.write("min_label_distance_to_edge = 0p\n")
    TICKS.write("label_separation = 5p\n")
    TICKS.write("label_offset     = 25p\n")
    TICKS.write("label_size = 8p\n")
    TICKS.write("multiplier = "+str(ticksMultiplier)+"\n")
    TICKS.write("color = black\n")
    # Major ticks
    TICKS.write("<tick>\n")
    TICKS.write("spacing        = "+majorTicksSpacing+"\n")
    TICKS.write("color          = black\n")
    TICKS.write("show_label     = yes\n")
    TICKS.write("suffix = \" "+majorTicksSuffix+"\"\n")
    TICKS.write("label_size     = 36p\n")
    TICKS.write("format         = %s\n")
    TICKS.write("grid           = yes\n")
    TICKS.write("grid_color     = dgrey\n")
    TICKS.write("grid_thickness = 1p\n")
    TICKS.write("thickness = 5p\n")
    TICKS.write("size      = 30p\n")
    TICKS.write("</tick>\n")
    # Minor ticks
    TICKS.write("<tick>\n")
    TICKS.write("spacing        = "+minorTicksSpacing+"\n")
    TICKS.write("color          = grey\n")
    TICKS.write("show_label     = no\n")
    TICKS.write("format         = %s\n")
    TICKS.write("grid           = yes\n")
    TICKS.write("grid_color     = dgrey\n")
    TICKS.write("grid_thickness = 1p\n")
    TICKS.write("thickness = 5p\n")
    TICKS.write("size      = 20p\n")
    TICKS.write("</tick>\n")
    TICKS.write("</ticks>\n")
    TICKS.close()

    # ***** COLORS ***** #
    COLORS = open(geminiset.pathTMP+"/colors.conf", 'w')
    COLORS.write("generev = 77,77,77\n")
    COLORS.write("genefor = 179,179,179\n")
    COLORS.write("generna = 0,0,0\n")
    # Gradient
    HEXlist, RBGlist = linear_gradient("#db073d", "#ffe9e7", 100)
    cpt = 1
    strLstGradient = ""
    for rgb in RBGlist[::-1]:
        COLORS.write("gradient"+str(cpt)+" = "+str(rgb[0])+","+str(rgb[1])+","+str(rgb[2])+"\n")
        strLstGradient += "gradient"+str(cpt)+","
        cpt += 1
    COLORS.close()
    strLstGradient = strLstGradient[:-1]

    # ***** HIGHLIGHTS ***** #
    HIGHLIGHTSFOR = open(geminiset.pathTMP+"/highlightsFOR.txt", 'w')
    totalSize = 0
    for contig in dicoGBK:
        ctgLen = len(dicoGBK[contig]['seq'])
        for lt in dicoGBK[contig]['dicoLT']:
            if dicoGBK[contig]['dicoLT'][lt]['type'] == "CDS" and dicoGBK[contig]['dicoLT'][lt]['strand'] == 1:
                HIGHLIGHTSFOR.write("chr1 "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['start'])+" "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['end'])+" fill_color=genefor\n")
        totalSize += ctgLen
    HIGHLIGHTSFOR.close()
    HIGHLIGHTSREV = open(geminiset.pathTMP+"/highlightsREV.txt", 'w')
    totalSize = 0
    for contig in dicoGBK:
        ctgLen = len(dicoGBK[contig]['seq'])
        for lt in dicoGBK[contig]['dicoLT']:
            if dicoGBK[contig]['dicoLT'][lt]['type'] == "CDS" and dicoGBK[contig]['dicoLT'][lt]['strand'] == -1:
                HIGHLIGHTSREV.write("chr1 "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['start'])+" "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['end'])+" fill_color=generev\n")
        totalSize += ctgLen
    HIGHLIGHTSREV.close()
    HIGHLIGHTSRNA = open(geminiset.pathTMP+"/highlightsRNA.txt", 'w')
    totalSize = 0
    for contig in dicoGBK:
        ctgLen = len(dicoGBK[contig]['seq'])
        for lt in dicoGBK[contig]['dicoLT']:
            if "RNA" in dicoGBK[contig]['dicoLT'][lt]['type']:
                HIGHLIGHTSRNA.write("chr1 "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['start'])+" "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['end'])+" fill_color=generna\n")
        totalSize += ctgLen
    HIGHLIGHTSRNA.close()

    # ***** ORTHOLOGOUS HEATMAP ***** #
    totalSize = 0
    if pathIN2 != "None":
        for contig in dicoGBK:
            ctgLen = len(dicoGBK[contig]['seq'])
            for lt in dicoGBK[contig]['dicoLT']:
                for orgQuery in dicoRBH[lt]:
                    ORTHOORG = open(geminiset.pathTMP+"/"+orgQuery+"_heatmap.txt", 'a')
                    ORTHOORG.write("chr1 "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['start'])+" "+str(totalSize+dicoGBK[contig]['dicoLT'][lt]['end'])+" "+str(dicoRBH[lt][orgQuery])+"\n")
                    ORTHOORG.close()
            totalSize += ctgLen

    # ***** CIRCOSCONF ***** #
    CIRCOSCONF = open(geminiset.pathTMP+"/circos.conf", 'w')
    CIRCOSCONF.write("karyotype = "+geminiset.pathTMP+"/karyotype.txt\n")
    CIRCOSCONF.write("<<include "+geminiset.pathTMP+"/ideogram.conf>>\n")
    CIRCOSCONF.write("<<include "+geminiset.pathTMP+"/ticks.conf>>\n")
    CIRCOSCONF.write("<image>\n")
    CIRCOSCONF.write("dir = "+os.path.realpath(os.path.dirname(pathOUT))+"\n")
    CIRCOSCONF.write("file = "+os.path.basename(pathOUT)+"\n")
    CIRCOSCONF.write("png = yes\n")
    CIRCOSCONF.write("svg = no\n")
    CIRCOSCONF.write("radius = 1500p\n")
    CIRCOSCONF.write("angle_offset = -90\n")
    CIRCOSCONF.write("auto_alpha_colors = yes\n")
    CIRCOSCONF.write("auto_alpha_steps  = 5\n")
    CIRCOSCONF.write("</image>\n")
    CIRCOSCONF.write("<colors>\n")
    CIRCOSCONF.write("<<include "+geminiset.pathTMP+"/colors.conf>>\n")
    CIRCOSCONF.write("</colors>\n")
    CIRCOSCONF.write("<highlights>\n")
    CIRCOSCONF.write("<highlight>\n")
    CIRCOSCONF.write("file = "+geminiset.pathTMP+"/highlightsFOR.txt\n")
    CIRCOSCONF.write("r0   = 0.94r\n")
    CIRCOSCONF.write("r1   = 0.99r\n")
    CIRCOSCONF.write("</highlight>\n")
    CIRCOSCONF.write("<highlight>\n")
    CIRCOSCONF.write("file = "+geminiset.pathTMP+"/highlightsREV.txt\n")
    CIRCOSCONF.write("r0   = 0.89r\n")
    CIRCOSCONF.write("r1   = 0.94r\n")
    CIRCOSCONF.write("</highlight>\n")
    CIRCOSCONF.write("<highlight>\n")
    CIRCOSCONF.write("file = "+geminiset.pathTMP+"/highlightsRNA.txt\n")
    CIRCOSCONF.write("r0   = 0.84r\n")
    CIRCOSCONF.write("r1   = 0.89r\n")
    CIRCOSCONF.write("</highlight>\n")
    CIRCOSCONF.write("</highlights>\n")
    if pathIN2 != "None":
        CIRCOSCONF.write("<plots>\n")
        rStart = 0.83
        for file in os.listdir(geminiset.pathTMP):
            if "_heatmap.txt" in file:
                CIRCOSCONF.write("<plot>\n")
                CIRCOSCONF.write("type    = heatmap\n")
                CIRCOSCONF.write("file    = "+geminiset.pathTMP+"/"+file+"\n")
                CIRCOSCONF.write("color   = "+strLstGradient+"\n")
                CIRCOSCONF.write("r0      = "+str(round(rStart-0.02, 2))+"r\n")
                CIRCOSCONF.write("r1      = "+str(rStart)+"r\n")
                CIRCOSCONF.write("min     = "+str(pident)+"\n")
                CIRCOSCONF.write("max     = 100\n")
                CIRCOSCONF.write("</plot>\n")
                rStart = round(rStart-0.03, 2)
        CIRCOSCONF.write("</plots>\n")
    CIRCOSCONF.write("<<include /etc/circos/colors_fonts_patterns.conf>>\n")
    CIRCOSCONF.write("<<include /etc/circos/housekeeping.conf>>\n")
    CIRCOSCONF.close()

    # ***** MAKE PLOT *****#
    print("circos -conf "+geminiset.pathTMP+"/circos.conf")
    os.system("circos -conf "+geminiset.pathTMP+"/circos.conf")

    exit()


@fct_checker
def circos_align(pathIN: str, pathOUT: str) -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |                   CIRCOS ALIGNMENT PLOT                    |
    |------------------------------------------------------------|
    |             Circular circos plot for alignment             |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of reference genbank file (required)      |
    |    pathOUT : path of output file (required)                |
     ------------------------------------------------------------
    '''
    dicoGeminiPath = get_gemini_path()
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "circos_align", [".fna", ".fasta"])

    # ***** CAT files ***** #
    cat_lstfiles(lstFiles, geminiset.pathTMP+"/catseq.fasta")

    # ***** LAUNCH MAFFT ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ MAFFT", side="right")
    spinner.start()
    title("MAFFT", None)
    cmdMAFFT = dicoGeminiPath['mafft']+" --quiet --thread "+str(geminiset.cpu)+" "+geminiset.pathTMP+"/catseq.fasta > "+geminiset.pathTMP+"/mafft.fasta"
    os.system(cmdMAFFT)
    spinner.stop()
    printcolor("♊ MAFFT"+"\n")

    # ***** PARSE MAFFT ***** #
    dicoFASTA = make_fasta_dict(geminiset.pathTMP+"/mafft.fasta")
    alignLen = len(list(dicoFASTA.values())[0])
    dicoMAFFT = {}
    for align in dicoFASTA:
        for i in range(0, alignLen, 1):
            if dicoFASTA[align][i] != "N":
                try:
                    dicoMAFFT[i].add(dicoFASTA[align][i])
                except KeyError:
                    dicoMAFFT[i] = set([dicoFASTA[align][i]])

    # ***** KARYOTYPE ***** #
    KARYOTYPE = open(geminiset.pathTMP+"/karyotype.txt", 'w')
    KARYOTYPE.write("chr - chr1 1 0 "+str(alignLen)+" black\n")
    KARYOTYPE.write("chr - chr2 1 "+str(alignLen+1)+" "+str(alignLen+2)+" white\n")
    KARYOTYPE.close()
    # Define ticks correspunding to total size
    if alignLen >= 4000000:
        ticksMultiplier = 0.000001
        majorTicksSpacing = "1000000u"
        minorTicksSpacing = "100000u"
        majorTicksSuffix = "Mb"
    elif alignLen >= 100000:
        ticksMultiplier = 0.00001
        majorTicksSpacing = "100000u"
        minorTicksSpacing = "10000u"
        majorTicksSuffix = "Kb"

    # ****** COLORS *****#
    COLORS = open(geminiset.pathTMP+"/colors.conf", 'w')
    COLORS.write("cons = 175,221,233\n")
    COLORS.write("diff = 255,85,85\n")
    COLORS.write("gap  = 255,255,255\n")
    COLORS.close()

    # ***** IDEOGRAM ***** #
    IDEOGRAM = open(geminiset.pathTMP+"/ideogram.conf", 'w')
    IDEOGRAM.write("<ideogram>\n")
    IDEOGRAM.write("<spacing>\n")
    IDEOGRAM.write("default = 0u\n")
    IDEOGRAM.write("break   = 0u\n")
    IDEOGRAM.write("<pairwise chr2 chr1>\n")
    IDEOGRAM.write("spacing = "+str(int(alignLen/33))+"u\n")
    IDEOGRAM.write("</pairwise>\n")
    IDEOGRAM.write("</spacing>\n")
    IDEOGRAM.write("thickness             = 20p\n")  # thickness (px) of chromosome ideogram
    IDEOGRAM.write("fill                  = yes\n")
    IDEOGRAM.write("fill_color            = black\n")  # the default chromosome color is set here and any value, defined in the karyotype file overrides it
    IDEOGRAM.write("radius                = 0.85r\n")  # fractional radius position of chromosome ideogram within image
    IDEOGRAM.write("show_label            = no\n")
    IDEOGRAM.write("label_font            = default\n")
    IDEOGRAM.write("label_radius          = dims(ideogram,radius) + 0.05r\n")
    IDEOGRAM.write("label_size            = 36\n")
    IDEOGRAM.write("label_parallel        = yes\n")
    IDEOGRAM.write("label_case            = upper\n")
    IDEOGRAM.write("band_stroke_thickness = 2\n")  # cytogenetic bands
    IDEOGRAM.write("show_bands            = yes\n")  # show_bands determines whether the outline of cytogenetic bands will be seen
    IDEOGRAM.write("fill_bands            = yes\n")  # in order to fill the bands with the color defined in the karyotype file you must set fill_bands
    IDEOGRAM.write("</ideogram>\n")
    IDEOGRAM.close()

    # ***** TICKS ***** #
    TICKS = open(geminiset.pathTMP+"/ticks.conf", 'w')
    TICKS.write("show_ticks          = yes\n")
    TICKS.write("show_tick_labels    = yes\n")
    TICKS.write("show_grid          = no\n")
    TICKS.write("grid_start         = dims(ideogram,radius_inner)-0.5r\n")
    TICKS.write("grid_end           = dims(ideogram,radius_inner)\n")
    TICKS.write("<ticks>\n")
    TICKS.write("skip_first_label     = yes\n")
    TICKS.write("skip_last_label      = no\n")
    TICKS.write("radius               = dims(ideogram,radius_outer)\n")
    TICKS.write("tick_separation      = 2p\n")
    TICKS.write("min_label_distance_to_edge = 0p\n")
    TICKS.write("label_separation = 5p\n")
    TICKS.write("label_offset     = 25p\n")
    TICKS.write("label_size = 8p\n")
    TICKS.write("multiplier = "+str(ticksMultiplier)+"\n")
    TICKS.write("color = black\n")
    # Major ticks
    TICKS.write("<tick>\n")
    TICKS.write("spacing        = "+majorTicksSpacing+"\n")
    TICKS.write("color          = black\n")
    TICKS.write("show_label     = yes\n")
    TICKS.write("suffix = \" "+majorTicksSuffix+"\"\n")
    TICKS.write("label_size     = 36p\n")
    TICKS.write("format         = %s\n")
    TICKS.write("grid           = yes\n")
    TICKS.write("grid_color     = dgrey\n")
    TICKS.write("grid_thickness = 1p\n")
    TICKS.write("thickness = 5p\n")
    TICKS.write("size      = 30p\n")
    TICKS.write("chromosomes = -chr2\n")
    TICKS.write("</tick>\n")
    # Minor ticks
    TICKS.write("<tick>\n")
    TICKS.write("spacing        = "+minorTicksSpacing+"\n")
    TICKS.write("color          = grey\n")
    TICKS.write("show_label     = no\n")
    TICKS.write("format         = %s\n")
    TICKS.write("grid           = yes\n")
    TICKS.write("grid_color     = dgrey\n")
    TICKS.write("grid_thickness = 1p\n")
    TICKS.write("thickness = 5p\n")
    TICKS.write("size      = 20p\n")
    TICKS.write("chromosomes = -chr2\n")
    TICKS.write("</tick>\n")
    TICKS.write("</ticks>\n")
    TICKS.close()

    # ***** HIGHLIGHTS ***** #
    cpt = 1
    lstHighlights = []
    dicoColor = {"X": "diff", "-": "gap", "*": "cons"}
    for align in dicoFASTA:
        dicoGroupByColor = {"*": [], "X": [], "-": []}
        lstHighlights.append(geminiset.pathTMP+"/highlights"+str(cpt)+".txt")
        for i in range(0, alignLen, 1):
            nucl = dicoFASTA[align][i].upper().replace("A", "X").replace("T", "X").replace("G", "X").replace("C", "X").replace("N", "X")
            if i not in dicoMAFFT or len(dicoMAFFT[i]) == 1:
                dicoGroupByColor["*"].append(i)
            else:
                dicoGroupByColor[nucl].append(i)
        # Range per organism and by color
        HIGHLIGHTSALIGN = open(geminiset.pathTMP+"/highlights"+str(cpt)+".txt", 'w')
        for nuclType in dicoGroupByColor:
            lstGrpByType = list(to_ranges(dicoGroupByColor[nuclType]))
            for grpByType in lstGrpByType:
                HIGHLIGHTSALIGN.write("chr1 "+str(grpByType[0])+" "+str(grpByType[1])+" fill_color="+dicoColor[nuclType]+"\n")
        HIGHLIGHTSALIGN.close()
        cpt += 1

    # ***** CIRCOSCONF ***** #
    CIRCOSCONF = open(geminiset.pathTMP+"/circos.conf", 'w')
    CIRCOSCONF.write("karyotype = "+geminiset.pathTMP+"/karyotype.txt\n")
    CIRCOSCONF.write("<<include "+geminiset.pathTMP+"/ideogram.conf>>\n")
    CIRCOSCONF.write("<<include "+geminiset.pathTMP+"/ticks.conf>>\n")
    CIRCOSCONF.write("<image>\n")
    CIRCOSCONF.write("dir = "+os.path.realpath(os.path.dirname(pathOUT))+"\n")
    CIRCOSCONF.write("file = "+os.path.basename(pathOUT)+"\n")
    CIRCOSCONF.write("png = yes\n")
    CIRCOSCONF.write("svg = no\n")
    CIRCOSCONF.write("radius = 1500p\n")
    CIRCOSCONF.write("angle_offset = -80\n")
    CIRCOSCONF.write("auto_alpha_colors = yes\n")
    CIRCOSCONF.write("auto_alpha_steps  = 5\n")
    CIRCOSCONF.write("</image>\n")
    CIRCOSCONF.write("<colors>\n")
    CIRCOSCONF.write("<<include "+geminiset.pathTMP+"/colors.conf>>\n")
    CIRCOSCONF.write("</colors>\n")
    CIRCOSCONF.write("<highlights>\n")
    rstart = 0.99
    for highlights in lstHighlights:
        CIRCOSCONF.write("<highlight>\n")
        CIRCOSCONF.write("file = "+highlights+"\n")
        CIRCOSCONF.write("r0   = "+str(round(rstart-0.04, 2))+"r\n")
        CIRCOSCONF.write("r1   = "+str(rstart)+"r\n")
        CIRCOSCONF.write("</highlight>\n")
        rstart = round(rstart-0.05, 2)
    CIRCOSCONF.write("</highlights>\n")
    CIRCOSCONF.write("<<include /etc/circos/colors_fonts_patterns.conf>>\n")
    CIRCOSCONF.write("<<include /etc/circos/housekeeping.conf>>\n")
    CIRCOSCONF.close()

    # ***** MAKE PLOT *****#
    os.system("echo \"circos -conf "+geminiset.pathTMP+"/circos.conf\" > "+geminiset.pathTMP+"/circos.log")
    os.system("circos -conf "+geminiset.pathTMP+"/circos.conf >> "+geminiset.pathTMP+"/circos.log 2>&1")
    # Display sequence name order
    print("Order from outside to inside: ")
    for align in dicoFASTA:
        print(align)
    exit()

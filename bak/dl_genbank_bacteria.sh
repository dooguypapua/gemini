#!/bin/bash
function usage(){
    printf "Usage:\n"
    printf "    $0 -o/--out <folder> -f/--filter <mode>\n\n"
    printf "Filtering modes:\n"
    printf "    all      : Any filter\n"
    printf "    complete : Filter uncomplete genome\n"
    printf "    refseq   : Filter refseq excluded\n\n"
    exit 2
}
function prog() {
    local w=50 p=$1;  shift
    # create a string of spaces, then change them to dots
    printf -v dots "%*s" "$(( $p*$w/100 ))" ""; dots=${dots// /.};
    # print those dots on a fixed-width space plus the percentage etc.
    printf "\r\e[K|%-*s| %3d%% %s" "$w" "$dots" "$p" "$*";
}
filterAssembly=/mnt/c/dgoudenege/Dev_prog/gemini/filter_assembly.py
delayProgress=10
# Stop if used disk space > limit%
maxPercentUsedDisk=95
# Banner
toilet -f term -F border " BACTERIA GENBANK DOWNLOAD "
echo ""
# Retrieve getopt arguments
inputArgs=$(getopt -n $(basename "$0") --options "-o","-f","-h" --longoptions out,filter,help -- "$@")
validArgs=$?
if [ "$validArgs" != "0" ] || [ "$validArgs" == "--" ]; then usage ; fi
eval set -- "$inputArgs"
while :
do
  case "$1" in
    --out) outDir="$2" ; shift 2 ;;
    -o) outDir="$2" ; shift 2 ;;
    --filter) filter="$2" ; shift 2 ;;
    -f) filter="$2" ; shift 2 ;;
    --help)  help="y" ; usage ;;
    -h)  help="y" ; usage ;;
    --) shift; break ;;
    *) printf "\033[31m[ERROR] Missing filtering option\033[0m\n\n"
       usage ;;
  esac
done
# Check missing arguments
missingArgs=""
if [ -z "${outDir}" ] ; then missingArgs="${missingArgs}-o/--out\n" ; fi
if [ -z "${filter}" ] ; then missingArgs="${missingArgs}-f/--filter \n" ; fi
if [ "$missingArgs" != "" ] ; then printf "\033[31m[ERROR] Missing arguments ${missingArgs}\033[0m\n" ; usage ; fi
# Check output directory
realOutDir=$(realpath -s ${outDir})
if ! [ -d ${realOutDir} ]; then printf "\033[31m[ERROR] Invalid output directory \"${realOutDir}\"\033[0m\n\n" ; usage ; fi
# Check bad filtering option
if [ "$filter" != "all" ] && [ "$filter" != "complete" ] && [ "$filter" != "refseq" ] ; then printf "\033[31m[ERROR] Invalid filtering option \"${filter}\"\033[0m\n\n" ; usage ; fi
# Informations
echo -e "❱ Start    : $(date "+%D %T")"
dldir=${realOutDir}/genbank_${filter}_$(date +%d_%b_%Y)
echo -e "❱ Output   : ${dldir}"
if [ "$filter" == "all" ]; then echo -e "❱ Filter   : all\n"-
elif [ "$filter" == "complete" ]; then echo -e "❱ Filter   : complete genome\n"
elif [ "$filter" == "refseq" ]; then echo -e "❱ Filter   : refseq\n"
fi
# Make output folder
mkdir -p $dldir
# Download assembly summary
echo -e "❱ Download assembly_summary"
if [ "$filter" == "all" ]; then wget -q --no-clobber -P $dldir ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
elif [ "$filter" == "complete" ]; then wget -q --no-clobber -P $dldir ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
elif [ "$filter" == "refseq" ]; then wget -q --no-clobber -P $dldir ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
fi
# Apply filtering
echo -e "\n❱ Filtering assembly_summary"
python3 ${filterAssembly} ${dldir}/assembly_summary.txt ${filter} > ${dldir}/assembly_summary_url.txt
nbassembly=$(wc -l ${dldir}/assembly_summary_url.txt | cut -d " " -f 1)
echo -e "==> ${nbassembly} assemblies"
# Download assembly file
echo -e "\n❱ Downloading"
aria2c --log=${dldir}/log.out --dir=${dldir} --input-file=${dldir}/assembly_summary_url.txt --max-concurrent-downloads=16 --max-connection-per-server=16 --quiet &
PIDaria=$!
while kill -0 $PIDaria >/dev/null 2>&1
    do
    nbDone=$(find db/genbank_refseq_18_Mar_2022 -name \*\.gz | wc -l)
    usedPercent=$(df -h --output=pcent / | tail -n 1 | xargs | sed s/"%"//)
    if (( $usedPercent > $maxPercentUsedDisk )) ; then printf "\033[31m[ERROR] Disk used % exceed limit \"${maxPercentUsedDisk}\"\033[0m\n\n" ; exit ; fi
    percent=$(( $nbDone * 100 / $nbassembly ))
    prog "$percent"
    sleep $delayProgress
done
prog "100"
echo -e "\n\n❱ Complete : $(date "+%D %T")\n"
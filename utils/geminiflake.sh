gemini_path=$(dirname $0 | sed s/"utils"/""/) #"/mnt/c/Users/dgoudenege/Dev_prog/gemini"
files=( $(find ${gemini_path}/*.py -type f | awk '{ print length, $0 }' | sort -n -s | cut -d' ' -f2-) )
longest_filename=$((${#files[-1]}-${#gemini_path}))
unvalid_file_count=0

for pyfile in ${gemini_path}/*.py
    do
    filename=$(basename ${pyfile})
    rjust=$(($longest_filename-${#filename}))
    # Count flake8 error/warning
    if ! [ "$filename" = "geminitmp.py" ] && ! [ "$filename" = "filter_assembly.py" ] && ! [ "$filename" = "gemini.py" ] ; then
        count=$(flake8 --max-line-length 500 ${pyfile} | wc -l)
        # Display "OK" or detailled error/warning
        if [ "$count" = "0" ]; then
            if [ "${1}" != "summary" ]
                then
                printf "  ${filename}%s%${rjust}s\e[32m✅\e[0m\n"
            fi
        else
            unvalid_file_count=$((unvalid_file_count+1))
            if [ "${1}" != "summary" ]
                then
                printf "  ${filename}%s%${rjust}s\e[31m⛔ (${count})\e[0m\n"
                flake8 --max-line-length 500 ${pyfile} | sed "s,${pyfile}:,    line,g"
            fi
        fi
    fi
    # for gemini exclude import error
    if [ "$filename" = "gemini.py" ] ; then
        count=$(flake8 --max-line-length 500 --extend-ignore=F403,F401 ${pyfile} | wc -l)
        # Display "OK" or detailled error/warning
        if [ "$count" = "0" ]; then
            if [ "${1}" != "summary" ]
                then
                printf "  ${filename}%s%${rjust}s\e[32m✅\e[0m\n"
            fi
        else
            unvalid_file_count=$((unvalid_file_count+1))
            if [ "${1}" != "summary" ]
                then
                printf "  ${filename}%s%${rjust}s\e[31m⛔ (${count})\e[0m\n"
                flake8 --max-line-length 500 --extend-ignore=F403,F401 ${pyfile} | sed "s,${pyfile}:,    line,g"
            fi
        fi
    fi
done

if [ "${1}" == "summary" ]
    then
        if [ $unvalid_file_count -ge 1 ]
            then echo "⛔(${unvalid_file_count})"
            else echo "✅"
        fi
fi
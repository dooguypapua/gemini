# Clean gemini temporary folders

# Get list of gemini temporary folders
tmp_dir="/tmp"
lst_dir_toclean=()
for gemini_tmp_dir in ${tmp_dir}/gemini*
    do
        if [[ -d "$gemini_tmp_dir" ]]
        then
            dirname=$(basename ${gemini_tmp_dir})
            lst_dir_toclean+=($dirname)
        fi
done

# get longest dirname
longest_dirname=0
for dir_toclean in "${lst_dir_toclean[@]}"
do
    if [ ${#dir_toclean} -gt ${longest_dirname} ]
        then longest_dirname=${#dir_toclean}
    fi
done

# Display number of folder to clean
if [ ${#lst_dir_toclean[@]} = 0 ]
    then
    echo "Any temporary gemini found"
elif [ ${#lst_dir_toclean[@]} = 1 ]
    then
    echo "Cleaning ${#lst_dir_toclean[@]} gemini folder"
else
    echo "Cleaning ${#lst_dir_toclean[@]} gemini folders"
fi

# Delete gemini temporary folders
for dir_toclean in "${lst_dir_toclean[@]}"
do
    rjust=$(($longest_dirname-${#dir_toclean}))
    printf "  ${dir_toclean}%s%${rjust}s"
    rm -r "${tmp_dir}/${dir_toclean}" 2> /dev/null
    if [ $? -eq 0 ]
        then
        printf "\e[32m ✅\e[0m\n"
    else
        printf "\e[32m ⛔\e[0m\n"
    fi
done

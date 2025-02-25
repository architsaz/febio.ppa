table_name="von_mises.txt"
study_names=("pst.1" "pst.2")
output_filename="combin_von_mises.txt"
if [ -f "$output_filename" ]; then
    echo "ERROR: file $output_filename exist."
    exit 1
fi
# List of cases
if [ -f /dagon1/achitsaz/runfebio/successful_cases.txt ];then 
    cases=($(cat /dagon1/achitsaz/runfebio/successful_cases.txt))
else
    dir_root=/dagon1/achitsaz/runfebio
    echo "ERROR: the successful_cases.txt does not exist in this path: $dir_root"  
    exit
fi
completed_cases=0
write_header_line=0

total_cases=${#cases[@]}
echo "Total # cases in the list: $total_cases"

run_case () {
    local case_dir=/dagon1/achitsaz/runfebio
    if [ -f $case_dir/$1/$2/$table_name ];then
        if [ $write_header_line -lt "1" ];then
            grep "Casename" $case_dir/$1/$2/$table_name >> $output_filename
            write_header_line=1
        fi
        grep "$1" $case_dir/$1/$2/$table_name >> $output_filename
        # echo "-> found a $table_name for case $1 "
    else    
        echo "ERROR: Can not found a $table_name for case $1 !" | tee -a ERROR.txt
    fi
}

while [ $completed_cases -lt $total_cases ]; do
    #echo "* starting working on case : ${cases[0]}"
    case=${cases[0]}
    cases=("${cases[@]:1}")
    for study in "${study_names[@]}"; do
        run_case "$case" "$study"
    done
    if [ ${#cases[@]} -eq 0 ]; then
        echo "All cases completed!"
        exit
    fi 
done













# found_tables=$(find /dagon1/achitsaz/runfebio/*/$study_name/* -type f -name $table_name)
# # Count how many files are found
# file_count=$(echo "$found_files" | wc -l)

# first_file=true

# if [ -z "$found_files" ]; then
#     echo "No file named $file_name found."
# else
#     echo "File(s) found: $file_count"
#     #echo "$found_files"
#     # Append the contents of each found file to the output file
#     for file in $found_files; do
#         if [ "$first_file" = true ]; then
#             cat "$file" >> "$output_file"
#             first_file=false
#         else
#             tail -n +2 "$file" >> "$output_file"
#         fi
#     done
# fi

# mean_meanSmax_aneu=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $4; count++ } } END { if (count > 0) print sum / count; }')
# mean_meanSmax_red=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $5; count++ } } END { if (count > 0) print sum / count; }')
# mean_meanSmax_wht=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $7; count++ } } END { if (count > 0) print sum / count; }')
# mean_meanSmax_dom=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $8; count++ } } END { if (count > 0) print sum / count; }')
# mean_meanSmax_bod=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $9; count++ } } END { if (count > 0) print sum / count; }')
# mean_meanSmax_nek=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $10; count++ } } END { if (count > 0) print sum / count; }')
# mean_meanSmax_part=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $11; count++ } } END { if (count > 0) print sum / count; }')
# mean_meanSmax_pres=$(grep "mean" $output_file | awk '{ if ($4 != 0) { sum += $12; count++ } } END { if (count > 0) print sum / count; }')

# mean_stddevSmax_aneu=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $4; count++ } } END { if (count > 0) print sum / count; }')
# mean_stddevSmax_red=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $5; count++ } } END { if (count > 0) print sum / count; }')
# mean_stddevSmax_wht=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $7; count++ } } END { if (count > 0) print sum / count; }')
# mean_stddevSmax_dom=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $8; count++ } } END { if (count > 0) print sum / count; }')
# mean_stddevSmax_bod=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $9; count++ } } END { if (count > 0) print sum / count; }')
# mean_stddevSmax_nek=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $10; count++ } } END { if (count > 0) print sum / count; }')
# mean_stddevSmax_part=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $11; count++ } } END { if (count > 0) print sum / count; }')
# mean_stddevSmax_pres=$(grep "stddev" $output_file | awk '{ if ($4 != 0) { sum += $12; count++ } } END { if (count > 0) print sum / count; }')


# echo -e "Allcases\tstudy\tmean\t$mean_meanSmax_aneu\t$mean_meanSmax_red\t0\t$mean_meanSmax_wht\t$mean_meanSmax_dom\t$mean_meanSmax_bod\t$mean_meanSmax_nek\t$mean_meanSmax_part\t$mean_meanSmax_pres" >> "$output_file"
# echo -e "Allcases\tstudy\tstddev\t$mean_stddevSmax_aneu\t$mean_stddevSmax_red\t0\t$mean_stddevSmax_wht\t$mean_stddevSmax_dom\t$mean_stddevSmax_bod\t$mean_stddevSmax_nek\t$mean_stddevSmax_part\t$mean_stddevSmax_pres" >> "$output_file"
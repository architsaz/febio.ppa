#!/bin/bash
# List of cases
if [ -f /dagon1/achitsaz/runfebio/all_cases.txt ];then 
    cases=($(cat /dagon1/achitsaz/runfebio/all_cases.txt))
else
    dir_root=/dagon1/achitsaz/runfebio
    echo "ERROR: the all_cases.txt does not exist in this path: $dir_root"  
    exit
fi
 
completed_cases=0

total_cases=${#cases[@]}
echo "Total # cases in the list: $total_cases"

run_case () {
    local case_dir=/dagon1/achitsaz/runfebio
    if [ -f /dagon1/jcebral/region/R01/wall/$1.bleb ];then
        cp /dagon1/jcebral/region/R01/wall/$1.bleb  $case_dir/$1/data/
        echo "-> found a mask for case $1 "
    fi    
}

while [ $completed_cases -lt $total_cases ]; do
    echo "* starting working on case : ${cases[0]}"
    case=${cases[0]}
    cases=("${cases[@]:1}")
    run_case "$case"
    if [ ${#cases[@]} -eq 0 ]; then
        echo "All cases completed!"
        exit
    fi 
done


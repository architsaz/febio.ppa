#!/bin/bash
code_name=EXEC_ppa
# List of machines
# machines=("ishtar" "loki" "hades" "attila" "marduk" "heise")
machines=("ishtar")
# List of cases
list_dir=/dagon1/achitsaz/runfebio/accept_runs_cases.txt
if [ -f "$list_dir" ]; then 
    cases=($(cat "$list_dir"))
else
    echo "ERROR: the list of case files does not exist in this path: $list_dir"  
    exit 1
fi
 
completed_cases=0

total_cases=${#cases[@]}
echo "Total # cases in the list: $total_cases"
echo "Total # of available mashines: ${#machines[@]}"

run_case_on_machine () {
    local case_dir=/dagon1/achitsaz/runfebio
    local code_dir=/dagon1/achitsaz/mylib/EXECs
    echo "-> Running $2 on machine $1"
    if [ $1 == "ishtar" ]; then
        #mkdir $code_dir/$dir_name/$2/comp.1.2
        # cd $code_dir/$dir_name/$2/msa.2/ 
        # cp /dagon1/achitsaz/FEBio/scripts/input_homo.txt  input.txt
        # cd $code_dir/$dir_name/$2/msa.2/ 
        # cp /dagon1/achitsaz/FEBio/scripts/input_hete.txt  input.txt
       cd $case_dir/$2/pst.1/
       nohup $code_dir/$code_name -c $2 -s msa.1 > run.log 2>&1 &
    else
        #ssh $1 "mkdir $code_dir/$dir_name/$2/comp.1.2 && cd $code_dir/$dir_name/$2/comp.1.2/ && nohup $code_dir/scripts/build/ppa_exec -c $2 -s msa.1 msa.2 -i 0 0 > run.log 2>&1 &" &
        ssh $1 "cd $case_dir/$2/pst.1/ && nohup $code_dir/$code_name -c $2 -s msa.1 > run.log 2>&1 &" &
    fi 
}

while [ $completed_cases -lt $total_cases ]; do
    echo "* starting to find a available machine for case : ${cases[0]}"
    find_machine=0
    for machine in "${machines[@]}"; do
        #check if the machine run a task
        if [ $machine == "ishtar" ]; then
            running_task=$(pgrep EXEC_ppa | wc -l)
        else
            running_task=$(ssh $machine "pgrep EXEC_ppa | wc -l")
        fi
        if [ "$running_task" -eq 0 ] && [ ${#cases[@]} -gt 0 ]; then
            find_machine=1
            case=${cases[0]}
            cases=("${cases[@]:1}")
            run_case_on_machine "$machine" "$case"
            if [ ${#cases[@]} -eq 0 ]; then
                echo "All cases completed!"
                exit
            fi 
        fi
    done    
    completed_cases=$(($total_cases-${#cases[@]}))
    if [ $find_machine -eq 0 ]; then
        sleep 1m
    else
        sleep 10
    fi    
done


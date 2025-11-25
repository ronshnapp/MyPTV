#!/bin/bash  
start=$(date +'%r')
python workflow.py params_file.yml matching
end=$(date +'%r')
echo "$start"
echo "$end"
start=$(date +'%r')
python workflow.py params_file.yml tracking
end=$(date +'%r')
echo "$start"
echo "$end"
start=$(date +'%r')
python workflow.py params_file.yml smoothing
end=$(date +'%r')
echo "$start"
echo "$end"

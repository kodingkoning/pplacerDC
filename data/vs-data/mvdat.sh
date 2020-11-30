#!/bin/bash
for size in 500 1000 5000 10000; do
  cd $size
  cp /home/malachi2/cs581-project/runs/variable-size/data/$size/delta_error_approach1.txt .
  cp /home/malachi2/cs581-project/runs/variable-size/data/$size/time_approach1.txt .
  cp /home/malachi2/cs581-project/runs/variable-size/data/$size/run_output.log .

  # stuff for approach 2
  cp /home/malachi2/cs581-project/runs/variable-size/data/$size/delta_error_approach2.txt .
  cp /home/malachi2/cs581-project/runs/variable-size/data/$size/delta_error_apples.txt .
  cp /home/malachi2/cs581-project/runs/variable-size/data/$size/time_approach2.txt .
  cp /home/malachi2/cs581-project/runs/variable-size/data/$size/run_output2.log .
  cd ../
done

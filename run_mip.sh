#!/bin/sh

instances="a1_1 a1_2 a1_3 a1_4 a1_5"

timelimit=21600  # 6h

models="mip1 lp1 mip2 lp2 mip3 lp3"

for instance in ${instances}; do
    for model in ${models}; do
        bin/${model} -t ${timelimit} -p data/model_${instance}.txt -i data/assignment_${instance}.txt -o output/out_${model}_${instance}.txt
    done
done
         

#!/bin/bash

nextflow run main.nf --input=test/data/samplesheet.csv --outdir=results -resume -w /home/sturm/scratch/projects/2020/single-cell-analysis-nf/work -profile icbi

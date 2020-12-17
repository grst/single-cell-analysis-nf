#!/bin/bash

nextflow run main.nf --input=test/data/samplesheet.csv --outdir=results -ansi-log false -resume

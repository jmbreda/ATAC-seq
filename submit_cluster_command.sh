#!/bin/bash
snakemake -s workflow/atac_seq_preprocessing.smk -j 999 --cluster-config config/cluster.json --cluster "sbatch  --job-name {cluster.name} --output {cluster.stdout} --error {cluster.stderr} --qos {cluster.qos} --time {cluster.time} --mem {cluster.mem} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus-per-task}"

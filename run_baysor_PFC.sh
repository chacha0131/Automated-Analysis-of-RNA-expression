#!/bin/bash
path="/home/bagotlab/ugrad/Emily/PFC_witness/cohort3/" 
for file in ${path}out/nuclei_masked/*.tif; do
    name=${file##*/}
    id=${name%.tif}
    echo ${id}
    mkdir "${path}out/quant/${id}"
    nuclei="${path}out/nuclei_masked/${id}.tif"
    rna="${path}out/rna/rna_clean/${id}.csv"
    config="${path}config.toml"
    out="${path}out/quant/${id}/${id}.csv"
    /home/bagotlab/ugrad/Emily/PFC_witness/cohort2/baysor_ubuntu-latest_x64_build/Baysor/bin/Baysor run --prior-segmentation-confidence 0.95 \
        --n-clusters 2 \
        --no-ncv-estimation \
        --iters 500 -o $out -c $config $rna $nuclei
done
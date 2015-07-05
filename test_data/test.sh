#!/bin/bash

mkdir 10genes
mkdir sequences
for j in {1..5}
do
    mkdir 10genes/t${j}
    mkdir 10genes/t${j}/seqs
    python ../genetree.py -n 5 -sp test_tree --out_dir 10genes/t${j}/
    python ../phase2.py -sp 10genes/t${j}/test_tree -n 10 -r roots/root_seqs.txt --out_dir 10genes/t${j}/seqs/
    python ../agalma_format.py -i 10genes/t${j}/seqs/tree -d sequences -n 10
done
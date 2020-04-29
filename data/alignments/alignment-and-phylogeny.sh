#!/bin/bash

#script to align Spike and RBD sequences using mafft, and infer phylogeny with RAxML

mafft --reorder --op 4.5 ./unaligned-sequences/Spikes.fasta > ./Spikes_aligned.fasta

mafft --reorder --op 4.5 ./unaligned-sequences/RBDs.fasta > ./RBDs_aligned.fasta

mafft --reorder --op 3.0 ./unaligned-sequences/Spikes_nt.fasta > ./Spikes_nt_aligned.fasta

mafft --reorder --op 3.0 ./unaligned-sequences/RBDs_nt.fasta > ./RBDs_nt_aligned.fasta

mkdir RBD_tree
cd ./RBD_tree
raxmlHPC -s ../RBDs_aligned.fasta -n RBD_tree.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -x 10

mkdir ../Spike_tree
cd ../Spike_tree
raxmlHPC -s ../Spikes_aligned.fasta -n Spike_tree.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -x 10

mkdir ../RBD_nt_tree
cd ../RBD_nt_tree
raxmlHPC -s ../RBDs_nt_aligned.fasta -n RBD_nt_tree.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10

mkdir ../Spike_nt_tree
cd ../Spike_nt_tree
raxmlHPC -s ../Spikes_nt_aligned.fasta -n Spike_nt_tree.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10

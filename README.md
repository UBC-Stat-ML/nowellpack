Summary
-------

[![Build Status](https://travis-ci.org/UBC-Stat-ML/nowellpack.png?branch=master)](https://travis-ci.org/UBC-Stat-ML/nowellpack)

``nowellpack`` is a Blang library for cancer genomics. The main current features focus on Bayesian phylogenetic tree inference from single cell data. We make the approximation that change points in the copy number profiles are *perfect phylogeny markers*. However those markers are observed in a noisy fashion hence the name we give to this model: **corrupt phylogenies**. 


Tree inference
--------

Run an end-to-end pipeline (sampling, summarizing posterior, viz) with:

```
git clone --depth=1 https://github.com/UBC-Stat-ML/corrupt-nextflow.git
cd corrupt-nextflow
./nextflow run main.nf -resume --tipInclusionProbabilities DATA.csv
```

Change ``DATA.csv`` into the data you are interested in (more below). Provided you have Oracle Java 8 in your ``PATH`` variable this will create a directory called ``deliverables/main`` containing an inferred tree and a bunch of other outputs.

### Input format 

The input for all included phylogenetic inference methods is a csv file with the following header:

```
cells,loci,tipInclusionProbabilities
```

followed by lines of the form

```
myfirstcell,somelocus,1
...
``` 

Note that for the binary model, tipInclusionProbabilities are set to 0 and 1, and a false positive/false negative rates model is inferred jointly.

The loci should follow the format ``1_100500001_101000000`` where 
``1`` here is the chromosome index, ``100500001`` is the left boundary of the bin, and ``101000000`` is the right boundary, both inclusive.

Cells can use arbitrary unique identifiers.



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


Tree growing
------------

Let us say you have a pre-computed tree inferred from a cell/loci matrix M and you would like to either:

- add more cells that were measured with the same set of loci as M
- add more loci that were measured with the same set of cells as M

In such cases you can use the ``corrupt-grow`` utility to quickly grow the tree by using the extra data. This is done using maximum a posteriori placement which has the advantage of being very fast but does not provide measures of uncertainty. 

To use the ``corrupt-grow`` utility:

1. Clone the repo
2. Build using ``./setup-cli.sh``
3. Add ``build/install/nowellpack/bin`` to your path
4. Invoke with ``corrupt-grow --matrix ReadOnlyCLMatrix TIDY_FILE --phylo NEWICK_FILE`` where ``NEWICK_FILE`` is the pre-computed phylogeny in newick format, and ``TIDY_FILE`` is a matrix in the format

```
cells,loci,tipInclusionProbabilities
```

followed by lines of the form

```
myfirstcell,somelocus,0.123
...
``` 

where the entries are the "local posterior distributions" described in the preprint. Alternatively, you can use ``corrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix TIDY_FILE --matrix.fpr 0.1 --matrix.fnr 0.2 --phylo NEWICK_FILE`` for a binary matrix ``TIDY_FILE`` with prescribed false positive and false negative rates.




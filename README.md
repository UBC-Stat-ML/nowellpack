Summary
-------

<!-- [![Build Status](https://travis-ci.org/alexandrebouchard/nowellpack.png?branch=master)](https://travis-ci.org/alexandrebouchard/nowellpack) -->

``nowellpack`` is a Blang library for cancer genomics. The main current features focus on Bayesian phylogenetic tree inference from single cell data. We make the approximation that change points in the copy number profiles are *perfect phylogeny markers*. However those markers are observed in a noisy fashion hence the name we give to this model: **corrupt phylogenies**. 


Find me a tree&mdash;quick!!!
--------

Run an end-to-end pipeline (sampling, summarizing posterior, viz) with:

```
git clone https://github.com/UBC-Stat-ML/corrupt-nextflow.git
cd corrupt-nextflow
./nextflow run corrupt-infer-pipeline.nf -resume --tipInclusionProbabilities demo-data.csv
```

Provided you have Oracle Java 8 in your ``PATH`` variable this will create a directory called ``deliverables/corrupt-infer-pipeline`` containing an inferred tree. Change ``demo-data.csv`` into the data you are interested in (more below). For the example tree, you will get this (left, consensus tree based on cell-locus 0-1 loss; right binary matrix where x axis indexes loci, y axis, cells; the black-and-white columns are indicators that in the consensus tree the x-axis cell has the y-axis locus' trait, the paired greyscale are uncertainties obtained from the posterior for the same quantities): 

![https://www.stat.ubc.ca/~bouchard/pub/output.png](https://www.stat.ubc.ca/~bouchard/pub/output.png)

Here is a more interesting example from a 1000 cell by 50 locus synthetic datas, showing the same information (rotated) with also a third column in red showing the indicators for the tree used to generated the dataset. 

![https://www.stat.ubc.ca/~bouchard/pub/output-eg-large.png](https://www.stat.ubc.ca/~bouchard/pub/output-eg-large.png)

Computing a tree of that scale takes less than one hour. If more loci are available, we recommend using fewer, but more reliable ones. See [synthetic benchmarking workflow](https://github.com/UBC-Stat-ML/corrupt-nextflow/blob/master/synthetic-benchmark.nf) for details on how the above example was computed.



Installation for more advanced usage
------------

You will need a unix environment with Oracle Java 8.

- Check out the source ``git clone https://github.com/UBC-Stat-ML/nowellpack.git``
- Compile using ``./gradlew installDist``
- Add the directory ``build/install/nowellpack/bin`` to your ``$PATH`` variable.


User guide
-----

### Input format 

The input for all included phylogenetic inference methods is a csv file with the following header:

```
cells,loci,tipInclusionProbabilities
```

followed by lines of the form

```
myfirstcell,somelocus,0.6
...
``` 

containing unique string identifier for the cell and locus followed by an estimated probability, here `0.6` that the given cell `myfirstcell` has the marker specified by locus `somelocus`.

From this file several phylogenetic inference methods are available. For example, to setup runs based on parallel tempering, you can look at the options with 

```
corrupt-infer --model.tipInclusionProbabilities data.csv --engine PT --help
```

For instead going a an adaptive annealing sequential change of measure route, use:

```
corrupt-infer --model.tipInclusionProbabilities data.csv --engine SCM --help
```

### Output format for posterior distributions

The command ``corrupt-infer`` creates a symlink ``results/latest`` to a unique execution folder in ``results/all``.


In ``results/latest/samples/phylo.csv``, sampled tree are stored as follows using the *newick* format:

```
sample,value
0,((cell_myfirstcell,cell_another)locus_somelocus,locus_anotherlocus);
...
```

### Summarizing the posterior distribution

To summarize the samples, we first compute for each cell *c* and locus *l*, the fraction of the sampled tree reconstructions where the cell *c* is deemed to have the marker at the locus *l*. For this use:

```
corrupt-average --csvFile path/to/samples/phylo.csv
```

which create a new exec foler with averages in ``results/latest/average.csv`` These averages are directly interpretable as reconstruction uncertainties. However they do not provide a single-tree summary. 

To do so, next, we compute the minimum Bayes risk tree reconstruction, which is defined as the tree that is closest to the posterior averages under a reasonable notion of distance (we use L1 distance on the induced cell-locus indicator functions).  

For a quick greedy approximation of the distance minimization problem, use 

```
corrupt-greedy --tipInclusionProbabilities path/to/average.csv
```

The file ``results/latest/trees.csv`` contains the list of trees visited in the greedy process of adding at each step the placement of the locus inducing the largest reduction in L1 distance to the posterior averages. The last one, also available in ``results/latest/tree.newick`` for convenience, is typically the most accurate however early stopping can perform in difficult cases.

For a less greedy alternative, instead of ``corrupt-greedy``, you can use again the adaptive annealing sequential change of measure algorithm but this time pushed to temperatures greater than one to achieve global optimization. For this use:

```
corrupt-infer --model.tipInclusionProbabilities path/to/average.csv --engine.maxAnnealingParameter 1000 --engine.nFinalRejuvenations 0
``` 

the file ``results/latest/samples/phylo.csv`` should now contain equivalent trees (otherwise further increase the maximal temperature beyond 1000.



### Alternative inference methods

The greedy and annealed procedures can also be applied directly to the data (tips inclusion probabilities), in which case we get an approximation to the maximum likelihood reconstruction instead of the minimum Bayes risk reconstruction. 

Posterior samples can also be computed using different sampling methods. For example, straight MCMC can be used via

```
corrupt-infer --model.tipInclusionProbabilities data.csv --stripped true --engine PT --engine.nChains 1
```

Here ``--stripped true`` gets rid of the random initialization and annealing. To enable several parallel tempering annealed chains, use ``--engine.nChains``. 

Use ``--help`` and see [Blang's documentation](https://www.stat.ubc.ca/~bouchard/blang/) for more information on the many configurations available. 
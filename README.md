# SPIn: model Selection in Phylogenetics based on linear INvariants 

> In phylogenetic inference an evolutionary model describes the substitution processes  along each edge of a phylogenetic tree. Misspecification of the model has important implications for the analysis of phylogenetic data.  
> Conventionally, however, the selection of a suitable evolutionary model is based on heuristics or relies on the choice of an approximate input tree. SPIn uses linear invariants to characterize a model of nucleotide evolution for phylogenetic mixtures on any number of components.   Linear invariants are constraints among the joint probabilities of the bases in the operational taxonomic units that hold irrespective of the tree topologies appearing in the mixtures.

> SPIn requires NO input tree and is designed to deal with non-homogeneous phylogenetic data consisting of multiple sequence alignments showing different patterns of evolution, e.g.  concatenated genes, exons and/or introns.  


### Online version of the package: [SPIn]

### Compilation
```
$ g++ -Wno-deprecated SPIn_main.cpp -o SPIn
```

To run, simply call the executable *./SPIn* and follow given instructions.

*Input*: multiple sequence alignment in [fasta format]. 
### Notes
Currently, maximum number of OUTs is 21 and sequence length < 10^6 bp, and uses the Akaike Information Criterion (AICc) to score the models. The output contains the weights of support for each of the model and an upper bound on the number of mixtures for which the parameters (both continuous and discrete) are identifiable. Currently supported evolutionary models are the non-homogeneous Kimura 2-paramater (K80*),  Kimura 3-parameter (K81*), Jukes-Cantor (JC69*) and the Strand Symmetric Model (SMM)

### Data from the paper: [grid search]

### Publications
>A. M. Kedzierska, M. Drton, R. Guigo and M. Casanellas, "SPIn: model selection for phylogenetic mixtures via linear invariants." (Mol. Biol. Evol., 29(3): 929-937, 2012)

For theory on this and related algebraic methods see:
* Casanellas M., Fernandez-Sanchez J. and Kedzierska A. M.: "The space of phylogenetic mixtures of equivariant models" (Algorithms for Molecular Biology,
2012, 7:33)
* Kedzierska A.M. and Casanellas M.: “ GenNon-h: Simulating multiple sequence alignments under the nonhomogeneous dna models." (BMC Bioinformatics 2012, 13:216)
*  Kedzierska A.M. and Casanellas M.: “EM for parameter estimation of Markov models on trees” (revision).

[SPIn]: <http://genome.crg.es/cgi-bin/phylo_mod_sel/AlgModelSelection.pl>
[fasta format]: <http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml>
[grid search]: <http://genome.crg.es/phylo_mod/DATA_SPIn/>

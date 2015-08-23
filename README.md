# SPIn 
 
* To compile 
 g++ -Wno-deprecated SPIn_main.cpp -o SPIn



A. M. Kedzierska, M. Drton, R. Guigo and M. Casanellas, "SPIn: model selection for phylogenetic mixtures via linear invariants." (Mol. Biol. Evol., 29(3): 929-937, 2012)

Currently supported evolutionary models are non-homogeneous the Kimura 2-paramater (K80*),  Kimura 3-parameter (K81*), Jukes-Cantor (JC69*) 
and the Strand Symmetric Model (SMM)

For theory behind the method check:
M. Casanellas, J. Fernandez-Sanchez and A. M. Kedzierska, "The space of phylogenetic mixtures of equivariant models",
    

Input format to SPIn is a fasta file. Current maximum number of operational taxonomic units is 21 and sequence length of 1 million bases.
This release of the software uses the Akaike Information Criterion (AICc) to score among the candidate non-homogeneous classes of models. 
The best-fit model minimizes the AICc score. 

In addition, the output reports the weights of support for each of the model and an upper bound on the number of mixtures, above which the non-identifiability of the parameters (both continuous and discrete) holds.


Online version: http://genome.crg.es/cgi-bin/phylo_mod_sel/AlgModelSelection.pl
Data used for tests: http://genome.crg.es/phylo_mod/DATA_SPIn/

###NOTE: the code is undergoing modularity revamp and will be available very soon!

FLU.py

reads in Fasta sequences of flu amino acid sequences and organizes them into a data structure called seqlist
reads in Q matrix of amino acid substitution rates and implements a function to exponentiate it i.e. exp(Qt)

MAKEPLOTS.py

imports phylogenetic tree and calculates mean evolutionary distance between timepoints (years)
creates an array of frequencies at each site
At each site, calculates relative entropy and probability of data fitting the neutral model for each year to year transition
Makes a bunch of plots:
	relative entropy vs. position
	sites under selection
	entropy of sequences chunks (windows) over time
	number of variant sites in chunks (windows) of the sequences over time


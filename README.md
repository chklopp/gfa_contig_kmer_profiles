# Generate kmer profile from a gfa or a fasta file

Given a genome assembly in fasta or gfa format, our python script will generate a kmer count table which can be uploaded in our R script to produce graphical visualisation of the contig kmer profiles

Dependencies 
- jellyfish (https://github.com/gmarcais/Jellyfish)

Running the script 

<pre>
python gkp.py --fasta assembly.fasta --output profiles.tsv
    
</pre>

Script help 

<pre>
usage: gkp.py [-h] --fasta ASSEMBLY.fasta --gfa ASSEMBLY.gfa [--output OUTPUT] [--mpthreads MPTHREADS]

generate kmer profile.

options:
  -h, --help            show this help message and exit
  --fasta ASSEMBLY_FASTA   input assembly fasta file
  --fga ASSEMBLY_GFA   input assembly gfa file
  --output OUTPUT   output tsv file

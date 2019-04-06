Usage

Simple alignment:

- `clustalo --infile=myccseqs.fasta --hmm-in=PfamGH5.hmm --full --full-iter --percent-id --distmat-out=myccpidm.txt --guidetree-out=mycctree.nw --outfile=myccaln.clustal --outfmt=clu`

Go to each fullseqs folders from rounds 1,2,4 folders and run hmmsearch on proteinseqs.gbk, e.g.:

`hmmsearch --noali -E 1e-10 --incE 1e-10 ../../HMMprofiles/PfamGHs.hmm proteinseqs.gbk > pfamres.txt`

- Currently use `-E 1e-1 --incE 1e-1` in CAZYfullseqs folder 
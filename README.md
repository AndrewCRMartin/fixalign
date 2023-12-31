FixAlign
========

(c) 2010 UCL, Andrew C.R. Martin
--------------------------------

This is very much a work in progress that hasn't been worked on since 2010!

`FixAlign` is designed to move gaps in loops to a site which causes
minimal disruption to the loop conformation.

It takes a sequence alignment together with secondary structure
assignments and a PDB file and writes a new, modified sequence
alignment.

If not available locally, Get a DSSP file using 

```
wget ftp://ftp.ebi.ac.uk/pub/databases/dssp/1mcp.dssp.Z
```

If not available locally, Get a PDB file using 
```
wget ftp://ftp.ebi.ac.uk/pub/databases/rcsb/pdb/data/structures/all/pdb/pdb1crn.ent.gz
```
Program `dssp2ss` extracts the secondary structure from a DSSP file

Program fixalign takes as parameters:

1. A secondary structure file as produced by dssp2ss
2. A PDB file
3. Optionally an alignment file in PIR format (or stdin if not specified)
4. Sub-optionally an output alignment file in PIR format (or stdout if not specified)
5. Optionally, a -c flag specifies the chain to work on - if not specified, it assumes the first chain

Quick test...
-------------

Obtained
```
   12e8H1.1neu00.seqaln.pir
   12e8H1.1neu00.strucaln.pir
```
Assumed 1neu00 was the parent and 12e8H1 the target, so also obtained
```
   pdb1neu.ent
   1neu.dssp
```
Using dssp2ss, created
```
   1neu.ss
```
Used fixalign to fix the sequence alignment, creating:
```
   12e8H1.1neu00.fixalign.pir
```

Built models based on all 3 alignments, using SwissModel and calculated
CA-RMS with ProFit:
```
   strucaln 3.250  (i.e. the best we can hope for)
   seqaln   6.150  (Needleman & Wunsch)
   fixaln   6.204  (N&W after running fixalign)
```

Oh well!!!

Note... This was a bad example where most of the alignment is wrong -
only 44% correct - and <18% sequence identity whereas this method is
really aimed just at refining the precise position of a deletion once
a reasonable alignment has been established


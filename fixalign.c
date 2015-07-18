/*************************************************************************

   Program:    FixAlign
   File:       fixalign.c
   
   Version:    V1.0
   Date:       25.02.10
   Function:   Modify alignments to minimize structural disruption in
               loop regions
   
   Copyright:  (c) UCL / Dr. Andrew C.R. Martin, 2010
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   FixAlign is a program to take an alignment of a target and template
   sequence together with a structure and secondary structure assignments.
   It identifies gaps in the target sequence which occur in, or near to
   loops and shifts the gaps to positions that will cause the least
   structural disruption when the appropriate residues are deleted from
   the template structure

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  25.02.10 Original

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 800
#define MAXSEQS 3
#define FLEX    3
#define OPTCACADIST 3.8
#define CADISTSQ (OPTCACADIST*OPTCACADIST)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *ssfile, char *pdbfile, 
                  char *infile, char *outfile, char *chain, 
                  BOOL *verbose, int *flex);
void Usage(void);
PDB *GetChain(PDB *pdbin, char chain);
BOOL ReadSS(FILE *fp, char chain, char **ss, char **seq);
BOOL CheckData(FILE *fp, char *seq, char *ssSeq, PDB *pdb);
char *CleanSeqMalloc(char *inseq);
BOOL FixAlignment(char *seq1, char *seq2, char *ss, PDB *pdb, 
                  BOOL verbose, int flex);
BOOL InLoop(int pos, int gaplen, char *ss, int flex, 
            int *start, int *stop);
int FindBestGapPosition(PDB **pdbindex, int pos, int start, int stop, 
                        int gaplen, BOOL verbose);
BOOL MoveGap(char *seq, int oldpos, int gaplen, int bestpos);
void WritePIR(FILE *out, char *seq, SEQINFO seqinfo);
int FindGapLen(int i, char *seq);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------

   25.02.10  Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char ssfile[MAXBUFF],
        pdbfile[MAXBUFF],
        infile[MAXBUFF],
        outfile[MAXBUFF];
   char chain[8];
   char *ss, *ssSeq;
   FILE    *in     = stdin,
           *out    = stdout,
           *fpSS   = NULL,
           *fpPDB  = NULL;
   PDB     *pdb    = NULL;
   int     natoms  = 0,
           nchain  = 0,
           flex    = FLEX;
   BOOL    ssError = FALSE,
           verbose = FALSE;
   SEQINFO seqinfo[2];
   char    *seqs[MAXSEQS];

   /* Parse the command line                                            */
   if(!ParseCmdLine(argc, argv, ssfile, pdbfile, infile, outfile, 
                    chain, &verbose, &flex))
   {
      Usage();
      return(0);
   }

   /* Open the alignment files if specified or use stdin/out            */
   if(!OpenStdFiles(infile, outfile, &in, &out))
   {
      fprintf(stderr,"Unable to open alignment files\n");
      return(1);
   }

   /* Open the secondary structure file                                 */
   if((fpSS = fopen(ssfile, "r")) == NULL)
   {
      fprintf(stderr,"Unable to open SS file: %s\n", ssfile);
      return(1);
   }

   /* Open the PDB file                                                 */
   if((fpPDB = fopen(pdbfile, "r")) == NULL)
   {
      fprintf(stderr,"Unable to open PDB file: %s\n", pdbfile);
      return(1);
   }

   /* Read the coordinates                                              */
   if((pdb = ReadPDBAtoms(fpPDB, &natoms))==NULL)
   {
      fprintf(stderr,"Unable to read PDB coordinates from %s\n", pdbfile);
      return(1);
   }

   /* Get the chain of interest                                         */
   if((pdb = GetChain(pdb, chain[0]))==NULL)
   {
      fprintf(stderr,"Unable to extract chain %s from PDB file %s\n", 
              chain, pdbfile);
      return(1);
   }

   /* Get only the CA atoms                                             */
   pdb = SelectCaPDB(pdb);

   /* Read the SS information for the specified chain                   */
   if(!ReadSS(fpSS, chain[0], &ss, &ssSeq))
   {
      fprintf(stderr,"Unable to read SS data from %s for %s chain\n", 
              ssfile, (chain[0]?chain:"first"));
      return(1);
   }

   /* Read the first sequence from the alignment file                  */
   nchain = ReadRawPIR(in, seqs, MAXSEQS, TRUE,
                       &(seqinfo[0]), &ssError);
   if((nchain != 1) || ssError)
   {
      fprintf(stderr,"Unable to read first sequence from input \
alignment file\n");
      if(nchain > 1)
      {
         fprintf(stderr,"There were multiple chains\n");
      }
      return(1);
   }

   /* Read the second sequence from the alignment file                  */
   nchain = ReadRawPIR(in, seqs+1, MAXSEQS, TRUE, 
                       &(seqinfo[1]), &ssError);
   if((nchain != 1) || ssError)
   {
      fprintf(stderr,"Unable to read second sequence from input \
alignment file\n");
      if(nchain > 1)
      {
         fprintf(stderr,"There were multiple chains\n");
      }
      return(1);
   }

   /* Check that the sequence data from the 3 sources agrees            */
   if(!CheckData(stderr, seqs[0], ssSeq, pdb))
   {
      return(1);
   }

   /* Do the actual work of fixing the alignment                        */
   FixAlignment(seqs[0], seqs[1], ss, pdb, verbose, flex);

   /* Print the resulting alignment                                     */
   WritePIR(out, seqs[0], seqinfo[0]);
   WritePIR(out, seqs[1], seqinfo[1]);
   
   /* Free memory allocated by ReadSS()                                 */
   if(ss!=NULL) free(ss);
   if(ssSeq!=NULL) free(ssSeq);
   
   return(0);
}


/************************************************************************/
/*>void WritePIR(FILE *out, char *seq, SEQINFO seqinfo)
   ----------------------------------------------------
   Input:   FILE    *out      Output file pointer
   Output:  char    *seq      Sequence to be written
   Returns: SEQINFO seqinfo   A sequence information structure

   A rather simple PIR writer - needs improvements!

   25.02.10  Original   By: ACRM
*/
void WritePIR(FILE *out, char *seq, SEQINFO seqinfo)
{
   fprintf(out, ">P1;%s\n", seqinfo.code);
   fprintf(out, "%s",       seqinfo.name);
   if(seqinfo.source[0])
   {
      fprintf(out, " - %s\n", (seqinfo.source));
   }
   else
   {
      fprintf(out, "\n");
   }
   fprintf(out, "%s*\n", seq);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *ssfile, char *pdbfile, 
                     char *infile, char *outfile, char *chain, 
                     BOOL *verbose, int *flex)
   ----------------------------------------------------------------------
   Input:   int   argc     Command line argc
            char  **argv   Command line argv
   Output:  char  *ssfile  SS filename
            char  *pdbfile PDB filename
            char  *infile  Input alignment filename (or blank string)
            char  *outfile Output alignment filename (or blank string)
            char  *chain   Chain name (1-char) - maybe null
            BOOL  *verbose Verbose mode
            int   *flex    Extend loop boundaries by this number of 
                           residues
   Returns: BOOL           Success

   Parse the command line

   25.02.10  Original   By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *ssfile, char *pdbfile, 
                  char *infile, char *outfile, char *chain, BOOL *verbose,
                  int *flex)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   ssfile[0]  = pdbfile[0]   = '\0';
   chain[0] = '\0';

   if(!argc) return(FALSE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         case 'c':
            argc--;
            argv++;
            strcpy(chain, argv[0]);
            break;
         case 'e':
            argc--;
            argv++;
            sscanf(argv[0], "%d", flex);
            break;
         case 'v':
            *verbose = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2, 3 or 4 arguments left               */
         if(argc < 2 || argc > 4)
            return(FALSE);
         
         /* Copy the first to ssfile and second to pdbfile              */
         strcpy(ssfile, argv[0]);
         argc--;
         argv++;
         strcpy(pdbfile, argv[0]);
         argc--;
         argv++;
        
         /* If another, Copy to infile                                  */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         }
         
         /* If there's another, copy it to outfile                      */
         if(argc)
         {
            strcpy(outfile, argv[0]);
            argc--;
            argv++;
         }
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Print usage message

   25.02.10  Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr, "\nfixalign V1.0 (C) 2010, Dr. Andrew C.R. Martin, \
UCL\n\n");
   fprintf(stderr, "Usage: fixalign [-v][-c chain][-e extend] ssfile \
pdbfile \n");
   fprintf(stderr, "                [align.in [align.out]]\n");
   fprintf(stderr, "       -v Verbose\n");
   fprintf(stderr, "       -c Specify chain to handle (default is first \
chain)\n");
   fprintf(stderr, "       -e Specify amount to extend loop boundaries \
to consider a gap as \n");
   fprintf(stderr, "          in a loop (default %d)\n", FLEX);
   fprintf(stderr, "       ssfile    - Secondary structure file from \
dssp2ss for template\n");
   fprintf(stderr, "       pdbfile   - Coordinates of template\n");
   fprintf(stderr, "       align.in  - PIR alignment file - template \
first\n");
   fprintf(stderr, "       align.out - revised PIR alignment file - \
template first\n");

   fprintf(stderr, "\nFixAlign is a program to take an alignment of a \
target and template\n");
   fprintf(stderr, "sequence together with a structure and secondary \
structure assignments.\n");
   fprintf(stderr, "It identifies gaps in the target sequence which \
occur in, or near to\n");
   fprintf(stderr, "loops and shifts the gaps to positions that will \
cause the least\n");
   fprintf(stderr, "structural disruption when the appropriate residues \
are deleted from\n");
   fprintf(stderr, "the template structure.\n\n");
}


/************************************************************************/
/*>PDB *GetChain(PDB *pdbin, char chain)
   -------------------------------------
   Input:   PDB   *pdbin     PDB linked list
            char  chain      Chain name to extract or null character to 
                             extract the first chain
   Returns: PDB*             PDB linked list

   Extracts the specified chain from a PDB linked list. The original list
   is modified to free any other chains. Designed to be called with
   something like:
      pdb = GetChain(pdb, 'A')
   to extract chain A, or:
      pdb = GetChain(pdb, '\0')
   to extract the first chain.

   25.02.10  Original   By: ACRM
*/
PDB *GetChain(PDB *pdbin, char chain)
{
   PDB *p, *start, *prev, *stop;
   
   /* If chain is a null character then set it to the first chain       */
   if(!chain)
   {
      chain = pdbin->chain[0];
   }

   /* Search for the start of the specified chain                       */
   prev = NULL;
   start = NULL;
   for(p=pdbin; p!=NULL; NEXT(p))
   {
      if(p->chain[0] == chain)
      {
         start = p;
         break;
      }
      prev = p;
   }

   /* If not found, then free the PDB linked list and return NULL       */
   if(!start)
   {
      FREELIST(pdbin, PDB);
      return(NULL);
   }
   
   /* If the chain we wanted isn't the first one, free the region 
      before it 
   */
   if(prev)
   {
      prev->next = NULL;
      FREELIST(pdbin, PDB);
      prev = NULL;
      pdbin = start;
   }
   
   /* Find the end of this chain                                        */
   prev = NULL;
   stop = NULL;
   for(p=pdbin; p!=NULL; NEXT(p))
   {
      if(p->chain[0] != chain)
      {
         stop = p;
         break;
      }
      prev = p;
   }

   /* If there was another chain, terminate our chain, and free 
      anything that follows 
   */
   if(stop)
   {
      prev->next = NULL;
      FREELIST(stop, PDB);
   }
   
   /* Finally return our new PDB linked list                            */
   return(pdbin);
}


/************************************************************************/
/*>BOOL ReadSS(FILE *fp, char chain, char **ss, char **seq)
   --------------------------------------------------------
   Input:   FILE  *fp    File pointer for SS file
            char  chain  Chain to extract (or null character to extract
                         first chain)
   Output:  char  **ss   Malloc'd secondary structure string
            char  **seq  Malloc's sequence string
   Returns: BOOL         Success

   Reads a simple secondary structure file of the form:
   >Chain_L
   DIVMTQSPSSLSVSAGERVTMS...
   CCCEEEECCEEEECCCCCEEEE...
   >Chain_H
   EVKLVESGGGLVQPGGSLRLSC...
   CCEEEEECCCEECCCCEEEEEE...

   Extracts the sequence and SS assignments for the specified chain,
   allocating memory as required. This will need to be freed by the
   calling routine.

   25.02.10  Original   By: ACRM
*/
BOOL ReadSS(FILE *fp, char chain, char **ss, char **seq)
{
   char *ptr = NULL;
   BOOL inChain = 0;
   int  ok = 0;
   
   while((ptr=fgetsany(fp))!=NULL)
   {
      TERMINATE(ptr);
      
      if(ptr[0] == '>')
      {
         if((ptr[7] == chain) || (chain == '\0'))
         {
            inChain = 1;
         }
         else
         {
            inChain = 0;
         }
         free(ptr);
      }
      else if(inChain)
      {
         if(inChain == 1)
         {
            *seq = ptr;
            ok++;
         }
         else if(inChain == 2)
         {
            *ss = ptr;
            ok++;
         }
         else
         {
            return(TRUE);
         }
         inChain++;
      }
      else
      {
         free(ptr);
      }
      if(ok == 2)
      {
         return(TRUE);
      }
   }
   if(ok == 2)
   {
      return(TRUE);
   }
   
   return(FALSE);
}

/************************************************************************/
/*>char *CleanSeqMalloc(char *inseq)
   ---------------------------------
   Input:   char  *inseq    Input sequence
   Returns: char  *         Cleaned sequence - malloc'd

   Takes a sequence string, malloc's space for a copy and copies only 
   alphabetical characters to that copy, upcasing any lower case
   characters as it goes.

   The calling routine will need to free any space allocated.

   25.02.10  Original   By: ACRM
*/
char *CleanSeqMalloc(char *inseq)
{
   char *outseq = NULL;
   int  i, j;
   
   if((outseq = (char *)malloc((strlen(inseq)+1) * sizeof(char))) == NULL)
      return(NULL);

   for(i=0, j=0; i<strlen(inseq); i++)
   {
      if(isalpha(inseq[i]))
      {
         outseq[j++] = islower(inseq[i]) ? toupper(inseq[i]) : inseq[i];
      }
   }
   outseq[j] = '\0';
   return(outseq);
}


/************************************************************************/
/*>BOOL CheckData(FILE *fp, char *seq, char *ssSeq, PDB *pdb)
   ----------------------------------------------------------
   Input:   FILE *fp    Output filehandle for printing error message
                        (or NULL for no messages)
            char *seq   First sequence
            char *ssSeq Second sequence
            PDB  *pdb   PDB linked list
   Returns: BOOL        TRUE if the 3 sequences are the same

   Takes 2 sequences and a PDB linked list. Extracts the sequence from 
   the PDB linked list ATOM records and then compares the 3 sequences.
   If they are the same returns TRUE, else FALSE

   25.02.10  Original   By: ACRM
*/
BOOL CheckData(FILE *fp, char *seq, char *ssSeq, PDB *pdb)
{
   char *cleanSeq,
        *cleanSSSeq,
        *cleanPDBSeq;
   BOOL retval = TRUE;

   cleanSeq    = CleanSeqMalloc(seq);
   cleanSSSeq  = CleanSeqMalloc(ssSeq);
   cleanPDBSeq = PDB2Seq(pdb);

   if((cleanSeq    != NULL) &&
      (cleanSSSeq  != NULL) &&
      (cleanPDBSeq != NULL))
   {
      if(strcmp(cleanSeq, cleanSSSeq) ||
         strcmp(cleanSeq, cleanPDBSeq))
      {
         if(fp!=NULL)
         {
            fprintf(fp,"Inconsistent sequence data...\n");
            fprintf(fp,"Alignment file:\n%s\n", cleanSeq);
            fprintf(fp,"Secondary structure file:\n%s\n", cleanSSSeq);
            fprintf(fp,"PDB file (ATOM sequence):\n%s\n", cleanPDBSeq);
         }
         retval = FALSE;
      }
   }
   
   if(cleanSeq    != NULL) free(cleanSeq);
   if(cleanSSSeq  != NULL) free(cleanSSSeq);
   if(cleanPDBSeq != NULL) free(cleanPDBSeq);

   return(retval);
}


/************************************************************************/
/*>BOOL FixAlignment(char *seq1, char *seq2, char *ss, PDB *pdb, 
                     BOOL verbose, int flex)
   --------------------------------------------------------------
   Input:   char  *seq1   First aligned sequence (template)
            char  *seq2   Second aligned sequence (target)
            char  *ss     Secondary structure assignments for template
            PDB   *pdb    PDB linked list for template. 
                          NOTE! Must be C-alphas only
            BOOL  verbose Verbose printing flag
            int   flex    Extend loop boundaries by this number of 
                          residues
   Returns: BOOL          Success?

   Does the actual work of fixing an alignment. Identifies gaps in the
   target sequence which occur in (or near to) loops. The optimal site
   within the loop for closing the gap is then identified and the
   alignment is modified.

   25.02.10  Original   By: ACRM
*/
BOOL FixAlignment(char *seq1, char *seq2, char *ss, PDB *pdb, 
                  BOOL verbose, int flex)
{
   PDB **pdbindex;
   int natoms, alnlength, i, start, stop, gaplen, bestpos;

   /* Index the PDB C-alpha linked list                                 */
   if((pdbindex = IndexPDB(pdb, &natoms)) == NULL)
      return(FALSE);

   /* Work through the target sequence                                  */
   alnlength = strlen(seq2);
   for(i=0; i<alnlength; i++)
   {
      /* If we have a gap in the target sequence                        */
      if(seq2[i] == '-')
      {
         /* See how long the gap is                                     */
         gaplen = FindGapLen(i, seq2);

         /* See if it's in a loop, or within flex residues of a loop. If
            it is, then start and stop will be the boundaries of the loop
            (extended to include this position if necessary - i.e. if it
            was within flex residues of the loop)
         */
         if(InLoop(i, gaplen, ss, flex, &start, &stop))
         {
            fprintf(stderr, "Gap of length %d at position %d found in \
loop. Bounds %d to %d\n", 
                   gaplen, i, start, stop);
            
            /* Find the best place within the loop to place the gap     */
            if((bestpos = FindBestGapPosition(pdbindex, i, start, stop, 
                                              gaplen, verbose))==(-1))
            {
               fprintf(stderr,"   Gap not processed\n");
            }
            else
            {
               fprintf(stderr,"   Best position for this gap is at %d\n", 
                       bestpos);

               /* Move the gap to the right place                       */
               MoveGap(seq2, i, gaplen, bestpos);
            }

            /* Skip to end of this gap                                  */
            i+=gaplen;
         }
      }
   }
   free(pdbindex);
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL InLoop(int pos, int gaplen, char *ss, int flex, 
               int *start, int *stop)
   ----------------------------------------------------
   Input:   int  pos    Position within the sequence
            int  gaplen Length of gap
            char *ss    Secondary structure sequence string
            int  flex   How many residues away from a loop we count as
                        being allowable for moving a gap into the loop
   Output:  int  *start Start of region defined as loop
            int  *stop  End of region defined as loop
   Returns: BOOL        Was pos in a loop?

   Examines a secondary structure string and determines whether the
   specified position (pos) is in a loop or within flex residues of a 
   loop. If so, then the boundaries of the loop (extended if necessary
   because the residue was within flex residues of the loop) are
   returned.

   25.02.10  Original   By: ACRM
*/
BOOL InLoop(int pos, int gaplen, char *ss, int flex, int *start, 
            int *stop)
{
   BOOL inLoop = FALSE;
   int  i;
   
   /* Define boundaries of the region in which we will look             */
   *start = pos - flex;
   *stop  = pos + flex;
   if(*start < 0)          *start = 0;
   if(*stop  > strlen(ss)) *stop = strlen(ss)-1;

   /* Search this region to see if it contains a loop character         */
   for(i=*start; i<*stop; i++)
   {
      if(ss[i] == 'C')
      {
         inLoop = TRUE;
         break;
      }
   }

   /* If it does, find the boundaries of the loop                       */
   if(inLoop)
   {
      *start = *stop = i;
      while((ss[*start] == 'C') && (*start > 0)) (*start)--;
      while((ss[*stop]  == 'C') && (*stop < strlen(ss))) (*stop)++;
      
      /* If the gap was out of the loop itself then extend the region to 
         include the gap
      */
      if(*start > pos) *start = pos;
      if(*stop  < pos) *stop  = pos+gaplen;
   }

   return(inLoop);
}


/************************************************************************/
/*>int FindGapLen(int i, char *seq)
   --------------------------------
   Input:   int  i     Position within sequence of start of a gap
            char *seq  Sequence (containing gaps)
   Returns: int        Length of gap starting at position i

   Calculates the length of a gap in a sequence (i.e. counts '-' 
   characters)

   25.02.10  Original   By: ACRM
*/
int FindGapLen(int i, char *seq)
{
   int gaplen = 0;
   while(seq[i++] == '-') gaplen++;
   return(gaplen);
}


/************************************************************************/
/*>int FindBestGapPosition(PDB **pdbindex, int pos, int start, int stop, 
                           int gaplen, BOOL verbose)
   ------------------------------------------------------------
   Input:   PDB  **pdbindex   Array of PDB pointers
            int  pos          Old position of gap
            int  start        Offset into array of first residue of loop
            int  stop         Offset into array of last residue of loop
            int  gaplen       Length of gap
            BOOL verbose      Print verbose messages
   Returns: int               Best position for gap (-1 if error)

   Looks for the best position to place a gap in a loop by looking for 
   the closure which will result in minimal structural disruption - 
   i.e. the gap to be closed is as close to the standard CA-CA distance 
   as possible.

   25.02.10  Original   By: ACRM
*/
int FindBestGapPosition(PDB **pdbindex, int pos, int start, int stop, 
                        int gaplen, BOOL verbose)
{
   REAL gapdist;
   REAL bestDist;
   int  i, bestI, count;
   
   /* First check integrity of loop - i.e. are all residues present     */
   for(i=start; i<stop-1; i++)
   {
      if(DISTSQ(pdbindex[i], pdbindex[i+1]) > (CADISTSQ + 1))
      {
         fprintf(stderr,"   PDB file has missing residues between %c%d%c \
and %c%d%c\n",
                 (pdbindex[i])->chain[0],
                 (pdbindex[i])->resnum,
                 (pdbindex[i])->insert[0],
                 (pdbindex[i+1])->chain[0],
                 (pdbindex[i+1])->resnum,
                 (pdbindex[i+1])->insert[0]);
         return(-1);
      }
   }

   /* Now find the distances between residues separated by gaplen       */
   bestDist = 1000.0;
   bestI    = pos; /* Assume initial position is the best               */
   count    = 0;
   
   /* Run through the loop looking at the CA-CA separation that a gap of
      the required length would introduce. Identify the one closes to
      the standard CA-CA distance
   */
   for(i=start; i<stop-gaplen; i++)
   {
      gapdist = DIST(pdbindex[i], pdbindex[i+gaplen+1]) - OPTCACADIST;
      gapdist = ABS(gapdist);

      if(gapdist < bestDist)
      {
         bestDist = gapdist;
         bestI = i+1; /* Note that we add 1 as the gap must be 
                         introduced *after* this residue 
                      */
      }

      if(verbose)
      {
         fprintf(stderr, "Gap at %d - penalty %f\n", i+1, gapdist);
      }
   }
   
   return(bestI);
}


/************************************************************************/
/*>BOOL MoveGap(char *seq, int oldpos, int gaplen, int newpos)
   -----------------------------------------------------------
   I/O:     char  *seq      Sequence (including gaps) - modified on
                            completion
   Input:   int   oldpos    Old gap position
            int   gaplen    Length of gap
            int   newpos    New (best) position for gap
   Returns: BOOL            success

   Takes a sequence containing a gap of length gaplen at position oldpos
   and moves the gap to position bestpos. The sequence as passed into
   the routine is modified in the process.

   25.02.10  Original   By: ACRM
*/
BOOL MoveGap(char *seq, int oldpos, int gaplen, int newpos)
{
   char *newseq;
   int  i,j,k,alnlen;
   
   if(oldpos == newpos)
   {
      return(TRUE);
   }
   
   /* Allocate memory for a copy of the sequence                        */
   alnlen=strlen(seq);
   if((newseq = (char *)malloc((alnlen+1) * sizeof(char)))==NULL)
   {
      return(FALSE);
   }

   /* Copy the sequence removing the old gap                           */
   for(i=0, j=0; i<alnlen; i++)
   {
      if((i<oldpos) || (i>=(oldpos+gaplen)))
      {
         newseq[j++] = seq[i];
      }
   }
   newseq[j] ='\0';
   
   /* Copy the sequence back again, inserting the new gap               */
   alnlen=strlen(newseq);
   for(i=0, j=0; i<alnlen; i++)
   {
      if(i==newpos)
      {
         for(k=0; k<gaplen; k++)
         {
            seq[j+k] = '-';
         }
         j=j+k;
      }
      seq[j++] = newseq[i];
   }
   seq[j] = '\0';

   return(TRUE);
}


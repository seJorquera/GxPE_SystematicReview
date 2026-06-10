# GxPE frequency calculation in GnomAD v4.0
**Author**: Samuel Jorquera

In this repository you can find the reproducible script for absolute frequency calculation across population using GnomAD v4.0 (```gnomad_r4```). Briefly, we extracted the data using graphiQL API. Then, we obtained the total allele number (AN) and alternative allele (AC) number for each population, and calculated absolute frequency (AF) as AC/AN.

This is a reproducible pipeline, however does not receive an input of variants for study. ```GnomAD_API_AF_Calculation.R``` script can be easily modified for the study of other variants using gnomAD SNP ID structure: CHROM-POS-Ref-Alt (e.g. 4-23795377-C-T for rs6821591).
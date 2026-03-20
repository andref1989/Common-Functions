
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin22.4.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> setwd('/Users/forbesa/Andre_F_functions_git/')
> RYR2_data <- read.table("/Users/forbesa/Downloads/clinvar_result.txt")
Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  line 2 did not have 52 elements
> RYR2_data <- read.table("/Users/forbesa/Downloads/clinvar_result.txt",header=T, sep='\t')
> str(RYR2_data)
'data.frame':	156 obs. of  25 variables:
 $ Name                                       : chr  "NM_001035.3(RYR2):c.28G>T (p.Glu10Ter)" "NM_001035.3(RYR2):c.44G>C (p.Arg15Pro)" "NM_001035.3(RYR2):c.139T>C (p.Cys47Arg)" "NM_001035.3(RYR2):c.230C>T (p.Ala77Val)" ...
 $ Gene.s.                                    : chr  "RYR2" "RYR2" "RYR2" "RYR2" ...
 $ Protein.change                             : chr  "E10*" "R15P" "C47R" "A77V" ...
 $ Condition.s.                               : chr  "not provided" "Catecholaminergic polymorphic ventricular tachycardia 1" "Catecholaminergic polymorphic ventricular tachycardia 1" "not provided|Catecholaminergic polymorphic ventricular tachycardia 1|Cardiovascular phenotype" ...
 $ Accession                                  : chr  "VCV001678118" "VCV000463600" "VCV001805396" "VCV000404190" ...
 $ GRCh37Chromosome                           : int  1 1 1 1 1 1 1 1 1 1 ...
 $ GRCh37Location                             : chr  "237205849" "237205865" "237433887" "237494239" ...
 $ GRCh38Chromosome                           : int  1 1 1 1 1 1 1 1 1 1 ...
 $ GRCh38Location                             : chr  "237042549" "237042565" "237270587" "237330939" ...
 $ VariationID                                : int  1678118 463600 1805396 404190 684806 572162 3776659 201193 427184 1745086 ...
 $ AlleleID.s.                                : int  1669807 447603 1862408 391039 672379 556831 3892806 196511 414779 1805932 ...
 $ dbSNP.ID                                   : chr  "rs2148061633" "rs865784613" "rs2528051422" "rs1060500142" ...
 $ Canonical.SPDI                             : chr  "NC_000001.11:237042548:G:T" "NC_000001.11:237042564:G:C" "NC_000001.11:237270586:T:C" "NC_000001.11:237330938:C:T" ...
 $ Variant.type                               : chr  "single nucleotide variant" "single nucleotide variant" "single nucleotide variant" "single nucleotide variant" ...
 $ Molecular.consequence                      : chr  "nonsense" "missense variant" "missense variant" "missense variant" ...
 $ Germline.classification                    : chr  "Likely pathogenic" "Likely pathogenic" "Likely pathogenic" "Pathogenic/Likely pathogenic" ...
 $ Germline.date.last.evaluated               : chr  "Oct 12, 2021" "Nov 14, 2022" "Mar 31, 2022" "Nov 21, 2025" ...
 $ Germline.review.status                     : chr  "criteria provided, single submitter" "criteria provided, single submitter" "criteria provided, single submitter" "criteria provided, multiple submitters, no conflicts" ...
 $ Somatic.clinical.impact                    : logi  NA NA NA NA NA NA ...
 $ Somatic.clinical.impact.date.last.evaluated: logi  NA NA NA NA NA NA ...
 $ Somatic.clinical.impact.review.status      : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.classification                : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.date.last.evaluated           : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.review.status                 : logi  NA NA NA NA NA NA ...
 $ X                                          : logi  NA NA NA NA NA NA ...
> RYR2_hg19_bed <- dplyr::select(RYR2_data,chr=GRCh37Chromosome,start=GRCh37Location,end=GRCh37Location,functional_impact=Germline.classification,molecular_classification=Molecular.consequence, AA_change=Protein.change,Clinical_Outcome=Condition.s.,dbSNP=dbSNP.ID)
> RYR2_hg38_bed <- dplyr::select(RYR2_data,chr=GRCh38Chromosome,start=GRCh38Location,end=GRCh378Location,functional_impact=Germline.classification,molecular_classification=Molecular.consequence, AA_change=Protein.change,Clinical_Outcome=Condition.s.,dbSNP=dbSNP.ID)
Error in `dplyr::select()`:
! Can't select columns that don't exist.
✖ Column `GRCh378Location` doesn't exist.
> RYR2_hg38_bed <- dplyr::select(RYR2_data,chr=GRCh38Chromosome,start=GRCh38Location,end=GRCh38Location,functional_impact=Germline.classification,molecular_classification=Molecular.consequence, AA_change=Protein.change,Clinical_Outcome=Condition.s.,dbSNP=dbSNP.ID)
> ANF_genome <- read.table("/Users/forbesa/Downloads/genome_Andre_Forbes_v4_Full_20170131080757.txt", header=T, sep='\t')
> str(ANF_genome)
'data.frame':	610525 obs. of  4 variables:
 $ rs12564807: chr  "rs3131972" "rs148828841" "rs12124819" "rs115093905" ...
 $ X1        : chr  "1" "1" "1" "1" ...
 $ X734462   : int  752721 760998 776546 787173 798959 824398 838555 846808 854250 861808 ...
 $ AA        : chr  "AA" "CC" "AA" "--" ...
> ANF_genome <- read.table("/Users/forbesa/Downloads/genome_Andre_Forbes_v4_Full_20170131080757.txt", header=F, sep='\t') 
> ANF_genome <- ANF_genome[c(1,2,3,3,4)]
> colnames(ANF_genome) <- c("dnSNP","chr","start","end","Genotype")
> head(ANF_genome)
        dnSNP chr  start    end Genotype
1  rs12564807   1 734462 734462       AA
2   rs3131972   1 752721 752721       AA
3 rs148828841   1 760998 760998       CC
4  rs12124819   1 776546 776546       AA
5 rs115093905   1 787173 787173       --
6  rs11240777   1 798959 798959       AA
> ANF_genome <- ANF_genome[c(2:5,1)]
> head(ANF_genome)
  chr  start    end Genotype       dnSNP
1   1 734462 734462       AA  rs12564807
2   1 752721 752721       AA   rs3131972
3   1 760998 760998       CC rs148828841
4   1 776546 776546       AA  rs12124819
5   1 787173 787173       -- rs115093905
6   1 798959 798959       AA  rs11240777
> source("~/Andre_F_functions_git/Andre_F_functions.R")
> ANF_genome_gr <- bed_to_granges_dynamic(ANF_genome)
Loading required package: GenomicRanges
Loading required package: stats4
Loading required package: BiocGenerics

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    setdiff, sort, table, tapply, union, unique, unsplit,
    which.max, which.min

Loading required package: S4Vectors

Attaching package: ‘S4Vectors’

The following object is masked from ‘package:utils’:

    findMatches

The following objects are masked from ‘package:base’:

    expand.grid, I, unname

Loading required package: IRanges
Loading required package: GenomeInfoDb
Error in read.table(bed, stringsAsFactors = F, sep = "\t", header = header,  : 
  'file' must be a character string or connection
> ANF_genome_gr <- table_to_granges(ANF_genome)
> ANF_genome_gr
GRanges object with 610526 ranges and 2 metadata columns:
           seqnames    ranges strand |    Genotype       dnSNP
              <Rle> <IRanges>  <Rle> | <character> <character>
       [1]        1    734462      * |          AA  rs12564807
       [2]        1    752721      * |          AA   rs3131972
       [3]        1    760998      * |          CC rs148828841
       [4]        1    776546      * |          AA  rs12124819
       [5]        1    787173      * |          -- rs115093905
       ...      ...       ...    ... .         ...         ...
  [610522]       MT     16524      * |           A    i4000693
  [610523]       MT     16526      * |           G    i4000757
  [610524]       MT     16527      * |           T    i4990307
  [610525]       MT     16540      * |           C    i4000756
  [610526]       MT     16547      * |           C    i3001931
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> RYR2_hg19_gr <- table_to_granges(RYR2_hg19_bed)
Error in .new_IRanges_from_start_end(start, end) : 
  'start' or 'end' cannot contain NAs
In addition: Warning messages:
1: In IRanges(start = as.numeric(table[, 2]), end = as.numeric(table[,  :
  NAs introduced by coercion
2: In IRanges(start = as.numeric(table[, 2]), end = as.numeric(table[,  :
  NAs introduced by coercion
> str(RYR2_hg19_bed)
'data.frame':	156 obs. of  8 variables:
 $ chr                     : int  1 1 1 1 1 1 1 1 1 1 ...
 $ start                   : chr  "237205849" "237205865" "237433887" "237494239" ...
 $ end                     : chr  "237205849" "237205865" "237433887" "237494239" ...
 $ functional_impact       : chr  "Likely pathogenic" "Likely pathogenic" "Likely pathogenic" "Pathogenic/Likely pathogenic" ...
 $ molecular_classification: chr  "nonsense" "missense variant" "missense variant" "missense variant" ...
 $ AA_change               : chr  "E10*" "R15P" "C47R" "A77V" ...
 $ Clinical_Outcome        : chr  "not provided" "Catecholaminergic polymorphic ventricular tachycardia 1" "Catecholaminergic polymorphic ventricular tachycardia 1" "not provided|Catecholaminergic polymorphic ventricular tachycardia 1|Cardiovascular phenotype" ...
 $ dbSNP                   : chr  "rs2148061633" "rs865784613" "rs2528051422" "rs1060500142" ...
> RYR2_hg19_bed$start <- as.integer(RYR2_hg19_bed$start)
Warning message:
NAs introduced by coercion 
> RYR2_hg19_bed$end <- as.integer(RYR2_hg19_bed$end)
Warning message:
NAs introduced by coercion 
> RYR2_hg19_gr <- table_to_granges(dplyr::filter(RYR2_hg19_bed, !is.na(start)))
> RYR2_hg19_gr
GRanges object with 154 ranges and 5 metadata columns:
        seqnames    ranges strand |      functional_impact
           <Rle> <IRanges>  <Rle> |            <character>
    [1]        1 237205849      * |      Likely pathogenic
    [2]        1 237205865      * |      Likely pathogenic
    [3]        1 237433887      * |      Likely pathogenic
    [4]        1 237494239      * | Pathogenic/Likely pa..
    [5]        1 237494252      * |      Likely pathogenic
    ...      ...       ...    ... .                    ...
  [150]        1 237995888      * |      Likely pathogenic
  [151]        1 237995892      * |      Likely pathogenic
  [152]        1 237995904      * |      Likely pathogenic
  [153]        1 237995919      * | Pathogenic/Likely pa..
  [154]        1 237995928      * | Pathogenic/Likely pa..
        molecular_classification   AA_change       Clinical_Outcome
                     <character> <character>            <character>
    [1]                 nonsense        E10*           not provided
    [2]         missense variant        R15P Catecholaminergic po..
    [3]         missense variant        C47R Catecholaminergic po..
    [4]         missense variant        A77V not provided|Catecho..
    [5]         missense variant        M81I Conduction disorder ..
    ...                      ...         ...                    ...
  [150]         missense variant      W4949R Catecholaminergic po..
  [151]         missense variant      E4950G           not provided
  [152]         missense variant      A4954G Cardiovascular pheno..
  [153]         missense variant      R4959Q Cardiovascular pheno..
  [154]         missense variant      Y4962C Cardiovascular pheno..
               dbSNP
         <character>
    [1] rs2148061633
    [2]  rs865784613
    [3] rs2528051422
    [4] rs1060500142
    [5] rs1572627115
    ...          ...
  [150]  rs794728810
  [151] rs1057517873
  [152] rs1663930286
  [153]  rs794728811
  [154]  rs794728832
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
> intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr)
GRanges object with 25 ranges and 2 metadata columns:
       seqnames    ranges strand |    Genotype       dnSNP
          <Rle> <IRanges>  <Rle> | <character> <character>
   [1]        1 237540665      * |          GG    i6017887
   [2]        1 237540665      * |          GG    i6017887
   [3]        1 237540686      * |          GG    i6017888
   [4]        1 237551437      * |          GG    i6017892
   [5]        1 237608774      * |          CC    i6017854
   ...      ...       ...    ... .         ...         ...
  [21]        1 237954780      * |          GG    i6017900
  [22]        1 237954780      * |          GG    i6017900
  [23]        1 237972213      * |          GG    i6017915
  [24]        1 237993857      * |          AA    i6017871
  [25]        1 237995919      * |          GG    i6017906
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> RYR2_hg38_bed$end <- as.integer(RYR2_hg38_bed$end);RYR2_hg38_bed$start <- as.integer(RYR2_hg38_bed$start);RYR2_hg38_gr <- table_to_granges(dplyr::filter(RYR2_hg38_bed, !is.na(start)))
Warning message:
NAs introduced by coercion 
Warning message:
NAs introduced by coercion 
> intersect_with_metadata(ANF_genome_gr,RYR2_hg38_gr)
GRanges object with 1 range and 2 metadata columns:
      seqnames    ranges strand |    Genotype       dnSNP
         <Rle> <IRanges>  <Rle> | <character> <character>
  [1]        1 237270587      * |          CT   rs2490346
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr)
GRanges object with 25 ranges and 2 metadata columns:
       seqnames    ranges strand |    Genotype       dnSNP
          <Rle> <IRanges>  <Rle> | <character> <character>
   [1]        1 237540665      * |          GG    i6017887
   [2]        1 237540665      * |          GG    i6017887
   [3]        1 237540686      * |          GG    i6017888
   [4]        1 237551437      * |          GG    i6017892
   [5]        1 237608774      * |          CC    i6017854
   ...      ...       ...    ... .         ...         ...
  [21]        1 237954780      * |          GG    i6017900
  [22]        1 237954780      * |          GG    i6017900
  [23]        1 237972213      * |          GG    i6017915
  [24]        1 237993857      * |          AA    i6017871
  [25]        1 237995919      * |          GG    i6017906
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> unique(intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr))
GRanges object with 20 ranges and 2 metadata columns:
       seqnames    ranges strand |    Genotype       dnSNP
          <Rle> <IRanges>  <Rle> | <character> <character>
   [1]        1 237540665      * |          GG    i6017887
   [2]        1 237540686      * |          GG    i6017888
   [3]        1 237551437      * |          GG    i6017892
   [4]        1 237608774      * |          CC    i6017854
   [5]        1 237632425      * |          CC    i6017934
   ...      ...       ...    ... .         ...         ...
  [16]        1 237947545      * |          AA    i6017849
  [17]        1 237954780      * |          GG    i6017900
  [18]        1 237972213      * |          GG    i6017915
  [19]        1 237993857      * |          AA    i6017871
  [20]        1 237995919      * |          GG    i6017906
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> unique(intersect_with_metadata(ANF_genome_gr,gr1=RYR2_hg19_gr))
GRanges object with 20 ranges and 5 metadata columns:
       seqnames    ranges strand |      functional_impact
          <Rle> <IRanges>  <Rle> |            <character>
   [1]        1 237540665      * | Pathogenic/Likely pa..
   [2]        1 237540686      * |      Likely pathogenic
   [3]        1 237551437      * |      Likely pathogenic
   [4]        1 237608774      * |      Likely pathogenic
   [5]        1 237632425      * |      Likely pathogenic
   ...      ...       ...    ... .                    ...
  [16]        1 237947545      * | Pathogenic/Likely pa..
  [17]        1 237954780      * | Pathogenic/Likely pa..
  [18]        1 237972213      * |      Likely pathogenic
  [19]        1 237993857      * |      Likely pathogenic
  [20]        1 237995919      * | Pathogenic/Likely pa..
       molecular_classification   AA_change       Clinical_Outcome
                    <character> <character>            <character>
   [1]         missense variant       R169L not provided|Cardiov..
   [2]         missense variant       R176L Cardiovascular pheno..
   [3]                 nonsense       E243*  RYR2-related disorder
   [4]         missense variant       T415I Catecholaminergic po..
   [5]         missense variant       A549V Arrhythmogenic right..
   ...                      ...         ...                    ...
  [16]         missense variant      N4178S Cardiovascular pheno..
  [17]         missense variant      A4510P not provided|Catecho..
  [18]         missense variant      V4771F Cardiovascular pheno..
  [19]         missense variant      N4895H Catecholaminergic po..
  [20]         missense variant      R4959Q Cardiovascular pheno..
              dbSNP
        <character>
   [1]  rs397516539
   [2]  rs794728708
   [3]  rs794728712
   [4] rs1288202574
   [5]  rs794728817
   ...          ...
  [16]  rs794728787
  [17]  rs397516510
  [18]  rs794728804
  [19] rs1185619003
  [20]  rs794728811
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
> gr_to_bed(unique(intersect_with_metadata(ANF_genome_gr,gr1=RYR2_hg19_gr)),metadata=T)
   Chrom     Start       End Strand            functional_impact
1      1 237540665 237540666      * Pathogenic/Likely pathogenic
2      1 237540686 237540687      *            Likely pathogenic
3      1 237551437 237551438      *            Likely pathogenic
4      1 237608774 237608775      *            Likely pathogenic
5      1 237632425 237632426      *            Likely pathogenic
6      1 237656273 237656274      *            Likely pathogenic
7      1 237777598 237777599      * Pathogenic/Likely pathogenic
8      1 237796969 237796970      *            Likely pathogenic
9      1 237801780 237801781      *            Likely pathogenic
10     1 237804240 237804241      * Pathogenic/Likely pathogenic
11     1 237804256 237804257      *            Likely pathogenic
12     1 237804283 237804284      * Pathogenic/Likely pathogenic
13     1 237942027 237942028      *            Likely pathogenic
14     1 237947458 237947459      * Pathogenic/Likely pathogenic
15     1 237947482 237947483      * Pathogenic/Likely pathogenic
16     1 237947545 237947546      * Pathogenic/Likely pathogenic
17     1 237954780 237954781      * Pathogenic/Likely pathogenic
18     1 237972213 237972214      *            Likely pathogenic
19     1 237993857 237993858      *            Likely pathogenic
20     1 237995919 237995920      * Pathogenic/Likely pathogenic
   molecular_classification AA_change
1          missense variant     R169L
2          missense variant     R176L
3                  nonsense     E243*
4          missense variant     T415I
5          missense variant     A549V
6          missense variant     S616L
7          missense variant    E1724K
8          missense variant    D2216G
9          missense variant    V2306F
10         missense variant    A2387T
11         missense variant    Y2392F
12         missense variant    R2401L
13         missense variant    G3946V
14         missense variant    Y4149C
15         missense variant    R4157Q
16         missense variant    N4178S
17         missense variant    A4510P
18         missense variant    V4771F
19         missense variant    N4895H
20         missense variant    R4959Q
                                                                                                                                      Clinical_Outcome
1                                                                                                                not provided|Cardiovascular phenotype
2                                                                                                                Cardiovascular phenotype|not provided
3                                                                                                                                RYR2-related disorder
4                                                                                              Catecholaminergic polymorphic ventricular tachycardia 1
5                                                 Arrhythmogenic right ventricular dysplasia 2|Catecholaminergic polymorphic ventricular tachycardia 1
6                                                                                 not provided|Catecholaminergic polymorphic ventricular tachycardia 1
7                                       Cardiovascular phenotype|not provided|Long QT syndrome|Catecholaminergic polymorphic ventricular tachycardia 1
8                                                                                              Catecholaminergic polymorphic ventricular tachycardia 1
9                                                                                              Catecholaminergic polymorphic ventricular tachycardia 1
10 Cardiovascular phenotype|Catecholaminergic polymorphic ventricular tachycardia|Catecholaminergic polymorphic ventricular tachycardia 1|not provided
11                                                                                               Catecholaminergic polymorphic ventricular tachycardia
12                                                                                not provided|Catecholaminergic polymorphic ventricular tachycardia 1
13                                                                                                                            Cardiovascular phenotype
14                                                                    Cardiovascular phenotype|Catecholaminergic polymorphic ventricular tachycardia 1
15                                                                                not provided|Catecholaminergic polymorphic ventricular tachycardia 1
16                                                       Cardiovascular phenotype|not provided|Catecholaminergic polymorphic ventricular tachycardia 1
17                                                                                not provided|Catecholaminergic polymorphic ventricular tachycardia 1
18                                                                    Cardiovascular phenotype|Catecholaminergic polymorphic ventricular tachycardia 1
19                                                                                             Catecholaminergic polymorphic ventricular tachycardia 1
20                                                       Cardiovascular phenotype|not provided|Catecholaminergic polymorphic ventricular tachycardia 1
          dbSNP
1   rs397516539
2   rs794728708
3   rs794728712
4  rs1288202574
5   rs794728817
6   rs730880187
7   rs794728740
8  rs1328318082
9   rs794728746
10  rs794728753
11             
12  rs794728756
13             
14 rs1234449785
15  rs794728786
16  rs794728787
17  rs397516510
18  rs794728804
19 rs1185619003
20  rs794728811
> RYR2_hg19_bed$end <- as.integer(RYR2_hg19_bed$end)  C-c C-c
> colnames(RYR2_hg19_bed)
[1] "chr"                      "start"                   
[3] "end"                      "functional_impact"       
[5] "molecular_classification" "AA_change"               
[7] "Clinical_Outcome"         "dbSNP"                   
> RYR2_hg19_gr <- table_to_granges(dplyr::filter(RYR2_hg19_bed, !is.na(start)))  C-c C-c
>RYR2_hg19_bed <- dplyr::select(RYR2_data,chr=GRCh37Chromosome,start=GRCh37Location,end=GRCh37Location,functional_impact=Germline.classification,molecular_classification=Molecular.consequence, AA_change=Protein.change,Clinical_Outcome=Condition.s.,dbSNP=dbSNP.ID)  C-c C-c
> str(RYR2_data)
'data.frame':	156 obs. of  25 variables:
 $ Name                                       : chr  "NM_001035.3(RYR2):c.28G>T (p.Glu10Ter)" "NM_001035.3(RYR2):c.44G>C (p.Arg15Pro)" "NM_001035.3(RYR2):c.139T>C (p.Cys47Arg)" "NM_001035.3(RYR2):c.230C>T (p.Ala77Val)" ...
 $ Gene.s.                                    : chr  "RYR2" "RYR2" "RYR2" "RYR2" ...
 $ Protein.change                             : chr  "E10*" "R15P" "C47R" "A77V" ...
 $ Condition.s.                               : chr  "not provided" "Catecholaminergic polymorphic ventricular tachycardia 1" "Catecholaminergic polymorphic ventricular tachycardia 1" "not provided|Catecholaminergic polymorphic ventricular tachycardia 1|Cardiovascular phenotype" ...
 $ Accession                                  : chr  "VCV001678118" "VCV000463600" "VCV001805396" "VCV000404190" ...
 $ GRCh37Chromosome                           : int  1 1 1 1 1 1 1 1 1 1 ...
 $ GRCh37Location                             : chr  "237205849" "237205865" "237433887" "237494239" ...
 $ GRCh38Chromosome                           : int  1 1 1 1 1 1 1 1 1 1 ...
 $ GRCh38Location                             : chr  "237042549" "237042565" "237270587" "237330939" ...
 $ VariationID                                : int  1678118 463600 1805396 404190 684806 572162 3776659 201193 427184 1745086 ...
 $ AlleleID.s.                                : int  1669807 447603 1862408 391039 672379 556831 3892806 196511 414779 1805932 ...
 $ dbSNP.ID                                   : chr  "rs2148061633" "rs865784613" "rs2528051422" "rs1060500142" ...
 $ Canonical.SPDI                             : chr  "NC_000001.11:237042548:G:T" "NC_000001.11:237042564:G:C" "NC_000001.11:237270586:T:C" "NC_000001.11:237330938:C:T" ...
 $ Variant.type                               : chr  "single nucleotide variant" "single nucleotide variant" "single nucleotide variant" "single nucleotide variant" ...
 $ Molecular.consequence                      : chr  "nonsense" "missense variant" "missense variant" "missense variant" ...
 $ Germline.classification                    : chr  "Likely pathogenic" "Likely pathogenic" "Likely pathogenic" "Pathogenic/Likely pathogenic" ...
 $ Germline.date.last.evaluated               : chr  "Oct 12, 2021" "Nov 14, 2022" "Mar 31, 2022" "Nov 21, 2025" ...
 $ Germline.review.status                     : chr  "criteria provided, single submitter" "criteria provided, single submitter" "criteria provided, single submitter" "criteria provided, multiple submitters, no conflicts" ...
 $ Somatic.clinical.impact                    : logi  NA NA NA NA NA NA ...
 $ Somatic.clinical.impact.date.last.evaluated: logi  NA NA NA NA NA NA ...
 $ Somatic.clinical.impact.review.status      : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.classification                : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.date.last.evaluated           : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.review.status                 : logi  NA NA NA NA NA NA ...
 $ X                                          : logi  NA NA NA NA NA NA ...
> RYR2_hg19_bed <- dplyr::select(RYR2_data,chr=GRCh37Chromosome,start=GRCh37Location,end=GRCh37Location,functional_impact=Germline.classification,molecular_classification=Molecular.consequence, AA_change=Protein.change,Clinical_Outcome=Condition.s.,dbSNP=dbSNP.ID, Nucleotide_Change=Name)
> RYR2_hg19_bed$end <- as.integer(RYR2_hg19_bed$end);RYR2_hg19_bed$start <- as.integer(RYR2_hg19_bed$start);RYR2_hg19_gr <- table_to_granges(dplyr::filter(RYR2_hg19_bed, !is.na(start)))
Warning message:
NAs introduced by coercion 
Warning message:
NAs introduced by coercion 
> intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr)
GRanges object with 25 ranges and 2 metadata columns:
       seqnames    ranges strand |    Genotype       dnSNP
          <Rle> <IRanges>  <Rle> | <character> <character>
   [1]        1 237540665      * |          GG    i6017887
   [2]        1 237540665      * |          GG    i6017887
   [3]        1 237540686      * |          GG    i6017888
   [4]        1 237551437      * |          GG    i6017892
   [5]        1 237608774      * |          CC    i6017854
   ...      ...       ...    ... .         ...         ...
  [21]        1 237954780      * |          GG    i6017900
  [22]        1 237954780      * |          GG    i6017900
  [23]        1 237972213      * |          GG    i6017915
  [24]        1 237993857      * |          AA    i6017871
  [25]        1 237995919      * |          GG    i6017906
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> intersect_with_metadata(ANF_genome_gr,gr1=RYR2_hg19_gr)
GRanges object with 25 ranges and 6 metadata columns:
       seqnames    ranges strand |      functional_impact
          <Rle> <IRanges>  <Rle> |            <character>
   [1]        1 237540665      * | Pathogenic/Likely pa..
   [2]        1 237540665      * | Pathogenic/Likely pa..
   [3]        1 237540686      * |      Likely pathogenic
   [4]        1 237551437      * |      Likely pathogenic
   [5]        1 237608774      * |      Likely pathogenic
   ...      ...       ...    ... .                    ...
  [21]        1 237954780      * | Pathogenic/Likely pa..
  [22]        1 237954780      * |      Likely pathogenic
  [23]        1 237972213      * |      Likely pathogenic
  [24]        1 237993857      * |      Likely pathogenic
  [25]        1 237995919      * | Pathogenic/Likely pa..
       molecular_classification   AA_change       Clinical_Outcome
                    <character> <character>            <character>
   [1]         missense variant       R169L not provided|Cardiov..
   [2]         missense variant       R169Q not provided|Catecho..
   [3]         missense variant       R176L Cardiovascular pheno..
   [4]                 nonsense       E243*  RYR2-related disorder
   [5]         missense variant       T415I Catecholaminergic po..
   ...                      ...         ...                    ...
  [21]         missense variant      A4510P not provided|Catecho..
  [22]         missense variant      A4510S not provided|Arrhyth..
  [23]         missense variant      V4771F Cardiovascular pheno..
  [24]         missense variant      N4895H Catecholaminergic po..
  [25]         missense variant      R4959Q Cardiovascular pheno..
              dbSNP      Nucleotide_Change
        <character>            <character>
   [1]  rs397516539 NM_001035.3(RYR2):c...
   [2]  rs397516539 NM_001035.3(RYR2):c...
   [3]  rs794728708 NM_001035.3(RYR2):c...
   [4]  rs794728712 NM_001035.3(RYR2):c...
   [5] rs1288202574 NM_001035.3(RYR2):c...
   ...          ...                    ...
  [21]  rs397516510 NM_001035.3(RYR2):c...
  [22]  rs397516510 NM_001035.3(RYR2):c...
  [23]  rs794728804 NM_001035.3(RYR2):c...
  [24] rs1185619003 NM_001035.3(RYR2):c...
  [25]  rs794728811 NM_001035.3(RYR2):c...
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
> overlaps <- findOverlaps(ANF_genome_gr, RYR2_hg19_gr)
> head(overlaps)
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]     44721          11
  [2]     44721          12
  [3]     44722          14
  [4]     44731          18
  [5]     44766          26
  [6]     44766          27
  -------
  queryLength: 610526 / subjectLength: 154
> intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr)
GRanges object with 25 ranges and 2 metadata columns:
       seqnames    ranges strand |    Genotype       dnSNP
          <Rle> <IRanges>  <Rle> | <character> <character>
   [1]        1 237540665      * |          GG    i6017887
   [2]        1 237540665      * |          GG    i6017887
   [3]        1 237540686      * |          GG    i6017888
   [4]        1 237551437      * |          GG    i6017892
   [5]        1 237608774      * |          CC    i6017854
   ...      ...       ...    ... .         ...         ...
  [21]        1 237954780      * |          GG    i6017900
  [22]        1 237954780      * |          GG    i6017900
  [23]        1 237972213      * |          GG    i6017915
  [24]        1 237993857      * |          AA    i6017871
  [25]        1 237995919      * |          GG    i6017906
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> unique(intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr))
GRanges object with 20 ranges and 2 metadata columns:
       seqnames    ranges strand |    Genotype       dnSNP
          <Rle> <IRanges>  <Rle> | <character> <character>
   [1]        1 237540665      * |          GG    i6017887
   [2]        1 237540686      * |          GG    i6017888
   [3]        1 237551437      * |          GG    i6017892
   [4]        1 237608774      * |          CC    i6017854
   [5]        1 237632425      * |          CC    i6017934
   ...      ...       ...    ... .         ...         ...
  [16]        1 237947545      * |          AA    i6017849
  [17]        1 237954780      * |          GG    i6017900
  [18]        1 237972213      * |          GG    i6017915
  [19]        1 237993857      * |          AA    i6017871
  [20]        1 237995919      * |          GG    i6017906
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> reduce(intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr))
GRanges object with 20 ranges and 0 metadata columns:
       seqnames    ranges strand
          <Rle> <IRanges>  <Rle>
   [1]        1 237540665      *
   [2]        1 237540686      *
   [3]        1 237551437      *
   [4]        1 237608774      *
   [5]        1 237632425      *
   ...      ...       ...    ...
  [16]        1 237947545      *
  [17]        1 237954780      *
  [18]        1 237972213      *
  [19]        1 237993857      *
  [20]        1 237995919      *
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> intersect_with_metadata(ANF_genome_gr,RYR2_hg19_gr)
GRanges object with 25 ranges and 2 metadata columns:
       seqnames    ranges strand |    Genotype       dnSNP
          <Rle> <IRanges>  <Rle> | <character> <character>
   [1]        1 237540665      * |          GG    i6017887
   [2]        1 237540665      * |          GG    i6017887
   [3]        1 237540686      * |          GG    i6017888
   [4]        1 237551437      * |          GG    i6017892
   [5]        1 237608774      * |          CC    i6017854
   ...      ...       ...    ... .         ...         ...
  [21]        1 237954780      * |          GG    i6017900
  [22]        1 237954780      * |          GG    i6017900
  [23]        1 237972213      * |          GG    i6017915
  [24]        1 237993857      * |          AA    i6017871
  [25]        1 237995919      * |          GG    i6017906
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths
> intersect_with_metadata(ANF_genome_gr,gr1=RYR2_hg19_gr)
GRanges object with 25 ranges and 6 metadata columns:
       seqnames    ranges strand |      functional_impact
          <Rle> <IRanges>  <Rle> |            <character>
   [1]        1 237540665      * | Pathogenic/Likely pa..
   [2]        1 237540665      * | Pathogenic/Likely pa..
   [3]        1 237540686      * |      Likely pathogenic
   [4]        1 237551437      * |      Likely pathogenic
   [5]        1 237608774      * |      Likely pathogenic
   ...      ...       ...    ... .                    ...
  [21]        1 237954780      * | Pathogenic/Likely pa..
  [22]        1 237954780      * |      Likely pathogenic
  [23]        1 237972213      * |      Likely pathogenic
  [24]        1 237993857      * |      Likely pathogenic
  [25]        1 237995919      * | Pathogenic/Likely pa..
       molecular_classification   AA_change       Clinical_Outcome
                    <character> <character>            <character>
   [1]         missense variant       R169L not provided|Cardiov..
   [2]         missense variant       R169Q not provided|Catecho..
   [3]         missense variant       R176L Cardiovascular pheno..
   [4]                 nonsense       E243*  RYR2-related disorder
   [5]         missense variant       T415I Catecholaminergic po..
   ...                      ...         ...                    ...
  [21]         missense variant      A4510P not provided|Catecho..
  [22]         missense variant      A4510S not provided|Arrhyth..
  [23]         missense variant      V4771F Cardiovascular pheno..
  [24]         missense variant      N4895H Catecholaminergic po..
  [25]         missense variant      R4959Q Cardiovascular pheno..
              dbSNP      Nucleotide_Change
        <character>            <character>
   [1]  rs397516539 NM_001035.3(RYR2):c...
   [2]  rs397516539 NM_001035.3(RYR2):c...
   [3]  rs794728708 NM_001035.3(RYR2):c...
   [4]  rs794728712 NM_001035.3(RYR2):c...
   [5] rs1288202574 NM_001035.3(RYR2):c...
   ...          ...                    ...
  [21]  rs397516510 NM_001035.3(RYR2):c...
  [22]  rs397516510 NM_001035.3(RYR2):c...
  [23]  rs794728804 NM_001035.3(RYR2):c...
  [24] rs1185619003 NM_001035.3(RYR2):c...
  [25]  rs794728811 NM_001035.3(RYR2):c...
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
> overlaps
Hits object with 25 hits and 0 metadata columns:
       queryHits subjectHits
       <integer>   <integer>
   [1]     44721          11
   [2]     44721          12
   [3]     44722          14
   [4]     44731          18
   [5]     44766          26
   ...       ...         ...
  [21]     44939         118
  [22]     44939         119
  [23]     44956         137
  [24]     44972         145
  [25]     44976         153
  -------
  queryLength: 610526 / subjectLength: 154
> str(overlaps)
Formal class 'SortedByQueryHits' [package "S4Vectors"] with 6 slots
  ..@ from           : int [1:25] 44721 44721 44722 44731 44766 44766 44776 44789 44829 44845 ...
  ..@ to             : int [1:25] 11 12 14 18 26 27 28 29 32 35 ...
  ..@ nLnode         : int 610526
  ..@ nRnode         : int 154
  ..@ elementMetadata: NULL
  ..@ metadata       : list()
> overlaps <- findOverlaps(ANF_genome_gr, RYR2_hg19_gr)  C-c C-c
> sapply(overlaps, function(x) data.frame(ANF_genome[x@from,],RYR2_hg19_gr[x@to,]))
Error in x@from : 
  no applicable method for `@` applied to an object of class "integer"
> do.call("rbind",lapply(1:length(overlaps), function(x) data.frame(ANF_genome[x@from,],RYR2_hg19_gr[x@to,])))
Error in h(simpleError(msg, call)) : 
  error in evaluating the argument 'args' in selecting a method for function 'do.call': no applicable method for `@` applied to an object of class "integer"
> lapply(overlaps, function(x) length(x))


> length(overlaps)
[1] 25
> do.call("rbind",lapply(1:length(overlaps), function(x) data.frame(ANF_genome[overlaps@from[x],],RYR2_hg19_gr[overlaps@to[x],])))
       chr     start       end Genotype    dnSNP seqnames   start.1
44721    1 237540665 237540665       GG i6017887        1 237540665
447211   1 237540665 237540665       GG i6017887        1 237540665
44722    1 237540686 237540686       GG i6017888        1 237540686
44731    1 237551437 237551437       GG i6017892        1 237551437
44766    1 237608774 237608774       CC i6017854        1 237608774
447661   1 237608774 237608774       CC i6017854        1 237608774
44776    1 237632425 237632425       CC i6017934        1 237632425
44789    1 237656273 237656273       CC i6017954        1 237656273
44829    1 237777598 237777598       GG i6017893        1 237777598
44845    1 237796969 237796969       AA i6017925        1 237796969
44850    1 237801780 237801780       GG i6017966        1 237801780
448501   1 237801780 237801780       GG i6017966        1 237801780
44858    1 237804240 237804240       GG i6017902        1 237804240
44859    1 237804256 237804256       AA i6017951        1 237804256
44861    1 237804283 237804283       GG i6017899        1 237804283
448611   1 237804283 237804283       GG i6017899        1 237804283
44907    1 237942027 237942027       GG i6017940        1 237942027
44925    1 237947458 237947458       AA i6017952        1 237947458
44926    1 237947482 237947482       GG i6017890        1 237947482
44929    1 237947545 237947545       AA i6017849        1 237947545
44939    1 237954780 237954780       GG i6017900        1 237954780
449391   1 237954780 237954780       GG i6017900        1 237954780
44956    1 237972213 237972213       GG i6017915        1 237972213
44972    1 237993857 237993857       AA i6017871        1 237993857
44976    1 237995919 237995919       GG i6017906        1 237995919
           end.1 width strand            functional_impact
44721  237540665     1      * Pathogenic/Likely pathogenic
447211 237540665     1      * Pathogenic/Likely pathogenic
44722  237540686     1      *            Likely pathogenic
44731  237551437     1      *            Likely pathogenic
44766  237608774     1      *            Likely pathogenic
447661 237608774     1      *            Likely pathogenic
44776  237632425     1      *            Likely pathogenic
44789  237656273     1      *            Likely pathogenic
44829  237777598     1      * Pathogenic/Likely pathogenic
44845  237796969     1      *            Likely pathogenic
44850  237801780     1      *            Likely pathogenic
448501 237801780     1      * Pathogenic/Likely pathogenic
44858  237804240     1      * Pathogenic/Likely pathogenic
44859  237804256     1      *            Likely pathogenic
44861  237804283     1      * Pathogenic/Likely pathogenic
448611 237804283     1      * Pathogenic/Likely pathogenic
44907  237942027     1      *            Likely pathogenic
44925  237947458     1      * Pathogenic/Likely pathogenic
44926  237947482     1      * Pathogenic/Likely pathogenic
44929  237947545     1      * Pathogenic/Likely pathogenic
44939  237954780     1      * Pathogenic/Likely pathogenic
449391 237954780     1      *            Likely pathogenic
44956  237972213     1      *            Likely pathogenic
44972  237993857     1      *            Likely pathogenic
44976  237995919     1      * Pathogenic/Likely pathogenic
       molecular_classification AA_change
44721          missense variant     R169L
447211         missense variant     R169Q
44722          missense variant     R176L
44731                  nonsense     E243*
44766          missense variant     T415I
447661         missense variant     T415K
44776          missense variant     A549V
44789          missense variant     S616L
44829          missense variant    E1724K
44845          missense variant    D2216G
44850          missense variant    V2306F
448501         missense variant    V2306I
44858          missense variant    A2387T
44859          missense variant    Y2392F
44861          missense variant    R2401L
448611         missense variant    R2401H
44907          missense variant    G3946V
44925          missense variant    Y4149C
44926          missense variant    R4157Q
44929          missense variant    N4178S
44939          missense variant    A4510P
449391         missense variant    A4510S
44956          missense variant    V4771F
44972          missense variant    N4895H
44976          missense variant    R4959Q
                                                                                                                                          Clinical_Outcome
44721                                                                                                                not provided|Cardiovascular phenotype
447211 not provided|Catecholaminergic polymorphic ventricular tachycardia 1|Cardiovascular phenotype|Catecholaminergic polymorphic ventricular tachycardia
44722                                                                                                                Cardiovascular phenotype|not provided
44731                                                                                                                                RYR2-related disorder
44766                                                                                              Catecholaminergic polymorphic ventricular tachycardia 1
447661                                                                                             Catecholaminergic polymorphic ventricular tachycardia 1
44776                                                 Arrhythmogenic right ventricular dysplasia 2|Catecholaminergic polymorphic ventricular tachycardia 1
44789                                                                                 not provided|Catecholaminergic polymorphic ventricular tachycardia 1
44829                                       Cardiovascular phenotype|not provided|Long QT syndrome|Catecholaminergic polymorphic ventricular tachycardia 1
44845                                                                                              Catecholaminergic polymorphic ventricular tachycardia 1
44850                                                                                              Catecholaminergic polymorphic ventricular tachycardia 1
448501                                                                                not provided|Catecholaminergic polymorphic ventricular tachycardia 1
44858  Cardiovascular phenotype|Catecholaminergic polymorphic ventricular tachycardia|Catecholaminergic polymorphic ventricular tachycardia 1|not provided
44859                                                                                                Catecholaminergic polymorphic ventricular tachycardia
44861                                                                                 not provided|Catecholaminergic polymorphic ventricular tachycardia 1
448611                                                       Cardiovascular phenotype|not provided|Catecholaminergic polymorphic ventricular tachycardia 1
44907                                                                                                                             Cardiovascular phenotype
44925                                                                     Cardiovascular phenotype|Catecholaminergic polymorphic ventricular tachycardia 1
44926                                                                                 not provided|Catecholaminergic polymorphic ventricular tachycardia 1
44929                                                        Cardiovascular phenotype|not provided|Catecholaminergic polymorphic ventricular tachycardia 1
44939                                                                                 not provided|Catecholaminergic polymorphic ventricular tachycardia 1
449391                                                                                           not provided|Arrhythmogenic right ventricular dysplasia 2
44956                                                                     Cardiovascular phenotype|Catecholaminergic polymorphic ventricular tachycardia 1
44972                                                                                              Catecholaminergic polymorphic ventricular tachycardia 1
44976                                                        Cardiovascular phenotype|not provided|Catecholaminergic polymorphic ventricular tachycardia 1
              dbSNP                           Nucleotide_Change
44721   rs397516539    NM_001035.3(RYR2):c.506G>T (p.Arg169Leu)
447211  rs397516539    NM_001035.3(RYR2):c.506G>A (p.Arg169Gln)
44722   rs794728708    NM_001035.3(RYR2):c.527G>T (p.Arg176Leu)
44731   rs794728712    NM_001035.3(RYR2):c.727G>T (p.Glu243Ter)
44766  rs1288202574   NM_001035.3(RYR2):c.1244C>T (p.Thr415Ile)
447661 rs1288202574   NM_001035.3(RYR2):c.1244C>A (p.Thr415Lys)
44776   rs794728817   NM_001035.3(RYR2):c.1646C>T (p.Ala549Val)
44789   rs730880187   NM_001035.3(RYR2):c.1847C>T (p.Ser616Leu)
44829   rs794728740  NM_001035.3(RYR2):c.5170G>A (p.Glu1724Lys)
44845  rs1328318082  NM_001035.3(RYR2):c.6647A>G (p.Asp2216Gly)
44850   rs794728746  NM_001035.3(RYR2):c.6916G>T (p.Val2306Phe)
448501  rs794728746  NM_001035.3(RYR2):c.6916G>A (p.Val2306Ile)
44858   rs794728753  NM_001035.3(RYR2):c.7159G>A (p.Ala2387Thr)
44859                NM_001035.3(RYR2):c.7175A>T (p.Tyr2392Phe)
44861   rs794728756  NM_001035.3(RYR2):c.7202G>T (p.Arg2401Leu)
448611  rs794728756  NM_001035.3(RYR2):c.7202G>A (p.Arg2401His)
44907               NM_001035.3(RYR2):c.11837G>T (p.Gly3946Val)
44925  rs1234449785 NM_001035.3(RYR2):c.12446A>G (p.Tyr4149Cys)
44926   rs794728786 NM_001035.3(RYR2):c.12470G>A (p.Arg4157Gln)
44929   rs794728787 NM_001035.3(RYR2):c.12533A>G (p.Asn4178Ser)
44939   rs397516510 NM_001035.3(RYR2):c.13528G>C (p.Ala4510Pro)
449391  rs397516510 NM_001035.3(RYR2):c.13528G>T (p.Ala4510Ser)
44956   rs794728804 NM_001035.3(RYR2):c.14311G>T (p.Val4771Phe)
44972  rs1185619003 NM_001035.3(RYR2):c.14683A>C (p.Asn4895His)
44976   rs794728811 NM_001035.3(RYR2):c.14876G>A (p.Arg4959Gln)
> ANF_RYR2_status <- do.call("rbind",lapply(1:length(overlaps), function(x) data.frame(ANF_genome[overlaps@from[x],],RYR2_hg19_gr[overlaps@to[x],])))
> ANF_RYR2_status <- ANF_RYR2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"))
Error in ANF_RYR2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean = stringr::str_extract(Nucleotide_Change,  : 
  could not find function "%>%"
> library(dplyr)

Attaching package: ‘dplyr’

The following objects are masked from ‘package:GenomicRanges’:

    intersect, setdiff, union

The following object is masked from ‘package:GenomeInfoDb’:

    intersect

The following objects are masked from ‘package:IRanges’:

    collapse, desc, intersect, setdiff, slice, union

The following objects are masked from ‘package:S4Vectors’:

    first, intersect, rename, setdiff, setequal, union

The following objects are masked from ‘package:BiocGenerics’:

    combine, intersect, setdiff, union

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

> ANF_RYR2_status <- ANF_RYR2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"))
> ANF_RYR2_status <- ANF_RYR2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"), Mut_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[2]) %>% ungroup
> str(ANF_RYR2_status)
tibble [25 × 18] (S3: tbl_df/tbl/data.frame)
 $ chr                     : chr [1:25] "1" "1" "1" "1" ...
 $ start                   : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ end                     : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ Genotype                : chr [1:25] "GG" "GG" "GG" "GG" ...
 $ dnSNP                   : chr [1:25] "i6017887" "i6017887" "i6017888" "i6017892" ...
 $ seqnames                : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
 $ start.1                 : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ end.1                   : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ width                   : int [1:25] 1 1 1 1 1 1 1 1 1 1 ...
 $ strand                  : Factor w/ 3 levels "+","-","*": 3 3 3 3 3 3 3 3 3 3 ...
 $ functional_impact       : chr [1:25] "Pathogenic/Likely pathogenic" "Pathogenic/Likely pathogenic" "Likely pathogenic" "Likely pathogenic" ...
 $ molecular_classification: chr [1:25] "missense variant" "missense variant" "missense variant" "nonsense" ...
 $ AA_change               : chr [1:25] "R169L" "R169Q" "R176L" "E243*" ...
 $ Clinical_Outcome        : chr [1:25] "not provided|Cardiovascular phenotype" "not provided|Catecholaminergic polymorphic ventricular tachycardia 1|Cardiovascular phenotype|Catecholaminergic"| __truncated__ "Cardiovascular phenotype|not provided" "RYR2-related disorder" ...
 $ dbSNP                   : chr [1:25] "rs397516539" "rs397516539" "rs794728708" "rs794728712" ...
 $ Nucleotide_Change       : chr [1:25] "NM_001035.3(RYR2):c.506G>T (p.Arg169Leu)" "NM_001035.3(RYR2):c.506G>A (p.Arg169Gln)" "NM_001035.3(RYR2):c.527G>T (p.Arg176Leu)" "NM_001035.3(RYR2):c.727G>T (p.Glu243Ter)" ...
 $ Nucleotide_Change_Clean : chr [1:25] "G>T" "G>A" "G>T" "G>T" ...
 $ Mut_Allele              : chr [1:25] "T" "A" "T" "T" ...
> ANF_RYR2_status <- ANF_RYR2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"), Mut_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[2],Ref_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[1]) %>% ungroup
> str(ANF_RYR2_status)
tibble [25 × 19] (S3: tbl_df/tbl/data.frame)
 $ chr                     : chr [1:25] "1" "1" "1" "1" ...
 $ start                   : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ end                     : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ Genotype                : chr [1:25] "GG" "GG" "GG" "GG" ...
 $ dnSNP                   : chr [1:25] "i6017887" "i6017887" "i6017888" "i6017892" ...
 $ seqnames                : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
 $ start.1                 : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ end.1                   : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ width                   : int [1:25] 1 1 1 1 1 1 1 1 1 1 ...
 $ strand                  : Factor w/ 3 levels "+","-","*": 3 3 3 3 3 3 3 3 3 3 ...
 $ functional_impact       : chr [1:25] "Pathogenic/Likely pathogenic" "Pathogenic/Likely pathogenic" "Likely pathogenic" "Likely pathogenic" ...
 $ molecular_classification: chr [1:25] "missense variant" "missense variant" "missense variant" "nonsense" ...
 $ AA_change               : chr [1:25] "R169L" "R169Q" "R176L" "E243*" ...
 $ Clinical_Outcome        : chr [1:25] "not provided|Cardiovascular phenotype" "not provided|Catecholaminergic polymorphic ventricular tachycardia 1|Cardiovascular phenotype|Catecholaminergic"| __truncated__ "Cardiovascular phenotype|not provided" "RYR2-related disorder" ...
 $ dbSNP                   : chr [1:25] "rs397516539" "rs397516539" "rs794728708" "rs794728712" ...
 $ Nucleotide_Change       : chr [1:25] "NM_001035.3(RYR2):c.506G>T (p.Arg169Leu)" "NM_001035.3(RYR2):c.506G>A (p.Arg169Gln)" "NM_001035.3(RYR2):c.527G>T (p.Arg176Leu)" "NM_001035.3(RYR2):c.727G>T (p.Glu243Ter)" ...
 $ Nucleotide_Change_Clean : chr [1:25] "G>T" "G>A" "G>T" "G>T" ...
 $ Mut_Allele              : chr [1:25] "T" "A" "T" "T" ...
 $ Ref_Allele              : chr [1:25] "G" "G" "G" "G" ...
> ANF_RYR2_status <- ANF_RYR2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"), Mut_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[2],Ref_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[1], Flag=!grepl(Ref_Allele, Genotype)) %>% ungroup
> str(ANF_RYR2_status)
tibble [25 × 20] (S3: tbl_df/tbl/data.frame)
 $ chr                     : chr [1:25] "1" "1" "1" "1" ...
 $ start                   : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ end                     : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ Genotype                : chr [1:25] "GG" "GG" "GG" "GG" ...
 $ dnSNP                   : chr [1:25] "i6017887" "i6017887" "i6017888" "i6017892" ...
 $ seqnames                : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
 $ start.1                 : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ end.1                   : int [1:25] 237540665 237540665 237540686 237551437 237608774 237608774 237632425 237656273 237777598 237796969 ...
 $ width                   : int [1:25] 1 1 1 1 1 1 1 1 1 1 ...
 $ strand                  : Factor w/ 3 levels "+","-","*": 3 3 3 3 3 3 3 3 3 3 ...
 $ functional_impact       : chr [1:25] "Pathogenic/Likely pathogenic" "Pathogenic/Likely pathogenic" "Likely pathogenic" "Likely pathogenic" ...
 $ molecular_classification: chr [1:25] "missense variant" "missense variant" "missense variant" "nonsense" ...
 $ AA_change               : chr [1:25] "R169L" "R169Q" "R176L" "E243*" ...
 $ Clinical_Outcome        : chr [1:25] "not provided|Cardiovascular phenotype" "not provided|Catecholaminergic polymorphic ventricular tachycardia 1|Cardiovascular phenotype|Catecholaminergic"| __truncated__ "Cardiovascular phenotype|not provided" "RYR2-related disorder" ...
 $ dbSNP                   : chr [1:25] "rs397516539" "rs397516539" "rs794728708" "rs794728712" ...
 $ Nucleotide_Change       : chr [1:25] "NM_001035.3(RYR2):c.506G>T (p.Arg169Leu)" "NM_001035.3(RYR2):c.506G>A (p.Arg169Gln)" "NM_001035.3(RYR2):c.527G>T (p.Arg176Leu)" "NM_001035.3(RYR2):c.727G>T (p.Glu243Ter)" ...
 $ Nucleotide_Change_Clean : chr [1:25] "G>T" "G>A" "G>T" "G>T" ...
 $ Mut_Allele              : chr [1:25] "T" "A" "T" "T" ...
 $ Ref_Allele              : chr [1:25] "G" "G" "G" "G" ...
 $ Flag                    : logi [1:25] FALSE FALSE FALSE FALSE FALSE FALSE ...
> ANF_RYR2_status <- ANF_RYR2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"), Mut_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[2],Ref_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[1], Flag=!grepl(Ref_Allele, Genotype)) %>% ungroup
> dplyr::filter(ANF_RYR2_status, Flag==TRUE)
# A tibble: 0 × 20
# ℹ 20 variables: chr <chr>, start <int>, end <int>, Genotype <chr>,
#   dnSNP <chr>, seqnames <fct>, start.1 <int>, end.1 <int>, width <int>,
#   strand <fct>, functional_impact <chr>, molecular_classification <chr>,
#   AA_change <chr>, Clinical_Outcome <chr>, dbSNP <chr>,
#   Nucleotide_Change <chr>, Nucleotide_Change_Clean <chr>,
#   Mut_Allele <chr>, Ref_Allele <chr>, Flag <lgl>
> dplyr::select(ANF_RYR2_status, Mut_Allele, Ref_Allele, Genotype)
# A tibble: 25 × 3
   Mut_Allele Ref_Allele Genotype
   <chr>      <chr>      <chr>   
 1 T          G          GG      
 2 A          G          GG      
 3 T          G          GG      
 4 T          G          GG      
 5 T          C          CC      
 6 A          C          CC      
 7 T          C          CC      
 8 T          C          CC      
 9 A          G          GG      
10 G          A          AA      
# ℹ 15 more rows
# ℹ Use `print(n = ...)` to see more rows
> dplyr::select(ANF_RYR2_status, Ref_Allele, Mut_Allele, Genotype) %>% data.frame)
Error: unexpected ')' in "dplyr::select(ANF_RYR2_status, Ref_Allele, Mut_Allele, Genotype) %>% data.frame)"
> dplyr::select(ANF_RYR2_status, Ref_Allele, Mut_Allele, Genotype) %>% data.frame
   Ref_Allele Mut_Allele Genotype
1           G          T       GG
2           G          A       GG
3           G          T       GG
4           G          T       GG
5           C          T       CC
6           C          A       CC
7           C          T       CC
8           C          T       CC
9           G          A       GG
10          A          G       AA
11          G          T       GG
12          G          A       GG
13          G          A       GG
14          A          T       AA
15          G          T       GG
16          G          A       GG
17          G          T       GG
18          A          G       AA
19          G          A       GG
20          A          G       AA
21          G          C       GG
22          G          T       GG
23          G          T       GG
24          A          C       AA
25          G          A       GG
> PKP2_data <- read.table("/Users/forbesa/Downloads/clinvar_result (1).txt")
Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  line 2 did not have 52 elements
> PKP2_data <- read.table("/Users/forbesa/Downloads/clinvar_result (1).txt",header=T, sep='\t')
> str(PKP2_data)
'data.frame':	210 obs. of  25 variables:
 $ Name                                       : chr  "NM_001005242.3(PKP2):c.2446-1G>A" "NM_001005242.3(PKP2):c.2411G>A (p.Trp804Ter)" "NM_001005242.3(PKP2):c.2390C>T (p.Ser797Phe)" "NM_001005242.3(PKP2):c.2361T>A (p.Tyr787Ter)" ...
 $ Gene.s.                                    : chr  "PKP2" "PKP2" "PKP2" "PKP2" ...
 $ Protein.change                             : chr  "" "G632R, G741R, W695*, W749*, W804*, W848*" "S841F, S797F" "M615K, M724K, Y831*, Y678*, Y732*, Y787*" ...
 $ Condition.s.                               : chr  "Arrhythmogenic right ventricular dysplasia 9" "not provided" "not provided" "Arrhythmogenic right ventricular dysplasia 9|not provided|Cardiovascular phenotype" ...
 $ Accession                                  : chr  "VCV003065037" "VCV003338738" "VCV000201966" "VCV002735817" ...
 $ GRCh37Chromosome                           : int  12 12 12 12 12 12 12 12 12 12 ...
 $ GRCh37Location                             : chr  "32945427" "32945612" "32945633" "32945662" ...
 $ GRCh38Chromosome                           : int  12 12 12 12 12 12 12 12 12 12 ...
 $ GRCh38Location                             : chr  "32792493" "32792678" "32792699" "32792728" ...
 $ VariationID                                : int  3065037 3338738 201966 2735817 936764 978281 419485 45071 6757 188663 ...
 $ AlleleID.s.                                : int  3225158 3498008 198388 2898827 940264 966505 408638 54238 21796 186419 ...
 $ dbSNP.ID                                   : chr  "rs2541184345" "" "rs794729099" "rs1956082308" ...
 $ Canonical.SPDI                             : chr  "NC_000012.12:32792492:C:T" "NC_000012.12:32792677:C:T" "NC_000012.12:32792698:G:A" "NC_000012.12:32792727:A:T" ...
 $ Variant.type                               : chr  "single nucleotide variant" "single nucleotide variant" "single nucleotide variant" "single nucleotide variant" ...
 $ Molecular.consequence                      : chr  "splice acceptor variant" "nonsense|missense variant" "missense variant" "nonsense|missense variant" ...
 $ Germline.classification                    : chr  "Pathogenic" "Likely pathogenic" "Likely pathogenic" "Pathogenic/Likely pathogenic" ...
 $ Germline.date.last.evaluated               : chr  "Mar 25, 2024" "Oct 18, 2023" "Mar 26, 2014" "Nov 15, 2024" ...
 $ Germline.review.status                     : chr  "criteria provided, single submitter" "criteria provided, single submitter" "criteria provided, single submitter" "criteria provided, multiple submitters, no conflicts" ...
 $ Somatic.clinical.impact                    : logi  NA NA NA NA NA NA ...
 $ Somatic.clinical.impact.date.last.evaluated: logi  NA NA NA NA NA NA ...
 $ Somatic.clinical.impact.review.status      : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.classification                : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.date.last.evaluated           : logi  NA NA NA NA NA NA ...
 $ Oncogenicity.review.status                 : logi  NA NA NA NA NA NA ...
 $ X                                          : logi  NA NA NA NA NA NA ...
> PKP2_hg19_bed <- dplyr::select(PKP2_data,chr=GRCh37Chromosome,start=GRCh37Location,end=GRCh37Location,functional_impact=Germline.classification,molecular_classification=Molecular.consequence, AA_change=Protein.change,Clinical_Outcome=Condition.s.,dbSNP=dbSNP.ID, Nucleotide_Change=Name); PKP2_hg19_bed$end <- as.integer(PKP2_hg19_bed$end);PKP2_hg19_bed$start <- as.integer(PKP2_hg19_bed$start);PKP2_hg19_gr <- table_to_granges(dplyr::filter(PKP2_hg19_bed, !is.na(start))); PKP2_hg19_gr <- table_to_granges(PKP2_hg19_bed)
Warning message:
NAs introduced by coercion 
Warning message:
NAs introduced by coercion 
Error in .new_IRanges_from_start_end(start, end) : 
  'start' or 'end' cannot contain NAs
> str(PKP2_hg19_bed)
'data.frame':	210 obs. of  9 variables:
 $ chr                     : int  12 12 12 12 12 12 12 12 12 12 ...
 $ start                   : int  32945427 32945612 32945633 32945662 NA 32945667 32949039 32949042 32949042 NA ...
 $ end                     : int  32945427 32945612 32945633 32945662 NA 32945667 32949039 32949042 32949042 NA ...
 $ functional_impact       : chr  "Pathogenic" "Likely pathogenic" "Likely pathogenic" "Pathogenic/Likely pathogenic" ...
 $ molecular_classification: chr  "splice acceptor variant" "nonsense|missense variant" "missense variant" "nonsense|missense variant" ...
 $ AA_change               : chr  "" "G632R, G741R, W695*, W749*, W804*, W848*" "S841F, S797F" "M615K, M724K, Y831*, Y678*, Y732*, Y787*" ...
 $ Clinical_Outcome        : chr  "Arrhythmogenic right ventricular dysplasia 9" "not provided" "not provided" "Arrhythmogenic right ventricular dysplasia 9|not provided|Cardiovascular phenotype" ...
 $ dbSNP                   : chr  "rs2541184345" "" "rs794729099" "rs1956082308" ...
 $ Nucleotide_Change       : chr  "NM_001005242.3(PKP2):c.2446-1G>A" "NM_001005242.3(PKP2):c.2411G>A (p.Trp804Ter)" "NM_001005242.3(PKP2):c.2390C>T (p.Ser797Phe)" "NM_001005242.3(PKP2):c.2361T>A (p.Tyr787Ter)" ...
> PKP2_hg19_gr
GRanges object with 135 ranges and 6 metadata columns:
        seqnames    ranges strand |      functional_impact
           <Rle> <IRanges>  <Rle> |            <character>
    [1]       12  32945427      * |             Pathogenic
    [2]       12  32945612      * |      Likely pathogenic
    [3]       12  32945633      * |      Likely pathogenic
    [4]       12  32945662      * | Pathogenic/Likely pa..
    [5]       12  32945667      * |      Likely pathogenic
    ...      ...       ...    ... .                    ...
  [131]       12  33049599      * |             Pathogenic
  [132]       12  33049664      * |      Likely pathogenic
  [133]       12  33049665      * |      Likely pathogenic
  [134]       12  33049665      * | Pathogenic/Likely pa..
  [135]       12  33049665      * |             Pathogenic
        molecular_classification              AA_change
                     <character>            <character>
    [1]   splice acceptor vari..                       
    [2]   nonsense|missense va.. G632R, G741R, W695*,..
    [3]         missense variant           S841F, S797F
    [4]   nonsense|missense va.. M615K, M724K, Y831*,..
    [5]   splice acceptor vari..                       
    ...                      ...                    ...
  [131]   nonsense|5 prime UTR..                   G23*
  [132]   missense variant|ini..                    M1T
  [133]   missense variant|ini..                    M1L
  [134]   missense variant|ini..                    M1V
  [135]   missense variant|ini..                    M1L
              Clinical_Outcome        dbSNP      Nucleotide_Change
                   <character>  <character>            <character>
    [1] Arrhythmogenic right.. rs2541184345 NM_001005242.3(PKP2)..
    [2]           not provided              NM_001005242.3(PKP2)..
    [3]           not provided  rs794729099 NM_001005242.3(PKP2)..
    [4] Arrhythmogenic right.. rs1956082308 NM_001005242.3(PKP2)..
    [5] Arrhythmogenic right.. rs1956082474 NM_001005242.3(PKP2)..
    ...                    ...          ...                    ...
  [131] Arrhythmogenic right..  rs770498155 NM_001005242.3(PKP2)..
  [132] Arrhythmogenic right.. rs1957129506 NM_001005242.3(PKP2)..
  [133] Arrhythmogenic right..  rs794729107 NM_001005242.3(PKP2)..
  [134] Cardiovascular pheno..  rs794729107 NM_001005242.3(PKP2)..
  [135]           not provided  rs794729107 NM_001005242.3(PKP2)..
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
> overlaps <- findOverlaps(ANF_genome_gr, PKP2_hg19_gr) ;ANF_PKP2_status <- do.call("rbind",lapply(1:length(overlaps), function(x) data.frame(ANF_genome[overlaps@from[x],],PKP2_hg19_gr[overlaps@to[x],])));  ANF_PKP2_status <- ANF_PKP2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"), Mut_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[2],Ref_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[1], Flag=!grepl(Ref_Allele, Genotype)) %>% ungroup
> dplyr::filter(ANF_RYR2_status, Flag==TRUE)
# A tibble: 0 × 20
# ℹ 20 variables: chr <chr>, start <int>, end <int>, Genotype <chr>,
#   dnSNP <chr>, seqnames <fct>, start.1 <int>, end.1 <int>, width <int>,
#   strand <fct>, functional_impact <chr>, molecular_classification <chr>,
#   AA_change <chr>, Clinical_Outcome <chr>, dbSNP <chr>,
#   Nucleotide_Change <chr>, Nucleotide_Change_Clean <chr>,
#   Mut_Allele <chr>, Ref_Allele <chr>, Flag <lgl>
> dplyr::filter(ANF_PKP2_status, Flag==TRUE)
# A tibble: 21 × 20
   chr     start    end Genotype dnSNP seqnames start.1  end.1 width strand
   <chr>   <int>  <int> <chr>    <chr> <fct>      <int>  <int> <int> <fct> 
 1 12     3.29e7 3.29e7 GG       i601… 12        3.29e7 3.29e7     1 *     
 2 12     3.29e7 3.29e7 GG       i601… 12        3.29e7 3.29e7     1 *     
 3 12     3.29e7 3.29e7 AA       i601… 12        3.29e7 3.29e7     1 *     
 4 12     3.30e7 3.30e7 GG       i500… 12        3.30e7 3.30e7     1 *     
 5 12     3.30e7 3.30e7 GG       i601… 12        3.30e7 3.30e7     1 *     
 6 12     3.30e7 3.30e7 CC       rs19… 12        3.30e7 3.30e7     1 *     
 7 12     3.30e7 3.30e7 GG       i601… 12        3.30e7 3.30e7     1 *     
 8 12     3.30e7 3.30e7 GG       i601… 12        3.30e7 3.30e7     1 *     
 9 12     3.30e7 3.30e7 CC       i601… 12        3.30e7 3.30e7     1 *     
10 12     3.30e7 3.30e7 GG       i601… 12        3.30e7 3.30e7     1 *     
# ℹ 11 more rows
# ℹ 10 more variables: functional_impact <chr>,
#   molecular_classification <chr>, AA_change <chr>,
#   Clinical_Outcome <chr>, dbSNP <chr>, Nucleotide_Change <chr>,
#   Nucleotide_Change_Clean <chr>, Mut_Allele <chr>, Ref_Allele <chr>,
#   Flag <lgl>
# ℹ Use `print(n = ...)` to see more rows
> dplyr::filter(ANF_PKP2_status, Flag==TRUE) %>% data.frame
   chr    start      end Genotype       dnSNP seqnames  start.1    end.1
1   12 32949111 32949111       GG    i6016890       12 32949111 32949111
2   12 32949111 32949111       GG    i6016890       12 32949111 32949111
3   12 32949146 32949146       AA    i6016877       12 32949146 32949146
4   12 32955433 32955433       GG    i5007217       12 32955433 32955433
5   12 32955460 32955460       GG    i6016870       12 32955460 32955460
6   12 32955491 32955491       CC rs193922674       12 32955491 32955491
7   12 32974316 32974316       GG    i6016881       12 32974316 32974316
8   12 32974340 32974340       GG    i6016863       12 32974340 32974340
9   12 32974407 32974407       CC    i6016897       12 32974407 32974407
10  12 32974457 32974457       GG    i6016869       12 32974457 32974457
11  12 32975460 32975460       GG    i6016868       12 32975460 32975460
12  12 32975524 32975524       GG    i6016889       12 32975524 32975524
13  12 32975524 32975524       GG    i6016889       12 32975524 32975524
14  12 32975528 32975528       GG    i6016893       12 32975528 32975528
15  12 32994037 32994037       CC    i6016896       12 32994037 32994037
16  12 33003841 33003841       GG    i6016883       12 33003841 33003841
17  12 33021869 33021869       GG    i6016895       12 33021869 33021869
18  12 33021899 33021899       GG    i6016864       12 33021899 33021899
19  12 33031183 33031183       GG    i6016880       12 33031183 33031183
20  12 33031417 33031417       GG    i6016867       12 33031417 33031417
21  12 33031955 33031955       GG    i6016871       12 33031955 33031955
   width strand            functional_impact molecular_classification
1      1      *                   Pathogenic  nonsense|intron variant
2      1      *                   Pathogenic  nonsense|intron variant
3      1      *                   Pathogenic         missense variant
4      1      *                   Pathogenic                 nonsense
5      1      *                   Pathogenic                 nonsense
6      1      *                   Pathogenic  splice acceptor variant
7      1      *            Likely pathogenic                 nonsense
8      1      *                   Pathogenic                 nonsense
9      1      *                   Pathogenic                 nonsense
10     1      *                   Pathogenic                 nonsense
11     1      *                   Pathogenic                 nonsense
12     1      *                   Pathogenic  nonsense|intron variant
13     1      *                   Pathogenic                 nonsense
14     1      *                   Pathogenic         missense variant
15     1      *                   Pathogenic                 nonsense
16     1      *                   Pathogenic                 nonsense
17     1      * Pathogenic/Likely pathogenic         missense variant
18     1      *                   Pathogenic                 nonsense
19     1      *                   Pathogenic                 nonsense
20     1      *                   Pathogenic                 nonsense
21     1      *                   Pathogenic                 nonsense
                    AA_change
1  Y654*, Y708*, Y763*, Y807*
2  Y654*, Y708*, Y763*, Y807*
3                C796R, C752R
4                R735*, R691*
5                Q682*, Q726*
6                            
7                Q707*, Q663*
8  Q655*, Q546*, Q600*, Q699*
9                W676*, W632*
10               Q660*, Q616*
11               Q638*, Q594*
12        Y572*, Y463*, Y616*
13               Y616*, Y572*
14               S615F, S571F
15               W538*, W494*
16                      R413*
17                      R388W
18                      Q378*
19               Q102*, Q211*
20                      Q133*
21                       R79*
                                                                                                                                                                                                 Clinical_Outcome
1                                                                                                                                           Arrhythmogenic right ventricular dysplasia 9|Cardiovascular phenotype
2                                                                                                                                        not provided|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
3                                                  Cardiovascular phenotype|Familial isolated arrhythmogenic right ventricular dysplasia|Arrhythmogenic right ventricular dysplasia 9|Cardiomyopathy|not provided
4                                                               Cardiovascular phenotype|Arrhythmogenic right ventricular cardiomyopathy|not provided|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
5                                                                                                       Familial isolated arrhythmogenic right ventricular dysplasia|Arrhythmogenic right ventricular dysplasia 9
6                                                               Cardiovascular phenotype|not provided|Arrhythmogenic right ventricular cardiomyopathy|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
7                                                                                                                                                    Arrhythmogenic right ventricular cardiomyopathy|not provided
8                                                                                                                                                       Arrhythmogenic right ventricular dysplasia 9|not provided
9                                                                                          Familial isolated arrhythmogenic right ventricular dysplasia|not provided|Arrhythmogenic right ventricular dysplasia 9
10                                                                                                    Arrhythmogenic right ventricular dysplasia 9|Arrhythmogenic right ventricular cardiomyopathy|Cardiomyopathy
11                          Familial isolated arrhythmogenic right ventricular dysplasia|not provided|Arrhythmogenic right ventricular cardiomyopathy|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
12                                                                                                                                                                   Arrhythmogenic right ventricular dysplasia 9
13                                                                                                                                                                   Arrhythmogenic right ventricular dysplasia 9
14                                                                                                                                                                   Arrhythmogenic right ventricular dysplasia 9
15 Cardiovascular phenotype|Familial isolated arrhythmogenic right ventricular dysplasia|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9|not provided|Arrhythmogenic right ventricular cardiomyopathy
16 Cardiovascular phenotype|Familial isolated arrhythmogenic right ventricular dysplasia|not provided|Cardiomyopathy|Arrhythmogenic right ventricular cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
17                                                                                       not provided|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9|Arrhythmogenic right ventricular cardiomyopathy
18 Cardiovascular phenotype|Familial isolated arrhythmogenic right ventricular dysplasia|not provided|Arrhythmogenic right ventricular cardiomyopathy|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
19                                                                                                                                                                                                   not provided
20                                                                not provided|Arrhythmogenic right ventricular dysplasia 9|Cardiovascular phenotype|Familial isolated arrhythmogenic right ventricular dysplasia
21                                           Cardiac arrhythmia|Cardiovascular phenotype|not provided|Arrhythmogenic right ventricular cardiomyopathy|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
          dbSNP                            Nucleotide_Change
1               NM_001005242.3(PKP2):c.2289C>G (p.Tyr763Ter)
2               NM_001005242.3(PKP2):c.2289C>A (p.Tyr763Ter)
3   rs794729098 NM_001005242.3(PKP2):c.2254T>C (p.Cys752Arg)
4   rs121434421 NM_001005242.3(PKP2):c.2071C>T (p.Arg691Ter)
5  rs1325285497 NM_001005242.3(PKP2):c.2044C>T (p.Gln682Ter)
6   rs193922674             NM_001005242.3(PKP2):c.2014-1G>C
7   rs397517017 NM_001005242.3(PKP2):c.1987C>T (p.Gln663Ter)
8  rs1214999277 NM_001005242.3(PKP2):c.1963C>T (p.Gln655Ter)
9   rs193922673 NM_001005242.3(PKP2):c.1896G>A (p.Trp632Ter)
10  rs762753884 NM_001005242.3(PKP2):c.1846C>T (p.Gln616Ter)
11  rs397517012 NM_001005242.3(PKP2):c.1780C>T (p.Gln594Ter)
12 rs1486464304 NM_001005242.3(PKP2):c.1716C>G (p.Tyr572Ter)
13 rs1486464304 NM_001005242.3(PKP2):c.1716C>A (p.Tyr572Ter)
14 rs1060501186 NM_001005242.3(PKP2):c.1712C>T (p.Ser571Phe)
15  rs193922672 NM_001005242.3(PKP2):c.1481G>A (p.Trp494Ter)
16  rs372827156 NM_001005242.3(PKP2):c.1237C>T (p.Arg413Ter)
17  rs766209297 NM_001005242.3(PKP2):c.1162C>T (p.Arg388Trp)
18  rs397516986 NM_001005242.3(PKP2):c.1132C>T (p.Gln378Ter)
19               NM_001005242.3(PKP2):c.631C>T (p.Gln211Ter)
20  rs794729132  NM_001005242.3(PKP2):c.397C>T (p.Gln133Ter)
21  rs121434420   NM_001005242.3(PKP2):c.235C>T (p.Arg79Ter)
   Nucleotide_Change_Clean Mut_Allele Ref_Allele Flag
1                      C>G          G          C TRUE
2                      C>A          A          C TRUE
3                      T>C          C          T TRUE
4                      C>T          T          C TRUE
5                      C>T          T          C TRUE
6                      G>C          C          G TRUE
7                      C>T          T          C TRUE
8                      C>T          T          C TRUE
9                      G>A          A          G TRUE
10                     C>T          T          C TRUE
11                     C>T          T          C TRUE
12                     C>G          G          C TRUE
13                     C>A          A          C TRUE
14                     C>T          T          C TRUE
15                     G>A          A          G TRUE
16                     C>T          T          C TRUE
17                     C>T          T          C TRUE
18                     C>T          T          C TRUE
19                     C>T          T          C TRUE
20                     C>T          T          C TRUE
21                     C>T          T          C TRUE
> dplyr::select(ANF_PKP2_status, Ref_Allele, Mut_Allele, Genotype) %>% data.frame
   Ref_Allele Mut_Allele Genotype
1           C          G       GG
2           C          A       GG
3           T          C       AA
4           C          T       GG
5           C          T       GG
6           G          C       CC
7           C          T       GG
8           C          T       GG
9           G          A       CC
10          C          T       GG
11          C          T       GG
12          C          G       GG
13          C          A       GG
14          C          T       GG
15          G          A       CC
16          C          T       GG
17          C          T       GG
18          C          T       GG
19          C          T       GG
20          C          T       GG
21          C          T       GG
> dplyr::filter(ANF_RYR2_status, Flag==TRUE) %>% dplyr::mutate(Flag  C-c C-c
> dplyr::select(ANF_PKP2_status, Ref_Allele, Mut_Allele, Genotype) %>% dplyr::mutate(Flag2 =grepl(Mut_Allele,Genotype)) %>% data.frame
   Ref_Allele Mut_Allele Genotype Flag2
1           C          G       GG  TRUE
2           C          A       GG  TRUE
3           T          C       AA FALSE
4           C          T       GG  TRUE
5           C          T       GG  TRUE
6           G          C       CC FALSE
7           C          T       GG  TRUE
8           C          T       GG  TRUE
9           G          A       CC FALSE
10          C          T       GG  TRUE
11          C          T       GG  TRUE
12          C          G       GG  TRUE
13          C          A       GG  TRUE
14          C          T       GG  TRUE
15          G          A       CC FALSE
16          C          T       GG  TRUE
17          C          T       GG  TRUE
18          C          T       GG  TRUE
19          C          T       GG  TRUE
20          C          T       GG  TRUE
21          C          T       GG  TRUE
Warning message:
There was 1 warning in `dplyr::mutate()`.
ℹ In argument: `Flag2 = grepl(Mut_Allele, Genotype)`.
Caused by warning in `grepl()`:
! argument 'pattern' has length > 1 and only the first element will be used 
> dplyr::select(ANF_PKP2_status, Ref_Allele, Mut_Allele, Genotype, Flag) %>% data.frame
   Ref_Allele Mut_Allele Genotype Flag
1           C          G       GG TRUE
2           C          A       GG TRUE
3           T          C       AA TRUE
4           C          T       GG TRUE
5           C          T       GG TRUE
6           G          C       CC TRUE
7           C          T       GG TRUE
8           C          T       GG TRUE
9           G          A       CC TRUE
10          C          T       GG TRUE
11          C          T       GG TRUE
12          C          G       GG TRUE
13          C          A       GG TRUE
14          C          T       GG TRUE
15          G          A       CC TRUE
16          C          T       GG TRUE
17          C          T       GG TRUE
18          C          T       GG TRUE
19          C          T       GG TRUE
20          C          T       GG TRUE
21          C          T       GG TRUE
> overlaps <- findOverlaps(ANF_genome_gr, PKP2_hg19_gr) ;ANF_PKP2_status <- do.call("rbind",lapply(1:length(overlaps), function(x) data.frame(ANF_genome[overlaps@from[x],],PKP2_hg19_gr[overlaps@to[x],])));  ANF_PKP2_status <- ANF_PKP2_status %>% group_by(Nucleotide_Change) %>% dplyr::mutate(Nucleotide_Change_Clean=stringr::str_extract(Nucleotide_Change,"[A-Z]>[A-Z]"), Mut_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[2],Ref_Allele=unlist(strsplit(Nucleotide_Change_Clean,">"))[1], Flag=grepl(Mut_Allele, Genotype)) %>% ungroup
> dplyr::select(ANF_PKP2_status, Ref_Allele, Mut_Allele, Genotype, Flag) %>% data.frame
   Ref_Allele Mut_Allele Genotype  Flag
1           C          G       GG  TRUE
2           C          A       GG FALSE
3           T          C       AA FALSE
4           C          T       GG FALSE
5           C          T       GG FALSE
6           G          C       CC  TRUE
7           C          T       GG FALSE
8           C          T       GG FALSE
9           G          A       CC FALSE
10          C          T       GG FALSE
11          C          T       GG FALSE
12          C          G       GG  TRUE
13          C          A       GG FALSE
14          C          T       GG FALSE
15          G          A       CC FALSE
16          C          T       GG FALSE
17          C          T       GG FALSE
18          C          T       GG FALSE
19          C          T       GG FALSE
20          C          T       GG FALSE
21          C          T       GG FALSE
> dplyr::select(ANF_PKP2_status, Ref_Allele, Mut_Allele, Genotype, Flag) %>% dplyr::filter(Flag==TRUE) %>% data.frame 
  Ref_Allele Mut_Allele Genotype Flag
1          C          G       GG TRUE
2          G          C       CC TRUE
3          C          G       GG TRUE
> dplyr::filter(ANF_PKP2_status, Flag==TRUE)
# A tibble: 3 × 20
  chr      start    end Genotype dnSNP seqnames start.1  end.1 width strand
  <chr>    <int>  <int> <chr>    <chr> <fct>      <int>  <int> <int> <fct> 
1 12    32949111 3.29e7 GG       i601… 12        3.29e7 3.29e7     1 *     
2 12    32955491 3.30e7 CC       rs19… 12        3.30e7 3.30e7     1 *     
3 12    32975524 3.30e7 GG       i601… 12        3.30e7 3.30e7     1 *     
# ℹ 10 more variables: functional_impact <chr>,
#   molecular_classification <chr>, AA_change <chr>,
#   Clinical_Outcome <chr>, dbSNP <chr>, Nucleotide_Change <chr>,
#   Nucleotide_Change_Clean <chr>, Mut_Allele <chr>, Ref_Allele <chr>,
#   Flag <lgl>
> dplyr::filter(ANF_PKP2_status, Flag==TRUE) %>% data.frame
  chr    start      end Genotype       dnSNP seqnames  start.1    end.1
1  12 32949111 32949111       GG    i6016890       12 32949111 32949111
2  12 32955491 32955491       CC rs193922674       12 32955491 32955491
3  12 32975524 32975524       GG    i6016889       12 32975524 32975524
  width strand functional_impact molecular_classification
1     1      *        Pathogenic  nonsense|intron variant
2     1      *        Pathogenic  splice acceptor variant
3     1      *        Pathogenic  nonsense|intron variant
                   AA_change
1 Y654*, Y708*, Y763*, Y807*
2                           
3        Y572*, Y463*, Y616*
                                                                                                                                   Clinical_Outcome
1                                                                             Arrhythmogenic right ventricular dysplasia 9|Cardiovascular phenotype
2 Cardiovascular phenotype|not provided|Arrhythmogenic right ventricular cardiomyopathy|Cardiomyopathy|Arrhythmogenic right ventricular dysplasia 9
3                                                                                                      Arrhythmogenic right ventricular dysplasia 9
         dbSNP                            Nucleotide_Change
1              NM_001005242.3(PKP2):c.2289C>G (p.Tyr763Ter)
2  rs193922674             NM_001005242.3(PKP2):c.2014-1G>C
3 rs1486464304 NM_001005242.3(PKP2):c.1716C>G (p.Tyr572Ter)
  Nucleotide_Change_Clean Mut_Allele Ref_Allele Flag
1                     C>G          G          C TRUE
2                     G>C          C          G TRUE
3                     C>G          G          C TRUE
> 
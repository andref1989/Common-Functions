table_to_granges <- function(table,strand=NULL){ require(GenomicRanges)
                                                 if(is.null(strand)){
                                                     gr <- GRanges(seqnames=Rle(table[,1]), ranges= IRanges(start=as.numeric(table[,2]), end=as.numeric(table[,3])),strand="*")
                                                     ##print(setdiff(colnames(table),x=colnames(table)[4:ncol(table)]))
                                                     if(ncol(table) >=4){
                                                     gr@elementMetadata@listData <- as.list(table[setdiff(colnames(table)[1:3],x=colnames(table)[4:ncol(table)])])}
                                                 } else {gr <- GRanges(seqnames=Rle(table[,1]), ranges= IRanges(start=table[,2], end=table[,3]),strand=table[,strand])
                                                                                                                                                               gr@elementMetadata@listData <- as.list(table[setdiff(grep("strand", colnames(table),ignore.case=T),x=4:ncol(table))])}
                        return(gr)}


bed_to_granges_dynamic <- function(bed,header=FALSE,nlines=-1){
    require(GenomicRanges)


    bed <- read.table(bed, stringsAsFactors=F, sep="\t",header=header,nrows=nlines)
##    str(bed)
    if(nrow(bed) <=1 && header==TRUE){ print("Doing Nothing")} else if(nrow(bed) <1 && header==FALSE){ print("Doing nothing, empty file")} else if (nrow(bed) ==1 && header==FALSE){
    gr <- GRanges(seqnames=Rle(bed[,1]), ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])))
    gr@elementMetadata@listData <- as.list(bed[4:ncol(bed)])

    return(gr)}
                                                                else{
                                                                    gr <- GRanges(seqnames=Rle(bed[,1]), ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])))
    gr@elementMetadata@listData <- as.list(bed[4:ncol(bed)])

                                                                    return(gr)}}


bedpe_to_granges <- function(bedpe,ID_col=NULL,ranges_1=1:3,ranges_2=4:6,strand=NULL,header=FALSE,filter_trans=FALSE,el_names="Link_"){

    if(is.character(bedpe)){
    bedpe <- read.delim(bedpe, stringsAsFactors=F, sep='\t',header=header)} else if(is.data.frame(bedpe)){ bedpe <- bedpe}

    if(filter_trans==TRUE){ bedpe <- bedpe[which(bedpe[,ranges_1[1]] ==bedpe[,ranges_2[1]]),]
                        print("Finished filtering trans interactions") } else{ bedpe <- bedpe}

##    str(bedpe)
    metadata_cols <- setdiff(1:ncol(bedpe),c(ranges_1,ranges_2))

    table1 <- bedpe[c(ranges_1,metadata_cols)]
    table2 <- bedpe[c(ranges_2,metadata_cols)]
   ## print(nrow(bedpe))
    link_names <- if(!is.null(ID_col)){
       link_names <- bedpe[,ID_col]} else { link_names <-  paste0(el_names,1:nrow(bedpe))}

   ## str(table1)
   ## str(table2)
   ## str(link_names)

    colnames(table1)[1:3] <- c("Chr","Start","End")
    colnames(table2)[1:3] <- c("Chr","Start","End")

    table1$Name <- link_names
    table2$Name <- link_names

    table_out <- rbind(table1, table2, stringsAsFactors=F)
##    str(table_out)

    gr_out <- table_to_granges(table_out)
    ##print(gr_out)
    out_list <- split(gr_out,gr_out$Name)


##    out_list <- list();
##    for(i in 1:length(link_names)){
##        if(!is.null(strand)){ gr1 <- table_to_granges(table1[i,],strand=last(ranges_1))
##                              gr2 <-table_to_granges(table2[i,],strand=last(ranges_1))} else{
##                                  gr1 <- table_to_granges(table1[i,])
##                                  gr2 <- table_to_granges(table2[i,])}

                                    ##print(gr1)
                                    ##print(gr2)
##                                    gr <- union(gr1,gr2)
##                                    gr$Name <- link_names[i]
                                    ##print(gr)
##                                    out_list[[link_names[[i]]]] <- gr
##                                }

##    out_list <- GRangesList(out_list)
    return(out_list)

}

sort_gr <- function(gr,...){
    require(GenomicRanges)
    gr <- sort(sortSeqlevels(gr),...)

    return(gr)}





gr_to_bed <- function(gr,outfile=NULL, metadata=FALSE,verbose=FALSE,header=FALSE,biostrings=FALSE,append=FALSE,drop_strand=TRUE){
    require(GenomicRanges)

if(length(gr) >0){
    Chrom <- as.character(gr@seqnames)
    Start <- gr@ranges@start
    End <- gr@ranges@start+gr@ranges@width
    Strand <- gr@strand

                                   if(metadata==FALSE){
                                       df <- data.frame(Chrom,Start,End)}
                                   else if(metadata==TRUE){
                                       metadata <- as.data.frame(gr@elementMetadata@listData,stringsAsFactors=F)
                                       if(biostrings==TRUE){
                                           metadata <- metadata[,setdiff(colnames(metadata),c("ALT.group","ALT.group_name"))]}
                                       df <- data.frame(Chrom,Start,End,Strand,metadata,stringsAsFactors=F)

                                       }
    if(is.null(outfile)){
        return(df)} else{
            if(drop_strand== TRUE){ df["Strand"] <- NULL} else{ df <- df}
            write.table(df, outfile,quote=F, sep='\t', row.names=F,col.names=header,append=append)}
                                   if(verbose== TRUE) {print(head(df))} else{}} else{ print("Zero length gr. Stopping")}}


collapse_gr_keep_metadata <- function(gr){
    require(GenomicRanges)
    reduced <- reduce(gr)
    final_gr <- reduced
    for (i in 1:length(reduced)){
        overlap_enhancer = findOverlaps(reduced[i], gr)
        row_names = vector()
        for( j in 1:length(overlap_enhancer)){
            names = gr[overlap_enhancer[j]@to]$Gene
            row_names = c(names, row_names)

        }

        final_gr$nearby[i] = row_names
    }
    return(final_gr)}


intersect_with_metadata <- function(gr1,gr2, ignore.strand=FALSE){
    require(GenomicRanges)
    if(ignore.strand == FALSE){ out <- gr1[findOverlaps(gr1,gr2)@from,]} else if( ignore.strand ==TRUE){ out <- gr1[findOverlaps(gr1,gr2,ignore.strand=TRUE)@from,]}
    ; return(out)}

setdiff_with_metadata <- function(gr1,gr2){ require(GenomicRanges)
                                            int <- findOverlaps(gr1,gr2); out <- gr1[setdiff(1:length(gr1),unique(int@from)),]; return(out)}

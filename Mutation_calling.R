call_mutations_from_bam <- function(bam, outfile, reference=NULL,regions=NULL,report_all_bases=FALSE){
### bam path to bam file
### where you want your output to go, output will be in vcf format (string)
### reference is a path to the reference fasta file that the bam was aligned to. (string)
### regions file is a path to a bed file (without headers) or a coordinate in chr:start-end format (string)
### report_all_bases will either return all bases in the regions of interest or only the bases that are different from the reference genome. (Boolean)

    if(!file.exists(reference)){ print("Not a valid path and/or reference");stop()}
    ##This requires a working installation of bcftools and samtools
    if(!is.null(regions)){
        if(is.character(regions)){
            file_check <- file.info(regions)
            if(is.na(file_check$size)){
                coord_check <- grepl("[0-9]:[0-9]{1,9}-[0-9]{1,9}|chr[a-zA-Z0-9_]{1,100}:[0-9]{1,9}-[0-9]{1,9}",regions)
                if(coord_check){ region_flag= "-r"} else{ print("Not a valid location. Region format is \"chr:start-end");stop()}} else{ region_flag="-R"}}
        if(!report_all_bases){
        cmd <- paste("bcftools mpileup",bam,region_flag,regions,"-f",reference, "|bcftools call -vmO v -o",outfile,sep=" ")} else { cmd <- paste("bcftools mpileup",bam,region_flag,regions,"-f",reference, "|bcftools call -mO v -o",outfile,sep=" ")

                                                                                                                            }
        system(cmd)
    } else if(!report_all_bases){
        cmd <- paste("bcftools mpileup",bam,"-f",reference, "|bcftools call -vmO v -o",outfile,sep=" ")} else {cmd <- paste("bcftools mpileup",bam,region_flag,regions,"-f",reference, "|bcftools call -mO v -o",outfile,sep=" ")}
    system(cmd)

##    if(!grepl(".gz", outfile)){ intfile <- paste0(outfile,".gz")
##    system(paste("mv",intfile, outfile))
##      system(paste0("gunzip ",outfile))
##      } else{ outfile <- outfile}

    print(outfile)
}

concatenate_vcfs <- function(input,sample_ID=NULL,verbose=F){
### input can be a path to a directory or a vector of paths
    ### sample id is whether to use the filenames as the sample identifier for each vcf or a single identifier for the provided vcfs
  out_vcf <- list()
  require(VariantAnnotation)
  if(is.directory(input)){
    file_list <- list.files(input,".vcf",full.names=TRUE)

    for(i in file_list){ in_vcf <- tryCatch({readVcfAsVRanges(i)}, error=function(e) { print("No variants")})
        if(!is.null(sample_ID) && sample_ID!="filename" && class(in_vcf) == "VRanges"){ sampleNames(in_vcf) <- sample_ID} else if(sample_ID=="filename" && class(in_vcf)=="VRanges"){ sampleNames(in_vcf) <- i}
        if(class(in_vcf)=="VRanges"){
    in_vcf$Sample <- unique(sampleNames(in_vcf))
    out_vcf[[as.character(in_vcf$Sample[1])]] <- in_vcf
        } else{ print("No variants")}
        if(verbose){ print(i)}}
  }
  else if( is.vector(input)){
    if(all(file.exists(input))){
      for(i in input){ in_vcf <- tryCatch({readVcfAsVRanges(i)}, error=function(e) { print("No variants")})
      if(!is.null(sample_ID) && sample_ID!="filename" && class(in_vcf) =="VRanges"){ sampleNames(in_vcf) <- sample_ID} else if(sample_ID=="filename" && class(in_vcf)=="VRanges"){ sampleNames(in_vcf) <- i}

      if(class(in_vcf)=="VRanges"){
      in_vcf$Sample <- unique(sampleNames(in_vcf))
      out_vcf[[as.character(in_vcf$Sample[1])]] <- in_vcf  } else{ print("No variants")}

    if(verbose){print(i)}}}}

  out <- collapse_granges_list(out_vcf)
  return(out)
}

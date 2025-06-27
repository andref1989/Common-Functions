call_mutations_from_bam <- function(bam, outfile, reference=NULL,regions=NULL,report_all_bases=FALSE){

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

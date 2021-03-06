library(CNprep)

#
# Retrieve information of cytobands (downloaded for genome hg19)
#
retrieveCytobands <- function(dir = "cytoBand.txt"){
  cytobands <- read.table(dir, header=F, sep = "\t", stringsAsFactors = F)
  names(cytobands) <- c("chrom", "start", "end", "cytoloc", "stain")
  return(cytobands)
}

#
# With the chromosome and chromosome location, retrieve the cytoband location
#
findCytolocation <- function(cytobands, chrom, chrom.position){
  chrom = if(chrom == 23) "X" else if (chrom == 24 ) "Y" else chrom
  row <- cytobands[cytobands$chrom == paste("chr", chrom, sep = "") & cytobands$start <= chrom.position & cytobands$end >= chrom.position, ]
  returnme <- data.frame(row)$cytoloc
  return(returnme)
}

#
# Retrieve the norminput argument for CNprep::CNpreprocessing()
#
retrieveNormInput <- function(normalSegments){
  norminput <- data.frame(stringsAsFactors = FALSE)
  for(normalSegments.index in seq(1, nrow(normalSegments))){
    norminput.entry <- data.frame(length = normalSegments[normalSegments.index, ]$end - normalSegments[normalSegments.index, ]$start, segmedian = normalSegments[normalSegments.index, ]$cnlr)
    norminput <- rbind(norminput, norminput.entry)
  }
  return(norminput)
}

#
# Filter norminput from artifacts
#
filterNormInput <- function(norminput, length_threshold=10000000){
  # Determine cutoff
  oNI <- norminput[order(norminput$length), ]
  #plot(oNI$length, oNI$segmedian, ylim = c(-2.25, 2.25), xlim = c(0, 300000000))
  #dev.off()
  nNI <- norminput[norminput$length > length_threshold & abs(norminput$segmedian) < 0.5,]
  #plot(nNI$length, nNI$segmedian, ylim = c(-2.25, 2.25), xlim = c(0, 300000000))
  return(nNI)
}

#
# Retrieve the seginput argument for CNprep::CNpreprocessing()
#
retrieveSegInput <- function(facets_segment_data, sample, chromosomeSizes, cytobands){
  seginput <- data.frame(stringsAsFactors = FALSE)
  # Iterate through each segment
  for(facets_segment_data.index in seq(1, nrow(facets_segment_data))){
    # Get absolute position of segment
    abs_position <- chromsomeToAbsoluteBPConversionForSingleEntry(facets_segment_data[facets_segment_data.index,]$X.chrom., facets_segment_data[facets_segment_data.index,]$X.start., facets_segment_data[facets_segment_data.index,]$X.end., chromosomeSizes)
    
    probes.start = 0
    probes.end = 0
    for(i in seq(1, facets_segment_data.index)){
      if(i != facets_segment_data.index){
        probes.start <- probes.start + facets_segment_data[i,]$X.num.mark.
      }
      probes.end <- probes.end + facets_segment_data[i,]$X.num.mark.
    }
    probes.start <- probes.start + 1
    
    
    cytoband.my.start <- findCytolocation(cytobands = cytobands, chrom = facets_segment_data[facets_segment_data.index,]$X.chrom., chrom.position = facets_segment_data[facets_segment_data.index,]$X.start.)
    cytoband.my.end <- findCytolocation(cytobands = cytobands, chrom = facets_segment_data[facets_segment_data.index,]$X.chrom., chrom.position = facets_segment_data[facets_segment_data.index,]$X.end.)
    
    seginput.entry <- data.frame(ID = sample, start = probes.start, end = probes.end, 
                                 num.probes = facets_segment_data[facets_segment_data.index,]$X.num.mark., seg.median = facets_segment_data[facets_segment_data.index,]$X.cnlr.median., 
                                 chrom = facets_segment_data[facets_segment_data.index,]$X.chrom., chrom.pos.start = facets_segment_data[facets_segment_data.index,]$X.start., 
                                 chrom.pos.end = facets_segment_data[facets_segment_data.index,]$X.end., cytoband.start = cytoband.my.start, 
                                 cytoband.end = cytoband.my.end, abs.pos.start = abs_position$start,
                                 abs.pos.end = abs_position$end)
    
    
    seginput <- rbind(seginput, seginput.entry)
  }
  return(seginput)
}

#
# Retrieve the ratinput argument for CNprep::CNpreprocessing()
#
retrieveRatInput <- function(facets_snp_data, sample){
  ratinput <- data.frame(facets_snp_data$X.cnlr.)
  names(ratinput) <- c(sample)
  return(ratinput)
}

#
# Wrapper function for CNprep::CNpreprocessing()
#
runCNpreprocessing <- function(seginput, ratinput, norminput,
                               blsize=50, minjoin=0.25, cweight=0.4, bstimes=50,
                               chromrange=1:22, distrib="Rparallel", njobs=4, modelNames="E", ntrial= 10){
  segtable<-CNpreprocessing(segall=seginput,ratall=ratinput,"ID","start","end",
                            chromcol="chrom",bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",blsize=blsize,
                            minjoin=minjoin,cweight=cweight,bstimes=bstimes,chromrange=chromrange,distrib=distrib,njobs=njobs,
                            modelNames=modelNames,normalength=norminput[,1],normalmedian=norminput[,2], ntrial=ntrial)
  return(segtable)
}

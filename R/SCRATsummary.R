#' SCRATsummary
#' 
#' Compile SCRAT summary table
#'
#' This function will compile a SCRAT summary table from bam files. The results should be the same as run on GUI.
#' @param dir The folder where the bam files are stored. If bamfile is NULL, all bam files within the folder will be analyzed.
#' @param genome The mapped genome of the bam files. Should be one the following: "hg19", "hg38", "mm9", "mm10"
#' @param bamfile A character vector of bam files. If NULL, all files in dir will be included.
#' @param singlepair Whether the original sequencing files are single-end or paired-end. Should be one of the following: "automated", "single", "pair". Default is "automated" where SCRAT will automatically determine the type.
#' @param removeblacklist Logical value indicating whether black list regions should be removed.
#' @param log2transform Logical value indicating whether the read counts should be log2 transformed (after adding pseudo-count of 1).
#' @param featurelist A character vector specifying what kind of features should be considered. Should be from the following: "GENE","ENCL","MOTIF_TRANSFAC","MOTIF_JASPAR","GSEA". By default all features are included. Note that "GSEA" features could be slow to run.
#' @param Genestarttype For "GENE" features, type of starting site. Should be one of the following: "TSSup", "TSSdown", "TESup", "TESdown". The four options stands for TSS upstream, TSS downstream, TES upstream and TES downstream.
#' @param Geneendtype For "GENE" features, type of ending site. Options same as Genestarttype
#' @param Genestartbp For "GENE" features, how many base pairs away from starting TSS/TES. For example, Genestarttype="TSSup" and Genestartbp=500 means 500 bp upstream of TSS.
#' @param Geneendbp For "GENE" features, how many base pairs away from ending TSS/TES.
#' @param ENCLclunum Number of clusters for ENCL features. Should be one of 1000, 2000 and 5000
#' @param Motifflank Defines the size of flanking region of motif sites in base pairs.
#' @param GSEAterm The GSEA terms included in the analysis. Only useful when "GSEA" is included in featurelist. Should be one of the following: "h.all","c1.all","c2.cgp","c2.cp","c3.mir","c3.tft","c4.cgn","c4.cm","c5.bp","c5.cc","c5.mf","c6.all","c7.all".
#' @param GSEAstarttype For "GSEA" features, type of starting site.
#' @param GSEAendtype For "GSEA" features, type of ending site.
#' @param GSEAstartbp For "GSEA" features, how many base pairs away from starting TSS/TES.
#' @param GSEAendbp For "GSEA" features, how many base pairs away from ending TSS/TES.
#' @export
#' @import GenomicAlignments
#' @author Zhicheng Ji, Weiqiang Zhou, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' \dontrun{
#'    SCRATsummary(dir="bamfiledir",genome="hg19")
#' }

SCRATsummary <- function(dir="",genome,bamfile=NULL,singlepair="automated",removeblacklist=T,log2transform=T,featurelist=c("GENE","ENCL","MOTIF_TRANSFAC","MOTIF_JASPAR","GSEA"),Genestarttype="TSSup",Geneendtype="TSSdown",Genestartbp=3000,Geneendbp=1000,ENCLclunum=2000,Motifflank=100,GSEAterm="c5.bp",GSEAstarttype="TSSup",GSEAendtype="TSSdown",GSEAstartbp=3000,GSEAendbp=1000) {
      if (is.null(bamfile)) {
            bamfile <- list.files(dir,pattern = ".bam$")
      }
      datapath <- system.file("extdata",package=paste0("SCRATdata",genome))
      bamdata <- list()
      
      for (i in bamfile) {
            filepath <- file.path(dir,i)
            if (singlepair=="automated") {
                  bamfile <- BamFile(filepath)
                  tmpsingle <- readGAlignments(bamfile)
                  tmppair <- readGAlignmentPairs(bamfile)
                  pairendtf <- testPairedEndBam(bamfile)
                  if (pairendtf) {
                        tmp <- tmppair
                        startpos <- pmin(start(first(tmp)),start(last(tmp)))
                        endpos <- pmax(end(first(tmp)),end(last(tmp)))
                        tmp <- GRanges(seqnames=seqnames(tmp),IRanges(start=startpos,end=endpos))
                  } else {
                        tmp <- GRanges(tmpsingle)            
                  }      
            } else if (singlepair=="single") {
                  tmp <- GRanges(readGAlignments(filepath))                  
            } else if (singlepair=="pair") {
                  tmp <- readGAlignmentPairs(filepath)
                  startpos <- pmin(start(first(tmp)),start(last(tmp)))
                  endpos <- pmax(end(first(tmp)),end(last(tmp)))
                  tmp <- GRanges(seqnames=seqnames(tmp),IRanges(start=startpos,end=endpos))
            }
            if (removeblacklist) {
                  load(paste0(datapath,"/gr/blacklist.rda"))
                  tmp <- tmp[-as.matrix(findOverlaps(tmp,gr))[,1],]                  
            }
            bamdata[[i]] <- tmp 
      }
      bamsummary <- sapply(bamdata,length)
      allres <- NULL
      datapath <- system.file("extdata",package=paste0("SCRATdata",genome))
      if ("GENE" %in% featurelist) {      
            print("Processing GENE features")
            load(paste0(datapath,"/gr/generegion.rda"))
            if (Genestarttype == "TSSup") {
                  grstart <- ifelse(as.character(strand(gr))=="+",start(gr)-as.numeric(Genestartbp),end(gr)+as.numeric(Genestartbp))
            } else if (Genestarttype == "TSSdown") {
                  grstart <- ifelse(as.character(strand(gr))=="+",start(gr)+as.numeric(Genestartbp),end(gr)-as.numeric(Genestartbp))
            } else if (Genestarttype == "TESup") {
                  grstart <- ifelse(as.character(strand(gr))=="+",end(gr)-as.numeric(Genestartbp),start(gr)+as.numeric(Genestartbp))
            } else if (Genestarttype == "TESdown") {
                  grstart <- ifelse(as.character(strand(gr))=="+",end(gr)+as.numeric(Genestartbp),start(gr)-as.numeric(Genestartbp))
            }
            if (Geneendtype == "TSSup") {
                  grend <- ifelse(as.character(strand(gr))=="+",start(gr)-as.numeric(Geneendbp),end(gr)+as.numeric(Geneendbp))
            } else if (Geneendtype == "TSSdown") {
                  grend <- ifelse(as.character(strand(gr))=="+",start(gr)+as.numeric(Geneendbp),end(gr)-as.numeric(Geneendbp))
            } else if (Geneendtype == "TESup") {
                  grend <- ifelse(as.character(strand(gr))=="+",end(gr)-as.numeric(Geneendbp),start(gr)+as.numeric(Geneendbp))
            } else if (Geneendtype == "TESdown") {
                  grend <- ifelse(as.character(strand(gr))=="+",end(gr)+as.numeric(Geneendbp),start(gr)-as.numeric(Geneendbp))
            }
            ngr <- names(gr)
            gr <- GRanges(seqnames=seqnames(gr),IRanges(start=pmin(grstart,grend),end=pmax(grstart,grend)))      
            names(gr) <- ngr
            tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
            tmp <- sweep(tmp,2,bamsummary,"/") * 10000
            if (log2transform) {
                  tmp <- log2(tmp + 1)
            }
            allres <- rbind(allres,tmp)      
      }
      if ("ENCL" %in% featurelist) {      
            print("Processing ENCL features")
            load(paste0(datapath,"/gr/ENCL",ENCLclunum,".rda"))
            tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
            tmp <- sweep(tmp,2,bamsummary,"/") * 10000
            if (log2transform) {
                  tmp <- log2(tmp + 1)
            }
            tmp <- tmp[rowSums(tmp) > 0,,drop=F]       
            allres <- rbind(allres,tmp)      
      }      
      if ("MOTIF_TRANSFAC" %in% featurelist) {  
            print("Processing MOTIF_TRANSFAC features")
            load(paste0(datapath,"/gr/transfac1.rda"))
            gr <- flank(gr,as.numeric(Motifflank),both = T)
            tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
            tmp <- sweep(tmp,2,bamsummary,"/") * 10000
            if (log2transform) {
                  tmp <- log2(tmp + 1)
            }
            allres <- rbind(allres,tmp)
            load(paste0(datapath,"/gr/transfac2.rda"))
            gr <- flank(gr,as.numeric(Motifflank),both = T)
            tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
            tmp <- sweep(tmp,2,bamsummary,"/") * 10000
            if (log2transform) {
                  tmp <- log2(tmp + 1)
            }
            allres <- rbind(allres,tmp)
            if (genome %in% c("hg19","hg38")) {
                  load(paste0(datapath,"/gr/transfac3.rda"))
                  gr <- flank(gr,as.numeric(Motifflank),both = T)
                  tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
                  tmp <- sweep(tmp,2,bamsummary,"/") * 10000
                  if (log2transform) {
                        tmp <- log2(tmp + 1)
                  }
                  allres <- rbind(allres,tmp)
            }
      }
      if ("MOTIF_JASPAR" %in% featurelist) {     
            print("Processing MOTIF_JASPAR features")
            load(paste0(datapath,"/gr/jaspar1.rda"))
            gr <- flank(gr,as.numeric(Motifflank),both = T)
            tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
            tmp <- sweep(tmp,2,bamsummary,"/") * 10000
            if (log2transform) {
                  tmp <- log2(tmp + 1)
            }
            allres <- rbind(allres,tmp)
            load(paste0(datapath,"/gr/jaspar2.rda"))
            gr <- flank(gr,as.numeric(Motifflank),both = T)
            tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
            tmp <- sweep(tmp,2,bamsummary,"/") * 10000
            if (log2transform) {
                  tmp <- log2(tmp + 1)
            }
            allres <- rbind(allres,tmp)
      }
      if ("GSEA" %in% featurelist) {
            print("Processing GSEA features")
            #for (i in c("h.all","c1.all","c2.cgp","c2.cp","c3.mir","c3.tft","c4.cgn","c4.cm","c5.bp","c5.cc","c5.mf","c6.all","c7.all")) {
            for (i in GSEAterm) {
                  load(paste0(datapath,"/gr/GSEA",i,".rda"))
                  allgr <- gr
                  for (sgrn in names(allgr)) {
                        gr <- allgr[[sgrn]]
                        if (GSEAstarttype == "TSSup") {
                              grstart <- ifelse(as.character(strand(gr))=="+",start(gr)-as.numeric(GSEAstartbp),end(gr)+as.numeric(GSEAstartbp))
                        } else if (GSEAstarttype == "TSSdown") {
                              grstart <- ifelse(as.character(strand(gr))=="+",start(gr)+as.numeric(GSEAstartbp),end(gr)-as.numeric(GSEAstartbp))
                        } else if (GSEAstarttype == "TESup") {
                              grstart <- ifelse(as.character(strand(gr))=="+",end(gr)-as.numeric(GSEAstartbp),start(gr)+as.numeric(GSEAstartbp))
                        } else if (GSEAstarttype == "TESdown") {
                              grstart <- ifelse(as.character(strand(gr))=="+",end(gr)+as.numeric(GSEAstartbp),start(gr)-as.numeric(GSEAstartbp))
                        }
                        if (GSEAendtype == "TSSup") {
                              grend <- ifelse(as.character(strand(gr))=="+",start(gr)-as.numeric(GSEAendbp),end(gr)+as.numeric(GSEAendbp))
                        } else if (GSEAendtype == "TSSdown") {
                              grend <- ifelse(as.character(strand(gr))=="+",start(gr)+as.numeric(GSEAendbp),end(gr)-as.numeric(GSEAendbp))
                        } else if (GSEAendtype == "TESup") {
                              grend <- ifelse(as.character(strand(gr))=="+",end(gr)-as.numeric(GSEAendbp),start(gr)+as.numeric(GSEAendbp))
                        } else if (GSEAendtype == "TESdown") {
                              grend <- ifelse(as.character(strand(gr))=="+",end(gr)+as.numeric(GSEAendbp),start(gr)-as.numeric(GSEAendbp))
                        }
                        ngr <- names(gr)
                        gr <- GRanges(seqnames=seqnames(gr),IRanges(start=pmin(grstart,grend),end=pmax(grstart,grend)))      
                        names(gr) <- ngr   
                        allgr[[sgrn]] <- gr
                  }
                  gr <- allgr
                  tmp <- sapply(bamdata,function(i) countOverlaps(gr,i))
                  tmp <- sweep(tmp,2,bamsummary,"/") * 10000
                  if (log2transform) {
                        tmp <- log2(tmp + 1)
                  }
                  allres <- rbind(allres,tmp)
            }
      }
      allres[rowSums(allres) > 0,]
}




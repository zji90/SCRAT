######################################################
##                      SCRAT                       ##
##             Interactive User Interface           ##
##                     Server File                  ##
##   Author:Zhicheng Ji, Weiqiang Zhou, Hongkai Ji  ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################
library(shiny)
library(GenomicAlignments)
library(ggplot2)
library(reshape2)
library(gplots)
library(pheatmap)
library(scatterD3)
library(DT)
library(mclust)
library(tsne)

options(shiny.maxRequestSize=10*1024^10)

shinyServer(function(input, output,session) {
      
      Maindata <- reactiveValues()
      
      output$InputGenomeui <- renderUI({
            genomelist <- list()
            if (require("SCRATdatahg19")) {
                  genomelist[["hg19 (Human)"]] <- c("hg19")
            }
            if (require("SCRATdatahg38")) {
                  genomelist[["hg38 (Human)"]] <- c("hg38")
            }
            if (require("SCRATdatamm9")) {
                  genomelist[["mm9 (Mouse)"]] <- c("mm9")
            }
            if (require("SCRATdatamm10")) {
                  genomelist[["mm10 (Mouse)"]] <- c("mm10")
            }
            
            selectInput("InputGenome","",genomelist)
            #list("hg19 (Human)"="hg19","hg38 (Human)"="hg38","mm9 (Mouse)"="mm9","mm10 (Mouse)"="mm10")
      })
      
      
      
      ### Input ###
      observe({
            if (input$Inputreadin > 0) {
                  isolate({
                        FileHandle <- input$InputFile
                        Maindata$Rawbamfile <- list()
                        Maindata$Rawbampairtf <- list()
                        filenum <- length(FileHandle$datapath)
                        if (!is.null(FileHandle)) {
                              withProgress(message = 'Reading in',detail="0%",{
                                    datapath <- system.file("extdata",package=paste0("SCRATdata",input$InputGenome))
                                    load(paste0(datapath,"/gr/blacklist.rda"))
                                    for (i in 1:filenum) {
                                          incProgress(1/filenum,detail=paste0(round(i/filenum*100),"%"))
                                          name <- FileHandle$name[i]
                                          bamfile <- BamFile(FileHandle$datapath[i])
                                          tmpsingle <- readGAlignments(bamfile)
                                          tmppair <- readGAlignmentPairs(bamfile)
                                          pairendtf <- testPairedEndBam(bamfile)
                                          if (pairendtf) {
                                                tmp <- tmppair
                                                startpos <- pmin(start(first(tmp)),start(last(tmp)))
                                                endpos <- pmax(end(first(tmp)),end(last(tmp)))
                                                tmp <- GRanges(seqnames=seqnames(tmp),IRanges(start=startpos,end=endpos))
                                                Maindata$Rawbampairtf[[name]] <- "paired-end"
                                          } else {
                                                tmp <- GRanges(tmpsingle)            
                                                Maindata$Rawbampairtf[[name]] <- "single-end"
                                          }
                                          if (input$Inputblacklist) {
                                                overres <- findOverlaps(tmp,gr)
                                                if ("from" %in% slotNames(overres)) {
                                                      tmp <- tmp[-overres@from,]                                                            
                                                } else {
                                                      tmp <- tmp[-overres@queryHits,]                                                            
                                                }
                                          }
                                          Maindata$Rawbamfile[[name]] <- tmp
                                    }
                              })
                              Maindata$Rawbamlength <- sapply(Maindata$Rawbamfile,length)
                        }
                  })
            }
      })
      
      observe({
            if (!is.null(Maindata$Rawbamfile)) {
                  selected <- Maindata$Rawbamlength[Maindata$Rawbamlength >= as.numeric(input$InputFilterreads)]
                  selected <- selected[!names(selected) %in% input$InputFiltersample]
                  Maindata$bamfile <- Maindata$Rawbamfile[names(Maindata$Rawbamfile) %in% names(selected)]
                  Maindata$bamsummary <- data.frame(BAM=names(selected),Reads=selected,Type=sapply(names(selected),function(i) Maindata$Rawbampairtf[[i]]),stringsAsFactors = F)
            }
      })
      
      output$Inputbamnumui <- renderUI(
            if (!is.null(Maindata$Rawbamfile)) {
                  tagList(
                        h5(paste(length(Maindata$Rawbamfile),"bam files read in")),
                        h5(paste(length(Maindata$bamfile),"bam files retained")),
                        h5("Reads for each bam file:")
                  )
            }
      )
      
      output$Inputbamsummary <- DT::renderDataTable(DT::datatable(Maindata$bamsummary,filter='top',rownames = F, options = list(columnDefs = list(list(className="dt-body-left","targets"="_all")))))
      
      output$InputFiltersampleui <- renderUI({
            selectInput("InputFiltersample","Exclude specific samples",names(Maindata$Rawbamfile),multiple = T)
      })
      
      observe({
            if (input$Inputnextstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 2: Feature summarization")
                  })
            }
      })
      
      ### Summarizing ###
      
      observe({input$SumuploadFile})
      
      observe({
            if (input$Sumuploadreadin > 0) {
                  isolate({
                        if (input$SumuploadFilemultipletf=="One") {
                              FileHandle <- input$SumuploadFile      
                              if (!is.null(FileHandle)) {
                                    tmp <- read.table(FileHandle$datapath,as.is=T,sep="\t")
                                    if (ncol(tmp)==3) {
                                          gr <- GRanges(seqnames=tmp[,1],IRanges(start=tmp[,2],end=tmp[,3]))
                                          names(gr) <- paste0("UPLOAD:",tmp[,1],":",tmp[,2],"-",tmp[,3])
                                    } else if (ncol(tmp)==4) {
                                          gr <- GRangesList()
                                          for (i in unique(tmp[,4])) {
                                                tmpgr <- GRanges(seqnames=tmp[tmp[,4]==i,1],IRanges(start=tmp[tmp[,4]==i,2],end=tmp[tmp[,4]==i,3]))
                                                gr[[paste0("MOTIF:",i)]] <- tmpgr
                                          }
                                    }
                                    Maindata$uploadgr <- gr      
                              }
                        } else {
                              FileHandle <- input$SumuploadFile2
                              if (!is.null(FileHandle)) {
                                    gr <- GRangesList()
                                    for (i in 1:length(FileHandle$datapath)) {
                                          tmp <- read.table(FileHandle$datapath[[i]],as.is=T,sep="\t")
                                          gr[[paste0("UPLOAD",i)]] <- GRanges(seqnames=tmp[,1],IRanges(start=tmp[,2],end=tmp[,3]))
                                    }
                                    Maindata$uploadgr <- gr
                              }
                        }
                  })
            }
      })
      
      output$Sumuploaddetail <- renderText({
            if (!is.null(Maindata$uploadgr)) {
                  if (input$SumuploadFilemultipletf=="One") {
                        paste("Uploaded bed file has",length(Maindata$uploadgr),"loci.")
                  } else {
                        paste("Uploaded",length(Maindata$uploadgr),"bed files.")     
                  }
            }
      })
      
      observe({
            if (input$Sumuploadsumtablereadin > 0) {
                  isolate({
                        withProgress(message = 'Reading in...',{
                              FileHandle <- input$SumuploadsumtableFile                        
                              if (!is.null(FileHandle)) {
                                    tmp <- read.table(FileHandle$datapath,sep="\t",header=T,as.is=T,row.names = 1,check.names=FALSE)                                                                        
                                    tmp <- tmp[,colnames(tmp) != "CV"]
                                    Maindata$allsumtable <- as.matrix(tmp)
                                    Maindata$sumtablenametype <- sapply(row.names(Maindata$allsumtable),function(i) strsplit(i,":")[[1]][1])                              
                              }
                        })
                  })
            }
      })
      
      observe({
            if (input$Sumrunbutton > 0)
                  isolate({
                        if (length(input$Sumselectmet) > 0) {
                              allres <- NULL
                              datapath <- system.file("extdata",package=paste0("SCRATdata",input$InputGenome))
                              # if ("TSS" %in% input$Sumselectmet) {
                              #       withProgress(message = 'Counting Overlaps for TSS',{
                              #             load(paste0(datapath,"/gr/generegion.rda"))
                              #             gr <- promoters(resize(gr,1),as.numeric(input$SumTSSupregionbp),as.numeric(input$SumTSSdownregionbp))
                              #             tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                              #             tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                              #             if (input$Sumlogtf) {
                              #                   tmp <- log2(tmp + 1)
                              #             }
                              #             if (input$Sumfiltertf) {
                              #                   tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]
                              #                   cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                              #                   tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                              #             }
                              #             allres <- rbind(allres,tmp)
                              #       })
                              # }
                              if ("ENCL" %in% input$Sumselectmet) {
                                    withProgress(message = 'Counting Overlaps for ENCODE Cluster',{
                                          load(paste0(datapath,"/gr/ENCL",input$SumENCLclunum,".rda"))
                                          tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                          tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                          if (input$Sumlogtf) {
                                                tmp <- log2(tmp + 1)
                                          }
                                          if (input$Sumfiltertf) {
                                                tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                          }                        
                                          allres <- rbind(allres,tmp)
                                    })
                              }
                              if ("MOTIF" %in% input$Sumselectmet) {
                                    if ("TRANSFAC" %in% input$SumMOTIFselect) {
                                          withProgress(message = 'Counting Overlaps for TRANSFAC',{
                                                load(paste0(datapath,"/gr/transfac1.rda"))
                                                gr <- flank(gr,as.numeric(input$SumMOTIFflank),both = T)
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                      cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                      tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                                }                        
                                                allres <- rbind(allres,tmp)
                                                load(paste0(datapath,"/gr/transfac2.rda"))
                                                gr <- flank(gr,as.numeric(input$SumMOTIFflank),both = T)
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                      cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                      tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                                }                        
                                                allres <- rbind(allres,tmp)
                                                if (input$InputGenome %in% c("hg19","hg38")) {
                                                      load(paste0(datapath,"/gr/transfac3.rda"))
                                                      gr <- flank(gr,as.numeric(input$SumMOTIFflank),both = T)
                                                      tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                      tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                                      if (input$Sumlogtf) {
                                                            tmp <- log2(tmp + 1)
                                                      }
                                                      if (input$Sumfiltertf) {
                                                            tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                            cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                            tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                                      }                        
                                                      allres <- rbind(allres,tmp)
                                                }
                                          })
                                    }
                                    if ("JASPAR" %in% input$SumMOTIFselect) {
                                          withProgress(message = 'Counting Overlaps for JASPAR',{
                                                load(paste0(datapath,"/gr/jaspar1.rda"))
                                                gr <- flank(gr,as.numeric(input$SumMOTIFflank),both = T)
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                      cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                      tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                                }                        
                                                allres <- rbind(allres,tmp)
                                                load(paste0(datapath,"/gr/jaspar2.rda"))
                                                gr <- flank(gr,as.numeric(input$SumMOTIFflank),both = T)
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                      cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                      tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                                }                        
                                                allres <- rbind(allres,tmp)
                                          })
                                    }
                              }
                              if ("GSEA" %in% input$Sumselectmet) {
                                    withProgress(message = 'Counting Overlaps for GSEA',{
                                          for (i in input$SumGSEAselect) {
                                                load(paste0(datapath,"/gr/GSEA",i,".rda"))
                                                gr <- promoters(resize(gr,1),as.numeric(input$SumGSEAupregionbp),as.numeric(input$SumGSEAdownregionbp))             
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                      cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                      tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                                }                        
                                                allres <- rbind(allres,tmp)
                                          }
                                    })
                              }
                              if ("Upload" %in% input$Sumselectmet & !is.null(Maindata$uploadgr)) {
                                    withProgress(message = 'Counting Overlaps for Upload bed',{
                                          tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(Maindata$uploadgr,i))
                                          tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * 10000
                                          if (input$Sumlogtf) {
                                                tmp <- log2(tmp + 1)
                                          }
                                          if (input$Sumfiltertf) {
                                                tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                cv <- apply(tmp,1,function(i) sd(i)/mean(i))
                                                tmp <- tmp[cv >= as.numeric(input$Sumfiltercv),,drop=F]
                                          }                        
                                          allres <- rbind(allres,tmp)
                                    })
                              }
                              Maindata$allsumtable <- allres
                              Maindata$sumtablenametype <- sapply(row.names(allres),function(i) strsplit(i,":")[[1]][1])
                        }
                  })
      })
      
      output$Sumtypesummarytable <- DT::renderDataTable({
            if (!is.null(Maindata$sumtablenametype)) {
                  tmp <- table(Maindata$sumtablenametype)
                  DT::datatable(data.frame(Type=names(tmp),Number=as.vector(tmp)),filter='top',rownames = F, options = list(columnDefs = list(list(className="dt-body-left","targets"="_all"))))
            }
      })
      
      output$Sumtableoutput <- DT::renderDataTable({
            if (!is.null(Maindata$allsumtable)) {
                  if (input$Sumaddcv) {
                        cv <- apply(Maindata$allsumtable,1,sd)/rowMeans(Maindata$allsumtable)                  
                        tmp <- cbind(row.names(Maindata$allsumtable),cv,Maindata$allsumtable)
                        colnames(tmp)[1:2] <- c("Feature","CV")
                  } else {                        
                        tmp <- cbind(row.names(Maindata$allsumtable),Maindata$allsumtable)
                        colnames(tmp)[1] <- "Feature"
                  }     
                  DT::datatable(tmp,filter="top",rownames = F, options = list(columnDefs = list(list(className="dt-body-left","targets"="_all"))))
            }            
      })
      
      output$Sumtabledownload <- downloadHandler(
            filename = function() { "Summarization.txt" },
            content = function(file) {       
                  if (input$Sumaddcv) {
                        cv <- apply(Maindata$allsumtable,1,sd)/rowMeans(Maindata$allsumtable)                  
                        tmp <- cbind(row.names(Maindata$allsumtable),cv,Maindata$allsumtable)
                        colnames(tmp)[1:2] <- c("Feature","CV")                  
                  } else {                        
                        tmp <- cbind(row.names(Maindata$allsumtable),Maindata$allsumtable)
                        colnames(tmp)[1] <- "Feature"                  
                  }
                  write.table(tmp,file,row.names=F,quote=F,sep="\t")     
            }
      )
      
      observe({
            if (input$Sumpreviousstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 1: Data input and preprocessing")
                  })
            }
      })
      
      observe({
            if (input$Sumnextstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 3: Cell heterogeneity analysis")
                  })
            }
      })
      
      # output$Sumheatmap <- renderD3heatmap({
      #       d3heatmap(Maindata$sumtable[Maindata$sumtablenametype == input$Sumheatmapselecttype,,drop=F],dendrogram="both")
      # })
      # 
      # output$Sumheatmapselecttypeui <- renderUI({
      #       selectInput("Sumheatmapselecttype","Select Feature Type to be displayed",unique(Maindata$sumtablenametype))
      # })
      
      ### Sample-level Analysis ###
      
      #Sample clustering
      
      output$Sampselectfeattypeui <- renderUI({
            tmp <- unique(Maindata$sumtablenametype)
            if (length(grep("ENCL",tmp)) > 0) {
                  checkboxGroupInput("Sampselectfeattype","Select Feature Type",tmp,selected=tmp[grep("ENCL",tmp)])       
            } else {
                  checkboxGroupInput("Sampselectfeattype","Select Feature Type",tmp,selected=tmp)       
            }
      })
      
      output$Sampselectfeatincludebulkui <- renderUI({
            datapath <- system.file("extdata",package=paste0("SCRATdata",input$InputGenome))
            load(paste0(datapath,"/ENCODE/ENCL1000.rda"))
            selectInput("Sampselectfeatincludebulk","Select samples",colnames(ENCODEcount),multiple = T)
      })
      
      observeEvent(input$Sampselectfeatincludeallbulktf,{
            datapath <- system.file("extdata",package=paste0("SCRATdata",input$InputGenome))
            load(paste0(datapath,"/ENCODE/ENCL1000.rda"))
            updateSelectInput(session,"Sampselectfeatincludebulk","Select samples",colnames(ENCODEcount),colnames(ENCODEcount))    
      })
      
      observeEvent(input$Sampselectfeatremoveallbulktf,{
            datapath <- system.file("extdata",package=paste0("SCRATdata",input$InputGenome))
            load(paste0(datapath,"/ENCODE/ENCL1000.rda"))
            updateSelectInput(session,"Sampselectfeatincludebulk","Select samples",colnames(ENCODEcount),NULL)    
      })
      
      observe({
            if (input$MainMenu == "Step 3: Cell heterogeneity analysis" && !is.null(Maindata$allsumtable)) {
                  tmp <- Maindata$allsumtable[Maindata$sumtablenametype %in% input$Sampselectfeattype,,drop=F]
                  if (nrow(tmp)==0) {
                        tmp <- Maindata$allsumtable
                  } 
                  Maindata$sumtable <- tmp
                  if (input$Sampselectfeatincludebulktf & length(input$Sampselectfeatincludebulk) > 0) {
                        ENCODEcounttable <- NULL
                        datapath <- system.file("extdata",package=paste0("SCRATdata",input$InputGenome))
                        cortype <- input$Sampselectfeattype
                        # if ("TSS" %in% cortype) {
                        #       load(paste0(datapath,"/ENCODE/generegion.rda"))
                        #       ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),input$Sampselectfeatincludebulk,drop=F])
                        # }
                        if (sum(grepl("ENCL",cortype))==1) {
                              load(paste0(datapath,"/ENCODE/",cortype[grep("ENCL",cortype)],".rda"))
                              ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),input$Sampselectfeatincludebulk,drop=F])
                        }
                        if ("MOTIF" %in% cortype) {
                              load(paste0(datapath,"/ENCODE/MOTIF.rda"))
                              ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),input$Sampselectfeatincludebulk,drop=F])
                        }
                        if ("GSEA" %in% cortype) {
                              load(paste0(datapath,"/ENCODE/GSEA.rda"))
                              ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),input$Sampselectfeatincludebulk,drop=F])
                        }
                        if (input$Sumlogtf) {
                              ENCODEcounttable <- log2(ENCODEcounttable + 1)
                        }
                        Maindata$ENCODEcounttable <- ENCODEcounttable[row.names(tmp),,drop=F]
                        if (input$Sampcludimredscale) {
                              tmptable <- apply(Maindata$ENCODEcounttable,2,scale)
                              dimnames(tmptable) <- dimnames(Maindata$ENCODEcounttable)
                        }
                        tmp <- predict(Maindata$pcares,t(tmptable))[,1:Maindata$pcadim,drop=F]
                        colnames(tmp) <- paste0("PCA",1:ncol(tmp))
                        Maindata$ENCODEreducedata <- tmp
                  }
            }
      })
      
      observe({
            if (input$Sampcluuploadbutton > 0) {
                  isolate({
                        FileHandle <- input$SampcluInputFile                        
                        if (!is.null(FileHandle)) {
                              Maindata$uploadclulist <- read.table(FileHandle$datapath,as.is=T,sep="\t")
                        }                        
                  })
            }
      })
      
      output$Sampcluuploadstateui <- renderUI({
            if (!is.null(Maindata$uploadclulist)) {
                  if (sum(!colnames(Maindata$allsumtable) %in% Maindata$uploadclulist[,1]) == 0) {
                        helpText("Upload Successful")
                  } else {
                        helpText("Something is wrong with the uploaded list. Please check!")
                  }
            }
      })
      
      
      
      output$Sampcluplotdim2optionui <- renderUI({
            if (!is.null(input$Sampcluplotselectfeat) && (length(input$Sampcluplotselectfeat)==2)) {
                  tagList(checkboxInput("Sampcluplotshowlabtf","Show label",value=T),
                          helpText("Use computer mouse to drag and zoom the plot. Move the curser to individual points to reveal details. Move the curser to cluster id on the right to highlight points belonging to each cluster."))      
            }
      })
      
      observe({
            if (input$Sampclurunbutton > 0) {
                  isolate({
                        if (!is.null(Maindata$sumtable)) {
                              withProgress(message = 'Performing Clustering',{
                                    if (input$Sampcludimredscale) {
                                          sumtable <- apply(Maindata$sumtable,2,scale)
                                          dimnames(sumtable) <- dimnames(Maindata$sumtable)
                                    } else {
                                          sumtable <- Maindata$sumtable
                                    }
                                    if (input$Sampcludimredmet=='PCA') {                              
                                          pcares <- prcomp(t(sumtable),scale. = T)
                                          if (input$Sampcluoptdimnum) {
                                                tmp <- pcares$sdev
                                                sdev <- tmp[1:min(20,length(tmp))]
                                                x <- 1:length(sdev)
                                                pcadim <- which.min(sapply(x, function(i) {
                                                      x2 <- pmax(0,x-i)
                                                      sum(lm(sdev~x+x2)$residuals^2)
                                                }))
                                          } else {
                                                pcadim <- as.numeric(input$Sampcluchoosedimnum)
                                          }
                                          tmp <- pcares$x[,1:pcadim]
                                          colnames(tmp) <- paste0("PCA",1:ncol(tmp))
                                          Maindata$reducedata <- tmp
                                          Maindata$pcares <- pcares
                                          Maindata$pcadim <- pcadim
                                    } else if (input$Sampcludimredmet=='tSNE') {
                                          set.seed(12345)
                                          tmp <- tsne(dist(t(sumtable)),k=as.numeric(input$Sampcluchoosedimnum))
                                          row.names(tmp) <- colnames(sumtable)
                                          colnames(tmp) <- paste0("tSNE",1:ncol(tmp))
                                          Maindata$reducedata <- tmp
                                    } else {
                                          Maindata$reducedata <- t(sumtable)
                                    }
                                    if (input$Sampcluclumet=="kmeans") {
                                          if (input$Sampcluoptclunum) {
                                                x <- 2:min(20,ncol(Maindata$sumtable))
                                                reducedata <- Maindata$reducedata
                                                y <- sapply(x, function(k) {
                                                      set.seed(12345)
                                                      clu <- kmeans(reducedata,k,iter.max = 100)$cluster
                                                      cluSS <- sum(sapply(unique(clu),function(i) {
                                                            sum(rowSums(sweep(reducedata[clu==i,,drop=F],2,colMeans(reducedata[clu==i,,drop=F]),"-")^2))
                                                      }))
                                                      1-cluSS/sum((sweep(reducedata,2,colMeans(reducedata),"-"))^2)
                                                })      
                                                y <- c(0,y)
                                                x <- c(1,x)
                                                clunum <- which.min(sapply(x, function(i) {
                                                      x2 <- pmax(0,x-i)
                                                      sum(lm(y~x+x2)$residuals^2)
                                                }))
                                          } else {
                                                clunum <- as.numeric(input$Sampcluchooseclunum)
                                          }
                                          set.seed(12345)
                                          Maindata$cluster <- kmeans(Maindata$reducedata,clunum,iter.max = 100)$cluster
                                    } else if (input$Sampcluclumet=="hclust") {
                                          datahclust <- hclust(dist(Maindata$reducedata))
                                          if (input$Sampcluoptclunum) {
                                                x <- 2:min(20,ncol(Maindata$sumtable))
                                                reducedata <- Maindata$reducedata
                                                y <- sapply(x, function(k) {
                                                      clu <- cutree(datahclust,k)
                                                      cluSS <- sum(sapply(unique(clu),function(i) {
                                                            sum(rowSums(sweep(reducedata[clu==i,,drop=F],2,colMeans(reducedata[clu==i,,drop=F]),"-")^2))
                                                      }))
                                                      1-cluSS/sum((sweep(reducedata,2,colMeans(reducedata),"-"))^2)
                                                })      
                                                y <- c(0,y)
                                                x <- c(1,x)
                                                clunum <- which.min(sapply(x, function(i) {
                                                      x2 <- pmax(0,x-i)
                                                      sum(lm(y~x+x2)$residuals^2)
                                                }))
                                          } else {
                                                clunum <- as.numeric(input$Sampcluchooseclunum)
                                          }
                                          Maindata$cluster <- cutree(datahclust,clunum)
                                    } else if (input$Sampcluclumet=="mclust") {
                                          if (input$Sampcluoptclunum) {
                                                x <- 2:min(20,ncol(Maindata$sumtable)-1)
                                                clunum <- tryCatch(Mclust(Maindata$reducedata,G=x,modelNames="VVV",priorControl(functionName="defaultPrior", shrinkage=0.1))$G,warning=function(w) {}, error=function(e) {})
                                          } else {
                                                clunum <- as.numeric(input$Sampcluchooseclunum)
                                          }
                                          res <- tryCatch(Mclust(Maindata$reducedata,G=clunum,modelNames="VVV",priorControl(functionName="defaultPrior", shrinkage=0.1)),warning=function(w) {}, error=function(e) {})
                                          Maindata$cluster <- tryCatch(apply(res$z,1,which.max),warning=function(w) {}, error=function(e) {})
                                    } else {
                                          Maindata$cluster <- Maindata$uploadclulist[match(colnames(Maindata$sumtable),Maindata$uploadclulist[,1]),2]
                                    }
                              })
                        }
                  })
            }
      })
      
      observe({
            if (input$Sampclucludiagrun > 0) {
                  isolate({
                        clu <- Maindata$cluster
                        reducedata <- Maindata$reducedata
                        cluSS <- sum(sapply(unique(clu),function(i) {
                              sum(rowSums(sweep(reducedata[clu==i,,drop=F],2,colMeans(reducedata[clu==i,,drop=F]),"-")^2))
                        }))
                        Maindata$trueSSper <- 1-cluSS/sum((sweep(reducedata,2,colMeans(reducedata),"-"))^2)
                        set.seed(12345)
                        Maindata$simuSSper <- sapply (1:as.numeric(input$Sampclucludiagsimnum), function(d) {
                              clu <- sample(Maindata$cluster)
                              cluSS <- sum(sapply(unique(clu),function(i) {
                                    sum(rowSums(sweep(reducedata[clu==i,,drop=F],2,colMeans(reducedata[clu==i,,drop=F]),"-")^2))
                              }))
                              1-cluSS/sum((sweep(reducedata,2,colMeans(reducedata),"-"))^2)      
                        })
                  })
            }
      })
      
      output$Sampclucludiagplot <- renderPlot({
            if (!is.null(Maindata$trueSSper)) {
                  hist(Maindata$simuSSper,xlim=c(0,Maindata$trueSSper),main="",xlab="Percentage of between cluster variance and total variance")
                  abline(v=Maindata$trueSSper,col="red",lwd=2)
            }
      })
      
      output$SampclucludiagtrueSS <- renderText({
            if (!is.null(Maindata$trueSSper)) {
                  paste("True percentage of between cluster variance and total variance:",Maindata$trueSSper)
            }
      })
      
      output$SampclucludiagsimuSS <- renderText({
            if (!is.null(Maindata$simuSSper)) {
                  paste("Averaged simulated percentage of between cluster variance and total variance:",mean(Maindata$simuSSper))
            }
      })
      
      output$Sampcluplotselectfeatui <- renderUI({
            if (!input$Sampcluplotorifeattf) {
                  selectInput("Sampcluplotselectfeat","Select plot features",colnames(Maindata$reducedata),multiple=T,selected = colnames(Maindata$reducedata)[1:2])
            } else {
                  selectInput("Sampcluplotselectfeat","Select plot features",row.names(Maindata$allsumtable),multiple=T,selected = row.names(Maindata$allsumtable)[1:2])
            }
      })
      
      output$Sampcluplotdim1 <- renderPlot({
            if (!is.null(Maindata$cluster) & length(input$Sampcluplotselectfeat)==1) {
                  if (!input$Sampcluplotorifeattf) {
                        drawdata <- data.frame(Feature=Maindata$reducedata[,input$Sampcluplotselectfeat],Cluster=paste0("Cluster",Maindata$cluster))
                  } else {
                        drawdata <- data.frame(Feature=Maindata$allsumtable[input$Sampcluplotselectfeat,],Cluster=paste0("Cluster",Maindata$cluster))
                  }
                  ggplot(drawdata, aes(Cluster, Feature)) + geom_boxplot() + ylab(input$Sampcluplotselectfeat) +
                        theme(axis.text.x = element_text(size=17,color="black"),
                              axis.text.y = element_text(size=17,color='black'),
                              axis.title.x = element_text(size=20,vjust=-1),
                              axis.title.y = element_text(size=20,vjust=1),
                              plot.margin = unit(c(1,1,1,1), "cm"),
                              legend.text = element_text(size=15),
                              legend.title = element_text(size=15),
                              legend.position = "right",
                              legend.key.size=unit(1,"cm")
                        )
            }
      })
      
      output$Sampcluplotdim2 <- renderScatterD3({
            if (!is.null(Maindata$cluster) & length(input$Sampcluplotselectfeat)==2) {
                  if (!input$Sampcluplotorifeattf) {
                        if (input$Sampselectfeatincludebulktf && length(input$Sampselectfeatincludebulk) > 0) {
                              drawdata <- data.frame(x=c(Maindata$reducedata[,input$Sampcluplotselectfeat[1]],Maindata$ENCODEreducedata[,input$Sampcluplotselectfeat[1]]),y=c(Maindata$reducedata[,input$Sampcluplotselectfeat[2]],Maindata$ENCODEreducedata[,input$Sampcluplotselectfeat[2]]),Cluster=c(Maindata$cluster,rep("Bulk",nrow(Maindata$ENCODEreducedata))),stringsAsFactors = F)
                              lab <- c(row.names(Maindata$reducedata),row.names(Maindata$ENCODEreducedata))      
                        } else {
                              drawdata <- data.frame(x=Maindata$reducedata[,input$Sampcluplotselectfeat[1]],y=Maindata$reducedata[,input$Sampcluplotselectfeat[2]],Cluster=Maindata$cluster,stringsAsFactors = F)
                              lab <- row.names(Maindata$reducedata)
                        }
                  } else {
                        if (input$Sampselectfeatincludebulktf && length(input$Sampselectfeatincludebulk) > 0) {
                              drawdata <- data.frame(x=c(Maindata$allsumtable[input$Sampcluplotselectfeat[1],],Maindata$ENCODEcounttable[input$Sampcluplotselectfeat[1],]),y=c(Maindata$allsumtable[input$Sampcluplotselectfeat[2],],Maindata$ENCODEcounttable[input$Sampcluplotselectfeat[2],]),Cluster=c(Maindata$cluster,rep("Bulk",ncol(Maindata$ENCODEcounttable))),stringsAsFactors = F)
                              lab <- c(colnames(Maindata$allsumtable),colnames(Maindata$ENCODEcounttable))
                        } else {
                              drawdata <- data.frame(x=Maindata$allsumtable[input$Sampcluplotselectfeat[1],],y=Maindata$allsumtable[input$Sampcluplotselectfeat[2],],Cluster=Maindata$cluster,stringsAsFactors = F)
                              lab <- colnames(Maindata$allsumtable)
                        }
                  }
                  if (!input$Sampcluplotshowlabtf) {
                        lab <- NULL
                  }
                  scatterD3(x = drawdata$x, y = drawdata$y, col_lab = "Cluster", col_var = drawdata$Cluster, lab = lab, xlab = input$Sampcluplotselectfeat[1], ylab = input$Sampcluplotselectfeat[2], transitions = TRUE)      
            }
      })
      
      output$Sampcluplotdim3 <- renderPlot({
            if (!is.null(Maindata$cluster) & length(input$Sampcluplotselectfeat)>2) {
                  if (!input$Sampcluplotorifeattf) {
                        data <- Maindata$reducedata[,input$Sampcluplotselectfeat]
                  } else {
                        data <- t(Maindata$allsumtable[input$Sampcluplotselectfeat,])
                  }
                  anno = data.frame(Cluster = as.factor(Maindata$cluster))
                  rownames(anno) = row.names(data)           
                  pheatmap(data, annotation_row = anno,cluster_rows=F,cluster_cols=F)
            }
      })
      
      output$Sampcluplotdownload <- downloadHandler(
            filename = function() { 'SampleCluster.pdf' },
            content = function(file) {                  
                  if (!is.null(Maindata$cluster)) {
                        if (length(input$Sampcluplotselectfeat)==1) {
                              if (!input$Sampcluplotorifeattf) {
                                    drawdata <- data.frame(Feature=Maindata$reducedata[,input$Sampcluplotselectfeat],Cluster=paste0("Cluster",Maindata$cluster))
                              } else {
                                    drawdata <- data.frame(Feature=Maindata$allsumtable[input$Sampcluplotselectfeat,],Cluster=paste0("Cluster",Maindata$cluster))
                              }
                              pdf(file,width=12,height=10)
                              p1 <- ggplot(drawdata, aes(Cluster, Feature)) + geom_boxplot() + ylab(input$Sampcluplotselectfeat) +
                                    theme(axis.text.x = element_text(size=17,color="black"),
                                          axis.text.y = element_text(size=17,color='black'),
                                          axis.title.x = element_text(size=20,vjust=-1),
                                          axis.title.y = element_text(size=20,vjust=1),
                                          plot.margin = unit(c(1,1,1,1), "cm"),
                                          legend.text = element_text(size=15),
                                          legend.title = element_text(size=15),
                                          legend.position = "right",
                                          legend.key.size=unit(1,"cm")
                                    )
                              print(p1)
                              dev.off()      
                        } else if (length(input$Sampcluplotselectfeat)==2) {
                              if (!input$Sampcluplotorifeattf) {
                                    drawdata <- data.frame(x=Maindata$reducedata[,input$Sampcluplotselectfeat[1]],y=Maindata$reducedata[,input$Sampcluplotselectfeat[2]],Cluster=as.character(Maindata$cluster),stringsAsFactors = F)
                              } else {
                                    drawdata <- data.frame(x=Maindata$allsumtable[input$Sampcluplotselectfeat[1],],y=Maindata$allsumtable[input$Sampcluplotselectfeat[2],],Cluster=as.character(Maindata$cluster),stringsAsFactors = F)
                              }
                              pdf(file,width=12,height=10)
                              p1 <- ggplot(aes(x = x, y = y, color = Cluster),data=drawdata) + geom_point(size=3) + xlab(input$Sampcluplotselectfeat[1]) + ylab(input$Sampcluplotselectfeat[2]) + 
                                    theme(axis.text.x = element_text(size=17,color="black"),
                                          axis.text.y = element_text(size=17,color='black'),
                                          axis.title.x = element_text(size=20,vjust=-1),
                                          axis.title.y = element_text(size=20,vjust=1),
                                          plot.margin = unit(c(1,1,1,1), "cm"),
                                          legend.text = element_text(size=15),
                                          legend.title = element_text(size=15),
                                          legend.position = "right",
                                          legend.key.size=unit(1,"cm")
                                    )
                              print(p1)
                              dev.off()
                        } else {
                              if (!input$Sampcluplotorifeattf) {
                                    data <- Maindata$reducedata[,input$Sampcluplotselectfeat]
                              } else {
                                    data <- t(Maindata$allsumtable[input$Sampcluplotselectfeat,])
                              }
                              anno = data.frame(Cluster = as.factor(Maindata$cluster))
                              rownames(anno) = row.names(data)           
                              pheatmap(data, annotation_row = anno,cluster_rows=F,cluster_cols=F,filename=file,width=10,height=40)
                        }
                  }
            }
      )
      
      output$Sampclutabledownload <- downloadHandler(
            filename = function() { "SampleCluster.txt" },
            content = function(file) {       
                  tmp <- data.frame(Cell=colnames(Maindata$sumtable),Cluster=Maindata$cluster)
                  write.table(tmp,file,row.names=F,quote=F,sep="\t")     
            }
      )
      
      #Bulk correlation
      
      observe({
            if (input$Sampbulkcorrunbutton > 0) {
                  isolate({
                        withProgress(message = 'Calculating correlation',{
                              ENCODEcounttable <- NULL
                              datapath <- system.file("extdata",package=paste0("SCRATdata",input$InputGenome))
                              cortype <- input$Sampselectfeattype
                              if ("TSS" %in% cortype) {
                                    load(paste0(datapath,"/ENCODE/generegion.rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
                              #                               if ("ENLO" %in% input$Sumselectmet) {
                              #                                     load(paste0(datapath,"/ENCODE/ENLO.rda"))
                              #                                     ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              #                               }
                              if (sum(grepl("ENCL",cortype))==1) {
                                    load(paste0(datapath,"/ENCODE/",cortype[grep("ENCL",cortype)],".rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
                              if ("MOTIF" %in% cortype) {
                                    load(paste0(datapath,"/ENCODE/MOTIF.rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
                              if ("GSEA" %in% cortype) {
                                    load(paste0(datapath,"/ENCODE/GSEA.rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
                              if (input$Sumlogtf) {
                                    ENCODEcounttable <- log2(ENCODEcounttable + 1)
                              }
                              sumtable <- Maindata$sumtable[row.names(Maindata$sumtable) %in% row.names(ENCODEcounttable),]
                              if (input$Sampbulkcombinetf) {
                                    uniname <- sub("AlnRep.*$","",colnames(ENCODEcounttable))
                                    ENCODEcounttable <- sapply(unique(uniname),function(uname) {
                                          rowMeans(ENCODEcounttable[,uniname==uname,drop=F])
                                    })
                              } else {
                                    colnames(ENCODEcounttable) <- sub("Aln","",colnames(ENCODEcounttable))
                              }
                              corres <- t(sapply(1:ncol(sumtable), function(sumtableid) {
                                    sapply(1:ncol(ENCODEcounttable), function(ENCODEcounttableid) {
                                          cor(sumtable[,sumtableid],ENCODEcounttable[,ENCODEcounttableid])
                                    })
                              }))
                              row.names(corres) <- colnames(sumtable)
                              colnames(corres) <- colnames(ENCODEcounttable)
                              Maindata$bulkcorres <- corres
                        })
                  })
            }
      })
      
      output$Sampbulkcorheatmap <- renderPlot({
            if (!is.null(Maindata$bulkcorres)) {
                  tmp <- Maindata$bulkcorres
                  if (input$Sampbulkcorplotstdtf) {
                        tmp <- t(apply(tmp,1,scale))
                  }
                  dimnames(tmp) <- dimnames(Maindata$bulkcorres)
                  heatmap.2(tmp,col=colorRampPalette(c("blue", "red"))(100),trace="none")
            }
      })
      
      output$Sampbulkcorplotdownload <- downloadHandler(
            filename = function() { 'Corheatmap.pdf' },
            content = function(file) {   
                  if (!is.null(Maindata$bulkcorres)) {
                        tmp <- Maindata$bulkcorres
                        if (input$Sampbulkcorplotstdtf) {
                              tmp <- t(apply(tmp,1,scale))
                        }
                        dimnames(tmp) <- dimnames(Maindata$bulkcorres)
                        pdf(file, width=20, height=20)
                        heatmap.2(tmp,col=colorRampPalette(c("blue", "red"))(100),trace="none",margins=c(10,10),lhei=c(1,8),lwid=c(1,8))
                        dev.off()
                  }
            }
      )
      
      output$Sampbulkcortabledownload <- downloadHandler(
            filename = function() { 'Cortable.csv' },
            content = function(file) {   
                  write.csv(Maindata$bulkcorres,file=file, quote=F)
            }
      )
      
      output$Sampbulkcortable <- DT::renderDataTable({
            if (!is.null(Maindata$bulkcorres)) {
                  data <- melt(Maindata$bulkcorres)
                  colnames(data) <- c("Single-cell Samples","ENCODE Bulk Samples","Corrlations")
                  DT::datatable(data,filter="top",rownames = F, options = list(columnDefs = list(list(className="dt-body-left","targets"="_all"))))
            }
      })
      
      observe({
            if (input$Samppreviousstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 2: Feature summarization")
                  })
            }
      })
      
      observe({
            if (input$Sampnextstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 4: Differential feature analysis")
                  })
            }
      })
      
      ### Feature-level Analysis ###
      
      output$Featselectfeattypeui <- renderUI({
            checkboxGroupInput("Featselectfeattype","Select Feature Type",unique(Maindata$sumtablenametype),selected=unique(Maindata$sumtablenametype))
      })
      
      output$Featselectclusterui <- renderUI({
            if (!is.null(Maindata$cluster)) {
                  selectInput("Featselectcluster","Select clusters where ANOVA or t test will be performed (at least two should be selected)",unique(Maindata$cluster),selected = sort(unique(Maindata$cluster))[1:2],multiple = T)
            }
      })
      
      output$Featttestaltui <- renderUI({
            if ((input$Featrunallclustertf && length(unique(Maindata$cluster)) == 2) || (!input$Featrunallclustertf && length(input$Featselectcluster) == 2)) {
                  radioButtons("Featttestalt","Select alternative hypothesis for t test",c("two.sided", "less", "greater"))      
            }
      })
      
      output$Featttestalttextui <- renderUI({
            if ((input$Featrunallclustertf && length(unique(Maindata$cluster)) == 2) || (!input$Featrunallclustertf && length(input$Featselectcluster) == 2)) {
                  clu <- Maindata$cluster
                  if (!input$Featrunallclustertf) {
                        clu <- clu[clu %in% as.numeric(input$Featselectcluster)]
                  }
                  uclu <- unique(clu)
                  if (input$Featttestalt=="two.sided") {
                        p(paste0("Cluster ",uclu[1]," will be compared with cluster ",uclu[2],". The alternative hypothesis is that Cluster ",uclu[1]," is not equal to cluster ",uclu[2],"."))
                  } else {
                        p(paste0("Cluster ",uclu[1]," will be compared with cluster ",uclu[2],". The alternative hypothesis is that Cluster ",uclu[1]," is ",input$Featttestalt," than cluster ",uclu[2],"."))
                  }
                  
                  
                  
            }
      })
      
      observe({
            if (input$Featrunbutton > 0) {
                  isolate({
                        if ((input$Featrunallclustertf && length(unique(Maindata$cluster)) > 2) || (!input$Featrunallclustertf && length(input$Featselectcluster) > 2)) {
                              withProgress(message = 'Performing ANOVA',{
                                    clu <- Maindata$cluster
                                    tmp <- Maindata$allsumtable[Maindata$sumtablenametype %in% input$Featselectfeattype,,drop=F]
                                    if (nrow(tmp)==0) {
                                          tmp <- Maindata$allsumtable
                                    }
                                    if (input$Featrunallclustertf) {
                                          data <- tmp
                                    } else {
                                          data <- tmp[,clu %in% as.numeric(input$Featselectcluster)]
                                          clu <- clu[clu %in% as.numeric(input$Featselectcluster)]
                                    }
                                    totalSS <- rowSums((sweep(data,1,rowMeans(data),"-"))^2)
                                    cluSS <- sapply(unique(clu),function(i) {
                                          rowSums((sweep(data[,clu==i,drop=F],1,rowMeans(data[,clu==i,drop=F]),"-"))^2)
                                    })
                                    Fstat <- ((totalSS-rowSums(cluSS))/(length(unique(clu)) - 1 ))/(rowSums(cluSS)/(ncol(data) - length(unique(clu))))          
                                    FDR <- p.adjust(pf(Fstat,(length(unique(clu)) - 1 ),(ncol(data) - length(unique(clu))),lower.tail = F),method="fdr")
                              })  
                              Maindata$Featres <- data.frame(Feature=row.names(data),Fstatistics=Fstat,FDR=FDR,stringsAsFactors = F)
                        } else if ((input$Featrunallclustertf && length(unique(Maindata$cluster)) == 2) || (!input$Featrunallclustertf && length(input$Featselectcluster) == 2)) {
                              withProgress(message = 'Performing t-test',{
                                    clu <- Maindata$cluster
                                    tmp <- Maindata$allsumtable[Maindata$sumtablenametype %in% input$Featselectfeattype,,drop=F]
                                    if (nrow(tmp)==0) {
                                          tmp <- Maindata$allsumtable
                                    }
                                    if (input$Featrunallclustertf) {
                                          data <- tmp
                                    } else {
                                          data <- tmp[,clu %in% as.numeric(input$Featselectcluster)]
                                          clu <- clu[clu %in% as.numeric(input$Featselectcluster)]
                                    }
                                    uclu <- unique(clu)
                                    res <- t(apply(data,1,function(i) {
                                          if (length(unique(i[clu==uclu[1]]))==1 & length(unique(i[clu==uclu[2]]))==1) {
                                                c(NA,NA)
                                          } else{
                                                tmptest <- t.test(i[clu==uclu[1]],i[clu==uclu[2]],alternative = input$Featttestalt,var.equal = T)
                                                c(tmptest$statistic,tmptest$p.value)      
                                          }
                                    }))
                                    FDR <- p.adjust(res[,2],method="fdr")
                              })  
                              Maindata$Featres <- data.frame(Feature=row.names(data),tstatistics=res[,1],FDR=FDR,stringsAsFactors = F)
                        }
                  })
            }
      })
      
      output$Featrestable <- DT::renderDataTable({
            DT::datatable(Maindata$Featres,filter="top",rownames = F, options = list(columnDefs = list(list(className="dt-body-left","targets"="_all",render = JS(
                  "function(data, type, row, meta) {",
                  "return data === null ? 'NA' : data;",
                  "}")))))
      })
      
      output$Featdownloadbutton <- downloadHandler(
            filename = function() { "KeyFeatures.txt" },
            content = function(file) {                  
                  write.table(Maindata$Featres,file,row.names=F,quote=F,sep="\t")     
            }
      )
      
      
      output$FeatSumtext <- renderText({
            if (!is.null(Maindata$Featres)) {
                  paste("There are ",sum(Maindata$Featres$FDR < 0.05),"significant features, which is",round(mean(Maindata$Featres$FDR < 0.05),digits=4),"percent of all features")
            }
      })
      
      output$FeatSumplot <- renderPlot({
            if (!is.null(Maindata$Featres)) {
                  if (input$Featviewtab == "statistics") {
                        hist(Maindata$Featres[,2],main="Histogram for statistics",xlab="statistics",ylab="Frequency")
                  } else {
                        hist(Maindata$Featres[,3],main="Histogram for FDR",xlab="FDR",ylab="Frequency")
                  }
                  
            }
      })
      
      observe({
            if (input$Featpreviousstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 3: Cell heterogeneity analysis")
                  })
            }
      })
      
})



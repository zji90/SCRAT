######################################################
##                      SCRAT                       ##
##             Interactive User Interface           ##
##                     UI File                      ##
##   Author:Zhicheng Ji, Weiqiang Zhou, Hongkai Ji  ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################
library(shiny)
library(GenomicAlignments)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(scatterD3)
library(DT)
library(SCRATdata)

options(shiny.maxRequestSize=10*1024^3)

shinyServer(function(input, output,session) {
      
      Maindata <- reactiveValues()
      
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
                                    datapath <- system.file("extdata",package="SCRATdata")
                                    load(paste0(datapath,"/gr/",input$InputGenome,"/blacklist.rda"))
                                    for (i in 1:filenum) {
                                          incProgress(1/filenum,detail=paste0(round(i/filenum*100),"%"))
                                          name <- FileHandle$name[i]
                                          if (input$Inputcleanname) {
                                                name <- strsplit(name,"\\.")[[1]][1]
                                          }
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
                                                tmp <- tmp[-findOverlaps(tmp,gr)@queryHits,]                                                      
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
                        updateTabsetPanel(session,"MainMenu",selected = "Step 2: Summarizing Loci")
                  })
            }
      })
      
      ### Summarizing ###
      
      observe({
            if (input$Sumuploadreadin > 0) {
                  isolate({
                        FileHandle <- input$SumuploadFile
                        if (!is.null(FileHandle)) {
                              tmp <- read.table(FileHandle$datapath,as.is=T,sep="\t")
                              gr <- GRanges(seqnames=tmp[,1],IRanges(start=tmp[,2],end=tmp[,3]))
                              names(gr) <- paste0("UPLOAD:",tmp[,1],":",tmp[,2],"-",tmp[,3])
                              Maindata$uploadgr <- gr
                        }
                  })
            }
      })
      
      output$Sumuploaddetail <- renderText({
            if (!is.null(Maindata$uploadgr)) {
                  paste("Uploaded bed file has",length(Maindata$uploadgr),"loci.")
            }
      })
      
      observe({
            if (input$Sumrunbutton > 0)
                  isolate({
                        if (length(input$Sumselectmet) > 0) {
                              allres <- NULL
                              datapath <- system.file("extdata",package="SCRATdata")
                              if ("TSS" %in% input$Sumselectmet) {
                                    withProgress(message = 'Counting Overlaps for TSS',{
                                          load(paste0(datapath,"/gr/",input$InputGenome,"/generegion.rda"))
                                          gr <- promoters(resize(gr,1),as.numeric(input$SumTSSupregionbp),as.numeric(input$SumTSSdownregionbp))
                                          tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                          tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                          if (input$Sumlogtf) {
                                                tmp <- log2(tmp + 1)
                                          }
                                          if (input$Sumfiltertf) {
                                                tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]
                                          }
                                          allres <- rbind(allres,tmp)
                                    })
                              }
#                               if ("ENLO" %in% input$Sumselectmet) {
#                                     withProgress(message = 'Counting Overlaps for ENCODE Loci',{
#                                           load(paste0(datapath,"/gr/",input$InputGenome,"/ENLO.rda"))
#                                           if (input$SumENLOreducetf) {
#                                                 gr <- reduce(gr)
#                                                 names(gr) <- paste0("ENLO:",seqnames(gr),":",start(gr),"-",end(gr))
#                                           }                                                
#                                           splitgrnum <- round(length(Maindata$bamfile)/10)
#                                           splitgrlength <- round(length(gr)/splitgrnum)
#                                           for (splitid in 1:splitgrnum) {
#                                                 if (splitid == splitgrnum) {
#                                                       splitgr <- gr[((splitid-1)*splitgrlength+1):length(gr)]
#                                                 } else {
#                                                       splitgr <- gr[((splitid-1)*splitgrlength+1):(splitid*splitgrlength)]
#                                                 }
#                                                 tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(splitgr,i))
#                                                 tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
#                                                 if (input$Sumlogtf) {
#                                                       tmp <- log2(tmp + 1)
#                                                 }
#                                                 if (input$Sumfiltertf) {
#                                                       tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
#                                                 }                        
#                                                 allres <- rbind(allres,tmp)      
#                                           }                                                
#                                     })
#                               }
                              if ("ENCL" %in% input$Sumselectmet) {
                                    withProgress(message = 'Counting Overlaps for ENCODE Cluster',{
                                          load(paste0(datapath,"/gr/",input$InputGenome,"/ENCL",input$SumENCLclunum,".rda"))
                                          tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                          tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                          if (input$Sumlogtf) {
                                                tmp <- log2(tmp + 1)
                                          }
                                          if (input$Sumfiltertf) {
                                                tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                          }                        
                                          allres <- rbind(allres,tmp)
                                    })
                              }
                              if ("MOTIF" %in% input$Sumselectmet) {
                                    if ("TRANSFAC" %in% input$SumMOTIFselect) {
                                          withProgress(message = 'Counting Overlaps for TRANSFAC',{
                                                load(paste0(datapath,"/gr/",input$InputGenome,"/transfac1.rda"))
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                }                        
                                                allres <- rbind(allres,tmp)
                                                load(paste0(datapath,"/gr/",input$InputGenome,"/transfac2.rda"))
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                }                        
                                                allres <- rbind(allres,tmp)
                                                if (input$InputGenome %in% c("hg19","hg38")) {
                                                      load(paste0(datapath,"/gr/",input$InputGenome,"/transfac3.rda"))
                                                      tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                      tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                                      if (input$Sumlogtf) {
                                                            tmp <- log2(tmp + 1)
                                                      }
                                                      if (input$Sumfiltertf) {
                                                            tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                      }                        
                                                      allres <- rbind(allres,tmp)
                                                }
                                          })
                                    }
                                    if ("JASPAR" %in% input$SumMOTIFselect) {
                                          withProgress(message = 'Counting Overlaps for JASPAR',{
                                                load(paste0(datapath,"/gr/",input$InputGenome,"/jaspar1.rda"))
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                }                        
                                                allres <- rbind(allres,tmp)
                                                load(paste0(datapath,"/gr/",input$InputGenome,"/jaspar2.rda"))
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                }                        
                                                allres <- rbind(allres,tmp)
                                          })
                                    }
                              }
                              if ("GSEA" %in% input$Sumselectmet) {
                                    withProgress(message = 'Counting Overlaps for GSEA',{
                                          for (i in input$SumGSEAselect) {
                                                load(paste0(datapath,"/gr/",input$InputGenome,"/GSEA",i,".rda"))
                                                gr <- promoters(resize(gr,1),as.numeric(input$SumGSEAupregionbp),as.numeric(input$SumGSEAdownregionbp))             
                                                tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(gr,i))
                                                tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                                if (input$Sumlogtf) {
                                                      tmp <- log2(tmp + 1)
                                                }
                                                if (input$Sumfiltertf) {
                                                      tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                                }                        
                                                allres <- rbind(allres,tmp)
                                          }
                                    })
                              }
                              if ("Upload" %in% input$Sumselectmet & !is.null(Maindata$uploadgr)) {
                                    withProgress(message = 'Counting Overlaps for Upload bed',{
                                          tmp <- sapply(Maindata$bamfile,function(i) countOverlaps(Maindata$uploadgr,i))
                                          tmp <- sweep(tmp,2,Maindata$bamsummary[,2],"/") * median(Maindata$bamsummary[,2])
                                          if (input$Sumlogtf) {
                                                tmp <- log2(tmp + 1)
                                          }
                                          if (input$Sumfiltertf) {
                                                tmp <- tmp[rowMeans(tmp < as.numeric(input$Sumfilterreads)) < as.numeric(input$Sumfilterpercen)/100,,drop=F]       
                                          }                        
                                          allres <- rbind(allres,tmp)
                                    })
                              }
                              Maindata$sumtable <- allres
                              Maindata$sumtablenametype <- sapply(row.names(Maindata$sumtable),function(i) strsplit(i,":")[[1]][1])
                              Maindata$allpcares <- prcomp(t(Maindata$sumtable),center = T, scale = T)$x
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
            if (!is.null(Maindata$sumtable)) {
                  if (input$Sumaddcv) {
                        cv <- apply(Maindata$sumtable,1,sd)/rowMeans(Maindata$sumtable)                  
                        tmp <- cbind(row.names(Maindata$sumtable),cv,Maindata$sumtable)
                        colnames(tmp)[1:2] <- c("Feature","CV")
                  } else {                        
                        tmp <- cbind(row.names(Maindata$sumtable),Maindata$sumtable)
                        colnames(tmp)[1] <- "Feature"
                  }     
                  DT::datatable(tmp,filter="top",rownames = F, options = list(columnDefs = list(list(className="dt-body-left","targets"="_all"))))
            }            
      })
      
      output$Sumtabledownload <- downloadHandler(
            filename = function() { "Summarization.txt" },
            content = function(file) {       
                  if (input$Sumaddcv) {
                        cv <- apply(Maindata$sumtable,1,sd)/rowMeans(Maindata$sumtable)                  
                        tmp <- cbind(row.names(Maindata$sumtable),cv,Maindata$sumtable)
                        colnames(tmp)[1:2] <- c("Feature","CV")                  
                  } else {                        
                        tmp <- cbind(row.names(Maindata$sumtable),Maindata$sumtable)
                        colnames(tmp)[1] <- "Feature"                  
                  }
                  write.table(tmp,file,row.names=F,quote=F,sep="\t")     
            }
      )
      
      observe({
            if (input$Sumpreviousstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 1: Input Bam Files")
                  })
            }
      })
      
      observe({
            if (input$Sumnextstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 3: Sample-level Analysis")
                  })
            }
      })
      
      ### Sample-level Analysis ###
      
      output$Sampclupcanumui <- renderUI({
            if (!is.null(Maindata$sumtable))
                  selectInput("Sampclupcanum","Select number of PCs",2:min(nrow(Maindata$sumtable),ncol(Maindata$sumtable)),2)
      })
      
      observe({
            if (input$Sampclupcaoptbutton > 0) {
                  isolate({
                        if (!is.null(Maindata$sumtable)) {
                              x <- 1:min(ncol(Maindata$sumtable),nrow(Maindata$sumtable),20)
                              sdev <- prcomp(t(Maindata$sumtable),center = T, scale = T)$sdev[x]
                              optpoint <- which.min(sapply(x, function(i) {
                                    x2 <- pmax(0,x-i)
                                    sum(lm(sdev~x+x2)$residuals^2)
                              }))
                              updateSelectInput(session,"Sampclupcanum","Select number of PCs",2:min(nrow(Maindata$sumtable),ncol(Maindata$sumtable)),optpoint)
                        }            
                  })
            }
      })
      
      observe({
            if (input$Sampclurunbutton > 0) {
                  isolate({
                        if (!is.null(Maindata$sumtable)) {
                              withProgress(message = 'Performing Clustering',{
                                    if (input$Sampclumet=='PCA') {                              
                                          Maindata$pcares <- Maindata$allpcares[,1:as.numeric(input$Sampclupcanum)]
                                          Maindata$samphclust <- hclust(dist(Maindata$pcares))
                                          data <- t(Maindata$pcares)
                                    } else {
                                          if (!input$SampcluspecificFeattf) {                                                
                                                group <- input$SampcluFeatselectgroup
                                                if ("ENCL" %in% group) 
                                                      group[group=="ENCL"] <- paste0("ENCL",input$SumENCLclunum)
                                                data <- Maindata$sumtable[Maindata$sumtablenametype %in% group,]
                                          } else if (input$SampcluspecificFeattf & length(input$SampcluFeatselectfeat) > 0) {
                                                data <- Maindata$sumtable[input$SampcluFeatselectfeat,,drop=F]
                                          }
                                          Maindata$samphclust <- hclust(dist(t(data)))
                                    }
                                    x <- 2:min(20,ncol(Maindata$sumtable)-1)
                                    y <- sapply(x, function(k) {
                                          clu <- cutree(Maindata$samphclust, k)
                                          cluSS <- sum(sapply(unique(clu),function(i) {
                                                rowSums((sweep(data[,clu==i,drop=F],1,rowMeans(data[,clu==i,drop=F]),"-"))^2)
                                          }))
                                          1-cluSS/sum((sweep(data,1,rowMeans(data),"-"))^2)            
                                    })
                                    if (input$Sampclusimutf) {
                                          Maindata$hclustvarpropsimu <- rowMeans(sapply(1:1000,function(d) {
                                                simudata <- t(apply(data,1,sample))
                                                clustres <- hclust(dist(t(simudata)))
                                                sapply(x, function(k) {
                                                      clu <- cutree(clustres, k)
                                                      cluSS <- sum(sapply(unique(clu),function(i) {
                                                            rowSums((sweep(simudata[,clu==i,drop=F],1,rowMeans(simudata[,clu==i,drop=F]),"-"))^2)
                                                      }))
                                                      1-cluSS/sum((sweep(simudata,1,rowMeans(simudata),"-"))^2)            
                                                })
                                          }))      
                                    } else {
                                          Maindata$hclustvarpropsimu <- NULL
                                    }
                                    Maindata$hclustvarprop <- data.frame(clunum=x,varprop=y)
                                    Maindata$optclustnum <- x[which.min(sapply(x, function(i) {                                          
                                          x2 <- pmax(0,x-i)
                                          sum(lm(y~x+x2)$residuals^2)
                                    }))]
                              })
                        }
                  })
            }
      })
      
      output$Sampoptcluplot <- renderPlot({
            if (!is.null(Maindata$hclustvarpropsimu)) {
                  plot(Maindata$hclustvarprop$clunum,Maindata$hclustvarprop$varprop,xlab="Cluster Number",ylab="Proportion of variance explained",pch=19,cex=1.2,cex.lab=1.2,ylim=c(0.95*min(Maindata$hclustvarprop$varprop,Maindata$hclustvarpropsimu),1.05*max(Maindata$hclustvarprop$varprop,Maindata$hclustvarpropsimu)))
                  points(Maindata$optclustnum,Maindata$hclustvarprop$varprop[which(Maindata$hclustvarprop$clunum==Maindata$optclustnum)],col="red",pch=19,cex=1.5)
                  points(Maindata$hclustvarprop$clunum,Maindata$hclustvarpropsimu,col="blue",cex=1.2,pch=19)
                  text(Maindata$optclustnum,Maindata$hclustvarprop$varprop[which(Maindata$hclustvarprop$clunum==Maindata$optclustnum)] - 0.1 * (max(Maindata$hclustvarprop$varprop)-min(Maindata$hclustvarprop$varprop)),"Optimal Cluster Number",col="red",cex=1.2)
                  legend("topleft",c("True data","Simulated data"),col=c("black","blue"),pch=c(19,19))      
            } else {
                  plot(Maindata$hclustvarprop$clunum,Maindata$hclustvarprop$varprop,xlab="Cluster Number",ylab="Proportion of variance explained",pch=19,cex=1.2,cex.lab=1.2,ylim=c(0.95*min(Maindata$hclustvarprop$varprop),1.05*max(Maindata$hclustvarprop$varprop)))
                  points(Maindata$optclustnum,Maindata$hclustvarprop$varprop[which(Maindata$hclustvarprop$clunum==Maindata$optclustnum)],col="red",pch=19,cex=1.5)
                  text(Maindata$optclustnum,Maindata$hclustvarprop$varprop[which(Maindata$hclustvarprop$clunum==Maindata$optclustnum)] - 0.1 * (max(Maindata$hclustvarprop$varprop)-min(Maindata$hclustvarprop$varprop)),"Optimal Cluster Number",col="red",cex=1.2)
            }
      })
      
      output$SampcluFeatselectfeatui <- renderUI({
            if (input$SampcluspecificFeattf) {
                  selectInput("SampcluFeatselectfeat","Select features (Initialization may take time)",row.names(Maindata$sumtable),multiple = T)      
            } else {
                  checkboxGroupInput("SampcluFeatselectgroup","Select types of features",input$Sumselectmet,selected=input$Sumselectmet)
            }           
      })
      
      output$Sampcludendrogram <- renderPlot({
            if (!is.null(Maindata$samphclust)) {
                  plot(Maindata$samphclust)
            }                       
      })
      
      output$Sampvisplotui <- renderUI({
            if (input$Sampvismet=="PCA") {
                  tagList(
                        helpText("Use computer mouse to drag and zoom the plot. Move the curser to individual points to reveal details"),
                        scatterD3::scatterD3Output("SampvisPCAplot",width="600px",height="500px")      
                  )                  
            } else if (input$Sampvismet=="MDS") {
                  tagList(
                        helpText("Use computer mouse to drag and zoom the plot. Move the curser to individual points to reveal details"),
                        scatterD3::scatterD3Output("SampvisMDSplot",width="600px",height="500px")      
                  )  
            } else if (input$Sampvismet=="Features") {
                  if (length(input$SampvisFeatselectfeat) == 1) {
                        plotOutput("SampvisFeatplotdim1",width="600px",height="500px")
                  } else if (length(input$SampvisFeatselectfeat) == 2) {
                        scatterD3::scatterD3Output("SampvisFeatplotdim2",width="600px",height="500px")      
                  } else if (length(input$SampvisFeatselectfeat) > 2) {
                        plotOutput("SampvisFeatplotdim3",width="600px",height="500px")
                  }
            }            
      })
      
      output$SampvisPCAplot <- renderScatterD3({
            drawdata <- data.frame(x=Maindata$allpcares[,as.numeric(input$SampvisPCAxpcid)],y=Maindata$allpcares[,as.numeric(input$SampvisPCAypcid)],Cluster=as.character(cutree(Maindata$samphclust,k=as.numeric(input$Sampvisclusternum))),stringsAsFactors = F)
            scatterD3(x = drawdata$x, y = drawdata$y, col_var = drawdata$Cluster, lab = row.names(Maindata$allpcares), xlab = paste0("PCA",as.numeric(input$SampPCAxpcid)), ylab = paste0("PCA",as.numeric(input$SampPCAypcid)), col_lab = "Cluster", transitions = TRUE)      
      },env=environment())
      
      output$SampvisMDSplot <- renderScatterD3({
            data <- t(Maindata$sumtable)
            d <- dist(data)
            fit <- cmdscale(d,eig=TRUE, k=min(20,ncol(Maindata$sumtable)-1))
            x <- fit$points[,as.numeric(input$SampvisMDSxdimid)]
            y <- fit$points[,as.numeric(input$SampvisMDSydimid)]
            cluster <- as.character(cutree(Maindata$samphclust,k=as.numeric(input$Sampvisclusternum)))
            scatterD3(x = x, y = y, col_var = cluster, lab = row.names(data), xlab = paste0("Dimension",as.numeric(input$SampvisMDSxdimid)), ylab = paste0("Dimension",as.numeric(input$SampvisMDSydimid)), col_lab = "Cluster", transitions = TRUE)      
      },env=environment())
      
      output$SampvisFeatplotdim2 <- renderScatterD3({
            if (length(input$SampvisFeatselectfeat) == 2) {
                  x <- Maindata$sumtable[input$SampvisFeatselectfeat[1],]
                  y <- Maindata$sumtable[input$SampvisFeatselectfeat[2],]
                  cluster <- as.character(cutree(Maindata$samphclust,k=as.numeric(input$Sampvisclusternum)))
                  scatterD3(x = x, y = y, col_var = cluster, lab = row.names(data), xlab = input$SampvisFeatselectfeat[1], ylab = input$SampvisFeatselectfeat[2], col_lab = "Cluster", transitions = TRUE)
            }
      },env=environment())
      
      output$SampvisFeatplotdim1 <- renderPlot({
            if (length(input$SampvisFeatselectfeat) == 1) {
                  data <- data.frame(Feature=Maindata$sumtable[input$SampvisFeatselectfeat,],Cluster=paste0("Cluster",cutree(Maindata$samphclust,k=as.numeric(input$Sampvisclusternum))))             
                  ggplot(data, aes(Cluster, Feature)) + geom_boxplot()      
            }            
      })
      
      output$SampvisFeatplotdim3 <- renderPlot({
            if (length(input$SampvisFeatselectfeat) > 2) {
                  data <- Maindata$sumtable[input$SampvisFeatselectfeat,]
                  Cluster=cutree(Maindata$samphclust,k=as.numeric(input$Sampvisclusternum))
                  annotation_col = data.frame(Cluster = as.factor(Cluster))
                  rownames(annotation_col) = colnames(data)           
                  pheatmap(data[,Maindata$samphclust$order], annotation_col = annotation_col[Maindata$samphclust$order,,drop=F],cluster_rows=F,cluster_cols=F)
            }
      })
      
      output$Sampvisoptionui1 <- renderUI({
            sliderInput("Sampvisclusternum","Select number of clusters",2,ncol(Maindata$sumtable)-1,2,step = 1)
      })
      
      observe({
            if (input$Sampvisclunumoptbutton > 0) {
                  isolate({
                        updateSliderInput(session,"Sampvisclusternum","Select number of clusters",Maindata$optclustnum,2,ncol(Maindata$sumtable)-1,step = 1)
                  })
            }
      })
      
      output$Sampvisoptionui2 <- renderUI({
            if (input$Sampvismet=="PCA") {
                  tagList(
                        selectInput("SampvisPCAxpcid","PC shown on x-axis",1:min(nrow(Maindata$sumtable),ncol(Maindata$sumtable),20),selected = 1),
                        selectInput("SampvisPCAypcid","PC shown on y-axis",1:min(nrow(Maindata$sumtable),ncol(Maindata$sumtable),20),selected = 2)
                  )
            } else if (input$Sampvismet=="MDS") {
                  tagList(
                        selectInput("SampvisMDSxdimid","Dimension shown on x-axis",1:min(20,ncol(Maindata$sumtable)-1),selected = 1),
                        selectInput("SampvisMDSydimid","Dimension shown on y-axis",1:min(20,ncol(Maindata$sumtable)-1),selected = 2)
                  )
            } else if (input$Sampvismet=="Features") {
                  selectInput("SampvisFeatselectfeat","Select features (Initialization may take time)",row.names(Maindata$sumtable),multiple = T)  
            }
      })
      
      observe({
            if (input$Sampbulkrunbutton > 0) {
                  isolate({
                        withProgress(message = 'Calculating correlation',{
                              ENCODEcounttable <- NULL
                              datapath <- system.file("extdata",package="SCRATdata")
                              if ("TSS" %in% input$Sumselectmet) {
                                    load(paste0(datapath,"/ENCODE/",input$InputGenome,"/generegion.rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
#                               if ("ENLO" %in% input$Sumselectmet) {
#                                     load(paste0(datapath,"/ENCODE/",input$InputGenome,"/ENLO.rda"))
#                                     ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
#                               }
                              if ("ENCL" %in% input$Sumselectmet) {
                                    load(paste0(datapath,"/ENCODE/",input$InputGenome,"/ENCL",input$SumENCLclunum,".rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
                              if ("MOTIF" %in% input$Sumselectmet) {
                                    load(paste0(datapath,"/ENCODE/",input$InputGenome,"/MOTIF.rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
                              if ("GSEA" %in% input$Sumselectmet) {
                                    load(paste0(datapath,"/ENCODE/",input$InputGenome,"/GSEA.rda"))
                                    ENCODEcounttable <- rbind(ENCODEcounttable,ENCODEcount[row.names(ENCODEcount) %in% row.names(Maindata$sumtable),])
                              }
                              if (input$Sumlogtf) {
                                    ENCODEcounttable <- log2(ENCODEcounttable + 1)
                              }
                              sumtable <- Maindata$sumtable[row.names(Maindata$sumtable) %in% row.names(ENCODEcounttable),]
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
                  pheatmap(Maindata$bulkcorres,cluster_rows = F,cluster_cols = F)
            }
      })
      
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
                        updateTabsetPanel(session,"MainMenu",selected = "Step 2: Summarizing Loci")
                  })
            }
      })
      
      observe({
            if (input$Sampnextstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 4: Feature-level Analysis")
                  })
            }
      })
      
      ### Feature-level Analysis ###
      
      observe({
            if (input$Featrunbutton > 0) {
                  isolate({
                        withProgress(message = 'Identifying feature statistics',{
                              data <- Maindata$sumtable
                              cluidlist <- as.numeric(input$Keyclunum[1]):as.numeric(input$Keyclunum[2])                              
                              FDR <- Fstat <- matrix(0,nrow=nrow(data),ncol=length(cluidlist))                              
                              for (cluid in 1:length(cluidlist)) {
                                    clu <- cutree(Maindata$samphclust,k=cluidlist[cluid])
                                    totalSS <- rowSums((sweep(data,1,rowMeans(data),"-"))^2)
                                    cluSS <- sapply(unique(clu),function(i) {
                                          rowSums((sweep(data[,clu==i,drop=F],1,rowMeans(data[,clu==i,drop=F]),"-"))^2)
                                    })
                                    Fstat[,cluid] <- ((totalSS-rowSums(cluSS))/(length(unique(clu)) - 1 ))/(rowSums(cluSS)/(ncol(data) - length(unique(clu))))                                    
                                    FDR[,cluid] <- p.adjust(pf(Fstat[,cluid],(length(unique(clu)) - 1 ),(ncol(data) - length(unique(clu))),lower.tail = F),method="fdr")
                              }
                        })  
                        row.names(Fstat) <- row.names(data)
                        row.names(FDR) <- row.names(data)
                        colnames(Fstat) <- paste0("Cluster_",cluidlist)
                        colnames(FDR) <- paste0("Cluster_",cluidlist)
                        Maindata$FeatFstat <- Fstat
                        Maindata$FeatFDR <- FDR
                  })
            }
      })
      
      output$Featsampclunumui <- renderUI({
            if (!is.null(Maindata$sumtable))
                  sliderInput("Keyclunum","",2,ncol(Maindata$sumtable)-1,c(2,ncol(Maindata$sumtable)-1),step=1)
      })
      
      output$Featrestable <- DT::renderDataTable({
            if (input$Featviewtab == "Fstat") {
                  tmp <- data.frame(Maindata$FeatFstat)      
            } else {
                  tmp <- data.frame(Maindata$FeatFDR)
            }
            tmp <- cbind(row.names(tmp),tmp)
            colnames(tmp)[1] <- "Feature"
            tmp[,1] <- as.character(tmp[,1])
            DT::datatable(tmp,filter="top",rownames = F, options = list(columnDefs = list(list(className="dt-body-left","targets"="_all"))))
      })
      
      output$Featdownloadbutton <- downloadHandler(
            filename = function() { "KeyFeatures.txt" },
            content = function(file) {                  
                  if (input$Featviewtab == "Fstat") {
                        tmp <- data.frame(Maindata$FeatFstat)      
                  } else {
                        tmp <- data.frame(Maindata$FeatFDR)
                  }
                  tmp <- cbind(row.names(tmp),tmp)
                  colnames(tmp)[1] <- "Feature"
                  write.table(tmp,file,row.names=F,quote=F,sep="\t")     
            }
      )
      
      observe({
            if (input$Featpreviousstepbutton) {
                  isolate({
                        updateTabsetPanel(session,"MainMenu",selected = "Step 3: Sample-level Analysis")
                  })
            }
      })
      
})



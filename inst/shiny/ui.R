######################################################
##                      SCRAT                       ##
##             Interactive User Interface           ##
##                     UI File                      ##
##   Author:Zhicheng Ji, Weiqiang Zhou, Hongkai Ji  ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################

shinyUI(      
      navbarPage("SCRAT",
                 
                 tabPanel("Step 1: Input Bam Files",
                          tags$head(
                                tags$style(HTML('#Inputnextstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Inputreadin{font-weight: bold;color:blue}'))
                          ),
                          hr(),
                          hr(),
                          sidebarPanel(
                                fluidRow(actionButton("Inputnextstepbutton","Next Step"),align="center"),
                                helpText("To start the analysis, users should have the already aligned bam files. Select the corresponding genome used for alignment and upload the bam files."),
                                br(),
                                wellPanel(
                                      h4("Select Genome"),
                                      uiOutput("InputGenomeui")
                                ),
                                wellPanel(
                                      h5("To upload a summary table from previous SCRAT session, skip this step and go directly to step 2."),
                                      h4("Input Bam Files"),
                                      fileInput('InputFile', 'Choose File', multiple = T, accept = ".bam"),
                                      checkboxInput("Inputblacklist","Filter black list",value = T),
                                      checkboxInput("Inputcleanname","Clean up sample name",value = F),
                                      p(actionButton("Inputreadin","Read in"))      
                                ),
                                wellPanel(
                                      h4("Filter Bam Files"),
                                      textInput("InputFilterreads","Exclude samples with reads less than","1000"),
                                      uiOutput("InputFiltersampleui")
                                ),
                                width=3),
                          
                          mainPanel(
                                uiOutput("Inputbamnumui"),
                                DT:: dataTableOutput("Inputbamsummary")
                          )),
                 tabPanel("Step 2: Summarizing Loci",
                          tags$head(
                                tags$style(HTML('#Sumpreviousstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sumnextstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sumrunbutton{font-weight: bold;color:blue}'))
                          ),
                          hr(),
                          hr(),
                          sidebarPanel(
                                fluidRow(actionButton("Sumpreviousstepbutton","Previous Step"),
                                         actionButton("Sumnextstepbutton","Next Step"),align="center"),
                                helpText("For each sample, count the number of reads overlapping the genomic loci of each feature."),
                                br(),
                                checkboxInput("Sumuploadsumtable","Upload Summarization Table",value=F),
                                conditionalPanel(condition="input.Sumuploadsumtable==1",
                                                 wellPanel(
                                                       h5("The table should be exactly the same saved from previous SCRAT session"),
                                                       h4("Input Summary Table"),
                                                       fileInput('SumuploadsumtableFile', 'Choose File', accept = ".txt"),                                                       
                                                       p(actionButton("Sumuploadsumtablereadin","Read in"))      
                                                 )
                                ),
                                conditionalPanel(condition="input.Sumuploadsumtable==0",
                                                 wellPanel(
                                                       checkboxGroupInput("Sumselectmet","Choose Summarizing Method",list("TSS Upstream"="TSS","ENCODE Cluster"="ENCL","Motif Sites"="MOTIF","GSEA Gene Sets"="GSEA","Upload BED"="Upload"),selected="TSS"),
                                                       hr(),
                                                       checkboxInput("Sumlogtf","Log2 transformation",value=T),
                                                       checkboxInput("Sumaddcv","Add coefficient of variation information",value=T),
                                                       checkboxInput("Sumfiltertf","Filter Features",value=T),
                                                       conditionalPanel(condition="input.Sumfiltertf==1",
                                                                        textInput("Sumfilterpercen","Exclude features having more than","90"),
                                                                        textInput("Sumfilterreads","percent of samples whose (normalized) reads are less than","0.01")
                                                       ),
                                                       hr(),
                                                       actionButton("Sumrunbutton","Run Summarization")
                                                 ),                        
                                                 wellPanel(
                                                       selectInput("Sumdetailchoose","Method Details",list("TSS Up/Downstream"="TSS","ENCODE Cluster"="ENCL","Motif Sites"="MOTIF","GSEA Gene Sets"="GSEA","Upload BED"="Upload")),
                                                       conditionalPanel(condition="input.Sumdetailchoose=='TSS'",
                                                                        helpText('For each gene, sum all reads overlapping TSS upstream and downstream region. For gene with ENTREZ id of 18999, the feature name will be GENE:18999'),
                                                                        textInput("SumTSSupregionbp","TSS upstream base pair","1000"),
                                                                        textInput("SumTSSdownregionbp","TSS downstream base pair","500")
                                                       ),
                                                       #                                       conditionalPanel(condition="input.Sumdetailchoose=='ENLO'",
                                                       #                                                        helpText('200 bp windows of genomic loci that are potential regulatory elements were precompiled based on ENCODE DNase-seq data.'),
                                                       #                                                        checkboxInput("SumENLOreducetf","Merge adjacent windows",value=F)
                                                       #                                       ),
                                                       conditionalPanel(condition="input.Sumdetailchoose=='ENCL'",
                                                                        helpText('Clusters of genomic regions (1000,2000 or 5000 clusters) were precompiled based on ENCODE DNase-seq data. For each cluster, sum all reads overlapping any of its genomic regions. For cluster id 1, the feature name will be ENCL1000:Cluster1'),
                                                                        radioButtons("SumENCLclunum","Choose number of clusters",c("1000","2000","5000"))
                                                       ),
                                                       conditionalPanel(condition="input.Sumdetailchoose=='MOTIF'",
                                                                        helpText('Motif sites were mapped to whole genome and filtered for potential TFBS using ENCODE DNase-seq data. For each motif, sum all reads overlapping any of its binding site. For Motif MA0040.1_Foxq1, the feature name will be MOTIF:MA0040.1_Foxq1'),
                                                                        textInput("SumMOTIFflank","Take flank region of base pair","100"),
                                                                        checkboxGroupInput("SumMOTIFselect","Include motifs sites from",c("TRANSFAC","JASPAR"),selected = c("TRANSFAC","JASPAR")),
                                                                        helpText(a("TRANSFAC web site",href="http://www.gene-regulation.com/pub/databases.html",target="_blank")),
                                                                        helpText(a("JASPAR web site",href="http://jaspar.genereg.net/",target="_blank"))
                                                       ),
                                                       conditionalPanel(condition="input.Sumdetailchoose=='GSEA'",
                                                                        helpText('For each GSEA gene set, sum all reads overlapping any gene TSS upstream region from the gene set. For GSEA gene set HALLMARK_HYPOXIA, the feature name will be GSEA:HALLMARK_HYPOXIA'),
                                                                        textInput("SumGSEAupregionbp","TSS upstream base pair","1000"),
                                                                        textInput("SumGSEAdownregionbp","TSS downstream base pair","500"),
                                                                        checkboxGroupInput("SumGSEAselect","Include gene sets from",list("Hallmark gene sets (50 gene sets)"="h.all","Positional gene sets (326 gene sets)"="c1.all","Curated gene sets: chemical and genetic perturbations (3395 gene sets)"="c2.cgp","Curated gene sets: canonical pathways (1330 gene sets)"="c2.cp","Motif gene sets: microRNA targets (221 gene sets)"="c3.mir","Motif gene sets: transcription factor targets (615 gene sets)"="c3.tft","Computational gene sets: cancer gene neighborhoods (427 gene sets)"="c4.cgn","Computational gene sets: cancer modules (431 gene sets)"="c4.cm","GO gene sets: biological process (825 gene sets)"="c5.bp","GO gene sets: cellular component (233 gene sets)"="c5.cc","GO gene sets: molecular function (396 gene sets)"="c5.mf","Oncogenic signatures (189 gene sets)"="c6.all","Immunologic signatures (1910 gene sets)"="c7.all"),selected = c("h.all","c1.all","c2.cgp","c2.cp","c3.mir","c3.tft","c4.cgn","c4.cm","c5.bp","c5.cc","c5.mf","c6.all","c7.all")),
                                                                        helpText(a("Browse GSEA gene sets",href="http://software.broadinstitute.org/gsea/msigdb/genesets.jsp",target="_blank"))
                                                       ),
                                                       conditionalPanel(condition="input.Sumdetailchoose=='Upload'",
                                                                        helpText('Upload BED file to define genomic regions'),
                                                                        fileInput('SumuploadFile', 'Choose File', accept = ".bed"),
                                                                        p(actionButton("Sumuploadreadin","Read in")),
                                                                        textOutput("Sumuploaddetail")
                                                       )
                                                 )),
                                width=3),
                          
                          mainPanel(
                                br(),
                                tabsetPanel(
                                      tabPanel("Results",                                               
                                               br(),
                                               p(downloadButton("Sumtabledownload")),
                                               DT:: dataTableOutput("Sumtableoutput")               
                                      ),
                                      tabPanel("Type Summary",
                                               br(),
                                               DT:: dataTableOutput("Sumtypesummarytable")
                                      )
                                      #                                       tabPanel("Interactive Heatmap",
                                      #                                                br(),
                                      #                                                uiOutput("Sumheatmapselecttypeui"),
                                      #                                                d3heatmap::d3heatmapOutput("Sumheatmap"))
                                )                                
                          )
                 ),
                 tabPanel("Step 3: Sample-level Analysis",
                          tags$head(
                                tags$style(HTML('#Samppreviousstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sampnextstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sampclurunbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sampbulkrunbutton{font-weight: bold;color:blue}'))
                          ),
                          hr(),
                          hr(),
                          sidebarPanel(
                                fluidRow(actionButton("Samppreviousstepbutton","Previous Step"),
                                         actionButton("Sampnextstepbutton","Next Step"),align="center"),
                                helpText("Samples will be clustered using hierarchical clustering using the features obtained in step 2. The clustering can be based on the PCA of the features or the original features."),
                                br(),
                                wellPanel(uiOutput("Sampselectfeattypeui"),
                                          helpText("Select feature type to be included in the sample-level analysis. After changing selected feature type, please rerun all analysis. If no feature type is selected, all feature types will be used in the analsysis.")),
                                wellPanel(
                                      radioButtons("Sampmainmet","",list("Sample Clustering"="Clustering","Visualization"="Visualization","Compare with Bulk"="Bulk"))
                                ),
                                wellPanel(
                                      conditionalPanel(condition="input.Sampmainmet=='Clustering'",
                                                       radioButtons("Sampclumet","Select clustering method",list("PCA"="PCA","Features"="Features")),
                                                       conditionalPanel(condition="input.Sampclumet=='PCA'",    
                                                                        uiOutput("Sampclupcanumui"),                                                                        
                                                                        actionButton("Sampclupcaoptbutton","Use optimal number of PCs.")
                                                       ),                                
                                                       hr(),
                                                       checkboxInput("Sampclusimutf","Perform simulation test"),
                                                       actionButton("Sampclurunbutton","Perform Clustering"),
                                                       hr()
                                                       
                                      ),
                                      conditionalPanel(condition="input.Sampmainmet=='Visualization'",
                                                       uiOutput("Sampvisoptionui1"),
                                                       actionButton("Sampvisclunumoptbutton","Use optimal number of clusters"),
                                                       hr(),
                                                       radioButtons("Sampvismet","Select visualization method",c("PCA","MDS","Features")),
                                                       uiOutput("Sampvisoptionui2")
                                      ),
                                      conditionalPanel(condition="input.Sampmainmet=='Bulk'",
                                                       checkboxInput("Sampbulkcombinetf","Combine replicates",value=F),
                                                       actionButton("Sampbulkrunbutton","Run comparison")              
                                      )
                                )
                                ,width=3),
                          mainPanel(
                                br(),
                                conditionalPanel(condition="input.Sampmainmet=='Clustering'",
                                                 tabsetPanel(
                                                       tabPanel("Cluster plot",
                                                                hr(),
                                                                br(),
                                                                downloadButton("SampCluresdownloadbutton","Save Clustering Table"),
                                                                downloadButton("SampCluresdenddownloadbutton","Save Dendrogram"),
                                                                uiOutput("SampCluresselectclunumui"),
                                                                plotOutput("Sampcludendrogram",height="600px",width="1000px")
                                                       ),
                                                       tabPanel("Optimal Cluster Number",
                                                                plotOutput("Sampoptcluplot",width = "600px",height="600px")                                               
                                                       )
                                                 )
                                ),
                                conditionalPanel(condition="input.Sampmainmet=='Visualization'",
                                                 downloadButton("Sampvisplotdownload"),
                                                 uiOutput("Sampvisplotui") 
                                ),
                                conditionalPanel(condition="input.Sampmainmet=='Bulk'",
                                                 tabsetPanel(
                                                       tabPanel("Heatmap",br(),downloadButton("Sampbulkcorplotdownload"),d3heatmap::d3heatmapOutput("Sampbulkcorheatmap")),
                                                       tabPanel("Table",br(),downloadButton("Sampbulkcortabledownload"),DT:: dataTableOutput("Sampbulkcortable"))
                                                 )
                                )
                          )
                 ),
                 tabPanel("Step 4: Feature-level Analysis",
                          tags$head(
                                tags$style(HTML('#Featpreviousstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Featrunbutton{font-weight: bold;color:blue}'))
                          ),
                          hr(),
                          hr(),
                          sidebarPanel(
                                fluidRow(actionButton("Featpreviousstepbutton","Previous Step"),align="center"),
                                helpText("For sample clustering, identify key features that mostly explains the between cluster variance. By default, the analysis will be performed on multiple cluster numbers to ensure that the results are robust to different cluster numbers"),
                                wellPanel(
                                      h5("Choose number of sample clusters"),
                                      uiOutput("Featsampclunumui"),                                      
                                      actionButton("Featrunbutton","Calculate Feature Scores")                          
                                ),
                                wellPanel(
                                      selectInput("Featviewtab","Select table to view",list("F statistic"="Fstat","FDR"="FDR")),
                                      downloadButton("Featdownloadbutton")
                                )
                                ,width=3),
                          mainPanel(
                                br(),
                                tabsetPanel(
                                      tabPanel("Results",br(),DT:: dataTableOutput("Featrestable")),
                                      tabPanel("Summary",
                                               br(),
                                               uiOutput("FeatSumselectcolnameui"),
                                               textOutput("FeatSumtext"),
                                               plotOutput("FeatSumplot")
                                               )                                      
                                )                                
                          )
                 ),
                 tabPanel("About",
                          mainPanel(br(),
                                    br(),
                                    hr(),
                                    h4("SCRAT: Single-Cell Regulome Analysis Tool"),
                                    h5("Author: Zhicheng Ji, Weiqiang Zhou, Hongkai Ji"),
                                    h5("Maintainer: Zhicheng Ji (zji4@jhu.edu)"),
                                    h5("Version: 0.99.1")
                          )
                          
                 )
                 ,id="MainMenu",position = "fixed-top",inverse=T)
)

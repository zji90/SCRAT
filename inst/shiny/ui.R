######################################################
##                      SCRAT                       ##
##             Interactive User Interface           ##
##                     UI File                      ##
##   Author:Zhicheng Ji, Weiqiang Zhou, Hongkai Ji  ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################

library(shiny)
library(scatterD3)

shinyUI(      
      navbarPage("SCRAT",
                 
                 tabPanel("Step 1: Data input and preprocessing",
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
                                      checkboxInput("Inputblacklist","Filter blacklist",value = T),
                                      p(actionButton("Inputreadin","Read in"))      
                                ),
                                wellPanel(
                                      h4("Filter Bam Files"),
                                      textInput("InputFilterreads","Exclude samples with reads less than","500"),
                                      uiOutput("InputFiltersampleui")
                                ),
                                width=3),
                          
                          mainPanel(
                                uiOutput("Inputbamnumui"),
                                DT:: dataTableOutput("Inputbamsummary")
                          )),
                 tabPanel("Step 2: Feature summarization",
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
                                                       checkboxGroupInput("Sumselectmet","Choose Summarizing Method",list("ENCODE Cluster"="ENCL","Motif Sites"="MOTIF","GSEA Gene Sets"="GSEA","Upload BED"="Upload"),selected=c("ENCL","MOTIF","GSEA")),
                                                       hr(),
                                                       checkboxInput("Sumlogtf","Log2 transformation",value=T),
                                                       checkboxInput("Sumaddcv","Add coefficient of variation (sd/mean) information",value=T),
                                                       checkboxInput("Sumfiltertf","Filter Features",value=T),
                                                       conditionalPanel(condition="input.Sumfiltertf==1",
                                                                        textInput("Sumfilterpercen","Exclude features having more than","90"),
                                                                        textInput("Sumfilterreads","percent of samples whose (normalized) reads are less than","0.01"),
                                                                        textInput("Sumfiltercv","Exclude features with coefficient of variation (sd/mean) less than","0.01")
                                                       ),
                                                       hr(),
                                                       actionButton("Sumrunbutton","Run Summarization")
                                                 ),                        
                                                 wellPanel(
                                                       selectInput("Sumdetailchoose","Method Details",list("ENCODE Cluster"="ENCL","Motif Sites"="MOTIF","GSEA Gene Sets"="GSEA","Upload BED"="Upload")),
                                                       # conditionalPanel(condition="input.Sumdetailchoose=='TSS'",
                                                       #                  helpText('For each gene, sum all reads overlapping TSS upstream and downstream region. For gene with ENTREZ id of 18999, the feature name will be GENE:18999'),
                                                       #                  textInput("SumTSSupregionbp","TSS upstream base pair","1000"),
                                                       #                  textInput("SumTSSdownregionbp","TSS downstream base pair","500")
                                                       # ),
                                                       #                                       conditionalPanel(condition="input.Sumdetailchoose=='ENLO'",
                                                       #                                                        helpText('200 bp windows of genomic loci that are potential regulatory elements were precompiled based on ENCODE DNase-seq data.'),
                                                       #                                                        checkboxInput("SumENLOreducetf","Merge adjacent windows",value=F)
                                                       #                                       ),
                                                       conditionalPanel(condition="input.Sumdetailchoose=='ENCL'",
                                                                        helpText('Clusters of genomic regions (1000,2000 or 5000 clusters) were precompiled based on ENCODE DNase-seq data. For each cluster, sum all reads overlapping any of its genomic regions. For cluster id 1, the feature name will be ENCL1000:Cluster1'),
                                                                        radioButtons("SumENCLclunum","Choose number of clusters",c("1000","2000","5000"),"2000")
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
                                                                        textInput("SumGSEAupregionbp","TSS upstream base pair","3000"),
                                                                        textInput("SumGSEAdownregionbp","TSS downstream base pair","1000"),
                                                                        checkboxGroupInput("SumGSEAselect","Include gene sets from",list("Hallmark gene sets (50 gene sets)"="h.all","Positional gene sets (326 gene sets)"="c1.all","Curated gene sets: chemical and genetic perturbations (3395 gene sets)"="c2.cgp","Curated gene sets: canonical pathways (1330 gene sets)"="c2.cp","Motif gene sets: microRNA targets (221 gene sets)"="c3.mir","Motif gene sets: transcription factor targets (615 gene sets)"="c3.tft","Computational gene sets: cancer gene neighborhoods (427 gene sets)"="c4.cgn","Computational gene sets: cancer modules (431 gene sets)"="c4.cm","GO gene sets: biological process (825 gene sets)"="c5.bp","GO gene sets: cellular component (233 gene sets)"="c5.cc","GO gene sets: molecular function (396 gene sets)"="c5.mf","Oncogenic signatures (189 gene sets)"="c6.all","Immunologic signatures (1910 gene sets)"="c7.all"),selected = c("c5.bp","c5.cc","c5.mf")),
                                                                        helpText(a("Browse GSEA gene sets",href="http://software.broadinstitute.org/gsea/msigdb/genesets.jsp",target="_blank"))
                                                       ),
                                                       conditionalPanel(condition="input.Sumdetailchoose=='Upload'",
                                                                        radioButtons("SumuploadFilemultipletf",NULL,list("Upload One File"="One","Upload Multiple Files"="Multiple")),
                                                                        conditionalPanel(condition="input.SumuploadFilemultipletf=='One'",helpText('Upload BED file to define genomic regions. BED files are tab-delimited files. First column: chromosome name; second column: start site; third column: end site; Optional fourth column: feature id. If the fourth column is missing, then each row will be treated as seperate features. Otherwise rows with same feature id will be pulled together as a single feature.'),
                                                                                         fileInput('SumuploadFile', 'Choose File', accept = ".bed")
                                                                                         ),
                                                                        conditionalPanel(condition="input.SumuploadFilemultipletf=='Multiple'",helpText('Upload BED file to define genomic regions. BED files are tab-delimited files. First column: chromosome name; second column: start site; third column: end site. Each BED file will be treated as a single feature: all genomic features in one BED file will be pooled together.'),
                                                                                         fileInput('SumuploadFile2', 'Choose Files', accept = ".bed",multiple = T)
                                                                                         
                                                                        ),
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
                 tabPanel("Step 3: Cell heterogeneity analysis",
                          tags$head(
                                tags$style(HTML('#Samppreviousstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sampnextstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sampclurunbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Sampbulkcorrunbutton{font-weight: bold;color:blue}'))
                          ),
                          hr(),
                          hr(),
                          sidebarPanel(
                                fluidRow(actionButton("Samppreviousstepbutton","Previous Step"),
                                         actionButton("Sampnextstepbutton","Next Step"),align="center"),
                                helpText("Samples will be clustered using features obtained in step 2."),
                                br(),
                                wellPanel(uiOutput("Sampselectfeattypeui"),
                                          helpText("Select feature type to be included in the sample-level analysis. If no feature type is selected, all feature types will be used in the analsysis.")
                                ),
                                wellPanel(
                                      radioButtons("Sampmainmet","",list("Sample Clustering"="Clustering","Similarity to existing cell types"="Bulk"))
                                ),
                                wellPanel(
                                      conditionalPanel(condition="input.Sampmainmet=='Clustering'",
                                                       checkboxInput("Sampcludimredscale","Scale all features of each cell before dimension reduction and clustering (zero mean and unit variance).",value=T),
                                                       radioButtons("Sampcludimredmet","Dimension reduction method",list("PCA"="PCA","t-SNE"="tSNE","No Reduction"="None")),
                                                       conditionalPanel(condition="input.Sampcludimredmet=='tSNE'",helpText("Note that t-SNE could take some time to run.")),
                                                       conditionalPanel(condition="input.Sampcludimredmet=='PCA'",
                                                                        checkboxInput("Sampcluoptdimnum","Automatically choose optimal number of dimensions",value=T)
                                                       ),
                                                       conditionalPanel(condition="(input.Sampcludimredmet=='PCA' && input.Sampcluoptdimnum==0)||input.Sampcludimredmet=='tSNE'",textInput("Sampcluchoosedimnum","Choose number of dimensions",2)),
                                                       
                                                       hr(),
                                                       radioButtons("Sampcluclumet","Clustering method",list("Hierarchical Clustering"="hclust","K-means"="kmeans","Model-Based Clustering (mclust)"="mclust","Upload"="Upload")),
                                                       conditionalPanel(condition="input.Sampcluclumet!='Sample'",
                                                                        checkboxInput("Sampcluoptclunum","Automatically choose optimal number of clusters",value = T),
                                                                        conditionalPanel(condition="input.Sampcluoptclunum==0",textInput("Sampcluchooseclunum","Choose number of clusters",2))
                                                       ),
                                                       conditionalPanel(condition="input.Sampcluclumet=='Upload'",
                                                                        helpText("Upload a file specifying the cluster for each sample."),
                                                                        helpText("First column: sample name; Second column: cluster ID"),
                                                                        helpText("The sample names should be exactly the same as the names of the bam files. All bam files should be included."),
                                                                        helpText("The two columns should be seperated by Tab"),
                                                                        fileInput('SampcluInputFile', 'Choose File', accept = ".txt"),
                                                                        actionButton('Sampcluuploadbutton',"Upload"),
                                                                        uiOutput("Sampcluuploadstateui")
                                                       ),
                                                       hr(),
                                                       actionButton("Sampclurunbutton","Perform Clustering"),
                                                       hr()
                                                       
                                      ),
                                      conditionalPanel(condition="input.Sampmainmet=='Bulk'",
                                                       checkboxInput("Sampbulkcombinetf","Combine replicates",value=F),
                                                       actionButton("Sampbulkcorrunbutton","Calculate Correlations")              
                                      )
                                )
                                ,width=3),
                          mainPanel(
                                br(),
                                conditionalPanel(condition="input.Sampmainmet=='Clustering'",
                                                 tabsetPanel(
                                                       tabPanel("Clustering Result",
                                                                br(),
                                                                downloadButton("Sampcluplotdownload","Download Plot"),
                                                                downloadButton("Sampclutabledownload","Download Clustering Table"),
                                                                br(),
                                                                wellPanel(
                                                                      checkboxInput("Sampcluplotorifeattf","Plot original features",value=F),
                                                                      uiOutput("Sampcluplotselectfeatui"),
                                                                      conditionalPanel(condition="input.Sampcludimredmet=='PCA'",checkboxInput("Sampselectfeatincludebulktf","Include existing cell types")),
                                                                      conditionalPanel(condition="input.Sampselectfeatincludebulktf==1",
                                                                                       uiOutput("Sampselectfeatincludebulkui")),
                                                                      uiOutput("Sampcluplotdim2optionui")
                                                                ),
                                                                conditionalPanel("input.Sampcluplotselectfeat.length==2",
                                                                                 scatterD3::scatterD3Output("Sampcluplotdim2",width="700px",height="500px")      
                                                                ),
                                                                conditionalPanel("input.Sampcluplotselectfeat.length==1",
                                                                                 plotOutput("Sampcluplotdim1",width="600px",height="500px")
                                                                ),
                                                                conditionalPanel("input.Sampcluplotselectfeat.length > 2",
                                                                                 plotOutput("Sampcluplotdim3",width="500px",height="2000px")
                                                                )),
                                                       tabPanel("Clustering Diagnosis",
                                                                textInput("Sampclucludiagsimnum","Number of simulations","1000"),
                                                                actionButton("Sampclucludiagrun","Run simulations to diagnose clustering results."),
                                                                plotOutput("Sampclucludiagplot"),
                                                                textOutput("SampclucludiagtrueSS"),
                                                                textOutput("SampclucludiagsimuSS")
                                                       )
                                                 )
                                ),
                                conditionalPanel(condition="input.Sampmainmet=='Bulk'",
                                                 tabsetPanel(
                                                       tabPanel("Heatmap",br(),checkboxInput("Sampbulkcorplotstdtf","Standardize each row to have zero mean and unit variance",value=T),downloadButton("Sampbulkcorplotdownload"),d3heatmap::d3heatmapOutput("Sampbulkcorheatmap"),width="500px",height="500px"),
                                                       tabPanel("Table",br(),downloadButton("Sampbulkcortabledownload"),DT::dataTableOutput("Sampbulkcortable"))
                                                 )
                                )
                          )
                 ),
                 tabPanel("Step 4: Differential feature analysis",
                          tags$head(
                                tags$style(HTML('#Featpreviousstepbutton{font-weight: bold;color:blue}')),
                                tags$style(HTML('#Featrunbutton{font-weight: bold;color:blue}'))
                          ),
                          hr(),
                          hr(),
                          sidebarPanel(
                                fluidRow(actionButton("Featpreviousstepbutton","Previous Step"),align="center"),
                                helpText("Perform ANOVA or t tests to identify key features that mostly explains the between cluster variance. The sample clusters are obtained in step 3."),
                                wellPanel(uiOutput("Featselectfeattypeui"),
                                          helpText("Select feature type to be included in the sample-level analysis. If no feature type is selected, all feature types will be used in the analsysis.")
                                ),    
                                wellPanel(
                                      checkboxInput("Featrunallclustertf","Perform tests for all clusters",value=T),
                                      conditionalPanel(condition="input.Featrunallclustertf==0",uiOutput("Featselectclusterui")),
                                      uiOutput("Featttestaltui"),
                                      uiOutput("Featttestalttextui"),
                                      actionButton("Featrunbutton","Perform Test")      
                                )
                                ,width=3),
                          mainPanel(
                                br(),
                                tabsetPanel(
                                      tabPanel("Results",br(),downloadButton("Featdownloadbutton"),br(),DT::dataTableOutput("Featrestable")),
                                      tabPanel("Summary",
                                               br(),
                                               textOutput("FeatSumtext"),
                                               radioButtons("Featviewtab","",c("statistics","FDR")),
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

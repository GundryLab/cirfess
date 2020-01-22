# SurfaceGenie_0.1/server.R
library(RSQLite)
library(shiny)
library(plotly)
library(RColorBrewer) 
library(stringr)
library(svglite)
library(dplyr)
library(DT)
#library(xlsx)
source("functions.R")

conn <- dbConnect(RSQLite::SQLite(), "matt.db")


function(input, output, session) {
  
  ##########    Single Protein Lookup   ##########
  
  ##### SPC Table ##### 
  
  SPCtable <- reactive({ 
    sql <- paste0("SELECT Accession, SPC, SPCstring AS SPC_DB FROM prot WHERE Accession = '", toupper(input$swissprtID), "';")
    dbGetQuery(conn, sql)
  })
  
  output$SPC <- renderTable({
    spc<-SPCtable()
    spc$SPC_DB<-translateSPC(spc$SPC_DB)
    return(spc)
  })
  
  ##### Signal Peptide Predictions ##### 
  
  SigPeptable <- reactive({ 
    sql <- paste0("SELECT SigPepPhobius, ScorePhobius, SigPepSignalP, ScoreSignalP, SigPepPrediSI, ScorePrediSI FROM prot WHERE Accession = '", toupper(input$swissprtID), "';")
    dbGetQuery(conn, sql)
  })
  
  output$SigPep <- renderTable({
    SigPeptable()
  })
  
  ##### Topology Table ##### 
  
  Topotable <- reactive({
    select <- 'SELECT numTMPhobius AS "#TMPhobius", numICPhobius AS "#ICPhobius", numECPhobius AS "#ECPhobius", numTMTMHMM AS "numTMTMHMM", numICTMHMM AS "#ICTMHMM", numECTMHMM AS "#ECTMHMM", StringOutPhobius, StringOutTMHMM'
    sql <- paste0(select, " FROM prot WHERE Accession = '", toupper(input$swissprtID), "';")
    dbGetQuery(conn, sql)
  })
  
  output$Topo <- renderTable({
    Topotable()
  })
  
  # ##### Motif Summary Table ##### 
  # 
  # MotifSummarytable <- reactive({ 
  #   sql <- paste0("SELECT * FROM pep WHERE Accession = '",input$swissprtID, "';")
  #   peps <- dbGetQuery(conn, sql)
  #   data.frame(summarise(group_by(peps, numMissedCleavages), NXS=sum(numMotifsPhobiusNXS), NXT=sum(numMotifsPhobiusNXT), NXC=sum(numMotifsPhobiusNXC), NXV=sum(numMotifsPhobiusNXV), C=sum(numMotifsPhobiusC), K=sum(numMotifsPhobiusK)))
  # })
  # 
  # output$MotifSummary <- renderTable({
  #   MotifSummarytable()
  # })
  
  ##### Peptide Summary Table ##### 
  
  PepSummarytable <- reactive({
    select <- 'SELECT numMissedCleavages, OKforMS, numMotifsNXS, numMotifsNXT, numMotifsNXC, numMotifsNXV, numMotifsC, numMotifsK'
#    sql <- paste0("SELECT * FROM pep WHERE Accession = '", toupper(input$swissprtID), "';")
    sql <- paste0(select, " FROM pep WHERE Accession = '", toupper(input$swissprtID), "';")
    peps <- dbGetQuery(conn, sql)
    data.frame(summarise(group_by(peps, numMissedCleavages, OKforMS), NXS=sum(numMotifsNXS), NXT=sum(numMotifsNXT), NXC=sum(numMotifsNXC), NXV=sum(numMotifsNXV), C=sum(numMotifsC), K=sum(numMotifsK)))
#    colnames(d) <- c("#MissedCleavages", "OKforMS", "NXS", "NXT", "NXC", "NXV", "C", "K")
  })
  
  output$PepSummary <- renderTable({
    PepSummarytable()
  })
  
  translateSPC <- function(s) {
    s<-gsub("BF","Surfy",s)
    s<-gsub("T","Town",s)
    s<-gsub("dC","da Cunha",s)
    s<-gsub("DR","Diaz Ramos",s)
  }
  
    ##########          Batch Lookup TextArea         ##########

  
  
  
  ##########          Bulk Lookup TextArea         ##########
  
  

  
  
  
  
    ##########          SurfaceGenie         ##########
  
  # Load and process data
  data_input <- reactive({
    withProgress(message = 'Reading Data', value = 0, {
      parts <- strsplit(input$file1$name, "\\.")[[1]]
      ext <- parts[length(parts)]
      if(ext == "csv"){
        df <- read.csv(input$file1$datapath, header=TRUE)
      } else if(ext=="txt" || ext=="tsv" || ext=="tab") {
        df <- read.delim(input$file1$datapath, header=TRUE)
      } else if(ext=="xlsx" || ext=="xls") {
        df <- read.xlsx(input$file1$datapath,  1)
      } else {
        validate(need(ext=="csv" || ext=="tsv" || ext=="xlsx" || ext=="xls" || ext=="tab" || ext=="txt", "You have the wrong file extension.  Please see the instructions for possible extensions and associated file types." ))
      }
      
      # validation
      validate(
        need(ncol(df)>1,"There is only one column of data.  This is probably caused by one of two things.  Either you used the wrong field delimiter - e.g. tabs instead of commas - or you need to enter one or more samples and their abundance data.")
      )

      # validation
      validate(
        need(colnames(df)[1]=="Accession",'The name of the first column must be "Accession"')
      )
      
      
      # validation
      validate(
        need(ncol(df)>1,"There is only one column of data.  This is probably caused by one of two things.  Either you used the wrong field delimiter - e.g. tabs instead of commas - or you need to enter one or more samples and their abundance data."),
        need(length(which(as.vector(sapply(df[2:3],is.numeric))==FALSE))==0,"There is a non-numeric value in the abundance data columns.  The first column should be accessions and each column after that should be numbers representing the abundance of the protein for that sample.  Please refer to the instructions for more help.")
      )

      ## look for duplicates and for NA values.  These errors will not stop the program, but we will inform the user.
      # look for duplicates
      warningMsg <- ""
      e <- as.vector(df[duplicated(df),][,1])
      if( length(e)>0 ){
        warningMsg <- "You have duplicates in the list of accessions.  You may wish to address that."
      }
      # look for NA values.  This is not necessarily wrong, but want to make the user informed.
      if(length(which(is.na(df)))>0) {
        warningMsg <- c(warningMsg, "There are missing values in your data.  You may wish to address that.")
      }
      output$txtWarning<-renderText(warningMsg)

      input_size <- c(nrow(df), ncol(df))
      list(df, input_size)
    })
  })
  

  # Load annotation file
  annotation <- reactive({
#    print("Annotation")
    withProgress(message = 'Reading Annotations', value = 0, {
      if(input$species == "Human") {
        read.delim("ref/anno.tsv")
      } else if(input$species == "Rat") {
        read.delim("ref/annotation.rat.tsv")
      } else if(input$species == "Mouse") {
        read.delim("ref/annotation.mouse.tsv")
      } else if(input$species == "Other/Ignore") {
        read.delim("ref/annotation.none.tsv")
      }
    })
  })
    
  
  
  ##########  SurfaceGenie: Data Grouping  ##########
  
  data_output <- reactive({

    # setting up the progress meter
    progress<-shiny::Progress$new()
    progress$set(message = "Calculating Scores", value=0)
    on.exit(progress$close())
    
    # creating the progress function to pass to the Surface Genie function
    updateProgress <- function(value=NULL,detail=NULL) {
      if(is.null(value)){
        value <- progress$getValue()
        value <- value + (progress$getMax()-value)/10
      }
      progress$set(value=value, detail=detail)
    }

        
    if("grouping" %in% input$processing_opts){
      gtags <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
      groupcols <- list()
      if(input$numgroups >= 2){
        validate(
          need(input$group1, "Please indicate columns in Group 1") %then%
            need(laply(strsplit(input$group1, ","), as.integer), 
                 "Group 1 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group1,",")[[1]]),
                 "Group 1 Error: Column 1 should contain accession numbers and can not be grouped" ),
          need(input$group2, "Please indicate columns in Group 2") %then%
            need(laply(strsplit(input$group2, ",")[[1]], as.integer), 
                 "Group 2 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group2,",")[[1]]),
                 "Group 2 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[1] <- input$group1
        groupcols[2] <- input$group2
      }
      if(input$numgroups >= 3){
        validate(
          need(input$group3, "Please indicate columns in Group 3") %then%
            need(laply(strsplit(input$group3, ",")[[1]], as.integer), 
                 "Group 3 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group3,",")[[1]]),
                 "Group 3 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[3] <- input$group3
      }
      if(input$numgroups >= 4){
        validate(
          need(input$group4, "Please indicate columns in Group 4") %then%
            need(laply(strsplit(input$group4, ",")[[1]], as.integer), 
                 "Group 4 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group4,",")[[1]]),
                 "Group 4 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[4] <- input$group4
      }
      if(input$numgroups >=5){
        validate(
          need(input$group5, "Please indicate columns in Group 5") %then%
            need(laply(strsplit(input$group5, ",")[[1]], as.integer), 
                 "Group 5 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group5,",")[[1]]),
                 "Group 5 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[5] <- input$group5
      }
      # call this if grouped
      sg<-SurfaceGenie(data_input()[[1]], input$processing_opts, 
                   input$groupmethod, input$numgroups, groupcols, annotation(), updateProgress)
      
    }
    else{
      # call this if not grouped
      sg<-SurfaceGenie(data_input()[[1]], input$processing_opts, 
                   groupmethod=NULL, numgroups=0, groupcols=NULL, annotation(), updateProgress)
    }
    if(length(which(as.vector(sapply(sg["SPC"],is.na))==TRUE))>0) {
      warningMsg <- "There are accessions that do not belong to the most recent swissprot version for the species you selected. They will be ignored for SurfaceGenie Score and IsoGenie Score. Download the csv file to find which ones."
      output$txtWarning<-renderText(warningMsg)
    }
    return(sg)
  })
  
  ##########  SurfaceGenie: Output Display  ##########
  
  # Apply export options
  data_export <- reactive({
    if(input$species == "Human") {
      df <- SG_export(data_output(), input$export_options1, input$export_options2h, input$scoring_opts)
    } else if(input$species == "Rat") {
      df <- SG_export(data_output(), input$export_options1, input$export_options2r, input$scoring_opts)
    } else if(input$species == "Mouse") {
      df <- SG_export(data_output(), input$export_options1, input$export_options2m, input$scoring_opts)
    } else if(input$species == "Other/Ignore") {
      df <- SG_export(data_output(), input$export_options1o, input$export_options2, input$scoring_opts_o)
    }
          
    output_size <- c(nrow(df), ncol(df))
    list(df, output_size)
  })
  
  # Example data to display
  example_data <- read.csv("ref/example_data.csv", header=TRUE)
  output$example_data <- renderTable(head(example_data, 5))
  
  # Display input data
  output$data_input <- renderTable({
    req(input$file1)
    head(data_input()[[1]], 10)
  })
  output$input_size <- renderText({
    req(input$file1)
    sprintf("[ %d rows x %d columns ]", data_input()[[2]][1], data_input()[[2]][2])
  })
  
  # Display output data
  output$data_output <- renderTable({
    req(input$file1)
    head(data_export()[[1]], 10)
  })
  output$output_size <- renderText({
    req(input$file1)
    sprintf("[ %d rows x %d columns ]", data_export()[[2]][1], data_export()[[2]][2])
  })

  # SG: SPC histogram
  output$SG_SPC_hist <- renderPlot({
    req(input$file1)
    SPC_hist(data_output())
  })
  
  output$SG_SPC_hist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_SPC_hist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=3*300,
          res=300)
      SPC_hist(data_output())
      dev.off()
    }
  )
  output$SG_SPC_hist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_SPC_hist_PNGdl", " .png", class="download_this")
  })
  output$SG_SPC_hist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_SPC_hist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=3,
          pointsize=10)
      SPC_hist(data_output())
      dev.off()
    }
  )
  output$SG_SPC_hist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_SPC_hist_SVGdl", " .svg", class="download_this")
  })

    # SG: Genie Score plot
  output$SG_dist <- renderPlotly({
    req(input$file1)
    SG_dist(data_output())
  })
  
  output$SG_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      SG_dist(data_output())
      SG_dist_export(data_output())
      dev.off()
    }
  )
  output$SG_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_dist_PNGdl", " .png", class="download_this")
  })
  output$SG_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      SG_dist_export(data_output())
      dev.off()
    }
  )
  output$SG_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_dist_SVGdl", " .svg", class="download_this")
  })
  
  # IsoGenie: reverse Genie Score plot
  output$IsoGenie_dist <- renderPlotly({
    req(input$file1)
    IsoGenie_dist(data_output())
  })
  
  output$IsoGenie_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoGenie_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      IsoGenie_dist(data_output())
      IsoGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoGenie_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoGenie_dist_PNGdl", " .png", class="download_this")
  })
  output$IsoGenie_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoGenie_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      IsoGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoGenie_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoGenie_dist_SVGdl", " .svg", class="download_this")
  })

  # OmniGenie: SurfaceGenie w/o SPC Score plot
  output$OmniGenie_dist <- renderPlotly({
    req(input$file1)
    OmniGenie_dist(data_output())
  })
  
  output$OmniGenie_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_OmniGenie_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      OmniGenie_dist(data_output())
      OmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$OmniGenie_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("OmniGenie_dist_PNGdl", " .png", class="download_this")
  })
  output$OmniGenie_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_OmniGenie_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      OmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$OmniGenie_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("OmniGenie_dist_SVGdl", " .svg", class="download_this")
  })
  
  
  # IsoOmniGenie: Similarity w/o SPC score plot
  output$IsoOmniGenie_dist <- renderPlotly({
    req(input$file1)
    IsoOmniGenie_dist(data_output())
  })
  
  output$IsoOmniGenie_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoOmniGenie_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      IsoOmniGenie_dist(data_output())
      IsoOmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoOmniGenie_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoOmniGenie_dist_PNGdl", " .png", class="download_this")
  })
  output$IsoOmniGenie_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoOmniGenie_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      IsoOmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoOmniGenie_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoOmniGenie_dist_SVGdl", " .svg", class="download_this")
  })
  

  # Downloadable csv of selected dataset
  output$csv_download <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, "_SurfaceGenie.csv", sep = "")
    },
    content = function(filename) {
      write.csv(data_export()[[1]], filename, row.names = FALSE)
    }
  )
  output$csv_dlbutton <- renderUI({
    req(input$file1)
    downloadButton("csv_download", ".csv", class="download_this")
  })
  
  # Downloadable tsv of selected dataset
  output$tsv_download <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, "_SurfaceGenie.tsv", sep = "")
    },
    content = function(filename) {
      write.table(data_export()[[1]], filename, sep="\t", row.names = FALSE)
    }
  )
  output$tsv_dlbutton <- renderUI({
    req(input$file1)
    downloadButton("tsv_download", ".tsv", class="download_this")
  })

  # Downloadable xlsx of selected dataset
  output$xlsx_download <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, "_SurfaceGenie.xlsx", sep = "")
    },
      content = function(filename) {
      write.xlsx(data_export()[[1]], filename, row.names = FALSE)
    }
  )
  output$xlsx_dlbutton <- renderUI({
    req(input$file1)
    downloadButton("xlsx_download", ".xlsx", class="download_this")
  })
  

  
  ##########          Batch Lookup TextArea         ##########
  
  
  pepSQL <- function(pepOptions, cleavages) {
    sql <- 'SELECT Accession, ID'
    if( 'LocProtein' %in% pepOptions ) {
      sql <- paste0( sql, ', range') 
    }
    if( 'Missed' %in% pepOptions ) {
      sql <- paste0( sql, ', numMissedCleavages AS misClev') 
    }
    if( 'OK' %in% pepOptions ) {
      sql <- paste0( sql, ', OKforMS') 
    }
    if( 'NumMotifs' %in% pepOptions ) {
      sql <- paste0( sql, ', numMotifsNXS, numMotifsNXT, numMotifsNXC, numMotifsNXV, numMotifsC, numMotifsK, numMotifsPhobiusNXS, numMotifsPhobiusNXT, numMotifsPhobiusNXC, numMotifsPhobiusNXV, numMotifsPhobiusC, numMotifsPhobiusK, numMotifsTMHMMNXS, numMotifsTMHMMNXT, numMotifsTMHMMNXC, numMotifsTMHMMNXV, numMotifsTMHMMC, numMotifsTMHMMK')
    }
    if( 'LocMotifs' %in% pepOptions ) {
      sql <- paste0( sql, ', motifLocNXS, motifLocNXT, motifLocNXC, motifLocNXV, motifLocC, motifLocK, motifLocPhobiusNXS, motifLocPhobiusNXT, motifLocPhobiusNXC, motifLocPhobiusNXV, motifLocPhobiusC, motifLocPhobiusK, motifLocTMHMMNXS, motifLocTMHMMNXT, motifLocTMHMMNXC, motifLocTMHMMNXV, motifLocTMHMMC, motifLocTMHMMK')
    }
    sql <- paste0(sql, ', pepSeq FROM pep')
    return(sql)
  }
  
  protSQL <- function(protOptions, cleavages) {
    sql <- 'SELECT Accession, Length'
    if( 'SPC' %in% protOptions ) {
      sql <- paste0( sql, ', SPC, SPCstring') 
    }
    if( 'Signal' %in% protOptions ) {
      sql <- paste0( sql, ', SigPepPhobius, ScorePhobius, SigPepSignalP, ScoreSignalP, SigPepPrediSI, ScorePrediSI, numSPpedictions') 
    }
    if( 'Topo' %in% protOptions ) {
      sql <- paste0( sql, ', numTMPhobius, numICPhobius, numECPhobius, numTMTMHMM, numICTMHMM, numECTMHMM, StringOutPhobius, StringOutTMHMM') 
    }
    if( 'ProtMotif' %in% protOptions ) {
      sql <- paste0( sql, ', numMotifsNXS, motifsLocNXS, numMotifsNXT, motifsLocNXT, numMotifsNXC, motifsLocNXC, numMotifsNXV, motifsLocNXV, numMotifsC, motifsLocC, numMotifsK, motifsLocK, numMotifsPhobiusNXS	motifsLocPhobiusNXS, numMotifsPhobiusNXT, motifsLocPhobiusNXT, numMotifsPhobiusNXC, motifsLocPhobiusNXC, numMotifsPhobiusNXV, motifsLocPhobiusNXV, numMotifsPhobiusC, motifsLocPhobiusC, numMotifsPhobiusK, motifsLocPhobiusK, numMotifsTMHMMNXS, motifsLocTMHMMNXS, numMotifsTMHMMNXT, motifsLocTMHMMNXT, numMotifsTMHMMNXC, motifsLocTMHMMNXC, numMotifsTMHMMNXV')
    }
    if( 'Pep' %in% protOptions ) {
      sql <- paste0( sql, ', numPepMC0, numPepMC1, numPepMC2, numPepOKMC0, numPepOKMC1, numPepOKMC2, numPepMC0NXS, numPepMC0NXT, numPepMC0NXC, numPepMC0NXV, numPepMC0C, numPepMC0K, numPepMC1NXS, numPepMC1NXT, numPepMC1NXC, numPepMC1NXV, numPepMC1C, numPepMC1K, numPepMC2NXS, numPepMC2NXT, numPepMC2NXC, numPepMC2NXV, numPepMC2C, numPepMC2K, numPepMC0PhobiusNXS, numPepMC0PhobiusNXT, numPepMC0PhobiusNXC, numPepMC0PhobiusNXV, numPepMC0PhobiusC, numPepMC0PhobiusK, numPepMC1PhobiusNXS, numPepMC1PhobiusNXT, numPepMC1PhobiusNXC, numPepMC1PhobiusNXV, numPepMC1PhobiusC, numPepMC1PhobiusK, numPepMC2PhobiusNXS, numPepMC2PhobiusNXT, numPepMC2PhobiusNXC, numPepMC2PhobiusNXV, numPepMC2PhobiusC, numPepMC2PhobiusK, numPepMC0TMHMMNXS, numPepMC0TMHMMNXT, numPepMC0TMHMMNXC, numPepMC0TMHMMNXV, numPepMC0TMHMMC, numPepMC0TMHMMK, numPepMC1TMHMMNXS, numPepMC1TMHMMNXT, numPepMC1TMHMMNXC, numPepMC1TMHMMNXV, numPepMC1TMHMMC, numPepMC1TMHMMK, numPepMC2TMHMMNXS, numPepMC2TMHMMNXT, numPepMC2TMHMMNXC, numPepMC2TMHMMNXV, numPepMC2TMHMMC, numPepMC2TMHMMK')
    }
    sql <- paste0( sql, ", Seq FROM prot")
    return(sql)
  }
  
  # Load and process data
  Batch_prot_lookup <- function(accs) {
    if(input$ProtLevelExp) {
      sql <- protSQL(input$export_prot_level, 0)
      sql <- paste0(sql, " WHERE Accession IN('" )
      l <- paste( unlist(accs), collapse="', '")
      sql <- paste0( sql, l)
      sql <- paste0(sql, "');")
#      print(sql)
      q <-dbGetQuery(conn, sql)
    }
    return(q)
  }

  # Load and process data
  Batch_pep_lookup <- function(accs) {
    if(input$PepLevelExp) {
      sql <- pepSQL(input$export_pep_level, 0)
      sql <- paste0(sql, " WHERE Accession IN('" )
      l <- paste( unlist(accs), collapse="', '")
      sql <- paste0( sql, l)
      sql <- paste0(sql, "');")
#      print(sql)
      q <-dbGetQuery(conn, sql)
    }
    return(q)
  }
  
  # get accessions from textarea box
  Batch_Text_input <- reactive({
#    print(input$quicklookup)
    u<-toupper(input$quicklookup)
#    accs <- strsplit(input$quicklookup, "\\n")
    accs <- strsplit(u, "\\n")
#    print(accs)
    return(accs)
#    return(toupper(accs))
  })

  Batch_Text_prot_output <- reactive({
    Batch_prot_lookup(Batch_Text_input())
  })

  Batch_Text_pep_output <- reactive({
    Batch_pep_lookup(Batch_Text_input())
  })
  
  # Display datarenderDataTable
#  output$Batch_Text_prot_output <- renderTable({
  output$Batch_Text_prot_output <- renderDataTable({
    req(input$quicklookup, input$ProtLevelExp)
    Batch_Text_prot_output()
  })
  
  
  output$Batch_Text_pep_output <- renderDataTable({
    req(input$quicklookup, input$PepLevelExp)
    Batch_Text_pep_output()
  })

  # Load and process data
  Batch_File_input <- reactive({
    withProgress(message = 'Reading Data', value = 0, {
      parts <- strsplit(input$file2$name, "\\.")[[1]]
      ext <- parts[length(parts)]
      if(ext == "csv"){
        df <- read.csv(input$file2$datapath)
      } else if(ext=="txt" || ext=="tsv" || ext=="tab") {
        df <- read.delim(input$file2$datapath)
      } else if(ext=="xlsx" || ext=="xls") {
        df <- read.xlsx(input$file2$datapath,  1)
      } else {
        validate(need(ext=="csv" || ext=="tsv" || ext=="xlsx" || ext=="xls" || ext=="tab" || ext=="txt", "You have the wrong file extension.  Please use .csv, .tsv, .txt, .tab, .xlsx, or .xls depending on your file type." ))
      }
      toupper(df[,1])
    })
  })

  Batch_File_prot_output <- reactive({
    Batch_prot_lookup(Batch_File_input())
  })
  
  Batch_File_pep_output <- reactive({
    Batch_pep_lookup(Batch_File_input())
  })
  
  
  output$Batch_File_prot_output <- renderDataTable({
    req(input$file2, input$ProtLevelExp)
    Batch_File_prot_output()
  })
  
  
  output$Batch_File_pep_output <- renderDataTable({
    req(input$file2, input$PepLevelExp)
    Batch_File_pep_output()
  })
  
  
  
#Prot level Option 1
  # Downloadable csv of selected dataset ----
  output$Batch_Text_prot_download_csv <- downloadHandler(
    filename = function() {
      paste0("ProteinLevelData-", Sys.Date(), ".csv")
    },
    content = function(filename) {
      write.csv(Batch_Text_prot_output(), filename, row.names = FALSE)
    }
  )
  output$Batch_Text_prot_dlbutton_csv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_Text_prot_download_csv", ".csv", class="download_this")
  })
  
  # Downloadable tsv of selected dataset ----
  output$Batch_Text_prot_download_tsv <- downloadHandler(
    filename = function() {
      paste0("ProteinLevelData-", Sys.Date(), ".tsv")
    },
    content = function(filename) {
#      print(Batch_Text_prot_output())
      write.table(Batch_Text_prot_output(), filename, sep="\t", row.names = FALSE)
    }
  )
  output$Batch_Text_prot_dlbutton_tsv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_Text_prot_download_tsv", ".tsv", class="download_this")
  })

  
  
#Pep level Option 1
  # Downloadable csv of selected dataset ----
  output$Batch_Text_pep_download_csv <- downloadHandler(
    filename = function() {
      paste0("ProteinLevelData-", Sys.Date(), ".csv")
    },
    content = function(filename) {
      write.csv(Batch_Text_pep_output(), filename, row.names = FALSE)
    }
  )
  output$Batch_Text_pep_dlbutton_csv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_Text_pep_download_csv", ".csv", class="download_this")
  })
  
  # Downloadable tsv of selected dataset ----
  output$Batch_Text_pep_download_tsv <- downloadHandler(
    filename = function() {
      paste0("PeptideLevelData-", Sys.Date(), ".tsv")
    },
    content = function(filename) {
      write.table(Batch_Text_pep_output(), filename, sep="\t", row.names = FALSE)
    }
  )
  output$Batch_Text_pep_dlbutton_tsv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_Text_pep_download_tsv", ".tsv", class="download_this")
  })
  

  
  
#Prot level Option 1
  # Downloadable csv of selected dataset ----
  output$Batch_File_prot_download_csv <- downloadHandler(
    filename = function() {
      paste0("ProteinLevelData-", Sys.Date(), ".csv")
    },
    content = function(filename) {
      write.csv(Batch_File_prot_output(), filename, row.names = FALSE)
    }
  )
  output$Batch_File_prot_dlbutton_csv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_File_prot_download_csv", ".csv", class="download_this")
  })
  
  # Downloadable tsv of selected dataset ----
  output$Batch_File_prot_download_tsv <- downloadHandler(
    filename = function() {
      paste0("ProteinLevelData-", Sys.Date(), ".tsv")
    },
    content = function(filename) {
#      print(Batch_File_prot_output())
      write.table(Batch_File_prot_output(), filename, sep="\t", row.names = FALSE)
    }
  )
  output$Batch_File_prot_dlbutton_tsv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_File_prot_download_tsv", ".tsv", class="download_this")
  })
  
  
  
  #Pep level Option 1
  # Downloadable csv of selected dataset ----
  output$Batch_File_pep_download_csv <- downloadHandler(
    filename = function() {
      paste0("ProteinLevelData-", Sys.Date(), ".csv")
    },
    content = function(filename) {
      write.csv(Batch_File_pep_output(), filename, row.names = FALSE)
    }
  )
  output$Batch_File_pep_dlbutton_csv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_File_pep_download_csv", ".csv", class="download_this")
  })
  
  # Downloadable tsv of selected dataset ----
  output$Batch_File_pep_download_tsv <- downloadHandler(
    filename = function() {
      paste0("PeptideLevelData-", Sys.Date(), ".tsv")
    },
    content = function(filename) {
      write.table(Batch_File_pep_output(), filename, sep="\t", row.names = FALSE)
    }
  )
  output$Batch_File_pep_dlbutton_tsv <- renderUI({
    req(input$quicklookup)
    downloadButton("Batch_File_pep_download_tsv", ".tsv", class="download_this")
  })
  
}




# output$text_prot_size <- renderText({
#   req(input$quicklookup, input$ProtLevelExp)
#   #    print("got here")
#   #    sprintf("[ %d rows x %d columns ]", Batch_Text_prot_output()[[2]][1], Batch_Text_prot_output()[[2]][2])
# })


# SPC_bulk_output <- reactive({
#   SPC_lookup(SPC_bulk_input(), input$species2)
# })
# 
# SPC_bulk_output_for_export <- reactive({
#   SPC_lookup_for_export(SPC_bulk_input(), input$species2)
# })
# 
# # Display data
# output$SPC_bulk_output <- renderTable({
#   req(input$file2)
#   head(SPC_bulk_output(), 10)
# })
# output$SPC_bulk_hist <- renderPlot({
#   req(input$file2)
#   SPC_hist(SPC_bulk_output())
# })

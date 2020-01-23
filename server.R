# SurfaceGenie_0.1/server.R
library(RSQLite)
library(shiny)
library(dplyr)
library(DT)
#library(xlsx)

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
  
  


  # # Downloadable xlsx of selected dataset
  # output$xlsx_download <- downloadHandler(
  #   filename = function() {
  #     fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
  #     paste(fname, "_SurfaceGenie.xlsx", sep = "")
  #   },
  #     content = function(filename) {
  #     write.xlsx(data_export()[[1]], filename, row.names = FALSE)
  #   }
  # )
  # output$xlsx_dlbutton <- renderUI({
  #   req(input$file1)
  #   downloadButton("xlsx_download", ".xlsx", class="download_this")
  # })
  

  
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
    u<-toupper(input$quicklookup)
    accs <- strsplit(u, "\\n")
    return(accs)
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





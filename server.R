#library(RSQLite)
library(RMariaDB)
library(pool)
library(shiny)
library(dplyr)
library(DT)
library(plotly)
library(xtable)
library(httr)
library(tidyr) #for pivot_longer
#library(xlsx)

if( file.exists("www/data/P14384.json") ) {
  print("already unzipped")
} else {
  print("unzipping")
  unzip('all.zip')
}

#conn <- dbConnect(RSQLite::SQLite(), "cirfess.db")


conn <- dbPool(
  drv = RMariaDB::MariaDB(),
  default.file="./ssl/my.cnf",
  group="cirfess",
  dbname='cirfess'
)

function(input, output, session) {

##################################################  
###########    Single Protein Lookup   ###########
##################################################  

  protData <- eventReactive( input$go, {
    sqlstmnt <- paste0("SELECT * FROM prot WHERE Accession = '", toupper(input$swissprtID), "';")
    result <- dbGetQuery(conn, sqlstmnt) 
    return( result )
  })
  
  pepData <- eventReactive( input$go, {
    sqlstmnt <- paste0("SELECT * FROM pep WHERE OKforMS != 0 AND Accession = '", toupper(input$swissprtID), "';")
    result <- dbGetQuery(conn, sqlstmnt)
    return( result )
  })

  
######################    Protein Level Stuff   ######################
    
      
  ##### SPC Table ##### 
  
  output$SPC <- renderTable({
    req(length(protData()[,1])>0)
    spc <- protData()[c("Accession", "SPC", "SPCstring")]
    names(spc)[names(spc) == 'SPCstring'] <- 'SPC_DB'
    spc$SPC_DB<-translateSPC(spc$SPC_DB)
    return(spc)
  }, striped = TRUE, hover=TRUE, bordered = TRUE)
  
  ##### Signal Peptide Predictions ##### 
  
  output$SigPep <- renderTable({
    req(length(protData()[,1])>0)
    r<-protData()
    Method<-c('Phobius', 'SignalP', 'PrediSI')
    Prediction<-c(r$SigPepPhobius, r$SigPepSignalP, r$SigPepPrediSI)
    Score<-c(r$ScorePhobius, r$ScoreSignalP, r$ScorePrediSI)
    df<-data.frame(Method, Prediction, Score)
  }, striped = TRUE, hover=TRUE, bordered = TRUE)
  
  ##### Topology Table ##### 
  
  output$Topo <- renderTable({
    req(length(protData()[,1])>0)
    r<-protData()
    Measure <- c('EC Residues', 'TM Residues', 'IC Residues', 'Short Format')
    Phobius <- c(r$numECPhobius, r$numTMPhobius, r$numICPhobius, r$StringOutPhobius)
    TMHMM <- c(r$numECTMHMM, r$numTMTMHMM, r$numICTMHMM, r$ StringOutTMHMM)
    df <- data.frame(Phobius, TMHMM)
    rownames(df) <- Measure
    return(df)
  }, rownames = TRUE, striped = TRUE, hover=TRUE, bordered = TRUE)

    
  ##### Phobius Table ##### 
  
  output$pepPhobius <-renderTable({
    req(length(protData()[,1])>0)
    filtered <- filter(pepData(), OKforMS != 0 )
    phobius <- filtered[c("numMissedCleavages", "numMotifsPhobiusNXS", "numMotifsPhobiusNXT", "numMotifsPhobiusNXC", "numMotifsPhobiusNXV", "numMotifsPhobiusC", "numMotifsPhobiusK")]
    longer <- phobius %>% pivot_longer(cols = !numMissedCleavages, names_to = "motif", values_to = "count")
    longer$count[longer$count>0]<-1  #because xtabs sums when it makes crosstabs and I want a count
    xt <- xtabs(count~motif+numMissedCleavages, longer)
    rownames(xt)<-c("C", "K", "NXC", "NXS", "NXT", "NXV")
    xtable(xt)
  },
  include.rownames=TRUE, #Don't print rownames
  caption.placement="top",
  include.colnames=FALSE,
  add.to.row = list(pos = list(0),
  command = "<tr><th colspan=1><B>Phobius<B></th><th colspan='3'>Number of Missed Cleavages</th></tr> <tr> <th> Motif </th> <th> 0 MC</th> <th> 1 MC</th> <th> 2 MC</th> </tr>"
  ),striped=TRUE, bordered=TRUE, digits=0)
  
  
  ##### TMHMM Table ##### 

  output$pepTMHMM <-renderTable({
    req(length(protData()[,1])>0)
    filtered <- filter(pepData(), OKforMS != 0 )
    tmhmm <- filtered[c("numMissedCleavages", "numMotifsTMHMMNXS", "numMotifsTMHMMNXT", "numMotifsTMHMMNXC", "numMotifsTMHMMNXV", "numMotifsTMHMMC", "numMotifsTMHMMK")]
    longer <- tmhmm %>% pivot_longer(cols = !numMissedCleavages, names_to = "motif", values_to = "count")
    longer$count[longer$count>0]<-1  #because xtabs sums when it makes crosstabs and I want a count
    xt <- xtabs(count~motif+numMissedCleavages, longer)
    rownames(xt)<-c("C", "K", "NXC", "NXS", "NXT", "NXV")
    xtable(xt)
  },
    #size="footnotesize", #Change size; useful for bigger tables
    include.rownames=TRUE, #Don't print rownames
    caption.placement="top",
    include.colnames=FALSE,
    add.to.row = list(pos = list(0),
      command = "<tr><th colspan=1><B>TMHMM<B></th><th colspan='3'>Number of Missed Cleavages</th></tr><tr> <th> Motif </th> <th> 0 MC</th> <th> 1 MC</th> <th> 2 MC</th> </tr>"
  ),striped=TRUE, bordered=TRUE, digits=0)
  

  ##### Motifs Plot ##### 

  getPepList <- function(a, motif, method) {
    column_name <- paste0( "numMotifs", method, motif )
    r <- filter(pepData(), numMissedCleavages == 0 & OKforMS != 0 & !!as.name(column_name) > 0 ) 
    return(r[[1]])
  }

  output$pepPlot <-renderPlotly({
    req(length(protData()[,1])>0)
    acc = toupper(input$swissprtID)
    motifs <- c("NXS", "NXT", "NXC", "NXV", "C", "K")
    M<-matrix(ncol=4)
    for (motif in motifs) { 
      M<-rbind(M, peptabulate(getPepList(toupper(input$swissprtID), motif, ''), getPepList(toupper(input$swissprtID), motif, 'Phobius'), getPepList(toupper(input$swissprtID), motif, 'TMHMM'))) 
    }
    M<-M[-1,]
    both <- M[,1]
    phobius <- M[,2]
    tmhmm <- M[,3]
    neither <- M[,4]
    data <- data.frame(motifs, both, phobius, tmhmm, neither)
    fig <- plot_ly(data, x = ~motifs, y = ~neither, type = 'bar', name = 'Not Predicted EC', marker=list(color=c('#9987C8')))
    fig <- fig %>% add_trace(y = ~tmhmm, name = 'TMHMM only', marker = list(color=c('#ADF5D1')))
    fig <- fig %>% add_trace(y = ~phobius, name = 'Phobius only', marker=list(color=c('#5CEBA3')))
    fig <- fig %>% add_trace(y = ~both, name = 'Phobius and TMHMM', marker=list(color=c('#4DB07E')))
    fig <- fig %>% layout(yaxis = list(title = 'Number of Peptides Containing Motif'), barmode = 'stack', xaxis=list(categoryorder = "array", categoryarray = ~motifs))  
  })  
  
  peptabulate <- function(all,phob,tmhmm) {
    ec <- union(phob, tmhmm)
    both <- intersect(phob, tmhmm)
    fphob <- setdiff(phob, tmhmm)
    ftmhmm <- setdiff(tmhmm, phob)
    nonec <- setdiff(all, ec)
    return(c(length(both), length(fphob), length(ftmhmm), length(nonec)))
  }  
  
  output$txtWarning <- renderText({
    if( length(protData()[,1]) == 0 ){
      paste0("The accession, ", toupper(input$swissprtID), ", is not in our database")
    }
  })
  

  translateSPC <- function(s) {
    s<-gsub("BF","Surfy",s)
    s<-gsub("T","Town",s)
    s<-gsub("dC","da Cunha",s)
    s<-gsub("DR","Diaz Ramos",s)
  }
  


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
  

##################################################  
###########       Visualizations        ###########


output$protvista <- renderUI({
  req(length(protData()[,1])>0)  
  h3( tags$a(href=paste0('https://www.uniprot.org/uniprot/', toupper(input$swissprtID), '/protvista'), paste0("Protvista: ", toupper(input$swissprtID)) ) )
})


output$protter <- renderUI({
  req(length(protData()[,1])>0)  
  h3(  tags$a(href=paste0('http://wlab.ethz.ch/protter/#up=', toupper(input$swissprtID), '&tm=auto&mc=lightsalmon&lc=blue&tml=numcount&numbers&legend&n:signal%20peptide,fc:red,bc:red=UP.SIGNAL&n:disulfide%20bonds,s:box,fc:greenyellow,bc:greenyellow=UP.DISULFID&n:variants,s:diamond,fc:orange,bc:orange=UP.VARIANT&n:PTMs,s:box,fc:forestgreen,bc:forestgreen=UP.CARBOHYD,UP.MOD_RES&format=svg'), paste0("Protter: ", toupper(input$swissprtID)) ))
})

generateProtter <- function () {
  

  df <- pepData()[c("range", "numMissedCleavages", "OKforMS", "numMotifsPhobiusK", "numMotifsPhobiusC", "numMotifsPhobiusNXS", "numMotifsPhobiusNXV", "numMotifsPhobiusNXT", "numMotifsPhobiusNXC")]
  df<-filter(df, numMissedCleavages == 0 & OKforMS != 0)
  df$K <- ifelse(df$numMotifsPhobiusK >0, 1, 0)
  df$C <- ifelse(df$numMotifsPhobiusC >0, 1, 0)
  df$N <- ifelse( (df$numMotifsPhobiusNXS + df$numMotifsPhobiusNXT + df$numMotifsPhobiusNXV + df$numMotifsPhobiusNXC) >0, 1, 0)

  n <- c()
  k <- c()
  c <- c()
  nc <- c()
  nk <-c ()
  ck <-c ()
  nck <- c()
  zero <- c()
  
  
  for (i in 1:nrow(df)) {
    if( df[i,]$N==1 & df[i,]$K==0 & df[i,]$C == 0 ) {
      n <- c( n, df[i,]$range )
    } else if( df[i,]$N==0 & df[i,]$K==1 & df[i,]$C == 0 ) {
      k <- c( k, df[i,]$range )
    } else if( df[i,]$N==0 & df[i,]$K==0 & df[i,]$C == 1 ) {
      c <- c( c, df[i,]$range )
    } else if( df[i,]$N==1 & df[i,]$K==1 & df[i,]$C == 0 ) {
      nk <- c( nk, df[i,]$range )
    } else if( df[i,]$N==1 & df[i,]$K==0 & df[i,]$C == 1 ) {
      nc <- c( nc, df[i,]$range )
    } else if( df[i,]$N==0 & df[i,]$K==1 & df[i,]$C == 1 ) {
      ck <- c( ck, df[i,]$range )
    } else if( df[i,]$N==1 & df[i,]$K==1 & df[i,]$C == 1 ) {
      nck <- c( nck, df[i,]$range )
    } else {
      zero <- c( zero, df[i,]$range )
    }
  }
  
  r <- protData()[c("motifsLocPhobiusNXS", "motifsLocPhobiusNXT", "motifsLocPhobiusNXC", "motifsLocPhobiusNXV", "motifsLocPhobiusC", "motifsLocPhobiusK")]
  locs <- c()
  for( motif in c('motifsLocPhobiusNXS', 'motifsLocPhobiusNXT', 'motifsLocPhobiusNXC', 'motifsLocPhobiusNXV', 'motifsLocPhobiusC', 'motifsLocPhobiusK') ){
    loc <- strsplit(as.character(r[motif]), ',')[[1]]
    loc<-loc[!loc=='n/a']
    locs <- union(locs, loc)
  } 

  len <- protData()["Length"]
  zero <- paste0('&n:Not%20Captured%20or%20Not%20OK,bc:FFFFFF,fc:D0D0D0=1-', len)
  if( !is.null(n) ){
    n <- paste0('&n:N%20Capture,bc:507F84,fc:507F84,cc:FFFFFF=', gsub(' ', '', toString(n)))
  } else {
    n<-c('')
  }
  if( !is.null(c) ){
    c <- paste0('&n:C%20Capture,bc:F9DDA2,fc:F9DDA2=', gsub(' ', '', toString(c)))
  }
  # else {
  #   c<-''
  # }
  if( !is.null(k) ){
    k <- paste0('&n:K%20Capture,bc:D8D8ED,fc:D8D8ED=', gsub(' ', '', toString(k)))    
  }
  # else {
  #   k<-''
  # }
  if( !is.null(nc) ){
    nc <- paste0('&n:N%20or%20C%20Capture,bc:A99F85,fc:A99F85=', gsub(' ', '', toString(nc)))
  }
  # else {
  #   nc<-''
  # }
  if( !is.null(nk) ){
    nk <- paste0('&n:N%20or%20K%20Capture,bc:5E6184,fc:5E6184,cc:FFFFFF=', gsub(' ', '', toString(nk)))    
  }
  # else {
  #   nk<-''
  # }
  if( !is.null(ck) ){
    ck <- paste0('&n:C%20or%20K%20Capture,bc:B6D1A8,fc:B6D1A8=', gsub(' ', '', toString(ck)))
  }
  # else {
  #   ck<-''
  # }
  if( !is.null(nck) ){
    nck <- paste0('&n:N%20or%20C%20or%20K%20Capture,bc:9E7EB9,fc:9E7EB9,cc:FFFFFF=', gsub(' ', '', toString(nck)))
  }
  m <- paste0('&n:Consensus%20Motif,s:diamond=', gsub(' ', '', toString(locs)))
  # else {
  #   nck<-''
  # }
  sp <- '&n:disulfide%20bonds,s:box,fc:greenyellow,bc:greenyellow=UP.DISULFID&n:signal%20peptide,fc:salmon,cc:white,bc:salmon=UP.SIGNAL'
  s <- paste0("http://wlab.ethz.ch/protter/create?up=", toupper(input$swissprtID), "&tm=auto&mc=lightgoldenrodyellow&lc=blue&tml=none&numbers&cutAt=peptidecutter.Tryps&legend")
  s <- paste0(s, zero, n, c, k, nc, nk, nck, m, sp, '&format=svg')
  im <-GET(s, timeout(30))
  return(im)
}

output$image1 <- renderImage({
  req(length(protData()[,1])>0)  

  image <- generateProtter()
  writeBin(image$content, 'test.svg')
  
  list(src = 'test.svg',
       contentType = "image/svg+xml",
       width = "120%",
       alt = "This is alternate text")
  
}, deleteFile = TRUE)

output$txtVWarning <- renderText({
  if( length(protData()[,1]) == 0 ){
    paste0("The accession, ", toupper(input$swissprtID), ", is not in our database")
  }
})


##################################################  
###########       Batch  Lookup        ###########
##################################################  
  
  pepSQL <- function(pepOptions, cleavages) {
    sql <- 'SELECT Accession, ID'
    if( 'LocProtein' %in% pepOptions ) {
      sql <- paste0( sql, ', range') 
    }
    if( 'Missed' %in% pepOptions ) {
      sql <- paste0( sql, ', numMissedCleavages AS MissedCleavages') 
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
      sql <- paste0( sql, ', SigPepPhobius, ScorePhobius, SigPepSignalP, ScoreSignalP, SigPepPrediSI, ScorePrediSI, numSPpredictions') 
    }
    if( 'Topo' %in% protOptions ) {
      sql <- paste0( sql, ', numTMPhobius, numICPhobius, numECPhobius, numTMTMHMM, numICTMHMM, numECTMHMM, StringOutPhobius, StringOutTMHMM') 
    }
    if( 'ProtMotif' %in% protOptions ) {
      sql <- paste0( sql, ', numMotifsNXS, motifsLocNXS, numMotifsNXT, motifsLocNXT, numMotifsNXC, motifsLocNXC, numMotifsNXV, motifsLocNXV, numMotifsC, motifsLocC, numMotifsK, motifsLocK, numMotifsPhobiusNXS, motifsLocPhobiusNXS, numMotifsPhobiusNXT, motifsLocPhobiusNXT, numMotifsPhobiusNXC, motifsLocPhobiusNXC, numMotifsPhobiusNXV, motifsLocPhobiusNXV, numMotifsPhobiusC, motifsLocPhobiusC, numMotifsPhobiusK, motifsLocPhobiusK, numMotifsTMHMMNXS, motifsLocTMHMMNXS, numMotifsTMHMMNXT, motifsLocTMHMMNXT, numMotifsTMHMMNXC, motifsLocTMHMMNXC, numMotifsTMHMMNXV, motifsLocTMHMMNXV, numMotifsTMHMMC, motifsLocTMHMMC, numMotifsTMHMMK, motifsLocTMHMMK')
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
      paste0("PeptideLevelData-", Sys.Date(), ".csv")
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
  

  
  
#Prot level Option 2
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
    req(input$file2)
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
    req(input$file2)
    downloadButton("Batch_File_prot_download_tsv", ".tsv", class="download_this")
  })
  
  
  
  #Pep level Option 2
  # Downloadable csv of selected dataset ----
  output$Batch_File_pep_download_csv <- downloadHandler(
    filename = function() {
      paste0("PeptideLevelData-", Sys.Date(), ".csv")
    },
    content = function(filename) {
      write.csv(Batch_File_pep_output(), filename, row.names = FALSE)
    }
  )
  output$Batch_File_pep_dlbutton_csv <- renderUI({
    req(input$file2)
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
    req(input$file2)
    downloadButton("Batch_File_pep_download_tsv", ".tsv", class="download_this")
  })

}





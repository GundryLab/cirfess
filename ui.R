# SurfaceGenie_0.1/ui.R
library(shiny)
library(plotly)
library(DT)

shinyUI(navbarPage("", theme = "bootstrap.css",
  
  ##########    Home    ##########
  
  tabPanel(
    "CIRFESS",
    # JavaScript to make links to the "Contact Us" tab
    # Because the anchors are created dynamically, there is no way to know what to use for the 
    # href part of an <a> tag that wants to link to the "Contacts" tab.  Luckily, anchors keep
    # a value called "data-value" that keeps the name of the tab.  So we use JavaScript in the 
    # browser to find all the a tags, look for which one has the data-value set to "contact" 
    # and grab the dynamically generated anchor and click it.
    tags$head(tags$script(HTML('
      var fakeClick = function(tabName) {
        var dropdownList = document.getElementsByTagName("a");
        for (var i = 0; i < dropdownList.length; i++) {
          var link = dropdownList[i];
          if(link.getAttribute("data-value") == tabName) {
            link.click();
          };
        }
      };
      '))),
    
    tags$img(src="website_homepage.svg",  width="25%",align="right"),
    h4("Welcome to ", span(class ="text-success", "CIRFESS"),"!"),
     p(tags$i("Helping Mass Spectrometrists find peptides since 2020")),
#    p(tags$u("C"), "ompiled ", tags$u("I"), "nteractive ", tags$u("R"), "esource ", tags$u("f"), "or ", tags$u("E"), "xtracellular and ", tags$u("S"), "urface ", tags$u("S"), "tudies (CIRFESS) integrates multiple 
     p("Compiled Interactive Resource for Extracellular and Surface Studies (CIRFESS) integrates multiple 
        prediction strategies and annotations into a single 
        interface for interrogating the human proteome. 
        Results from CIRFESS reveal the cell surface 
        proteome theoretically detectable by current 
        approaches and highlights where current prediction 
        strategies provide concordant and discordant 
        information. Overall, CIRFESS is designed to aid the 
        selection of available proteomic strategies and inform 
        the development of new strategies to enhance 
        coverage of the cell surface and extracellular 
        proteome."),

    
#     p("SurfaceGenie was written in R and the web application was developed using the Shiny library. Source code and all reference lookup tables are publicly available ", a(href="https://github.com/GundryLab/SurfaceGenie", "at GitHub"), "."),
     br(), 
#     div(style="width:20%;display:block;margin-left:auto;margin-right:auto",
     div(style="width:25%;display:block;margin-right:auto",
             tags$script(type="text/javascript", id="clustrmaps", src="https://cdn.clustrmaps.com/map_v2.js?d=oOpnxVR26WnLiF16Exi-XrxGb3rCX9xwJ4nBPUTKf1E&cl=ffffff&w=a")     )
  ),

  
  ##########  Instructions ##########

#   tabPanel(
#     "Instructions",
#     div(
#       p(style="font-size: 17px", tags$i("Before you begin:")),
# #      p(em("Before you begin:")),
#       p("SurfaceGenie contains two separate,  though related,  tools – ", span(class ="text-success", tags$b("GenieScore Calculator")), " and ", span(class ="text-success", tags$b("SPC Score Lookup")), ".  It is strongly recommended that all users read the ", a(href="UserGuide.pdf", "User Guide"), " which contains step-by-step tutorials for both the GenieScore Calculator and the SPC Score Lookup tools. The User Guide comprehensively defines all of the features available in the SurfaceGenie web application including some background on the theory and calculations."),
#       p(style="font-size: 17px", tags$i("Conversion to Uniprot Accession IDs:")),
#       p("SurfaceGenie operates with Uniprot Accession IDs only. Bulk conversion of alternate IDs to Uniprot IDs can be performed using the ‘Retrieve/ID mapping tool’ available on the Uniprot website, found here. Note that conversion between IDs is not always one-to-one. Manual curation of the results from the ID mapping is advisable."),
#       p(style="font-size: 17px", tags$i("Species availability:")),
#       p("Currently, most functions on SurfaceGenie are available only for human, mouse, and rat data. Calculation of some GenieScore permutations do not require Accession numbers, and will work on any type input data (see User Guide for more information). If you have requests for additional species, please ", a(href="#Contact", "contact us", onclick = "fakeClick('Contact')"), "."),
#       p(style="font-size: 17px", tags$i("Example files:")),
#       p("Examples of files formatted correctly for the GenieScore Calculator and the SPC Score Lookup tools can be downloaded using the links below. For more information, please refer to the ", a(href="UserGuide.pdf", "User Guide"), "." ),
#       p(a(href="ExampleDataForSurfaceGenie.csv", "GenieScore Calculator example file")),
#       p(a(href="ExampleDataForSPCdownload.csv", "SPC Score Lookup example file")),
#       br(),
#       p(style="font-size: 17px", tags$i("GenieScore Quick Overview:")),
#       tags$img(src="GSC_instructions.png", width="800px"), #, align="right"),
#       p(style="font-size: 17px", tags$i("SPC Score Lookup Quick Overview:")),
#       tags$img(src="SSL_instructions.png", width="800px") #, align="right")
#     )
#     ),
  
  ##########  SurfaceGenie ##########
  
  tabPanel(
    "Single Protein Lookup",
    sidebarPanel(
      textInput("swissprtID", "Enter a Human UniProt Accession", placeholder = "P14384"),
      em("Note: if no data appears after entering an accession, it is not one of the 20,405 human accessions in our database"  ),
      br(),
      br(),
      h5(class="text-info", "Key:"),
      p("SPC -  Surface Protein Consensus Score from SurfaceGenie"),
      p("SPC_DB - SPC Databases that predict this protein is cell surface"),
      p("SigPepPhobius - Signal Peptide predicted by Phobius"),
      p("ScorePhobius - Phobius' Signal Peptide prediction score"),
      p("SigPepSignalP - Signal Peptide predicted by SignalP"),
      p("ScoreSignalP - SignalP's Signal Peptide prediction score"),
      p("SigPepPredisi - Signal Peptide predicted by Predisi"),
      p("ScorePredisi - Predisi's Signal Peptide prediction score"),
      p("numTMPhobius - Phobius' number of transmembrane domains"),
      p("numICPhobius - Phobius' number of intracellular domains"),
      p("numECPhobius - Phobius' number of extracellular domains"),
      p("numTMTMHMM - TMHMM's number of transmembrane domains"),
      p("numICTMHMM - TMHMM's number of intracellular domains"),
      p("numECTMHMM - TMHMM's number of extracellular domains"),
      p("StringOutPhobius - raw output from Phobius"),
      p("StringOutTMHMM - raw output from TMHMM")
      
  #     h5(class="text-info", "Data Input"),
  #     fileInput("file1", "Choose Input File", multiple =FALSE, 
  #               accept=c(".csv", ".tsv", ".txt", ".tab", ".xls", ".xlsx"), 
  #               buttonLabel = "Browse...", placeholder = "No file selected"),
  # 
  #     h5(class="text-info", "Species"),
  #     radioButtons(
  #       "species", NULL,
  #       choices = list(
  #         "Human",
  #         "Rat",
  #         "Mouse",
  #         "Other/Ignore"),
  #       selected = list("Human")
  #     ),      
  #     
  #     h1(),
  #     h5(class="text-info", "Scoring Options"),
  #     
  #       conditionalPanel(
  #       condition = "input.species!='Other/Ignore'",
  #       checkboxGroupInput(
  #         "scoring_opts", label=NULL,
  #         choiceNames = mapply(scores, images, FUN=function(score, imgloc) {
  #           tagList(
  #             score,
  #             tags$img(src=imgloc, width=75)
  #           )
  #         }, SIMPLIFY = FALSE, USE.NAMES = FALSE),
  #       choiceValues = list(
  #         "GS", "IsoGenie", "OmniGenie", "IsoOmniGenie")
  #       )
  #     ),
  #     conditionalPanel(
  #       condition = "input.species=='Other/Ignore'",
  #       checkboxGroupInput(
  #         "scoring_opts_o", label=NULL,
  #         choiceNames = mapply(scores[c(3,4)], images[c(3,4)], FUN=function(score, imgloc) {
  #           tagList(
  #             score,
  #             tags$img(src=imgloc, width=75)
  #           )
  #         }, SIMPLIFY = FALSE, USE.NAMES = FALSE),
  #         choiceValues = list(
  #         "OmniGenie", "IsoOmniGenie")
  #     )
  #     ),
  #     h1(),
  #     h5(class="text-info", "Processing Option"),
  #     checkboxGroupInput(
  #       "processing_opts", NULL,
  #       choiceNames = list(
  #         "Group samples"),
  #       choiceValues = list(
  #         "grouping")
  #     ),
  #     conditionalPanel(
  #       condition = "input.processing_opts.indexOf('smarker') > -1",
  #       h5(class="text-info", "Markers for Specific Sample"),
  #       textInput(
  #         "markersample", "Enter sample name:", placeholder="i.e. 'd00' or 'Group 1'"
  #       )
  #     ),
  #     conditionalPanel(
  #       condition = "input.processing_opts.indexOf('grouping') > -1",
  #       h5(class="text-info", "Sample Grouping"),
  #       selectInput("groupmethod", "Grouping method",
  #                   choices = c(
  #                     "Mean" = "ave",
  #                     "Median" = "med"),
  #                   selected = "ave"
  #       ),
  #       p("*Please see Sample Grouping section on the Home page for instructions 
  #         on how to enter grouping information."),
  #       sliderInput("numgroups", "Number of groups",
  #                   min=2, max=5, value=2, step=1, ticks=FALSE),
  #       textInput("group1", "Group 1", placeholder="Columns in Group 1"),
  #       textInput("group2", "Group 2", placeholder="Columns in Group 2"),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 3",
  #         textInput("group3", "Group 3", placeholder="Columns in Group 3")
  #       ),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 4",
  #         textInput("group4", "Group 4", placeholder="Columns in Group 4")
  #       ),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 5",
  #         textInput("group5", "Group 5", placeholder="Columns in Group 5")
  #       ),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 6",
  #         textInput("group6", "Group 6", placeholder="Columns in Group 6")
  #       ),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 7",
  #         textInput("group7", "Group 7", placeholder="Columns in Group 7")
  #       ),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 8",
  #         textInput("group8", "Group 8", placeholder="Columns in Group 8")
  #       ),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 9",
  #         textInput("group9", "Group 9", placeholder="Columns in Group 9")
  #       ),
  #       conditionalPanel(
  #         condition = "input.numgroups >= 10",
  #         textInput("group10", "Group 10", placeholder="Columns in Group 10")
  #       )
  #     ),
  # 
  #     h1(),
  #     h5(class="text-info", "Export Options (CSV Download Tab)"),
  #     conditionalPanel(
  #       condition = "input.species!='Other/Ignore'",
  #       checkboxGroupInput(
  #       'export_options1', "SurfaceGenie Components:",
  #       choiceNames = list(
  #         "SPC score (SPC)",
  #         "Gini coefficient (Gini)",
  #         "Signal strength (SS)"
  #       ),
  #       choiceValues = list(
  #         "SPC", "Gini", "SS")
  #     )
  #     ),
  #     conditionalPanel(
  #       condition = "input.species=='Other/Ignore'",
  #       checkboxGroupInput(
  #         'export_options1o', "SurfaceGenie Components:",
  #         choiceNames = list(
  #           "Gini coefficient (Gini)",
  #           "Signal strength (SS)"
  #         ),
  #         choiceValues = list(
  #           "Gini", "SS")
  #       )
  #     ),
  #     conditionalPanel(
  #       condition = "input.species=='Human'",
  #       checkboxGroupInput(
  #       'export_options2h', "Annotations / Link outs:",
  #       choiceNames = list(
  #            "HLA molecules",
  #            "CD molecules",
  #            "Gene Name",
  #            "Number of CSPA experiments",
  #            "Transmembrane",
  #            "Subcellular Location",
  #            "UniProt Linkout"),
  #       choiceValues = list(
  #         "HLA", "CD", "geneName", "CSPA..e", "Transmembrane", "CC", "UniProt Linkout")
  #       )
  # ),
  # conditionalPanel(
  #   condition = "input.species=='Rat'",
  #   checkboxGroupInput(
  #     'export_options2r', "Annotations / Link outs:",
  #     choiceNames = list(
  #       "CD molecules",
  #       "Gene Name",
  #       "Transmembrane",
  #       "Subcellular Location",
  #       "UniProt Linkout"),
  #     choiceValues = list(
  #        "CD", "geneName", "Transmembrane", "CC", "UniProt Linkout")
  #   )
  # ),
  # conditionalPanel(
  #   condition = "input.species=='Mouse'",
  #   checkboxGroupInput(
  #     'export_options2m', "Annotations / Link outs:",
  #     choiceNames = list(
  #       "CD molecules",
  #       "Gene Name",
  #       "Number of CSPA experiments",
  #       "Transmembrane",
  #       "Subcellular Location",
  #       "UniProt Linkout"),
  #     choiceValues = list(
  #       "CD", "geneName", "CSPA..e", "Transmembrane", "CC", "UniProt Linkout")
  #   )
  # )
  # 
    ),
    
  
  ##########    Single Protein Lookup   ##########
  
    mainPanel(
      h5(class="text-info", "SPC Information"),
      tableOutput('SPC'),
      h5(class="text-info", "Signal  Peptide Predictions"),
      tableOutput('SigPep'),
      h5(class="text-info", "Topology Predictions"),
      tableOutput('Topo'),
      h5(class="text-info", "Peptide Summary"),
      tableOutput('PepSummary')
#      h5(class="text-info", "Motif Summary"),
#      tableOutput('MotifSummary')
      
    #   span(textOutput("txtWarning"), style="color:red"),
     )
  ),
  
  ##########  Surface Protein Concensus (SPC) Score  Lookup ##########
  
  tabPanel(
    "Batch Retreival",
    sidebarPanel(
      h5(class="text-info", "Input Option 1"),
      textAreaInput("quicklookup", "Uniprot accession number:", 
                placeholder="Enter accession numbers, each on a new line. For example:
                                                              
                                        A0AVT1-1 
                                        A0FGR8-6 
                                        A1L0T0 
                                        A1X283",
                rows=10),
      br(),
      h5(class="text-info", "Input Option 2"),
      fileInput("file2", "Choose Input File", multiple =FALSE, 
                accept=c(".csv", ".tsv", ".txt", ".tab", ".xls", ".xlsx"),
                buttonLabel = "Browse...", placeholder = "No file selected"),
      h1(),
       h5(class="text-info", "Protein Export Options"),
        checkboxInput("ProtLevelExp", "Protein Level Export", value = TRUE, width = NULL),
        conditionalPanel(
        condition = "input.ProtLevelExp",
        checkboxGroupInput(
          'export_prot_level', "Annotations/LinkOuts:",
          choiceNames = list(
            "SPC score (SPC)",
            "Signal Peptide Info",
            "Topology Info",
            "Motif Info",
            "Peptide Info"
          ),
          choiceValues = list(
            "SPC", "Signal", "Topo", "ProtMotif", "Pep")
        )
      ),
      h5(class="text-info", "Peptide Export Options"),
      checkboxInput("PepLevelExp", "Peptide Level Export", value = TRUE, width = NULL),
      conditionalPanel(
        condition = "input.PepLevelExp",
        checkboxGroupInput(
          'export_pep_level', "Annotations/LinkOuts:",
          choiceNames = list(
            "Location in Protein",
            "Missed Cleavages",
            "Ok for MS",
            "Number of Motifs",
            "Motif Locations"
          ),
          choiceValues = list(
            "LocProtein", "Missed", "OK", "NumMotifs", "LocMotifs")
        )
      )
      
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Instructions",
                    h4("Data Upload Instructions"),
                    h5(class="text-info", "Input Option 1"),
                    p("Enter a UniProt accession number(s) for your protein(s) of interest (e.g. P14384) into the textbox on the left, one per line."),
                    p("If your data are in a form other than UniProt (e.g. ENSEMBL gene, UniGene), a conversion 
                      tool is available", a(href="https://www.uniprot.org/uploadlists/", "here"), "Under 'Select options', 
                      select your ID type in the 'From' field and then 'UniProt KB' in the 'To' field. "),
                    h5(class="text-info", "Bulk Lookup"),
                    p("Upload a file containing a single column of UniProt accession numbers. "),
                    p("Bulk conversion from a different protein ID type to UniProt is available ",
                      a(href="https://www.uniprot.org/uploadlists/", "here"), 
                      ". Under 'Select options', select your ID type in the 'From' field and then 'UniProt KB'
                      in the 'To' field."),
                    p("With this method, the original upload file will be returned as a downloadable csv file which includes a column containing SPC Scores appended to the original input file.")
                  ),
                  tabPanel(
                    "Output Option 1",
                    h4(class="text-info", "Protein Data"),
                    uiOutput("Batch_Text_prot_dlbutton_csv", class='download_this'), uiOutput("Batch_Text_prot_dlbutton_tsv", class='download_this'), uiOutput("Batch_Text_prot_dlbutton_xlsx", class='download_this'), 
                    dataTableOutput("Batch_Text_prot_output"),
                    br(),

                    h4(class="text-info", "Peptide Data"),
                    uiOutput("Batch_Text_pep_dlbutton_csv", class='download_this'), uiOutput("Batch_Text_pep_dlbutton_tsv", class='download_this'), uiOutput("Batch_Text_pep_dlbutton_xlsx", class='download_this'),
                    dataTableOutput("Batch_Text_pep_output")
                  ),
                  tabPanel(
                    "Output Option 2",
                    h4(class="text-info", "Protein Data"),
                    uiOutput("Batch_File_prot_dlbutton_csv", class='download_this'), uiOutput("Batch_File_prot_dlbutton_tsv", class='download_this'), uiOutput("Batch_File_prot_dlbutton_xlsx", class='download_this'), 
                    dataTableOutput("Batch_File_prot_output"),
                    br(),
                    
                    h4(class="text-info", "Peptide Data"),
                    uiOutput("Batch_File_pep_dlbutton_csv", class='download_this'), uiOutput("Batch_File_pep_dlbutton_tsv", class='download_this'), uiOutput("Batch_File_pep_dlbutton_xlsx", class='download_this'),
                    dataTableOutput("Batch_File_pep_output")
                  )
      )
    )
  ),
  
##########  Reverse Lookup ##########

tabPanel(
  "Reverse Lookup",
  sidebarPanel(
    fluidRow(
      splitLayout(cellWidths = c("40%", "20%", "40%"),
        selectInput("inputA", "Property", choices= c("Length", "M/Z", "numMotif"), selected = NULL, multiple = FALSE, selectize = FALSE, width = NULL, size = NULL),
#        textInput("inputA", "Property"),
        selectInput("inputB", "Comparator", choices= c(">", "=", "<", "!="), selected = NULL, multiple = FALSE, selectize = FALSE, width = NULL, size = NULL),
#        textInput("inputB", "Comparator"),
        textInput("inputC", "Value")
      )
    )
  ),
  mainPanel(
    em("Under Construction")
  )
),

##########    References   ##########
  
  tabPanel(
    "References",
    div(
      h4("How to reference ", span(class ="text-success", "Cirfess") ),
      p("If you use CIRFESS in your work, please cite the original manuscript:"),
      p("Waas M, Littrell J, Gundry RL, Construction of an interactive resource for informing 
extracellular and surface proteomic studies (CIRFESS), submitted.") 
#        tags$a(href="https://doi.org/10.1101/575969", "https://doi.org/10.1101/575969"))
    ),
    br(),
    div(
      h4("Publications that cite ", span(class ="text-success", "CIRFESS") ),
      p("Coming Soon!")
    ),
    br(),
    div(
      h4("Publications that support ", span(class ="text-success", "CIRFESS") ),
      tags$ol(
        tags$li( "http://www.predisi.de/"),
        tags$li( "http://www.cbs.dtu.dk/services/SignalP/"),
        tags$li( "http://www.cbs.dtu.dk/services/TMHMM/" ),
        tags$li( "http://www.cellsurfer.net/surfacegenie")
      )
    ),
    br()
#    div(
#      h4("Users:"),
#      tags$script(type="text/javascript", id="clustrmaps", src="https://cdn.clustrmaps.com/map_v2.js?d=VJztTvZJUQlwpFCwOOYTSK6ktP0YBoNDEMPj1OS_ID0&cl=ffffff&w=a")
#      )
    ),

  ##########    Contact   ##########

  tabPanel(
  "Contact",
  div(
    h4(span(class ="text-success", "Contact "), "Us!"),
    p("If you have questions or suggestions for additional features, please contact us by email:"),
    p(class="text-info", style="text-indent:1.5em", "rebekah.gundry at unmc.edu"),
    p("Additional cell surface-related information and tools can be found at our growing website:"),
    p(class="text-info", style="text-indent:1.5em", "www.cellsurfer.net")
  ),
  br()
)


  
  ##########    Footer   ##########

#  div(
#    br(), br(),
#    tags$em(p(style="font-size:12px", "Publication Info [Gundry Lab 2018]"))
# )
))

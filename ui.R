# SurfaceGenie_0.1/ui.R
library(shiny)
library(DT)
library(rsconnect)
library(plotly)

shinyUI(navbarPage("", theme = "bootstrap.css",


  ####################          Home               #######################
  
  tabPanel(
    "CIRFESS",
    # JavaScript to make links to the "Contact Us" tab
    # Because the anchors are created dynamically, there is no way to know what to use for the 
    # href part of an <a> tag that wants to link to the "Contacts" tab.  Luckily, anchors keep
    # a value called "data-value" that keeps the name of the tab.  So we use JavaScript in the 
    # browser to find all the a tags, look for which one has the data-value set to "contact" 
    # and grab the dynamically generated anchor and click it.

    # requirements for protvista.  js and css in www/, json files in www/data/
    tags$head(tags$script(src="protvista.js")),
    tags$head(tags$link(href="main.css", rel="stylesheet")),
    
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

    # JavaScript to make it so when the enter button is clicked,  the Find button runs (on the single
    # protein lookup page).  It uses jquery functionality built into shiny - that the $() notation.
    # So with every keypress we look to see if the keypress was the enter key (keycode 13) AND if 
    # the textbox (swissprtID) had the focus at that time.  If so, same as clicking the Find button.
    tags$head(tags$script(HTML('
        $(document).keyup(function(event) {
        if ($("#swissprtID").is(":focus") && (event.keyCode == 13)) {
          $("#go").click();
        }
      });
      '))),
    
    
    tags$img(src="Fig1_CIRFESS.png",  width="33%",align="right"),
    h4("Welcome to ", span(class ="text-success", "CIRFESS"),"!"),
#     p(tags$i("Helping Mass Spectrometrists find peptides since 2020")),
#    p(tags$u("C"), "ompiled ", tags$u("I"), "nteractive ", tags$u("R"), "esource ", tags$u("f"), "or ", tags$u("E"), "xtracellular and ", tags$u("S"), "urface ", tags$u("S"), "tudies (CIRFESS) integrates multiple 
     p("Compiled Interactive Resource for Extracellular and Surface Studies (CIRFESS) integrates multiple 
        prediction strategies and annotations into a single 
        interface for interrogating the human, mouse, and rat proteomes. 
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

  

      ####################  Single Protein Lookup #######################

  tabPanel(
    "Single Protein Lookup",
    sidebarPanel( width = 3,
      textInput("swissprtID", "Enter a Human UniProt Accession", placeholder = "P14384"),
      actionButton("go", "Find"),
      br(),
#      em("Note: if no data appear after clicking Find, the accession is not one of the 20,405 human accessions in our database"  ),
      br(),
      br(),
      h5(class="text-info", "Key:"),
      p("SPC -  Surface Protein Consensus Score from SurfaceGenie"),
      p("SPC_DB - SPC Databases that predict this protein is cell surface", br(),
        "- Surfy = Bausch-Fluck et al., PNAS, 2018", br(),
        "- Town = Town et all, PNAS, 2016", br(),
        "- daCunha = da Cunha et al., PNAS, 2009", br(),
        "- DiazRamoz = Diaz-Ramos et al., Immol Lett, 2011")
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Information",
      fluidRow(
      column(width=3, h5(class="text-info", "SPC Information"),
      tableOutput('SPC')),
      column(width=4, h5(class="text-info", "Signal  Peptide Predictions"),
      tableOutput('SigPep')),
      h5(class="text-info", "Topology Predictions"),
      column(width=5,tableOutput('Topo')) ),
      h5(class="text-info", "Peptides Predicted to be Present in Extracellular MS Experiment"),
      h6(class="text-info", "These Peptides have consensus motifs, m/z is ok for MS, and predicted extracellular by one or more prediction methods"),
      
      dataTableOutput('PepSummary'),
      fluidRow(      column( width=5, tableOutput('pepPhobius') ),
        column( width=5, offset=1, tableOutput('pepTMHMM') )),
        textOutput('txtWarning'),
        plotlyOutput('pepPlot')
      ),
tabPanel(
  "Visualizations",
  uiOutput('protvista'),
  div(id="pvDiv"),
  textOutput('txtVWarning'),
  uiOutput('protter'),
  imageOutput("image1"),
  
  # This is the javascript to run protvista on the single protein lookup visualization page.  The
  # The full protvista js library is in a file in the www directory.  There is an associated css
  # file called main.css that goes with it.  
  tags$script("go.onclick = function() {
              var uniID = swissprtID.value;
              var uuniID = uniID.toUpperCase();
              var pvDiv = document.getElementById('pvDiv');
              var ProtVista = require('ProtVista');
              var instance = new ProtVista({
                el: pvDiv,
                uniprotacc: uuniID,
                defaultSources: false,
                customDataSource: {
                  url: './data/',
                  useExtension: true
                },
                overwritePredictions: true,
              });
          }"        
  )
  
)
)
    )
  ),
  


#   span(textOutput("txtWarning"), style="color:red"),

  tabPanel(
    "Batch Retrieval",
    sidebarPanel(
      h5(class="text-info", "Input Option 1"),
      textAreaInput("quicklookup", "Uniprot accession number:", 
        placeholder="Enter accession numbers, each on a new line. For example:
                                                              
                                        P23508 
                                        Q6ZV80 
                                        Q9NX46 
                                        P41240",
        rows=10),
      br(),

      h5(class="text-info", "Input Option 2"),
      fileInput("file2", "Choose Input File", multiple =FALSE, 
                accept=c(".csv", ".tsv", ".txt", ".tab", ".xls", ".xlsx"),
                buttonLabel = "Browse...", placeholder = "No file selected"),
      h1(),
       h5(class="text-info", "Select information to include in the output"),
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
                    h5(class="text-info", "Input Option 2"),
                    p("Upload a file containing a single column of UniProt accession numbers. "),
                    p("Bulk conversion from a different protein ID type to UniProt is available ",
                      a(href="https://www.uniprot.org/uploadlists/", "here"), 
                      ". Under 'Select options', select your ID type in the 'From' field and then 'UniProt KB'
                      in the 'To' field."),
                    br(),
                    p("The results will  be displayed on the screen and can be downloaded as either csv or tsv files."),
                    br(),
                    h4("Legends"),
                    p("There are two sets of results.  Each ", em("protein"), " will have data associated with it - signal peptide predictions, 
                      number of NXS motifs, SPC score, etc. Similiarly, each ", em("peptide"), " from those proteins will also have data assoctiated 
                      with them - length, location of each consensus motif, phobius and TMHMM predictions, etc.  These files will explain the data 
                      in the fields."),
                    a(href="Legend_CIRFESS_ProteinLevel.pdf", "Proteins"),
                    br(),
                    a(href="legend_CIRFESS_peptidelevel.pdf", "Peptides")
                    
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

# tabPanel(
#   "Reverse Lookup",
#   sidebarPanel(
#     fluidRow(
#       splitLayout(cellWidths = c("40%", "20%", "40%"),
#         selectInput("inputA", "Property", choices= c("Length", "M/Z", "numMotif"), selected = NULL, multiple = FALSE, selectize = FALSE, width = NULL, size = NULL),
# #        textInput("inputA", "Property"),
#         selectInput("inputB", "Comparator", choices= c(">", "=", "<", "!="), selected = NULL, multiple = FALSE, selectize = FALSE, width = NULL, size = NULL),
# #        textInput("inputB", "Comparator"),
#         textInput("inputC", "Value")
#       )
#     )
#   ),
#   mainPanel(
#     em("Under Construction")
#   )
#),

##########    References   ##########
  
  tabPanel(
    "References",
    div(
      h4("How to reference ", span(class ="text-success", "CIRFESS") ),
      p("If you use CIRFESS in your work, please cite the original manuscript:"),
      p("Waas M, Littrell J, Gundry RL, CIRFESS: An interactive resource for querying the set of 
        theoretically detectable peptides for cell surface and extracellular enrichment proteomic studies, ", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/32212654/", "https://pubmed.ncbi.nlm.nih.gov/32212654/.")) 
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
        tags$li( "http://www.cellsurfer.net/surfacegenie"),
        tags$li ("Protter: interactive protein feature visualization and integration with experimental proteomic data.
        Omasits U, Ahrens CH, Müller S, Wollscheid B.
        Bioinformatics. 2014 Mar 15;30(6):884-6."),
        #doi: 10.1093/bioinformatics/btt607")
        tags$li ("ProtVista: Visualization of Protein Sequence Annotations, Xavier Watkins, Leyla J Garcia, Sangya Pundir, Maria J Martin, UniProt Consortium. Bioinformatics. 2017 Jul 1;33(13):2040-2041.") 
        #                 (doi: 10.1093/bioinformatics/btx120) )
      )
    ),
    br(),
    div( tags$img(src="CIRFESS_COVER.jpg",  width="50%",align="center"))
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

))

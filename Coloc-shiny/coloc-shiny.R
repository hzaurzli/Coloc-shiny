library(shiny)

ui <- tagList(
  fluidPage(
    titlePanel("Coloc (post-GWAS and Colocalization)"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        # 对应后面的 output$file2
        uiOutput('file2'),
        
        actionButton('reset', 'RESET'),
        hr(),
        sliderInput("prob", "Step 3: Post-probability",
                    min = 0, max = 1,
                    value = 0.95),
        hr(),
        downloadButton("downloadData", "Download"),
        hr(),
        h5('Developer:'),
        h6('Small runze (shiny app)'),
        br(),
        h5('Github: '),
        h6('https://github.com/hzaurzli (Small runze)'),
        br(),
        h5('Cition:'),
        h6('https://github.com/chr1swallace/coloc')
      ),
      mainPanel(
        h4("Coloc table"),
        br(),
        br(),
        shinycssloaders::withSpinner(
          dataTableOutput("table")
        )
      )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  gwas <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T)
  })
  
  eqtl <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T)
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose GWAS csv",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "Step 2: Choose eQTL csv",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  

  output$table <- renderDataTable({
    library(dplyr)
    library(remotes)
    library(coloc)
    
    
    gwas <- gwas()
    eqtl <- eqtl()
    
    
    if(is.null(gwas) | is.null(eqtl)){
      warning("Please upload files!")
    } 
    else{
      input_data <- merge(eqtl, gwas, by="rs_id", all=FALSE, suffixes=c("_eqtl","_gwas"))
      colnames(input_data)[1] = 'snp'
      
      result <- coloc.abf(dataset1=list(pvalues=input_data$pval_nominal_gwas, type="cc", s=0.33, N=50000,snp = input_data$snp), 
                          dataset2=list(pvalues=input_data$pval_nominal_eqtl, type="quant", N=10000,snp=input_data$snp), 
                          MAF=input_data$maf)
      condition = as.numeric(input$prob)
      need_result = result$results
      need_result$SNP.PP.H4 = as.numeric(need_result$SNP.PP.H4)
      print(class(need_result$SNP.PP.H4))
      need_result = subset(need_result,need_result$SNP.PP.H4 > condition)
      need_result <<- need_result[,c(1,2,9,16,17)]
    }
  }, options = list(pageLength = 10))
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(need_result,file,row.names = T,quote = F)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)


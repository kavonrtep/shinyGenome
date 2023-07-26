#!/usr/bin/env Rscript
library(shiny)
library(colourpicker)
library(DT)
library(rtracklayer)

get_density <- function(x, chr_size=NULL, tw=1000000){
  cvg <- coverage(x)
  bins <- tileGenome(chr_size, tilewidth = tw)
  d <- binnedAverage(unlist(bins), cvg, "coverage")
  d
}
max_chr_length <-function(g){
  x <- split(g, ~seqnames)
  L <- sapply(x,function(x)max(end(x),na.rm=TRUE))
  L
}

numToHex <- function(num) {
  return(
    substr(rgb(num, num, num, maxColorValue = 1),
           2,3
    )
  )
}





# Define UI for data upload app ----
ui <- fluidPage(
  # App title ----
  titlePanel("Choose Genome Data:"),
  # Sidebar layout with input and output definitions
  fluidRow(
    column(2,
      wellPanel(
      # Input: Select a file ----
      fileInput("file1", "Choose Genome File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      p("Tab delimited file with columns: chromosome_name, chromosome size, centromere
      start, centromere end"),
      br(),

      # Button
      actionButton("submit", "Submit"),
      # Output: display error message
      verbatimTextOutput("error_message")
    ),
    wellPanel(
      # Input: Select a file ----
      fileInput("file2", "Choose Genome Annotation File",
                multiple = FALSE),
      # set the label of the file from text input
      textInput("file2_label", "Enter the label for data"),
      # color picker
      colourInput("file2_color", "Choose the color for data", value = "red"),
      # Button
      actionButton("submit_data", "Add Annotation Data"),
      # Output: display error message
      verbatimTextOutput("error_message2"),
      br(),
      br(),
      dataTableOutput("uploaded_files_table"),
      actionButton("submit_selected", "Show only selected data"),
    )
           ),
    column(10,
           mainPanel(
             plotOutput("hist", width = "1000px", height = "1000px")
           )
    )
  )
)


# server.R
# Define server logic to read selected file ----
server <- function(input, output) {

  uploaded_files_info <- reactiveValues(data = data.frame(file = character(),
                                                          color = character(),
                                                          label = character(),
                                                          stringsAsFactors = FALSE))
  uploaded_data <- reactiveValues(
    gff_list = list(),
    density_list = list()
  )

  chrom_sizes <- reactiveValues(
    sizes = NULL,
    ctr_start = NULL,
    ctr_end = NULL
  )

  df <- reactiveValues(
    df = NULL
  )


  observeEvent(input$submit, {
    req(input$file1)
    df$df <- read.table(input$file1$datapath, sep = "\t", header = FALSE)
    # check the number of columns in the data
    if (ncol(df$df) < 4) {
      output$error_message <- renderText({
        "Error: The uploaded file should have at least four columns."
      })
      return()
    } else {
      colnames(df$df)[1:4] <- c('seqname', 'length', 'ctr_start', 'ctr_end')
      output$hist <- renderPlot({
        # assume the histogram needs to be built from the second column
        plot_chromosomes(df$df)
      })
      uploaded_data$chrom_sizes <- df$df$length
      names(uploaded_data$chrom_sizes) <- df$df$seqname
      uploaded_data$ctr_start <- df$df$ctr_start
      uploaded_data$ctr_end <- df$df$ctr_end
    }
  })

  observeEvent(input$submit_data, {
    req(input$file2)
    req(input$file2_label)
    req(input$file2_color)
    req(input$file1)

    new_data <- data.frame(file = input$file2$name,
                           color = input$file2_color,
                           label = input$file2_label,
                           selected = TRUE,
                           stringsAsFactors = FALSE)
    uploaded_files_info$data <- rbind(uploaded_files_info$data, new_data)
    print(input$file2)
    print(input$file2$datapath)
    g <- import(input$file2$datapath)

    uploaded_data$gff_list[[input$file2_label]] <- g


    g <- keepSeqlevels(g, names(uploaded_data$chrom_sizes)[names(uploaded_data$chrom_sizes) %in% seqlevels(g)], pruning.mode='coarse')
    missing_seqlevels <- names(uploaded_data$chrom_sizes)[!names(uploaded_data$chrom_sizes) %in% seqlevels(g)]
    seqlevels(g) <- c(seqlevels(g), missing_seqlevels)
    d <- get_density(g, uploaded_data$chrom_sizes[seqlevels(g)])
    # remove zero coverage
    d <- d[d$coverage > 0,]
    dd <- data.frame(seqnames = seqnames(d), start = start(d),
                     end = end(d),
                     value = d$coverage,
                     stringsAsFactors = FALSE,
                     alpha = numToHex(d$coverage^(1/2)),
                     chr_index = match(seqnames(d),names(uploaded_data$chrom_sizes))
    )
    n <- length(uploaded_data$density_list) + 1
    uploaded_data$density_list[[n]] <- dd

    output$hist <- renderPlot({
        # assume the histogram needs to be built from the second column
      plot_chromosomes(df$df)
      print(uploaded_data$density_list)
      for (i in seq_along(uploaded_data$density_list)){

        if (!uploaded_files_info$data$selected[i]){
          next
        }
        col <- paste(uploaded_files_info$data$color[i], uploaded_data$density_list[[i]]$alpha, sep = "")
        rect(
          uploaded_data$density_list[[i]]$chr_index - 0.26,
          uploaded_data$density_list[[i]]$start,
          uploaded_data$density_list[[i]]$chr_index + 0.26,
          uploaded_data$density_list[[i]]$end,
          col = col,
          border = NA)
      }



    })


  })

  observeEvent(input$submit_selected,{
    selected <-  rep(FALSE, nrow(uploaded_files_info$data))
    selected[input$uploaded_files_table_rows_selected] <- TRUE
    uploaded_files_info$data$selected <- selected
    print(uploaded_files_info$data)
  })


  observeEvent(input$uploaded_files_table_cell_edit, {
    info <- input$uploaded_files_table_cell_edit
    str(info)
    i <- info$row
    j <- info$col
    v <- info$value
    uploaded_files_info$data[i, j] <<- DT::coerceValue(v, uploaded_files_info$data[i, j])
    print(uploaded_files_info$data)
  })

  output$uploaded_files_table <- renderDataTable({
    uploaded_files_info$data
  }, editable = TRUE)


}


plot_chromosomes <- function(g) {
  print('plotting chromosomes')
  N <- nrow(g)
  # Create a blank plot with appropriate limits and labels
  M <- max(g$length)
  plot(0, ylim = c(M,0), xlim = c(0, N), type = "n", xlab = "Chromosome Length", ylab
    = "", axes = FALSE)
  # plot chromosomes as rectangles
  W <- 0.5
  rect( (1:N) - W / 2,0, (1:N) + W / 2, g$length, col = "white", border = "grey")
  # add labels
  axis(3, at = 1:N, labels = g$seqname, tick = FALSE, las = 1)
  # add centromeres
  rect((1:N) - W / 2, g$ctr_start, (1:N) + W / 2, g$ctr_end,  col = "grey", border =
    "grey")

}

# Run the application
app <- shinyApp(ui = ui, server = server)
runApp(app, port = 3044)



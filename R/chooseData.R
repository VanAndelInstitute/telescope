#' STChooseData
#'
#' Select region of interest in ST spot matrix for display and analysis
#'
#' @param x an ST data object such as returned by STLoad
#' @param gene gene whose expression will be used to provide
#'             visual indication of where spots of interest are.
#' @import shiny
#' @import plotly
#' @export
STChooseData <- function(x, gene) {
    ui <- shiny::fluidPage(
      shiny::titlePanel("Select region of interest"),

      # Sidebar layout with input and output definitions ----
        shiny::mainPanel(width = 6,
          plotlyOutput("plot"),
          shiny::actionButton("done", "Done")
        )
      )

    server <- function(input, output, session) {
      output$plot <- renderPlotly({
        p <- suppressMessages(STPlotGene(x, genes=gene, size=2)) +
          theme(legend.position = "none")
        ggplotly(p, tooltip=NA) %>%
          layout(dragmode = "lasso") %>%
          config(displayModeBar = F)
      })

      shiny::observeEvent(input$done, {
        shiny::stopApp(jsonlite::fromJSON(input$`plotly_selected-A`))
      })

    }
    sel <- suppressMessages(shiny::runApp(shiny::shinyApp(ui, server)))
    return(STSubset(x, sel$x, sel$y))
}

#' STSubset
#'
#' Subset an ST data object base on x,y pairs
#'
#' @param d an ST data object such as returned by STLoad
#' @param x Vector of x coords for points to retain
#' @param y Vector of x coords for points to retain
#' @export
STSubset <- function(d, x, y) {
  crd.sel <- paste0(x, "x", y)
  crd.cur <- paste0(d$x, "x", d$y)
  ix <- which(crd.cur %in% crd.sel)
  d$x <- d$x[ix]
  d$y <- d$y[ix]
  d$exp <- d$exp[ix,]
  d
}

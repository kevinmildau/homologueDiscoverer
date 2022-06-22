#' plotPeakTableInteractive
#'
#' Opens s shiny web-application with an interactive plot depicting all the series in the provided series data base. Accepts the series
#'
#' @param seriesdb
#'
plotSeriesDBInteractive <- function(seriesdb, max_legend_size = 100) {
  seriesdb <- mutate(seriesdb, homologue_id = as.factor(homologue_id),
                     timestamp = as.character(timestamp))
  if (length(unique(seriesdb$homologue_id)) > max_legend_size) {
    legend_setting = "none"
  } else {
    legend_setting = "right"
  }
  colourCount = length(unique(seriesdb$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(seriesdb$homologue_id))
  seriesdb <- arrange(seriesdb, desc(homologue_id))
  ui <- fluidPage(
    plotOutput("plot", width = '100%', height = 600,
               dblclick = "plot_dblclick",
               brush = brushOpts(
                 id = "plot_brush",
                 resetOnNew = TRUE)),
    tableOutput("seriesdb")
  )

  # Set standard values for x and y limits of the plot if no reactive values provided
  xlim_lower <- min(seriesdb$rt) - 0.005 * median(seriesdb$rt)
  xlim_upper <- max(seriesdb$rt) + 0.005 * median(seriesdb$rt)
  ylim_lower <- min(seriesdb$mz) - 0.005 * median(seriesdb$mz)
  ylim_upper <- max(seriesdb$mz) + 0.005 * median(seriesdb$mz)

  server <- function(input, output) {
    ranges <- reactiveValues(x = c(xlim_lower, xlim_upper), y = c(ylim_lower, ylim_upper)) # defaults provided here
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$plot_dblclick, {
      brush <- input$plot_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- c(xlim_lower, xlim_upper) # defaults provided here
        ranges$y <- c(ylim_lower, ylim_upper) # defaults provided here
      }
    })
    output$seriesdb <- renderTable({
      brushedPoints(seriesdb, input$plot_brush)
    })

    output$plot <- renderPlot({
      ggplot(seriesdb, aes(group = homologue_id)) +
        geom_point(aes(x = rt, y = mz, color = homologue_id, size = homologue_id, shape = homologue_id),
                   alpha = 0.7) +
        geom_line(data = filter(seriesdb, homologue_id != 0), aes(x = rt, y = mz, group = homologue_id), alpha = 0.5) +
        ggtitle("Peak Table") +
        scale_colour_manual(values = c("black", getPalette(colourCount-1))) +
        scale_size_manual(values = c(rep(1.5, times = ncolor))) +
        scale_shape_manual(values = c(rep(19, times = ncolor))) +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
        theme(legend.position=legend_setting)
    })

  }
  shinyApp(ui, server)
}

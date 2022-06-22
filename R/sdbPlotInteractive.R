#' sdbPlotInteractive
#'
#' Functions creates interactive R Shiny server (local) for visualizing a seriesdb created using sdbCreate().
#'
#' Figure can be interacted with in the following ways:
#'
#' Left click and drag to select an area of the figure for which to show the corresponding data in the table below the figure.
#'
#' Left click and drag to select an area of the plot. Use double click to zoom into the area.
#'
#' Double click with no selection zooms out to original view.
#'
#' Select a polymer id using the input option. Click on the toggle plot button to change the plot type and table to highlight the selected homologue series.
#'
#' @param seriesdb seriesdb table as created by sdbCreate()
#' @param legend_setting A legend setting (position i.e. "bottom", "left", or default = "none")
#' @param plot_height Plot height in pixels. Defaults to 600. Plot always adjusts width to window size.
#'
#' @return None. An interactive plot is opened.
#' @export
sdbPlotInteractive <- function(seriesdb, legend_setting = "none", plot_height = 600) {
  # Set default values for x and y axis limits of the plot if no reactive values provided
  xlim_lower <- min(seriesdb$rt) - 0.005 * median(seriesdb$rt)
  xlim_upper <- max(seriesdb$rt) + 0.005 * median(seriesdb$rt)
  ylim_lower <- min(seriesdb$mz) - 0.005 * median(seriesdb$mz)
  ylim_upper <- max(seriesdb$mz) + 0.005 * median(seriesdb$mz)
  # Create vector of valid homologue_id choices for interactive selection
  choices <- sort(unique(seriesdb$homologue_id))
  # Modify timestamp column for shiny compatibility
  seriesdb <- mutate(seriesdb, timestamp = as.character(timestamp))
  # Define user interface elements
  ui <- fluidPage(
    plotOutput("plot", width = '100%', height = plot_height,
               dblclick = "plot_dblclick",
               brush = brushOpts(
                 id = "plot_brush",
                 resetOnNew = TRUE)),
    actionButton("button", "Toggle Highlight Single Homologue Series"),
    selectizeInput("selected_homologue", "Select Homologue:",
                   choices = choices, multiple = FALSE),
    tableOutput("seriesdb")
  )
  # Define active server
  server <- function(input, output, session) {
    # Set defaults for reactive ranges
    ranges <- reactiveValues(x = c(xlim_lower, xlim_upper),
                             y = c(ylim_lower, ylim_upper))
    # Define Double Click Reactive Event
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
    # Define Reactive table selection
    output$seriesdb <- renderTable({
      brushedPoints(seriesdb, input$plot_brush)
    })
    # Set default state for reactive plot type
    whichplot <- reactiveVal(TRUE)
    # Set homologue_selection observe event
    observeEvent(input$selected_homologue, {
      input$selected_homologue
    })
    # Set buttpn toggle for plot type change
    observeEvent(input$button, {
      whichplot(!whichplot())
    })
    which_graph <- reactive({
      if (whichplot()) {
        ggplotAllHomologuesSDB(seriesdb, {ranges$x}, {ranges$y}, legend_setting)
      } else {
        ggplotFocusHomologueSDB(seriesdb, {input$selected_homologue},
                                {ranges$x}, {ranges$y}, legend_setting)
      }
    })
    # Defining Table boundaries
    output$seriesdb <- renderTable({
      brushedPoints(seriesdb, input$plot_brush)
    })
    # Set table selection dependent on plot type (all homologues or focus homologue)
    which_table <- reactive({
      if (whichplot()){
        brushedPoints(seriesdb, input$plot_brush)
      } else {
        seriesdb %>% filter(., homologue_id == {input$selected_homologue})
        #brushedPoints(seriesdb %>% filter(., homologue_id == {input$selected_homologue}), input$plot_brush)
      }
    })
    output$seriesdb <- renderTable({
      which_table()
    })
    output$plot <- renderPlot({
      which_graph()
    })
  }
  shinyApp(ui, server)
}



#' ggplotAllHomologuesSDB
#'
#' Helper function that creates the ggplot for all homologue series in the seriesdb.
#'
#' @param seriesdb Series database as create by sdbCreate()
#' @param x Range vector for xlims.
#' @param y Range vector for ylims.
#' @param legend_setting Legend setting, defaults to "none" to avoid crowded legends.
#'
#' @return A ggplot object.
#' @keywords Internal
ggplotAllHomologuesSDB <- function(seriesdb, x, y, legend_setting){
  seriesdb <- mutate(seriesdb, homologue_id = as.factor(homologue_id))
  seriesdb <- mutate(seriesdb, timestamp = as.character(timestamp))
  colourCount = length(unique(seriesdb$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(seriesdb$homologue_id))
  seriesdb <- arrange(seriesdb, desc(homologue_id))
  plot <- ggplot(seriesdb, aes(group = homologue_id)) +
    geom_point(aes(x = rt, y = mz, color = homologue_id,
                   size = homologue_id, shape = homologue_id),
               alpha = 0.7) +
    geom_line(data = filter(seriesdb, homologue_id != 0),
              aes(x = rt, y = mz, group = homologue_id), alpha = 0.5) +
    ggtitle("Series DataBase - All Homologues") +
    scale_colour_manual(values = c("black", getPalette(colourCount-1))) +
    scale_size_manual(values = c(rep(1.5, times = ncolor))) +
    scale_shape_manual(values = c(rep(19, times = ncolor))) +
    coord_cartesian(xlim = x, ylim = y, expand = FALSE) +
    theme(legend.position=legend_setting)
  return(plot)
}



#' ggplotFocusHomologuesSDB
#'
#' Helper function that creates the focus homologue plot in sdbPlotInteractive.
#'
#' @param seriesdb Series database as create by sdbCreate()
#' @param selected_homologue_num homologue_id of the homologue to be highlighted.
#' @param x Range vector for xlims.
#' @param y Range vector for ylims.
#' @param legend_setting Legend setting, defaults to "none" to avoid crowded legends.
#'
#' @return A ggplot object.
#' @keywords Internal
ggplotFocusHomologueSDB <- function(seriesdb, selected_homologue_num, x, y, legend_setting){
  seriesdb <- mutate(seriesdb, homologue_id = as.factor(homologue_id))
  seriesdb <- mutate(seriesdb, timestamp = as.character(timestamp))
  seriesdb <- arrange(seriesdb, desc(homologue_id))
  selected_homologue_char <- paste("Homologue ID = ", selected_homologue_num)
  plot2db <- seriesdb %>%
    mutate(., selected_homologue =
             if_else(homologue_id == selected_homologue_num, 1L,0L),
           selected_homologue =
             factor(selected_homologue, levels = c(0L,1L),
                    labels = c("other", selected_homologue_char))) %>%
    arrange(., selected_homologue, homologue_id)
  plot <- ggplot(plot2db, aes(color = selected_homologue,
                              group = homologue_id,
                              lty = selected_homologue,
                              size = selected_homologue)) +
    geom_point(aes(x = rt, y = mz, shape = selected_homologue,
                   alpha = selected_homologue)) +
    geom_line(data =
                filter(plot2db, selected_homologue == selected_homologue_char),
              aes(x = rt, y = mz, group = homologue_id),
              alpha = 0.5, size = 0.4) +
    ggtitle("Series DataBase - Single Homologue Highlighted") +
    scale_colour_manual(values = c("black", "magenta")) +
    scale_size_manual(values = c(0.5, 4)) +
    scale_linetype_manual(values = c("blank", "solid")) +
    scale_shape_manual(values = c(1, 13)) +
    scale_alpha_manual(values = c(0.4, 0.75)) +
    coord_cartesian(xlim = x, ylim = y, expand = FALSE) +
    theme(legend.position=legend_setting)
  return(plot)
}




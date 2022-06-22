#' ptbPlotInteractive
#'
#' Functions creates interactive R Shiny server (local) for visualizing a peak_table created using detectHomologues().
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
#' @param peak_table A annotated peak table as produced by detectHomologues()
#' @param legend_setting A legend setting (position i.e. "bottom", "left", or default = "none")
#' @param plot_height Plot height in pixels. Defaults to 600. Plot always adjusts width to window size.
#'
#' @return None. An interactive plot is opened.
#' @export
ptbPlotInteractive <- function(peak_table, legend_setting = "none", plot_height = 600) {
  # Set default values for x and y axis limits of the plot if no reactive values provided
  xlim_lower <- min(peak_table$rt) - 0.005 * median(peak_table$rt)
  xlim_upper <- max(peak_table$rt) + 0.005 * median(peak_table$rt)
  ylim_lower <- min(peak_table$mz) - 0.005 * median(peak_table$mz)
  ylim_upper <- max(peak_table$mz) + 0.005 * median(peak_table$mz)
  # Create vector of valid homologue_id choices for interactive selection
  choices <- sort(unique(peak_table$homologue_id))

  # Define user interface elements #############################################
  ui <- fluidPage(
    plotOutput("plot", width = '100%', height = plot_height,
               dblclick = "plot_dblclick",
               brush = brushOpts(
                 id = "plot_brush",
                 resetOnNew = TRUE)),
    actionButton("button", "Toggle Highlight Single Homologue Series"),
    selectizeInput("selected_homologue", "Select Homologue:",
                   choices = choices, multiple = FALSE),
    tableOutput("peak_table")
  )

  # Define active server #######################################################
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
    output$peak_table <- renderTable({
      brushedPoints(peak_table, input$plot_brush)
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
        ggplotAllHomologuesPTB(peak_table, {ranges$x}, {ranges$y}, legend_setting)
      } else {
        ggplotFocusHomologuePTB(peak_table, {input$selected_homologue},
                                {ranges$x}, {ranges$y}, legend_setting)
      }
    })
    # Defining Table boundaries
    output$peak_table <- renderTable({
      brushedPoints(peak_table, input$plot_brush)
    })
    # Set table selection dependent on plot type (all homologues or focus homologue)
    which_table <- reactive({
      if (whichplot()){
        brushedPoints(peak_table, input$plot_brush)
      } else {
        peak_table %>% filter(., homologue_id == {input$selected_homologue})
        #brushedPoints(peak_table %>% filter(., homologue_id == {input$selected_homologue}), input$plot_brush)
      }
    })
    output$peak_table <- renderTable({
      which_table()
    })
    output$plot <- renderPlot({
      which_graph()
    })
  }
  shinyApp(ui, server)
}



#' ggplotAllHomologuesPTB
#'
#' Helper function that creates the ggplot for all homologue series in the peak_table.
#'
#' @param peak_table Series database as create by detectHomologues()
#' @param x Range vector for xlims.
#' @param y Range vector for ylims.
#' @param legend_setting Legend setting, defaults to "none" to avoid crowded legends.
#'
#' @return A ggplot object.
#' @keywords Internal
ggplotAllHomologuesPTB <- function(peak_table, x, y, legend_setting){
  peak_table <- mutate(peak_table,
                       homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                       homologue_id = as.factor(homologue_id))
  colourCount = length(unique(peak_table$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(peak_table$homologue_id))

  peak_table <- arrange(peak_table, desc(homologue_id))


  plot <- ggplot(peak_table, aes(group = homologue_id)) +
    geom_point(aes(x = rt, y = mz, color = homologue_id,
                   size = homologue_id, shape = homologue_id),
               alpha = 0.7) +
    geom_line(data = filter(peak_table, homologue_id != 0),
              aes(x = rt, y = mz, group = homologue_id), alpha = 0.5) +
    ggtitle("Annotated Peak Table - All Homologues") +
    scale_colour_manual(values = c("black", getPalette(colourCount-1))) +
    scale_size_manual(values = c(0.1, rep(4, times = ncolor))) +
    scale_shape_manual(values = c(1, rep(19, times = ncolor))) +
    coord_cartesian(xlim = x, ylim = y, expand = FALSE) +
    theme(legend.position=legend_setting)
  return(plot)
}



#' ggplotFocusHomologuePTB
#'
#' Helper function that generates a focus plot for a specific homologue series.
#'
#' @param peak_table Series database as create by detectHomologues()
#' @param selected_homologue_num homologue_id of the homologue to be highlighted.
#' @param x Range vector for xlims.
#' @param y Range vector for ylims.
#' @param legend_setting Legend setting, defaults to "none" to avoid crowded legends.
#'
#' @return A ggplot object.
#' @keywords Internal
ggplotFocusHomologuePTB <- function(peak_table, selected_homologue_num, x, y, legend_setting){
  peak_table <- mutate(peak_table,
                       homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                       homologue_id = as.factor(homologue_id))
  peak_table <- arrange(peak_table, desc(homologue_id))
  selected_homologue_char <- paste("Homologue ID = ", selected_homologue_num)
  plot2db <- peak_table %>%
    mutate(., selected_homologue =
             if_else(homologue_id == selected_homologue_num, 1L,0L),
           selected_homologue =
             factor(selected_homologue, levels = c(0L,1L),
                    labels = c("other homologues & non-homologue", selected_homologue_char))) %>%
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
    ggtitle("Annotated Peak Table - Single Homologue Series Highlighted") +
    scale_colour_manual(values = c("black", "magenta")) +
    scale_size_manual(values = c(0.5, 4)) +
    scale_linetype_manual(values = c("blank", "solid")) +
    scale_shape_manual(values = c(1, 13)) +
    scale_alpha_manual(values = c(0.4, 0.75)) +
    coord_cartesian(xlim = x, ylim = y, expand = FALSE) +
    theme(legend.position=legend_setting)
  return(plot)
}




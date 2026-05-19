#' Assembly Index
#'
#' Index and rank genome assemblies derived from the same set of reads using Shiny
#'
#' @details Different genome assemblers and purging programs can output assemblies of
#' varying qualities in terms of contiguity, gene completeness, duplication, and
#' missassemblies even when working with the same set of reads. Despite this, no
#' single metric adequately embodies a "high" quality genome. Here, we designed
#' a Shiny application that will index and rank assemblies by user-defined criteria
#' and optionally user-defined weights.
#'
#' The intuitive graphical user interface allows the user to upload a data frame
#' with an identification column and any number of genome attributes. The user can
#' select any combination of attributes that will either be feature normalized
#' such that either the attributes scores are maximized (i.e., attributes
#' wherein higher values are better) or minimized (i.e., attributes wherein lower values
#' are better). A sum from these normalized values will be calculated for each row
#' which is used to rank the assemblies. Each attribute can be differentially weighted
#' and the applied weights can be any number > 0.
#'
#' Rankings are visualized using a lollipop plot, where the top two assemblies are highlighted.
#' The normalized data and the generated plot can be exported.
#' @importFrom shiny fluidPage navbarPage tabPanel sidebarLayout sidebarPanel mainPanel fileInput checkboxInput radioButtons uiOutput numericInput actionButton plotOutput downloadButton renderPlot renderUI renderText eventReactive reactive req
#' @importFrom shinyWidgets pickerInput
#' @importFrom DT dataTableOutput renderDataTable datatable formatRound
#' @importFrom dplyr select rename mutate across all_of filter slice_max inner_join
#' @importFrom purrr reduce
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text scale_y_continuous ggtitle labs coord_flip theme ggplot element_text
#' @importFrom cowplot theme_minimal_vgrid panel_border
#' @export
asmidx<-function() {
  ui<-fluidPage(
    navbarPage(title = 'asmidx',
               tabPanel(title = 'Upload assembly attribute data',
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            fileInput(inputId = 'file', label = 'Upload file'),
                            checkboxInput(inputId = 'file_header', 'First row contains column names', value = T),
                            radioButtons(inputId = 'sep', 'Separator',
                                         choices= c(comma = ',',
                                                    semicolon = ';',
                                                    tab = '\t',
                                                    space = ' '),
                                         selected = ',')),
                          mainPanel(DT::dataTableOutput('uploaded_data'))
                          )
                        ),
               tabPanel(title = 'Choose columns',
                        sidebarLayout(
                          sidebarPanel(uiOutput('id_col'),
                                       uiOutput('max_col'),
                                       uiOutput('min_col'),
                                       uiOutput('assembly_size_out'),
                                       numericInput(inputId = 'known_size',
                                                    label = 'Known genome size (for calculating relative assembly size difference)',
                                                    value = NA),
                                       actionButton('select', 'Choose these columns')
                                       ),
                          mainPanel(DT::dataTableOutput('subset_data')
                                    )
                          )
                        ),
               tabPanel(title = 'Results: Rankings by Normalized Metrics',
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            div(style="display: inline-block; vertical-align:center", class = "row-fluid",
                                downloadButton(outputId = 'dl_norm_dat', 'Download Normalized Data'),
                                downloadButton(outputId = 'dl_norm_plot', 'Download Normalized Rankings Plot')
                                )
                            ),
                          mainPanel(
                            plotOutput('plot_norm'),
                            DT::dataTableOutput('norm_data')
                            )
                          )
                        ),
               tabPanel(title = 'Results: Rankings by User-specified Weights for Normalized Metrics',
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            div(style = "display: inline-block; vertical-align:center; horizontal-align:center", class = "row-fluid",
                                h5(strong('Weights can be any positive value')),
                                uiOutput('weight_ui'),
                                textOutput('w_sum'),
                                actionButton('w_select', 'Set Weights'),
                                downloadButton(outputId = 'dl_w_dat', 'Download Weighted Data and Weights'),
                                downloadButton(outputId = 'dl_w_plot', 'Download Weighted Rankings Plot')
                                )
                            ),
                          mainPanel(plotOutput('plot_weighted'),
                                    DT::dataTableOutput('weighted_data')
                                    )
                          )
                        )
               )
    )


  server<-function(input, output, session) {
    data<-reactive({
      dat<-input$file
      req(dat)
      read.csv(file = dat$datapath, sep = input$sep, header = input$file_header)
      })

    output$uploaded_data<-DT::renderDataTable(data(),
                                              server = T,
                                              options = list(scrollX = T)
                                              )

  output$id_col<-renderUI({
    pickerInput(inputId = 'id_selection',
                label = 'Choose ID column(s)',
                choices = colnames(data()),
                options = list(`action-box` = T,
                               maxOptions = 1),
                selected = NULL,
                multiple = T)
    })

  output$assembly_size_out<-renderUI({
    pickerInput(inputId = 'assembly_size_in',
                label = 'Select column that contains assembly size for calculating relative assembly size difference',
                choices = colnames(data()),
                options = list(`action-box` = T,
                               maxOptions = 1),
                selected = NULL,
                multiple = T)
    })

  output$max_col<-renderUI({
    pickerInput(inputId = 'max_selection',
                label = 'Choose columns where a higher value would be considered better',
                choices = colnames(data()),
                options = list(`action-box`= T),
                selected = NULL,
                multiple = T)
    })

  output$min_col<-renderUI({
    pickerInput(inputId = 'min_selection',
                label = 'Choose columns where a lower value would be considered better',
                choices = colnames(data()),
                options = list(`action-box` = T),
                selected = NULL,
                multiple = T)
    })

  size_reactive<-reactive({
    tryCatch({
      known_size<-input$known_size
      assembly_size_in<-input$assembly_size_in
      data.frame(data() %>% select(input$id_selection, all_of(assembly_size_in)), known_size = known_size) %>%
        rename(assembly_size = assembly_size_in) %>%
        mutate(r_size_diff = abs((assembly_size - known_size)/known_size),
               norm_r_size_diff = (r_size_diff - max(r_size_diff))/(min(r_size_diff) - max(r_size_diff))) %>%
        select(-known_size) %>%
        rename(assembly = input$id_selection)
      }, error = function(e) NA)
    })

  data_max_reactive<-reactive({
    tryCatch({
      max_selection<-input$max_selection
      data() %>%
        select(input$id_selection, all_of(max_selection)) %>%
        rename(assembly = input$id_selection) %>%
        mutate(across(all_of(max_selection), ~ (.x - min(.x))/(max(.x) - min(.x)), .names = 'norm_{.col}'))
      }, error = function(e) NA)
    })

  data_min_reactive<-reactive({
    tryCatch({
      min_selection<-input$min_selection
      data() %>%
        select(input$id_selection, all_of(min_selection)) %>%
        rename(assembly = input$id_selection) %>%
        mutate(across(all_of(min_selection), ~ (.x - max(.x))/(min(.x) - max(.x)), .names = 'norm_{.col}'))
      }, error = function(e) NA)
    })

  data_list<-reactive({
    data_list<-list(size_reactive(), data_max_reactive(), data_min_reactive())
    data_list<-Filter(Negate(anyNA), data_list)
    data_list<-data_list %>% reduce(inner_join, by = 'assembly')
    return(data_list)
    })

  data_norm<-eventReactive(input$select, {
    data_list() %>%
      select(assembly, contains('norm_')) %>%
      mutate(Score = rowMeans(pick(where(is.numeric))) * 100) %>%
      relocate(Score, .after = assembly)
    })

  data_raw<-eventReactive(input$select, {
    data_list() %>%
      select(assembly, !contains('norm_'))
    })

  top2<-reactive({
    top2<-data_norm() %>%
      slice_max(Score, n = 2)
    return(top2)
    })

  output$norm_data<-DT::renderDataTable({
    DT::datatable(data_norm()) %>%
      DT::formatRound(which(sapply(data_norm(), is.numeric)), digits = 3)
    },
    server = T,
    options = list(scrollX = T)
    )

  output$subset_data<-DT::renderDataTable(data_raw(),
                                          server = T,
                                          options = list(scrollX = T)
                                          )

  output$plot_norm<-renderPlot({
    data_norm() %>%
      ggplot(aes(x = assembly, y = Score)) +
      geom_segment(aes(x = assembly, xend = assembly, y = 0, yend = Score)) +
      geom_point(shape = 21, fill = 'red') +
      geom_point(data = top2(), shape = 21, fill = 'dodgerblue', size = 3) +
      geom_text(data = top2(), aes(x = assembly, y = Score + 8, label = round(Score, 2))) +
      # scale_y_continuous(expand = expansion(mult = c(0, 0.05)), ylim = c(0, 100)) +
      scale_y_continuous(expand = c(0, 0.05), limits = c(0, 100)) +
      ggtitle('Unweighted assembly ranking') +
      xlab('Assembly') +
      coord_flip() +
      theme_minimal_vgrid() +
      panel_border() +
      theme(aspect.ratio = 1)
    })

  #weighting scheme ui
  output$weight_ui<-renderUI({
    req(data_norm())
    var_name<-names(data_norm())[c(-1, -2)]
    var_name<-paste0('w_', var_name)
    lapply(var_name, function(x) {
      numericInput(inputId = x,
                   label = paste('Enter a weight for', x),
                   value = NA)
      }
      )
    })

  w_sum<-eventReactive(input$w_select, {
    weights<-sapply(paste0('w_', names(data_norm())[c(-1, -2)]), function(x) {
      input[[x]]
      }, simplify = F, USE.NAMES = T)
    weights_sum<-sum(unlist(weights))
    return(weights_sum)
    })


  data_w<-reactive({
    req(w_sum())
    weight_inputs<-sapply(names(data_norm())[c(-1, -2)], function(x) {
      input[[paste0('w_', x)]]
      }, simplify = T, USE.NAMES = T)
    weight_inputs<-weight_inputs/sum(weight_inputs)
    data_w<-data_norm()
    data_w[,c(-1, -2)]<-sweep(data_w[,c(-1, -2)], 2, weight_inputs, `*`)
    colnames(data_w)[-1]<-paste('w', colnames(data_w)[-1], sep = '_')
    data_w<-data_w %>%
      select(-w_Score) %>%
      mutate(wScore = rowSums(across(where(is.numeric))) * 100) %>%
      relocate(wScore, .after = assembly)
    return(data_w)
    })

  output$w_sum<-renderText({paste('Your weights sum to: ', w_sum(), '. This value will be used to divide each weight and the resultant value will be applied to the respective column.')})

  output$weighted_data<-DT::renderDataTable({
    DT::datatable(data_w()) %>%
      DT::formatRound(which(sapply(data_w(), is.numeric)), digits = 3)
    },
    server = T,
    options = list(scrollX = T))

    top2_w<-reactive({
      top2_w<-data_w() %>%
        slice_max(wScore, n = 2)
      return(top2_w)
      })

  output$plot_weighted<-renderPlot({
    data_w() %>%
      ggplot(aes(x = assembly, y = wScore)) +
      geom_segment(aes(x = assembly, xend = assembly, y = 0, yend = wScore)) +
      geom_point(shape = 21, fill = 'red') +
      geom_point(data = top2_w(), shape = 21, fill = 'dodgerblue', size = 3) +
      geom_text(data = top2_w(), aes(x = assembly, y = wScore + 8, label = round(wScore, 2))) +
      # scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
      scale_y_continuous(expand = c(0, 0.05), limits = c(0, 100)) +
      ggtitle('Weighted assembly ranking') +
      labs(x = 'Assembly', y = 'Weighted Score') +
      coord_flip() +
      theme_minimal_vgrid() +
      panel_border() +
      theme(aspect.ratio = 1)
    })

  output$dl_norm_dat<-downloadHandler(
    filename = function() {
      paste(Sys.Date(), '_normalized_metric_ranks', '.csv', sep = '')
      },
    content = function(file) {
      write.csv(data_norm(), file)
      })

  dl_plot_norm<-reactive({
    d<-data_norm()
    p<-ggplot(d, aes(x = assembly, y = Score)) +
      geom_segment(aes(x = assembly, xend = assembly, y = 0, yend = Score)) +
      geom_point(shape = 21, fill = 'red') +
      geom_point(data = top2(), shape = 21, fill = 'dodgerblue', size = 3) +
      geom_text(data = top2(), aes(x = assembly, y = Score + 8, label = round(Score, 2))) +
      # scale_y_continuous(expand = expansion(mult = c(0, 0.05)), ylim = c(0, 100)) +
      scale_y_continuous(expand = c(0, 0.05), limits = c(0, 100)) +
      ggtitle('Unweighted assembly ranking') +
      xlab('Assembly') +
      coord_flip() +
      theme_minimal_vgrid() +
      panel_border() +
      theme(aspect.ratio = 1,
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9))
    })

  output$dl_norm_plot<-downloadHandler(
    filename = function() {
      paste(Sys.Date(), '_normalized_metric_ranks', '.pdf', sep = '')
      },
    content = function(file) {
      ggsave(file, plot = dl_plot_norm(), device = 'pdf',
             width  = 7,
             height = 7,
             units = 'in')
      })


  output$dl_w_dat<-downloadHandler(
    filename = function() {
      paste(Sys.Date(), '_weighted_metric_ranks', '.csv', sep = '')
      },
    content = function(file) {
      x<-sapply(paste0('w_', names(data_norm())[c(-1, -2)]), function(x) {
        input[[x]]
        }, simplify = F, USE.NAMES = T)
      y<-paste(names(x), x, sep = '_')
      y<-c('assembly', 'wScore', y)
      z<-data_w() %>%
        rename_with(~ y)
      write.csv(z, file)
      })

  dl_plot_w<-reactive({
    d<-data_w()
    p<-ggplot(d, aes(x = assembly, y = wScore)) +
      geom_segment(aes(x = assembly, xend = assembly, y = 0, yend = wScore)) +
      geom_point(shape = 21, fill = 'red') +
      geom_point(data = top2_w(), shape = 21, fill = 'dodgerblue', size = 3) +
      geom_text(data = top2_w(), aes(x = assembly, y = wScore + 8, label = round(wScore, 2))) +
      # scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
      scale_y_continuous(expand = c(0, 0.05), limits = c(0, 100)) +
      ggtitle('Weighted assembly ranking') +
      labs(x = 'Assembly', y = 'Weighted Score') +
      coord_flip() +
      theme_minimal_vgrid() +
      panel_border() +
      theme(aspect.ratio = 1,
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9))
    })

  output$dl_w_plot<-downloadHandler(
    filename = function() {
      paste(Sys.Date(), '_weighted_normalized_metric_ranks', '.pdf', sep = '')
      },
    content = function(file) {
      ggsave(file, plot = dl_plot_w(), device = 'pdf',
             width  = 7,
             height = 7,
             units = 'in')
      }
    )}

  shinyApp(ui, server)
}

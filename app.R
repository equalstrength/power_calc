library(shiny)
library(tidyverse)
library(broom)
library(purrr)
library(plotly)
library(MASS)
library(shinythemes)
library(shinyWidgets)
library(bsplus)
library(InteractionPoweR)
library(htmltools)



# Define the UI ###############################################################

# UI settings ===================================================

ui <- fluidPage(theme = shinytheme("cosmo"),
                tags$head(tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
  titlePanel("EqualStrength Power Calculator"),
  sidebarLayout(
# Side panel simulation settings =================================
    sidebarPanel(
      h4("Simulation settings"),
      hr(),
      fluidRow(
        column(6, 
               numericInput("NumberSim",
                            "# of Simulations:",
                            2000, min = 1, step = 1000, max = 10000)),
        column(6, 
               numericInput("Alpha",
                            "Alpha:",
                            0.050, min = 0.001, step = 0.01, max = 0.100))),
      strong("Sample sizes"),
      fluidRow(
        column(6, 
               numericInput("SampMin",
                            "Min:",
                            200, min = 1, step = 100, max = 2000)),
        column(6, 
               numericInput("SampMax",
                            "Max:",
                            2000, min = 500, step = 100, max = 14000))),
# Side panel parameters =========================================
      h4("Parameters"),
      hr(),
      sliderInput("CallBack",
                  "Baseline call-back rate:",
                  value = 0.4,
                  min = 0.0, step = 0.05, max = 0.6,
                  round = FALSE),
      strong("Effect sizes"),
      numericInput("MainEffect", 
                  "Main treat. (X1):",
                  -0.2,
                  min = -2, step = 0.1, max = 2),
      numericInput("SecondEffect",
                  "Second treat. (X2):",
                  -0.2,
                  min = -2, step = 0.1, max = 2),
      numericInput("InteractEffect",
                  "Interaction:",
                  0.1,
                  step = 0.1, min = -2, max = 2),
# Side additional settings =========================================
      materialSwitch(inputId = "Additional", label = "Additional settings", 
                     value = FALSE, status = "primary"),
        conditionalPanel(condition = "input.Additional == true",
           numericInput("NumberSamp",
                        "# of Samples:",
                        50, min = 2, step = 10, max = 100)),
    actionButton("run_simulation", "Run simulation")
    ),
# Main panel  =====================================================
    mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Interaction", 
                 br(),
                 p("The plot below shows the expected power for the", strong("interaction"), 
                   "effect size indicated on the left according to each samples size"),
                 p("To show the plot or refresh with new parameters, click on the button 
                   'Run simulation' on the left"),
                 plotlyOutput("plot_interact")),
        tabPanel("Treatment", 
                 br(),
                 p("The plot below shows the expected power for the", strong("treatment"), 
                   "effect size indicated on the left according to each samples size"),
                 p("To show the plot or refresh with new parameters, click on the button 
                   'Run simulation' on the left"),
                 plotlyOutput("plot_treat"))
      )
    )
  )
)

# Define the Server ########################################################

server <- function(input, output) {
  observeEvent(input$run_simulation, {


# Read input values =====================================================
    in_method <- input$Method
    in_alpha <- input$Alpha
    in_n_treats <- input$NumberTreat
    in_base_cb <- input$CallBack
    in_treat1 <- input$MainEffect
    in_treat2 <- input$SecondEffect
    in_interact <- input$InteractEffect
    in_n_samples <- input$NumberSamp
    in_n_simul <- input$NumberSim

  # Sample sizes to include in the simulations
    range <- input$SampMax - input$SampMin
    samp_sizes <- seq(input$SampMin, input$SampMax, round(range/in_n_samples))
    in_reps <- in_n_simul/in_n_samples

    # Function for simulation ===============================================
    sim_data <- function(n, base_cb, treat1, treat2, interact){
      
      # a) Create the two variables, and assign treatment distribution (probabilities) in last part
      dat <- tibble(
        minority = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5)),
        second_treat = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
      )

      # b) Create the outcome based on probabilities
      dat$p <- base_cb +
        (treat1 * dat$minority) +
        (treat2 * dat$second_treat) +
        (interact * dat$minority * dat$second_treat)

      # c) To determine outcome: first create a random number between 0 and 1,
      # then, if that number is above "p", we assign a callback ("1")
      dat$random <- runif(nrow(dat), min = 0, max = 1)
      dat$outcome <- ifelse(dat$random < dat$p, 1, 0)

      # d) Estimate regression model, tidy it, save output with sample size.
      m1 <- lm(outcome ~ minority * second_treat, data = dat)  
      m1_out <- tidy(m1) 
      m1_out$n <- n

      return(m1_out)
    }

    # Define the parameters data frame
    params <- data.frame(
      expand.grid(
        n = samp_sizes,
        base_cb = in_base_cb,
        treat1 = in_treat1,
        treat2 = in_treat2,
        interact = in_interact)
    )

    # Duplicate the parameters data frame
    df_duplicated <- replicate(in_reps, params, simplify = FALSE)
    df_param <- do.call(rbind, df_duplicated)
    
    # Apply function to simulated data
    df_result <- purrr::pmap_dfr(df_param, sim_data, .progress = TRUE)
    
    # Define significance
    df_result$significant <- if_else(df_result$p.value < in_alpha, 1, 0)

    # Group and summarise the results
    df_grp <-
      df_result  |>
      filter(term %in% c("minority", "second_treat", "minority:second_treat")) |>
      group_by(n, term)  |>
      summarise(power = mean(significant), .groups = "drop") |> 
      mutate(n = n + (n/3)) # 33% larger sample for three ethnic groups


    # Plot the results
    output$plot_interact <- renderPlotly({
      p_int <- 
        df_grp |> 
        filter(term == "minority:second_treat") |> 
        ggplot(aes(y = power, x = n)) +
        geom_point(color = "#2C3E50") +
        geom_hline(yintercept = 0.8, color = "red", linetype = 2) +
        geom_smooth(method = "loess") +
        labs(
          x = "Sample size",
          y = "Power",
          title = "Power of interaction term across different sample sizes") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme_bw()

      print(p_int)
    })

    output$plot_treat <- renderPlotly({
      p_trt <- 
        df_grp |> 
        filter(term != "minority:second_treat") |> 
        ggplot(aes(y = power, x = n)) +
        geom_point(color = "#2C3E50") +
        geom_hline(yintercept = 0.8, color = "red", linetype = 2) +
        geom_smooth(method = "loess") +
        labs(
          x = "Sample size",
          y = "Power",
          title = "Power of treatments across different sample sizes") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        facet_wrap(~term, nrow = 1)+
        theme_bw()
      
      print(p_trt)
    })
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

library(shiny)
library(tidyverse)
library(broom)
library(purrr)
library(plotly)
library(MASS)

# Define the UI
ui <- fluidPage(
  titlePanel("EqualStrength Power Calculator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("callback",
            "Baseline call-back rate",
            min = 0.0,
            max = 0.6,
            step = 0.05,
            value = 0.45,
            round = FALSE),
      sliderInput("maineffect",
            "Main treatment effect (minority)",
            min = -0.2,
            max = 0.2,
            step = 0.05,
            value = -0.1),
      sliderInput("othereffect",
            "Other treatment effect",
            min = -0.2,
            max = 0.2,
            step = 0.05,
            value = 0),
      sliderInput("interacteffect",
            "Interaction effect",
            min = -0.2,
            max = 0.2,
            step = 0.05,
            value = 0.05),
      numericInput('number_reps', 'Number of repetitions', 20, min = 1, max = 100),
      actionButton("run_simulation", "Run simulation")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Interaction", plotlyOutput("plot_interact")),
        tabPanel("Treatment", plotlyOutput("plot_treat")),
        tabPanel(
          "About",
          br(),
          h4("What is this?"),
          p("This app calculates the expected power based on simulations with different sample sizes"),
          p("Power is estimated from the proportion of times in which the term was found significant (p < 0.05)"),
          h4("Precision"),
          p("The plot was generated after the following total number of simulations:"),
          verbatimTextOutput("about"),
          br(),
          p("The estimates are based on the following pre-determined specifications:"),
          p("- There are only two main treatment groups (ethnicity). For an additional group, we should expect a 33% larger sample for the same power"),
          p("- The main treatment and control group are of the same size (50% and 50%)"),
          br(),
          h4("Credits"),
          p("Code prepared by Edvard and adapted by the EqualStrength team ")
          )
      )
    )
  )
)

# Define the server
server <- function(input, output) {
  observeEvent(input$run_simulation, {
    # Read input values
    in_base_cb <- input$callback
    in_treat1 <- input$maineffect
    in_treat2 <- input$othereffect
    in_interact <- input$interacteffect
    in_reps <- input$number_reps

    # Create function to generate simulated data
    sim.data <- function(n, base_cb, treat1, treat2, interact){
      
      # b) Create an empty dataset
      dat <- data.frame(matrix(nrow = n, ncol = 0))
      
      # c) Create the two variables, and assign treatment distribution (probabilities) in last part
      dat$minority <- sample(c(0, 1), nrow(dat), replace = TRUE, prob = c(0.5, 0.5))
      dat$other_treat <- sample(c(0, 1), nrow(dat), replace = TRUE, prob = c(0.5, 0.5))
      
      # d) Create the outcome based on probabilities
      dat$p <- base_cb +
        (treat1 * dat$minority) +
        (treat2 * dat$other_treat) +
        (interact * dat$minority * dat$other_treat)
      
      # e) To determine outcome: first create a random number between 0 and 1,
      # then, if that number is above "p", we assign a callback ("1")
      
      dat$random <- runif(nrow(dat), min = 0, max = 1)
      dat$outcome <- ifelse(dat$random < dat$p, 1, 0)
      
      # f) Estimate regression model, tidy it, and fix names/save result + parameters.
      
      m1 <- lm(outcome ~ minority * other_treat, data = dat)  |>
        tidy() %>%
        mutate(n = n,
               baseline_cb = base_cb,
               minority_effect = treat1,
               other_treat_effect = treat2,
               interaction_effect = interact,
               significant = ifelse(.$p.value < 0.05, 1, 0))
      
      return(m1)
    }

    # Define the parameters data frame
    params <- data.frame(
      expand.grid(
        n = seq(2000, 12000, 100),
        base_cb = in_base_cb,
        treat1 = in_treat1,
        treat2 = in_treat2,
        interact = in_interact
      )
    )

    # Duplicate the parameters data frame
    df_duplicated <- replicate(in_reps, params, simplify = FALSE)
    df_param <- do.call(rbind, df_duplicated)

    # Apply simulation on each row of the parameter data frame
    df_result <- purrr::pmap_dfr(df_param, sim.data, .progress = TRUE)

    # Group and summarize the results
    df_grp_interact <- 
      df_result  |>
      filter(term == "minority:other_treat") |>
      dplyr::select(n, baseline_cb, minority_effect, interaction_effect, significant)  |>
      group_by(baseline_cb, minority_effect, interaction_effect, n)  |>
      summarise(power = mean(significant), .groups = "drop")
    
    df_grp_treat <- 
      df_result  |>
      filter(term %in% c("minority", "other_treat")) |>
      dplyr::select(n, baseline_cb, minority_effect, interaction_effect, significant, term)  |>
      group_by(baseline_cb, minority_effect, interaction_effect, n, term)  |>
      summarise(power = mean(significant), .groups = "drop")
    
    # Plot the results
    output$plot_interact <- renderPlotly({
      p_int <- ggplot(data = df_grp_interact, aes(y = power, x = n)) +
        geom_point(color = "darkblue") +
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
      p_trt <- ggplot(data = df_grp_treat, aes(y = power, x = n)) +
        geom_point(color = "darkblue") +
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
    # Information about number of repetitions

    output$about <- renderText({
      input$number_reps*nrow(params)
    })
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

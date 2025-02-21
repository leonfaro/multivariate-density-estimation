library(shiny)
library(copula)
library(benchden)
library(tram)
library(np)
library(ks)

dist_map <- c(
  "Normal (11)"        = "11",
  "Lognormal (12)"     = "12",
  "Logistic (5)"       = "5",
  "Cauchy (6)"         = "6",
  "Exponential (2)"    = "2",
  "Skewed Bimodal(22)" = "22",
  "Claw (23)"          = "23",
  "Beta(2,2) (17)"     = "17"
  
)


generate_data_copula_benchden <- function(n, rho, dist_indices) {
  d <- length(dist_indices)
  stopifnot(d %in% c(2,3))
  
  # Normal-Copula
  if (d == 2) {
    myCop <- normalCopula(param = rho, dim = 2)
  } else {
    corMat <- matrix(rho, 3, 3)
    diag(corMat) <- 1
    myCop <- normalCopula(param = corMat[upper.tri(corMat)], 
                          dim = 3, dispstr = "un")
  }
  
  set.seed(123)
  U <- rCopula(n, myCop)
  
  data_mat <- matrix(NA, nrow=n, ncol=d)
  for (j in seq_len(d)) {
    data_mat[, j] <- qberdev(U[, j], dnum = dist_indices[j])
  }
  df <- as.data.frame(data_mat)
  names(df) <- paste0("y", seq_len(d))
  df
}


fit_factorized_tram <- function(df) {
  d <- ncol(df)
  stopifnot(d %in% c(2,3))
  # p(y1)
  m1 <- BoxCox(y1 ~ 1, data = df)
  # p(y2|y1)
  m2 <- BoxCox(y2 ~ y1, data = df)
  if (d == 3) {
    m3 <- BoxCox(y3 ~ y1 + y2, data = df)
  } else {
    m3 <- NULL
  }
  list(m1=m1, m2=m2, m3=m3)
}


get_conditional_density <- function(model, newdata, grid_size=100) {
  response_var <- names(model$data)[1]
  y_range <- range(model$data[[response_var]])
  y_grid  <- seq(y_range[1], y_range[2], length.out=grid_size)
  
  newdata   <- as.data.frame(newdata)
  densities <- matrix(NA, nrow=nrow(newdata), ncol=length(y_grid))
  
  for (i in seq_len(nrow(newdata))) {
    extended_data <- newdata[rep(i, grid_size), , drop=FALSE]
    extended_data[[response_var]] <- y_grid
    dens_vals <- predict(model, newdata = extended_data, type="density")
    densities[i, ] <- dens_vals
  }
  rowMeans(densities)
}


evaluate_factorized_tram <- function(model_list, df) {
  d <- ncol(df)
  names(df) <- paste0("y", seq_len(d))
  
  f1_vals <- get_conditional_density(model_list$m1, df["y1"])
  f2_vals <- get_conditional_density(model_list$m2, df[c("y1")])
  
  if (d == 2) {
    joint_dens <- f1_vals * f2_vals
  } else {
    f3_vals <- get_conditional_density(model_list$m3, df[c("y1","y2")])
    joint_dens <- f1_vals * f2_vals * f3_vals
  }
  sum(log(joint_dens), na.rm=TRUE)
}


compare_kernel <- function(df) {
  d <- ncol(df)
  # np-Paket
  form_np <- as.formula(paste0("~", paste(names(df), collapse="+")))
  bw_np  <- npudensbw(form_np, data=df)
  np_est <- npudens(bw_np)
  ll_np  <- sum(log(fitted(np_est)), na.rm=TRUE)
  
  # ks-Paket
  H_ks    <- Hpi(x=as.matrix(df))
  ks_est  <- kde(x=as.matrix(df), H=H_ks)
  dens_ks <- predict(ks_est, x=as.matrix(df))
  ll_ks   <- sum(log(dens_ks), na.rm=TRUE)
  
  list(ll_np=ll_np, ll_ks=ll_ks)
}

performFactorizedTram <- function(n, rho, dist_indices) {
  df_data <- generate_data_copula_benchden(n, rho, dist_indices)
  mod_list <- fit_factorized_tram(df_data)
  ll_factor <- evaluate_factorized_tram(mod_list, df_data)
  
  # Kernel
  cmp <- compare_kernel(df_data)
  

  list(
    df = df_data,
    mod_list = mod_list,
    ll_factor = ll_factor,
    ll_np = cmp$ll_np,
    ll_ks = cmp$ll_ks
  )
}


ui <- fluidPage(
  titlePanel("Faktorisierte tram (BoxCox) mit NormalCopula + benchden"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        inputId = "dim",
        label   = "Dimension wählen",
        choices = c("2D"=2, "3D"=3),
        selected = 2
      ),
      selectizeInput(
        inputId = "dist_sel",
        label   = "Wähle Verteilungen (2 oder 3):",
        choices = dist_map,
        multiple=TRUE,
        selected = c("11", "12"), # default Normal, Lognormal
        options = list(maxItems = 3, # max 3
                       minItems = 2) # min 2
      ),
      numericInput("n", "Stichprobengröße (n)", value=100, min=10, step=50),
      numericInput("rho", "Korrelation (rho)", value=0.3, step=0.1),
      
      actionButton("applyBtn", "Apply / Run")
    ),
    
    mainPanel(
      # Text-Ausgabe (Summary)
      verbatimTextOutput("summary_text"),
      
      # Plots: 1) Scatter/pairs, 2) Barplot
      plotOutput("main_plot", height="600px")
    )
  )
)

server <- function(input, output, session) {
  
  myResults <- eventReactive(input$applyBtn, {
    req(input$dist_sel)
    # Dimension
    dim_val <- as.numeric(input$dim)

    
    chosen_dist <- as.numeric(input$dist_sel)
    
    if (dim_val == 2 && length(chosen_dist) > 2) {
      chosen_dist <- chosen_dist[1:2]
    }
    if (dim_val == 3 && length(chosen_dist) < 3) {

      
      chosen_dist <- c(chosen_dist, chosen_dist[1])[1:3]
    }
    

    
    performFactorizedTram(
      n           = input$n,
      rho         = input$rho,
      dist_indices= chosen_dist
    )
  })
  
 
  
  output$summary_text <- renderPrint({
    res <- myResults()
    cat("Summary der generierten Daten:\n")
    print(summary(res$df))
    cat("\nLogLik Factorized:", res$ll_factor,
        "\nLogLik np       :", res$ll_np,
        "\nLogLik ks       :", res$ll_ks, "\n")
  })
  

  
  output$main_plot <- renderPlot({
    res <- myResults()
    df  <- res$df
    d   <- ncol(df)
    
    # 2 Felder nebeneinander
    par(mfrow=c(1,2))
    
    # Plot links: Scatter oder pairs
    if (d == 2) {
      plot(df[,1], df[,2],
           xlab="Y1", ylab="Y2",
           main=sprintf("2D-Scatter (n=%d, rho=%.2f)", 
                        nrow(df), input$rho))
    } else {
      pairs(df, main=sprintf("3D-Pairs (n=%d, rho=%.2f)", 
                             nrow(df), input$rho))
    }
    
    # Plot rechts: Balkendiagramm Log-Likelihood
    barplot(
      c("Factorized"=res$ll_factor, 
        "np"=res$ll_np, 
        "ks"=res$ll_ks),
      main="LogLik Vergleich",
      ylab="Log-Likelihood",
      las=2
    )
    
  })
}


shinyApp(ui, server)

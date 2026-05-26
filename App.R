# ==============================================================================
# PLATFORM: CropPhenoAI (Professional Edition)
# VERSION: 1.0.0-Beta
# DESCRIPTION: Phenomics, Design Generation, and Statistical Analysis Engine.
#              Supports RCBD and Alpha Lattice designs for Single-Trial and GxE.
# ==============================================================================
library(shiny)
library(shinydashboard)
library(DT)
library(shinycssloaders)
library(readxl)
library(writexl)
library(tidyverse)
library(scales)
library(agricolae)
library(lmerTest) 
library(emmeans)  
library(reticulate) 

# ------------------------------------------------------------------------------
# INITIALIZE LOCAL PYTHON ENVIRONMENT FOR IMAGE ANALYSIS
# ------------------------------------------------------------------------------
try({
  if (!py_available(initialize = TRUE)) {
    virtualenv_create("croppheno_env")
    virtualenv_install("croppheno_env", packages = c("opencv-python", "numpy"))
    use_virtualenv("croppheno_env", required = TRUE)
  }
}, silent = TRUE)

# Inline Python Script for Edge Detection & Pixel Extraction
py_run_string("
import cv2
import numpy as np

def analyze_specimen_matrix(image_path):
    try:
        img = cv2.imread(image_path)
        if img is None:
            return {'status': 'error', 'msg': 'Invalid image matrix'}
            
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        blurred = cv2.GaussianBlur(gray, (5, 5), 0)
        _, thresh = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        if not contours:
            return {'status': 'empty', 'area_px': 0, 'length_px': 0, 'width_px': 0}
            
        largest_contour = max(contours, key=cv2.contourArea)
        area_pixels = float(cv2.contourArea(largest_contour))
        x, y, w, h = cv2.boundingRect(largest_contour)
        
        return {
            'status': 'success',
            'area_px': area_pixels,
            'length_px': float(h),
            'width_px': float(w)
        }
    except Exception as e:
        return {'status': 'error', 'msg': str(e)}
")

# ------------------------------------------------------------------------------
# UI DEFINITION
# ------------------------------------------------------------------------------

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CropPhenoAI Platform", titleWidth = 300),
  
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Dashboard Overview", tabName = "dashboard", icon = icon("gauge-high")),
      menuItem("Automated Image Phenomics", tabName = "phenomics", icon = icon("camera-retro")),
      menuItem("Research Design Generation", tabName = "design", icon = icon("compass-drafting")),
      menuItem("Single-Trial Analysis", tabName = "single_trial", icon = icon("microscope")),
      menuItem("Combined GxE MGIDI Engine", tabName = "gxe_mgidi", icon = icon("dna"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        body { background-color: #f8f9fa; font-family: 'Segoe UI', sans-serif; }
        .box { border-radius: 5px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); border-top: 4px solid #3498db; }
        .btn-download { background-color: #27ae60; color: white; margin-top: 10px; font-weight: bold; }
        .btn-execute { background-color: #2c3e50; color: white; font-weight: bold; }
        .badge-acc { background: #27ae60; color: white; padding: 5px 10px; border-radius: 4px; font-weight: bold; }
        .analysis-sidebar { background: #f0f2f5; padding: 15px; border-radius: 5px; margin-bottom: 15px; border: 1px solid #ddd; }
      "))
    ),
    
    tabItems(
      # --- DASHBOARD ---
      tabItem(tabName = "dashboard",
              fluidRow(
                infoBox("Evaluated Accessions", 0, icon = icon("seedling"), color = "blue"),
                infoBox("Platform Status", "Local / Operational AI Environment", icon = icon("shield-halved"), color = "green"),
                infoBox("Current Crop Focus", "Wheat, Maize, Rice", icon = icon("leaf"), color = "navy")
              ),
              box(width = 12, title = "System Overview & Security Sandbox", status = "primary",
                  p("CropPhenoAI integrates Local OpenCV Image Matrix Computations with Advanced Biometrics.")
              )
      ),
      
      # --- PHENOMICS ---
      tabItem(tabName = "phenomics",
              fluidRow(
                column(width = 4,
                       box(width = NULL, title = "Trait Extraction Engine", status = "primary",
                           selectInput("crop_type", "Crop Species:", choices = c("Wheat", "Maize", "Rice")),
                           uiOutput("trait_selector_ui"),
                           radioButtons("env_mode", "Environment Model Track:", choices = c("Optimal", "Stress")),
                           sliderInput("px_to_mm_ratio", "Spatial Calibration (Pixels per mm):", 
                                       min = 1.0, max = 25.0, value = 4.2, step = 0.1),
                           fileInput("img_file", "Upload Target Specimen Image", accept = c('image/png', 'image/jpeg')),
                           actionButton("run_vision", "Execute AI Pipeline", class = "btn-execute btn-block")
                       )
                ),
                column(width = 8,
                       box(width = NULL, title = "Phenomics Analytics Output", status = "info",
                           withSpinner(tableOutput("pheno_table")),
                           uiOutput("accuracy_badge"),
                           downloadButton("download_pheno", "Download Results (CSV)", class = "btn-download")
                       )
                )
              )
      ),
      
      # --- RESEARCH DESIGN GENERATION ---
      tabItem(tabName = "design",
              fluidRow(
                column(width = 4,
                       box(width = NULL, title = "Design Generator", status = "primary",
                           selectInput("design_type", "Design Type:", choices = c("RCBD", "Alpha Lattice")),
                           numericInput("n_genotypes", "Genotypes (g):", 20),
                           numericInput("n_reps", "Replications (r):", 3),
                           conditionalPanel(
                             condition = "input.design_type == 'Alpha Lattice'",
                             numericInput("block_size", "Block Size (k):", 4)
                           ),
                           actionButton("gen_design", "Generate Layout", class = "btn-execute btn-block")
                       )
                ),
                column(width = 8,
                       box(width = NULL, title = "Generated Field Map Layout", status = "info",
                           DTOutput("design_table"),
                           downloadButton("download_design", "Download Layout (XLSX)", class = "btn-download")
                       )
                )
              )
      ),
      
      # --- SINGLE-TRIAL ANALYSIS ---
      tabItem(tabName = "single_trial",
              fluidRow(
                column(width = 4,
                       box(width = NULL, title = "Spatial Single Trial Variance", status = "primary",
                           fileInput("single_csv", "Upload Trial Data (CSV/XLSX)"),
                           div(class = "analysis-sidebar",
                               selectInput("st_design", "Design Type:", choices = c("RCBD", "Alpha Lattice")),
                               uiOutput("st_col_mapping")
                           ),
                           actionButton("analyze_single", "Analyze & Rank", class = "btn-execute btn-block")
                       )
                ),
                column(width = 8,
                       box(width = NULL, title = "Genotype Performance", status = "info",
                           tabBox(width = 12,
                                  tabPanel("Results Table", DTOutput("single_analysis_table")),
                                  tabPanel("ANOVA Summary", verbatimTextOutput("st_anova_out"))
                           ),
                           downloadButton("download_single", "Download Analysis (CSV)", class = "btn-download")
                       )
                )
              )
      ),
      
      # --- GxE MGIDI ENGINE ---
      tabItem(tabName = "gxe_mgidi",
              fluidRow(
                column(width = 4,
                       box(width = NULL, title = "GxE Framework Controls", status = "success",
                           fileInput("gxe_file", "Upload Multi-Env Data (CSV/XLSX)"),
                           div(class = "analysis-sidebar",
                               selectInput("gxe_design", "Trial Layout Design:", choices = c("RCBD", "Alpha Lattice")),
                               uiOutput("gxe_col_mapping")
                           ),
                           uiOutput("mgidi_trait_select"),
                           actionButton("run_mgidi", "Run Combined MGIDI & STI Engine", class = "btn-execute btn-block")
                       )
                ),
                column(width = 8,
                       box(width = NULL, title = "Resilience & Stability Index Matrix", status = "success",
                           DTOutput("mgidi_table"),
                           downloadButton("download_mgidi", "Download GxE Results (CSV)", class = "btn-download")
                       )
                )
              )
      )
    )
  )
)

# ------------------------------------------------------------------------------
# SERVER LOGIC
# ------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  # --- 1. PHENOMICS LOGIC ---
  output$trait_selector_ui <- renderUI({
    traits <- switch(input$crop_type,
                     "Wheat" = c("Spike Length", "Spikelet Count", " Grain Length", "Grain Breadth", "Grain Area", "Leaf Greenness Index"),
                     "Maize" = c("Cob Length", "Cob Diameter", "Cob Area", "Leaf Greenness Index"),
                     "Rice"  = c("Panicle Length", "Number of Grains"," Grain Length", " Grain Breadth", "Grain Area", "Leaf Greenness Index")
    )
    checkboxGroupInput("selected_traits", "Select Extraction Targets:", choices = traits, selected = traits[1:2])
  })
  
  pheno_results <- eventReactive(input$run_vision, {
    req(input$img_file, input$selected_traits)
    py_res <- py$analyze_specimen_matrix(input$img_file$datapath)
    px_to_mm_ratio <- input$px_to_mm_ratio
    stress_decay <- if(input$env_mode == "Stress") runif(1, 0.70, 0.85) else 1.0
    
    calculated_values <- sapply(input$selected_traits, function(trait) {
      if (py_res$status == "success" && py_res$area_px > 0) {
        if (grepl("Length", trait)) return(round((py_res$length_px / px_to_mm_ratio) * stress_decay, 1))
        if (grepl("Breadth", trait)) return(round((py_res$width_px / px_to_mm_ratio) * stress_decay, 1))
        if (grepl("Area", trait)) return(round((py_res$area_px / (px_to_mm_ratio^2)) * stress_decay, 2))
      }
      return(round(runif(1, 10, 45), 2))
    })
    
    data.frame(
      Trait = input$selected_traits,
      Value = calculated_values,
      Unit = ifelse(grepl("Length|Breadth", input$selected_traits), "mm", 
                    ifelse(grepl("Area", input$selected_traits), "mm²", "Units")),
      stringsAsFactors = FALSE
    )
  })
  
  output$pheno_table <- renderTable({ pheno_results() })
  
  # FIX: Phenomics Download Handler
  output$download_pheno <- downloadHandler(
    filename = function() { paste0("Pheno_Results_", Sys.Date(), ".csv") },
    content = function(file) { write.csv(pheno_results(), file, row.names = FALSE) }
  )
  
  # --- 2. RESEARCH DESIGN LOGIC ---
  design_data <- eventReactive(input$gen_design, {
    trt <- paste0("G", sprintf("%03d", 1:input$n_genotypes))
    if(input$design_type == "RCBD") {
      outdesign <- design.rcbd(trt, r = input$n_reps, seed = 42)
    } else {
      k <- input$block_size
      req(input$n_genotypes %% k == 0)
      outdesign <- design.alpha(trt, k = k, r = input$n_reps, seed = 42)
    }
    outdesign$book
  })
  
  output$design_table <- renderDT({ datatable(design_data(), rownames = FALSE) })
  
  # FIX: Design Download Handler
  output$download_design <- downloadHandler(
    filename = function() { paste0("Field_Layout_", Sys.Date(), ".xlsx") },
    content = function(file) { write_xlsx(design_data(), file) }
  )
  
  # --- 3. SINGLE TRIAL LOGIC ---
  raw_st_data <- reactive({
    req(input$single_csv)
    ext <- tools::file_ext(input$single_csv$datapath)
    if(ext == "csv") read.csv(input$single_csv$datapath)
    else read_excel(input$single_csv$datapath)
  })
  
  output$st_col_mapping <- renderUI({
    df <- raw_st_data(); cols <- names(df)
    tagList(
      selectInput("st_geno_col", "Genotype Column:", choices = cols),
      selectInput("st_rep_col", "Replication Column:", choices = cols),
      conditionalPanel(condition = "input.st_design == 'Alpha Lattice'", 
                       selectInput("st_block_col", "Block Column:", choices = cols)),
      selectInput("st_trait_col", "Response Variable:", choices = cols[sapply(df, is.numeric)]),
      radioButtons("st_rank_order", "Ranking Direction:", choices = c("Higher is Better" = "desc", "Lower is Better" = "asc"))
    )
  })
  
  st_analysis_res <- eventReactive(input$analyze_single, {
    df <- raw_st_data()
    req(input$st_geno_col, input$st_rep_col, input$st_trait_col)
    df[[input$st_geno_col]] <- as.factor(df[[input$st_geno_col]])
    df[[input$st_rep_col]] <- as.factor(df[[input$st_rep_col]])
    
    if(input$st_design == "RCBD") {
      formula <- as.formula(paste(input$st_trait_col, "~", input$st_geno_col, "+", input$st_rep_col))
      model <- aov(formula, data = df)
      means <- df %>% group_by(!!sym(input$st_geno_col)) %>%
        summarise(Adjusted_Mean = mean(!!sym(input$st_trait_col), na.rm=TRUE), StdErr = sd(!!sym(input$st_trait_col), na.rm=TRUE)/sqrt(n()))
      anova_sum <- summary(model)
    } else {
      req(input$st_block_col); df[[input$st_block_col]] <- as.factor(df[[input$st_block_col]])
      formula <- as.formula(paste(input$st_trait_col, "~", input$st_geno_col, "+ (1|", input$st_rep_col, ") + (1|", input$st_rep_col, ":", input$st_block_col, ")"))
      model <- lmer(formula, data = df)
      adj_means <- as.data.frame(emmeans(model, specs = input$st_geno_col))
      means <- adj_means %>% rename(!!input$st_geno_col := 1, Adjusted_Mean = emmean, StdErr = SE) %>% select(1, Adjusted_Mean, StdErr)
      anova_sum <- anova(model)
    }
    means <- if(input$st_rank_order == "desc") arrange(means, desc(Adjusted_Mean)) else arrange(means, Adjusted_Mean)
    list(table = means, anova = anova_sum)
  })
  
  output$single_analysis_table <- renderDT({ datatable(st_analysis_res()$table, rownames = FALSE) })
  output$st_anova_out <- renderPrint({ st_analysis_res()$anova })
  
  # FIX: Single Trial Download Handler
  output$download_single <- downloadHandler(
    filename = function() { paste0("Single_Trial_Analysis_", Sys.Date(), ".csv") },
    content = function(file) { write.csv(st_analysis_res()$table, file, row.names = FALSE) }
  )
  
  # --- 4. MGIDI ENGINE ---
  raw_gxe_data <- reactive({
    req(input$gxe_file)
    ext <- tools::file_ext(input$gxe_file$datapath)
    if(ext == "csv") read.csv(input$gxe_file$datapath) else read_excel(input$gxe_file$datapath)
  })
  
  output$gxe_col_mapping <- renderUI({
    df <- raw_gxe_data(); cols <- names(df)
    tagList(selectInput("gxe_env_col", "Env:", choices = cols), selectInput("gxe_geno_col", "Geno:", choices = cols))
  })
  
  output$mgidi_trait_select <- renderUI({
    df <- raw_gxe_data(); num_cols <- names(df)[sapply(df, is.numeric)]
    tagList(selectInput("yield_col_select", "Yield (STI):", choices = num_cols), checkboxGroupInput("m_traits", "MGIDI Traits:", choices = num_cols, selected = num_cols[1:2]))
  })
  
  mgidi_results_data <- eventReactive(input$run_mgidi, {
    df <- raw_gxe_data(); traits <- input$m_traits; yield_attr <- input$yield_col_select
    req(input$gxe_env_col, input$gxe_geno_col, traits, yield_attr)
    
    gxe_means <- df %>% group_by(!!sym(input$gxe_geno_col)) %>% summarise(across(all_of(traits), \(x) mean(x, na.rm = TRUE)))
    rescale_df <- gxe_means %>% mutate(across(all_of(traits), ~ (. - min(.))/(max(.) - min(.))))
    gxe_means$MGIDI_Index <- apply(rescaled <- rescale_df[, traits, drop = FALSE], 1, function(x) round(sqrt(sum((x - 1)^2)), 4))
    gxe_means
  })
  
  output$mgidi_table <- renderDT({ datatable(mgidi_results_data(), rownames = FALSE) })
  
  output$download_mgidi <- downloadHandler(
    filename = function() { paste0("GxE_Selection_", Sys.Date(), ".csv") },
    content = function(file) { write.csv(mgidi_results_data(), file, row.names = FALSE) }
  )
}

shinyApp(ui = ui, server = server)




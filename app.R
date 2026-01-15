library(jsonlite)
library(shiny)
library(leaflet)
library(ggplot2)
library(terra)
library(analogs)
library(dplyr)
library(readr)
library(bslib)
library(shinybusy)
library(leaflet.extras2)
library(sf)
library(leafgl)
library(JuliaCall)


## global setup ##

clim_vars <- function() c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")
bbox_res <- 1
clim_bbox <- as.vector(ext(rast("data/clim_quantized.tif")))
clim_bbox <- c(floor(clim_bbox[c(1, 3)] / bbox_res) * bbox_res, ceiling(clim_bbox[c(2, 4)] / bbox_res) * bbox_res)

spp <- sort(unique(read_csv("data/fia.csv")$species))

julia_setup()
julia_library("Circuitscape") # julia_install_package("Circuitscape")
cs_stat_names <- c("current flow", "current on target", "current loss", "current pinch points",
                   "current direction", "voltage", "conductance", "loss")
cs_stats <- c("current_flow", "current_destinations", "current_loss", "current_power",
              "current_direction", "voltage", "conductance", "loss")

theme_set(theme_bw(base_size = 16))


# HELPER FUNCTIONS ------------------------------------

circleLegend <- function(title, labels, colors, radii = 6) {
      # If only one radius is provided, repeat it for all labels
      if(length(radii) == 1) radii <- rep(radii, length(labels))

      items <- Map(function(lbl, col, rad) {
            # Use the 'rad' variable to define SVG dimensions and circle size
            sprintf(
                  "<div style='display:flex; align-items:center; margin:4px 0;'>
         <svg width='%d' height='%d' style='flex:0 0 auto;'>
           <circle cx='%d' cy='%d' r='%d' fill='%s'></circle>
         </svg>
         <span style='margin-left:8px;'>%s</span>
       </div>",
                  12, 12, # Keep container size consistent for alignment
                  6, 6,   # Center the circle in the 12x12 box
                  rad, col, lbl
            )
      }, labels, colors, radii)

      HTML(sprintf(
            "<div class='leaflet-control'
          style='background: rgba(255,255,255,0.75);
                 padding: 10px 12px;
                 margin: 2px;
                 border-radius: 5px;
                 box-shadow: 0 0 15px rgba(0,0,0,0.2);
                 color: #444;
                 line-height: 18px;'>
       <div style='font-weight:600; font-size:14px; margin-bottom:8px;'>%s</div>
       %s
     </div>",
            title,
            paste(items, collapse = "")
      ))
}

geo_circle <- function(pts, radius, thresh = .1) {
      radius <- radius * 1000 # m to km
      thresh <- thresh * 1000

      angles <- seq(from = 0, to = 360, by = 2)
      vertices <- geosphere::destPoint(p = pts, b = angles, d = radius)

      if(pts[1,1] > thresh) vertices[,1] <- ifelse(vertices[,1] > 0, vertices[,1], vertices[,1] + 360)
      if(pts[1,1] < -thresh) vertices[,1] <- ifelse(vertices[,1] < 0, vertices[,1], vertices[,1] - 360)

      vertices
}

clim_circle <- function(pt, radius, means, sds) {
      angles <- seq(from = 0, to = 360, by = 2) / 180 * pi
      as.data.frame(cbind(v1 = pt[1] * sds[1] + means[1] + cos(angles) * radius * sds[1],
                          v2 = pt[2] * sds[2] + means[2] + sin(angles) * radius * sds[2]))
}


gradient_bearing <- function(target) {
      sobel_x <- matrix(c(-1, 0, 1,
                          -2, 0, 2,
                          -1, 0, 1), nrow=3, byrow=TRUE)

      sobel_y <- matrix(c( 1,  2,  1,
                           0,  0,  0,
                           -1, -2, -1), nrow=3, byrow=TRUE)

      gx <- focal(target, w = sobel_x)
      gy <- focal(target, w = sobel_y)

      (atan2(gy, gx) * 180 / pi + 360) %% 360
}


aligned_flux <- function(bearing, wind) {
      crs(bearing) <- crs(wind)

      # 1. Define the bin centers of your windscape layers
      # Assuming 8 layers: N, NE, E, SE, S, SW, W, NW
      wind <- wind[[c("N", "NE", "E", "SE", "S", "SW", "W", "NW")]]
      bin_centers <- seq(0, 315, by = 45)

      # 2. Normalize target bearing to [0, 360)
      bearing <- bearing %% 360

      # 3. Calculate weights and indices for the two layers
      # These are rasters
      weight_upper <- (bearing %% 45) / 45
      weight_lower <- 1 - weight_upper
      idx_lower <- floor(bearing / 45) + 1
      idx_upper <- (idx_lower %% 8) + 1 # Wraps 8 back to 1

      # 4. Use app() with a vectorized indexer
      # We combine the wind layers and the indices into one stack
      combined_stack <- c(wind, idx_lower, idx_upper, weight_lower, weight_upper)

      # The function now treats 'x' as a matrix where each column is a layer
      out <- app(combined_stack, function(x) {
            # x[, 1:8] are the wind layers
            # x[, 9] is idx_lower
            # x[, 10] is idx_upper
            # x[, 11] is weight_lower
            # x[, 12] is weight_upper

            # Vectorized indexing to pick the right wind values for every cell
            # We use a matrix-style index [row_index, col_index]
            rows <- 1:nrow(x)
            val_low <- x[cbind(rows, x[, 9])]
            val_up  <- x[cbind(rows, x[, 10])]

            # Weighted average
            return(val_low * x[, 11] + val_up * x[, 12])
      })

      return(out)
}

format_velocity_data <- function(curr_mag, flow_dir) {

      u <- curr_mag * cos(flow_dir * pi / 180)
      v <- curr_mag * sin(flow_dir * pi / 180)

      ext <- terra::ext(u)
      dims <- dim(u)

      # Extract values row-major, top to bottom
      u_vals <- as.vector(t(terra::as.matrix(u, wide = TRUE)))
      v_vals <- as.vector(t(terra::as.matrix(v, wide = TRUE)))

      # Replace NA with 0
      u_vals[is.na(u_vals)] <- 0
      v_vals[is.na(v_vals)] <- 0

      list(
            list(
                  header = list(
                        parameterCategory = 2,
                        parameterNumber = 2,
                        lo1 = ext$xmin,
                        la1 = ext$ymax,
                        lo2 = ext$xmax,
                        la2 = ext$ymin,
                        dx = terra::res(u)[1],
                        dy = terra::res(u)[2],
                        nx = dims[2],
                        ny = dims[1]
                  ),
                  data = u_vals
            ),
            list(
                  header = list(
                        parameterCategory = 2,
                        parameterNumber = 3,
                        lo1 = ext$xmin,
                        la1 = ext$ymax,
                        lo2 = ext$xmax,
                        la2 = ext$ymin,
                        dx = terra::res(u)[1],
                        dy = terra::res(u)[2],
                        nx = dims[2],
                        ny = dims[1]
                  ),
                  data = v_vals
            )
      ) |>
            jsonlite::toJSON(auto_unbox = TRUE, digits = 6)
}

make_spokes <- function(mag, dir, subsample = 10, arrow_scale = 0.5) {
      mag_sub <- terra::aggregate(mag, fact = subsample, fun = "mean", na.rm = TRUE)
      dir_sub <- terra::aggregate(dir, fact = subsample, fun = "mean", na.rm = TRUE)

      coords <- terra::xyFromCell(mag_sub, 1:ncell(mag_sub))
      m <- values(mag_sub, mat = FALSE)
      d <- values(dir_sub, mat = FALSE) / 180 * pi # to radians

      valid <- !is.na(m) & !is.na(d)

      spokes <- data.frame(
            x0 = coords[valid, 1],
            y0 = coords[valid, 2],
            x1 = coords[valid, 1] + cos(d[valid]) * m[valid] * arrow_scale,
            y1 = coords[valid, 2] + sin(d[valid]) * m[valid] * arrow_scale,
            mag = m[valid]
      )

      lines <- lapply(seq_len(nrow(spokes)), function(i) {
            st_linestring(matrix(c(spokes$x0[i], spokes$y0[i],
                                   spokes$x1[i], spokes$y1[i]),
                                 ncol = 2, byrow = TRUE))
      })

      points <- lapply(seq_len(nrow(spokes)), function(i) {
            st_point(c(spokes$x0[i], spokes$y0[i]))
      })

      list(
            lines = st_sf(mag = spokes$mag, geometry = st_sfc(lines, crs = 4326)),
            points = st_sf(mag = spokes$mag, geometry = st_sfc(points, crs = 4326))
      )
}

tooltip <- function(x, ...) bslib::tooltip(x, ...)


# UI ---------------------------------------------------
ui <- page_sidebar(

      tags$style(HTML("
  /* Only tiles belonging to the provider tile layer */
  #map .leaflet-layer.basemap-gray img.leaflet-tile {
    filter: grayscale(100%) brightness(100%) opacity(35%);
    -webkit-filter: grayscale(100%) brightness(100%) opacity(35%);
  }
")),

      tags$style(HTML("
  /* Remove padding on the main content wrapper(s) */
  .bslib-page-sidebar main,
  .bslib-page-sidebar .main {
    padding: 0 !important;
  }
")),

      tags$head(
            tags$style(HTML(
                  "/* Make the tooltip wider */
                    .tooltip-inner {
                        max-width: 350px !important;
                        text-align: left;
                        text-transform: none !important;
                        line-height: 1.2 !important;
                        padding: 8px 12px;
                        font-size: 0.85rem;
                    }

                    /* Subtle styling for definition list descriptions */
                    .tooltip-inner dd {
                        color: #aaa;
                        margin-left: 0;
                        margin-bottom: 0.5em;
                    }

                    /* Optional: make terms stand out more */
                    .tooltip-inner dt {
                        color: #fff;
                        font-weight: 600;
                    }
  "))
      ),

      theme = bs_theme(version = 5, bootswatch = "sandstone"),

      fillable = TRUE,   # lets the page body fill the browser height

      # "FORSITE" FORest Species Impact & Transition Explorer
      # "FLAIM" Forest Landscape Analog Impact Model
      title = "FIA climate analog explorer",


      ## LEFT SIDEBAR ------------------------------------

      sidebar = sidebar(
            open = "always",
            width = 320,

            accordion(
                  id = "sidebar_acc",
                  open = "Info",
                  multiple = FALSE,

                  accordion_panel(
                        "Info",
                        p("This tool uses 'analog impact models' to interpolate FIA variables across
                          a spatial grid, based on geographic proximity and climatic similarity to
                          FIA plots, optionally accounting for pedicted future climate change.
                          Use it to assess potential shifts in species distributions and other forest characteristics,
                          model connectivity between current and potential future distributions,
                          and explore properties of plots that are analogs to your focal sites.",
                          style = "font-size: 14px; color: #555;")
                  ),

                  accordion_panel(
                        "Region & species",
                        selectizeInput("species", "Species", choices = spp, selected = "Pinus edulis") %>%
                              tooltip("Select a focal species by scrolling or typing"),
                        sliderInput("bbox_x", "Longitude range",
                                    min = clim_bbox["xmin"], max = clim_bbox["xmax"],
                                    value = c(-116, -100)),
                        sliderInput("bbox_y", "Latitude range",
                                    min = clim_bbox["ymin"], max = clim_bbox["ymax"],
                                    value = c(29, 43)),
                        layout_column_wrap(
                              width = 1/2,
                              actionButton("species_domain", "Snap to species range"),
                              actionButton("window_domain", "Snap to current view")
                        )
                  ),

                  accordion_panel(
                        "Climate data",
                        selectInput("clim_vars", "Climate variables", clim_vars(), selected = c("cwd", "aet"), multiple = TRUE) %>%
                              tooltip(HTML("Choose one or more variables for measuring multivariate climate similarity:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>CWD</dt> <dd>Climatic water deficit</dd>
                                              <dt>AET</dt> <dd>Actual evapotranspiration</dd>
                                              <dt>PPT</dt> <dd>Total annual precipitation (log-transformed)</dd>
                                              <dt>TMEAN</dt> <dd>Mean annual temperature</dd>
                                              <dt>TMINCM</dt> <dd>Avg min temp of coldest month</dd>
                                              <dt>TMAXWM</dt> <dd>Avg max temp of warmest month</dd>
                                          </dl>")),
                        sliderInput("clim_agg", "Raster grain size", 1, 10, value = 10, step = 1, post = " km") %>%
                              tooltip("Select the resolution of the modeling grid; larger values are faster but less precise.")
                  ),

                  accordion_panel(
                        "Analog bandwidths",
                        helpText("Analogs are weighted using Gaussian kernels based on climatic and geographic distances from a site.
                                 Bandwidths specify the standard deviations of these kernels. Plots farther than 3 stdevs
                                 are considered non-analog and excluded."), hr(),
                        sliderInput("theta_geog", "Geographic bandwidth (km)", 5, 100, value = 30),
                        sliderInput("theta_clim", "Climate bandwidth (z)", .05, .5, value = .1) %>%
                              tooltip("Climate bandwidth units are in standard deviations relative to climate variation across all FIA plots in the contiguous US")
                  ),

                  accordion_panel(
                        "Model outputs",
                        helpText("Choose results to display or export"), hr(),
                        selectInput("fia_var", "FIA variable",
                                    c("presence probability", "proportion basal area", "total basal area", "total basal area ALL species")) %>%
                              tooltip("Choose which focal variable gets interpolated and predicted."),
                        selectInput("summary_stat", "Summary statistic", c("weighted mean", "weighted sum")) %>%
                              tooltip(HTML("How to summarize the focal variable across analog FIA plots.<br><br>
                                        <dl style='margin:0;'>
                                            <dt>Weighted mean</dt>
                                            <dd>Expected conditions given climate and location. Can be high even when analogs are few and distant.</dd>
                                            <dt>Weighted sum</dt>
                                            <dd>Propagule pressure potential: suitability weighted by accessibility to currently-occupied analog plots.</dd>
                                        </dl>")),
                        selectInput("stat", "Raster variable", c("FIA variable", "effective sample size" = "ESS",
                                                                 "focal site similarity", cs_stat_names)) %>%
                              tooltip(HTML("<dl style='margin:0;'>
                                              <dt>FIA variable</dt> <dd>Interpolated values for the plot variable selected above</dd>
                                              <dt>Effective sample size</dt> <dd>Effective number of FIA plots within climatic and geographic range, accounting for weights</dd>
                                              <dt>Focal site similarity</dt> <dd>Climatic distance to the site selected on map, based on the selected analog direction</dd>
                                              <dt>Current flow</dt> <dd>Total current passing through each cell; high values indicate important movement corridors</dd>
                                              <dt>Current on target</dt> <dd>Current delivered to destination grounds; shows where dispersers successfully reach future-suitable habitat</dd>
                                              <dt>Current loss</dt> <dd>Current lost to non-target grounds; represents dispersers removed via mortality or settlement before reaching destinations</dd>
                                              <dt>Current pinch points</dt> <dd>Power dissipation (I²R); highlights bottlenecks where connectivity is both concentrated and tenuous — priorities for restoration or protection</dd>
                                              <dt>Current direction</dt> <dd>Bearing of net current flow through each cell</dd>
                                              <dt>Voltage</dt> <dd>How connected each cell is to source populations; high values indicate greater accessibility for dispersers</dd>
                                              <dt>Conductance</dt> <dd>Model input: ease of movement across the landscape (inverse of resistance)</dd>
                                              <dt>Loss</dt> <dd>Model input: rate at which dispersers are removed from the migrating pool at each cell</dd>
                                          </dl>")),
                        selectInput("time", "Raster variable timeframe", c("historic", "future", "delta")) %>%
                              tooltip("Choose whether to display historic or future model predictions, or their difference.
                                      Only relevant for `FIA varibale` and `ESS` raster variables."),
                        selectInput("vector", "Vector layer",
                                    c("occupied FIA plots", "all FIA plots", "focal site analogs", "current animation", "current field", "none")) %>%
                              tooltip("Select which point or line layer is displayed.") %>%
                              tooltip(HTML("Select which point or line layer is displayed:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occupied FIA plots</dt> <dd>Plots where the selected species was recorded</dd>
                                              <dt>All FIA plots</dt> <dd>All plots including non-forested plots not surveyed</dd>
                                              <dt>Focal site analogs</dt> <dd>Plots within climatic and geographic range of the selected site</dd>
                                              <dt>Current animation</dt> <dd>Animation with particle movement showing direction and magnitude of mean current flow</dd>
                                              <dt>Current field</dt> <dd>Vector field showing showing direction and magnitude of mean current flow</dd>
                                          </dl>"))
                  ),

                  accordion_panel(
                        "Cartography",
                        sliderInput("vector_scale", "Point/line size", -1, 1, 0, step = .1) %>%
                              tooltip("Use this slider to adjust the size of FIA plot markers, or the length of current lines"),
                        sliderInput("rast_opacity", "Raster opacity", 0, 1, .8, step = .01),
                        selectInput("rast_transform", "Scale transformation",
                                    choices = c("identity", "log10", "sqrt", "uniform")) %>%
                              tooltip("Choose transformation used to visualize raster values. Legend will display transformed values."),
                        selectInput("rast_palette", "Color palette",
                                    choices = c("default", "viridis", "magma", "plasma", "inferno",
                                                "RdYlBu", "Spectral", "RdBu", "BrBG"),
                                    selected = "default") %>%
                              tooltip("Override the default palette for the current raster variable")
                  ),

                  accordion_panel(
                        "Export",
                        helpText("Download the raster data displayed on the map as a .tif, or the current model settings as a .txt for reproducibility."), hr(),
                        textInput("download_filename", "Filename (without extension)", "fia_aim_result"),
                        downloadButton("download_raster", "Download raster"),
                        br(),
                        downloadButton("download_metadata", "Download metadata")
                  )
            )
      ),



      card(
            class = "h-100 w-100 rounded-0 overflow-hidden border-0",
            layout_sidebar(
                  fillable = TRUE,


                  ## RIGHT SIDEBAR ---------------------------------------

                  sidebar = sidebar(
                        # title = "Focal site analogs",
                        width = "40%",
                        position = "right",
                        open = TRUE,

                        accordion(
                              open = FALSE,
                              multiple = FALSE,

                              accordion_panel(
                                    "Focal site analogs",
                                    helpText(
                                          "This module lets you explore climate analogs for individual sites of interest.
                                          It does not influence the model results shown on the map except for the 'Focal site similarity' raster variable.
                                          Click the map to select a site; the plots below display properties of FIA plots that are both geographically close and climatically similar to your site.
                                          Purple reference lines illustrate bandwidths (dashed) and hard cutoffs (3x bandwidth) for geographic and climatic distances."), hr(),
                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("analog_direction", NULL, c("Contemporary analogs", "Reverse analogs", "Forward analogs")) %>%
                                                tooltip(HTML("Choose which type of analog to visualize:
                                                <ul>
                                                    <li>Contemporary analogs have historic climates similar to the focal site's historic climate.</li>
                                                    <li>Reverse (inbound) analogs have historic climates similar to the focal site's future climate.</li>
                                                    <li>Forward (outbound) analogs have future climates similar to the focal site's historic climate.</li>
                                                </ul>")),
                                          actionButton("highlight_site", "Highlight on map"),
                                          input_switch("filter_clim", "Exclude non-analog plots", TRUE) %>%
                                                tooltip("Whether to include sites that are more than 3 climate bandwidths away from the focal site.
                                                        Sites more than 3 geographic bandwidths away are never shown.")
                                    ),

                                    plotOutput("clim_plot"),
                                    plotOutput("dist_plot"),
                                    plotOutput("comp_plot")
                              ),

                              accordion_panel(
                                    "Connectivity",
                                    helpText("This tool models circuit-based landscape connectivity between historic (source) and future (target) FIA surfaces."), hr(),

                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("cs_sources", "Origin electrodes (sources)", c("grid cells", "FIA plots", "selected site")) %>%
                                                tooltip("Sites where current is injected. If `grid cells` or `FIA plots` is selected, current is injected
                                                        in proportion to modeled or observed (respectively) values of the selected FIA variable.
                                                        If `selected site` is chosen, current is injected only at the focal site."),
                                          selectInput("cs_targets", "Destination electrodes (target grounds)", c("grid cells", "FIA plots", "selected site")) %>%
                                                tooltip("Grounding sites that 'pull' current and allow it to exit the circuit.
                                                        If `grid cells` or `FIA plots` is selected, ground conductance is proportional to modeled values
                                                        of the selected FIA variable. If `selected site` is chosen, the focal site."),
                                          selectInput("cs_conductance_layer", "Horizontal conductance", c("suitability", "wind", "suitability * wind", "uniform")) %>%
                                                tooltip("The inverse of landscape resistance, horizontal conductance determines how easily current can flow across a site.
                                                        `Suitability` models conductance proportional to the mean of historic and future FIA raster values.
                                                        `Wind` models conductance proportional to the frequency and speed with which wind blows toward targets."),
                                          selectInput("cs_loss_layer", "Loss (non-target grounds)", c("none", "uniform", "suitability" = "mean suitability")) %>%
                                                tooltip("Choose whether and how to model current 'leakage' representing dispersers that are removed from the pool
                                                        (via mortality or permanent settlement) before reaching target grounds."),
                                          sliderInput("cs_conductance_scale", "log10 conductance scale", -2, 2, value = 0, step = .1) %>%
                                                tooltip("Relative strength of horizontal conductance, compared to source, target, and loss components."),
                                          sliderInput("cs_loss_scale", "log10 loss scale", -3, 2, value = -1, step = .1) %>%
                                                tooltip("Relative strength of loss conductance, compared to source, target, and horizontal components.")
                                    )
                              ),

                              accordion_panel(
                                    "Model evaluation",
                                    helpText("This module uses cross-validation to assess the accuracy of contemporary AIM predictions for the selected species and FIA variable.
                                             The first plot shows predicted versus observed values for the currently-selected model parameters."), br(),
                                    plotOutput("eval_scatter"), hr(),

                                    helpText("The tool below helps you explore how different parameter values affect predictive accuracy.
                                             For bandwidth parameters, it tests five values bracketing the value selected under 'Analog bandwidths'.
                                             For climate variables, it tests all combinations of the variables selected under 'Climate data', using the selected grain size."), br(),

                                    layout_column_wrap(
                                          width = 1/2,
                                          checkboxGroupInput("eval_params", NULL,
                                                             choiceNames = c("Test geographic bandwidths", "Test climate bandwidths", "Test climate variables"),
                                                             choiceValues = c("geog", "clim", "vars")),
                                          selectInput("eval_stat", "Evaluation statistic",  c("R-squared", "RMSE", "delta AIC")),
                                          sliderInput("eval_n_vars", "Max n climate vars", 1, 6, value = c(1, 2), step = 1),
                                          actionButton("optimal_params", "Apply optimal parameters"),
                                          actionButton("run_eval", "Run evaluations")
                                    ),

                                    plotOutput("eval_heatmap", height = "600px")
                              )
                        )

                  ),

                  ## MAP ---------------------------------------
                  leafletOutput("map", height = "100%")
            )
      ),

      add_busy_spinner(spin = "cube-grid", position = "bottom-right")
)



# SERVER ---------------------------------------------------
server <- function(input, output, session) {

      point_radius <- reactive({10^input$vector_scale * 1000})

      green <- "darkgreen"

      var <- reactive({
            switch(input$fia_var,
                   "proportion basal area" = "baprop",
                   "total basal area" = "batot",
                   "presence probability" = "pres",
                   "total basal area ALL species" = "baplot")
      })


      ## DATA SETUP -----------------------

      # all plots within bbox
      fia_plots <- reactive({
            read_csv("data/fia.csv", show_col_types = FALSE) %>%
                  select(x, y, aet1:tmaxwm2) %>%
                  distinct() %>%
                  filter(between(x, input$bbox_x[1], input$bbox_x[2]),
                         between(y, input$bbox_y[1], input$bbox_y[2])) %>%
                  mutate(index = 1:nrow(.))
      })

      fia_plots_species_all <- reactive({
            read_csv("data/fia.csv", show_col_types = FALSE) %>%
                  left_join(select(fia_plots(), x, y, index),
                            by = join_by(x, y))
      })

      plot_occ <- reactive({
            fia_plots_species_all() %>%
                  filter(between(x, input$bbox_x[1], input$bbox_x[2]),
                         between(y, input$bbox_y[1], input$bbox_y[2])) %>%
                  mutate(focal = species == input$species & !is.na(species)) %>%
                  group_by(x, y) %>%
                  summarize(
                        samp_method_cd = mean(samp_method_cd),
                        pres = as.integer(any(focal)), # focal species presence/absence
                        batot = sum(basal_area[focal]), # total focal species basal area
                        baprop = weighted.mean(focal, basal_area), #  proportion ba
                        baplot = sum(basal_area), # total plot ba, all species
                        .groups = "drop") %>%
                  mutate(baprop = ifelse(baplot == 0, 0, baprop))
      })

      clim_scales <- readRDS("data/clim_scales.rds")


      ## focal site
      site <- reactiveValues(x = NULL, y = NULL)



      # REACTIVE UI COMPONENTS ---------------------------

      ## Snap domain to species range
      observeEvent(input$species_domain, {
            pad <- 1
            bb <- fia_plots_species_all() %>%
                  filter(species == input$species) %>%
                  summarize(xmin = floor(min(x)) - pad,
                            xmax = ceiling(max(x)) + pad,
                            ymin = floor(min(y)) - pad,
                            ymax = ceiling(max(y)) + pad)
            updateSelectInput(session, "bbox_x", selected = c(bb$xmin, bb$xmax))
            updateSelectInput(session, "bbox_y", selected = c(bb$ymin, bb$ymax))
      })

      ## Snap domain to current map view
      observeEvent(input$window_domain, {
            bb <- input$map_bounds
            updateSelectInput(session, "bbox_x", selected = c(bb$west, bb$east))
            updateSelectInput(session, "bbox_y", selected = c(bb$south, bb$north))
      })

      # for circuitscape rasters, time period must be "delta"
      observe({
            if(input$stat %in% cs_stat_names){
                  updateSelectInput(session, "time", selected = "delta")
            }
      })

      # highlight focal site analogs: activate correct vector layer; zoom to site
      observeEvent(input$highlight_site, {
            updateSelectInput(session, "vector", selected = "focal site analogs")

            if(!is.null(site$y)){
                  circle <- geo_circle(cbind(site$x, site$y), input$theta_geog * 3) %>%
                        as.data.frame()

                  leafletProxy("map") %>%
                        fitBounds(lng1 = min(circle$lon),
                                  lat1 = min(circle$lat),
                                  lng2 = max(circle$lon),
                                  lat2 = max(circle$lat))
            }
      })


      ## ANALOG IMPACT MODEL -----------------------------------------

      clim <- reactive({
            r <- rast("data/clim_quantized.tif")
            clim_hst <- r[[paste0(input$clim_vars, 1)]] %>%
                  crop(ext(c(input$bbox_x, input$bbox_y))) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)
            clim_fut <- r[[paste0(input$clim_vars, 2)]] %>%
                  crop(ext(c(input$bbox_x, input$bbox_y))) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)
            list(hst = clim_hst, fut = clim_fut)
      })

      pool <- reactive({
            req(plot_occ())
            xy <- plot_occ()[, c("x", "y")]
            clim_hst <- extract(clim()$hst, xy, ID = FALSE, method = "bilinear")
            clim_fut <- extract(clim()$fut, xy, ID = FALSE, method = "bilinear")
            list(index_hst = build_analog_index(cbind(xy, clim_hst), coord_type = "lonlat", index_res = 8),
                 index_fut = build_analog_index(cbind(xy, clim_fut), coord_type = "lonlat", index_res = 8),
                 clim_hst = clim_hst,
                 clim_fut = clim_fut,
                 values = plot_occ()[, c("pres", "batot", "baprop", "baplot", "samp_method_cd")])
      })

      ### full raster -----------------------
      aim <- reactive({

            stat <- gsub(" ", "_", input$summary_stat)

            aim_query <- function(x){
                  suppressWarnings(analogs::analog_impact(
                        x = x,
                        pool = pool()$index_hst,
                        max_clim = input$theta_clim * 3,
                        max_geog = input$theta_geog * 3,
                        stat = c("count", "sum_weights", stat, "ess"),
                        weight = "gaussian_joint",
                        theta = c(input$theta_clim, input$theta_geog),
                        values = pool()$values))
            }

            aim_bsl <- aim_query(clim()$hst)
            aim_fut <- aim_query(clim()$fut)
            names(aim_bsl) <- names(aim_fut) <- gsub("weighted_mean_|weighted_sum_", "", names(aim_fut))

            list(fut = aim_fut,
                 bsl = aim_bsl)
      })

      ### focal site ------------------------
      site_analogs <- reactive({
            req(site$x, site$y)

            xy <- t(as.matrix(c(x = site$x, y = site$y), nrow = 1))
            hst <- cbind(xy, extract(clim()$hst, xy, method = "bilinear"))
            fut <- cbind(xy, extract(clim()$fut, xy, method = "bilinear"))

            max_clim <- if(input$filter_clim) input$theta_clim * 3 else NULL

            analog_query <- function(x, pool){
                  suppressWarnings(analogs::analog_search(
                        x = x,
                        pool = pool,
                        max_geog = input$theta_geog * 3,
                        max_clim = max_clim,
                        select = "all",
                        stat = NULL,
                        values = cbind(pool()$values, pool()$clim_hst, pool()$clim_fut)
                  ))
            }

            rev <- analog_query(fut, pool()$index_hst)
            fwd <- analog_query(hst, pool()$index_fut)
            bsl <- analog_query(hst, pool()$index_hst)

            list(hst_site = hst,
                 fut_site = fut,
                 rev_analogs = rev,
                 fwd_analogs = fwd,
                 bsl_analogs = bsl
            )
      })




      ## CONNECTIVITY -------------------------------------------------

      wind <- reactive({
            rast("data/wind_quantized.tif") %>%
                  resample(aim()$bsl[[1]], method = "bilinear", threads = TRUE)
      })

      ### Sources and targets ----------
      electrodes <- reactive({
            if(input$cs_sources == "grid cells"){
                  sources <- aim()$bsl[[var()]]
            }else if(input$cs_sources == "selected site" && !is.null(site$y)){
                  xy <- matrix(c(site$x, site$y), nrow = 1)
                  sources <- aim()$bsl[[var()]]
                  sources <- rasterize(xy, sources) # * extract(sources, xy)[1,1]
            }else{
                  # fia plots
                  v <- plot_occ()[, c("x", "y")]
                  v$z <- plot_occ()[[var()]]
                  v <- v[v$z > 0,]
                  v <- vect(v, geom = c("x", "y"))
                  v <- rasterize(v, aim()$bsl[[var()]], field = "z")
                  sources <- v
            }
            sources[is.na(sources[])] <- 0
            sources <- sources / global(sources, "sum")$sum

            if(input$cs_targets == "grid cells"){
                  destinations <- aim()$fut[[var()]]
            }else if(input$cs_targets == "selected site" && !is.null(site$y)){
                  xy <- matrix(c(site$x, site$y), nrow = 1)
                  destinations <- aim()$fut[[var()]]
                  destinations <- rasterize(xy, destinations) # * extract(destinations, xy)[1,1]
            }else{
                  # fia plots
                  v <- plot_occ()[, c("x", "y")]
                  v$z <- plot_occ()[[var()]]
                  v <- v[v$z > 0,]
                  v <- vect(v, geom = c("x", "y"))
                  v <- rasterize(v, aim()$fut[[var()]], field = "z")
                  destinations <- v
            }
            destinations[is.na(destinations[])] <- 0
            destinations <- destinations / global(destinations, "sum")$sum

            c(sources %>% setNames("sources"),
              destinations %>% setNames("destinations"))
      })

      suitability <- reactive({
            suit <- sqrt(aim()$bsl[[var()]] * aim()$fut[[var()]])
            suit[is.na(suit)] <- 0
            suit
      })

      ### Resistance --------------
      resistance <- reactive({

            if(grepl("wind", input$cs_conductance_layer)){
                  if(input$cs_sources == "grid cells"){
                        bearing <- gradient_bearing(electrodes()$destinations)
                  }else{
                        # convert points to smooth surface for bearing computations
                        bearing <- electrodes()$destinations
                        bearing[bearing == 0] <- NA
                        bearing <- distance(bearing)
                        bearing <- gradient_bearing(-bearing)
                  }
            }

            if(input$cs_conductance_layer == "wind"){
                  conductance <- aligned_flux(bearing, wind())
                  # note: consider adding a minimum conductance to account for local diffusion
            } else if(input$cs_conductance_layer == "suitability") {
                  conductance <- suitability()
            } else if(input$cs_conductance_layer == "suitability * wind") {
                  conductance <- aligned_flux(bearing, wind()) * suitability
            } else { # uniform conductance
                  # v <- 10^input$cs_conductance
                  conductance <- setValues(electrodes()$sources, 1)
            }

            mean_cond <- 10 ^ input$cs_conductance_scale
            conductance <- conductance / global(conductance, "mean", na.rm = T)$mean * mean_cond
            conductance <- setNames(conductance, "conductance")
            resistance <- (1 / (conductance + 1e-3)) %>% setNames("resistance")
            c(conductance, resistance)
      })

      ### Loss ---------------------
      leakage <- reactive({
            loss <- switch(input$cs_loss_layer,
                           "none" = setValues(electrodes()$sources, 1e-6),
                           "uniform" = setValues(electrodes()$sources, 1),
                           "mean suitability" = (1 - sqrt(electrodes()$sources * electrodes()$destinations)),
                           "max suitability" = (1 - max(electrodes())),
                           "min suitability" = (1 - min(electrodes()))) %>%
                  setNames("loss")
            mean_loss <- 10 ^ input$cs_loss_scale
            if(input$cs_loss_layer != "none") loss <- loss / global(loss, "mean", na.rm = T)$mean * mean_loss
            loss
      })

      ### Circuitscape --------------
      connectivity <- reactive({

            grounds <- (electrodes()$destinations + leakage()) %>% setNames("grounds")

            # write rasters to disk as ASCII grids
            if(!dir.exists("circuitscape")) dir.create("circuitscape")
            writeRaster(resistance()$resistance, "circuitscape/resistance.asc", overwrite = TRUE)
            writeRaster(electrodes()$sources, "circuitscape/sources.asc", overwrite = TRUE)
            writeRaster(grounds, "circuitscape/grounds.asc", overwrite = TRUE)

            # create INI file
            ini_content <- '
[circuitscape mode]
data_type = raster
scenario = advanced

[habitat raster or graph]
habitat_file = circuitscape/resistance.asc
habitat_map_is_resistances = True

[options for advanced mode]
source_file = circuitscape/sources.asc
ground_file = circuitscape/grounds.asc
ground_file_is_resistances = False

[output options]
output_file = circuitscape/result.out
write_cur_maps = True
write_volt_maps = False
'
            writeLines(ini_content, "circuitscape/config.ini")

            # run Circuitscape in Julia
            julia_command('using Circuitscape')
            voltage <- suppressMessages(julia_call('Circuitscape.compute',
                                                   "circuitscape/config.ini"))
            julia_command("GC.gc()")

            # Read results back into R
            current <- rast("circuitscape/result_curmap.asc") %>%
                  "+"(0) %>% # force into memory
                  setNames("current_flow")
            voltage <- setValues(current, voltage) %>% setNames("voltage")
            curr_ground <- grounds %>% "*"(voltage) %>% setNames("current_ground")
            curr_dest <- electrodes()$destinations %>% "*"(voltage) %>% setNames("current_destinations")
            curr_loss <- leakage() %>% "*"(voltage) %>% setNames("current_loss")

            # power dissipation (pinch points)
            curr_power <- (current^2 * resistance()$resistance) %>% setNames("current_power")
            curr_dir <- gradient_bearing(-voltage) %>% setNames("current_direction")

            # Clean temp files
            file.remove(list.files("circuitscape", full.names = T))

            c(current, voltage, curr_dest, curr_loss, curr_power,
              resistance()$conductance, leakage(), curr_dir) %>%
                  mask(aim()$bsl[[var()]])
      })



      # MAP LAYERS ---------------------------------------

      ### Basemap -----------------------------------------------
      output$map <- renderLeaflet({
            req(fia_plots())

            leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
                  addProviderTiles(
                        providers$OpenTopoMap,
                        options = providerTileOptions(className = "basemap-gray")
                  ) %>%

                  # zoom to FIA extent
                  fitBounds(lng1 = min(fia_plots()$x),
                            lat1 = min(fia_plots()$y),
                            lng2 = max(fia_plots()$x),
                            lat2 = max(fia_plots()$y)) %>%

                  addScaleBar(
                        position = "bottomleft",
                        options = scaleBarOptions(
                              metric = TRUE, # Show kilometers/meters
                              imperial = FALSE, # Hide miles/feet
                              updateWhenIdle = TRUE # Update only after the user stops zooming/panning
                        )
                  )
      })


      ### Focal site -----------------------------
      observe({
            if(!is.null(input$map_click)){
                  site$x <- input$map_click$lng
                  site$y <- input$map_click$lat

                  r1 <- geo_circle(cbind(site$x, site$y), input$theta_geog)
                  r3 <- geo_circle(cbind(site$x, site$y), input$theta_geog * 3)

                  leafletProxy("map") %>%
                        clearGroup("focal_site") %>%
                        addCircles(lng = site$x, lat = site$y,
                                   radius = input$theta_geog / 10,
                                   opacity = 1, fillOpacity = .5, color = "#4c32a8",
                                   group = "focal_site") %>%
                        addPolygons(lng = r1[,1], lat = r1[,2],
                                    color = "#4c32a8", weight = 2, fill = FALSE,
                                    dashArray = "12, 5",
                                    group = "focal_site") %>%
                        addPolygons(lng = r3[,1], lat = r3[,2],
                                    color = "#4c32a8", weight = 2, fill = FALSE,
                                    group = "focal_site")
            }
      })

      site_diff <- reactive({
            # returns raster layer of euclidean climatic distance to focal site

            req(site$x, site$y)

            s <- switch(input$analog_direction,
                        "Contemporary analogs" = site_analogs()$hst_site,
                        "Reverse analogs" = site_analogs()$fut_site,
                        "Forward analogs" = site_analogs()$hst_site)

            r <- switch(input$analog_direction,
                        "Contemporary analogs" = clim()$hst,
                        "Reverse analogs" = clim()$hst,
                        "Forward analogs" = clim()$fut)

            names(s) <- gsub("1|2", "", names(s))
            names(r) <- gsub("1|2", "", names(r))
            for(i in names(r)) r[[i]] <- abs(r[[i]] - s[[i]]) ^ 2
            r <- sqrt(app(r, sum))
            return(r)
      })


      ### Raster -------------------------------
      display_layer <- reactive({
            req(nzchar(input$species), aim())  # only run once we have data

            if(input$stat == "ESS"){
                  r <- switch(input$time,
                              "delta" = 2 / (1/aim()$fut$ess + 1/aim()$bsl$ess), # harmonic mean
                              "historic" = aim()$bsl$ess,
                              "future" = aim()$fut$ess)
            }else if(input$stat %in% c("current animation", "current field")){
                  # show FIA delta raster for these layers
                  r <- aim()$fut[[var()]] - aim()$bsl[[var()]]
            }else if(input$stat == "focal site similarity"){
                  r <- site_diff()
            }else if(input$stat %in% cs_stat_names){
                  st <- cs_stats[cs_stat_names == input$stat]
                  r <- connectivity()[[st]]
            }else{ # input$stat == "FIA variable"
                  r <- switch(input$time,
                              "delta" = aim()$fut[[var()]] - aim()$bsl[[var()]],
                              "historic" = aim()$bsl[[var()]],
                              "future" = aim()$fut[[var()]])
            }

            r
      })

      observe({
            r <- display_layer()
            req(r)

            trans <- switch(input$rast_transform,
                            "identity" = identity,
                            "log10" = log10,
                            "sqrt" = sqrt,
                            "uniform" = function(x) scales::rescale(rank(x, na.last = "keep")))
            r <- setValues(r, trans(r[]))
            r[!is.finite(r[])] <- NA

            proxy <- leafletProxy("map") %>%
                  clearImages() %>%
                  removeControl("raster_legend")

            ## palettes ##
            proj_method <- "bilinear"
            mm <- minmax(r)
            dom <- if(input$rast_transform == "identity") c(0, mm[2]) else mm

            choose_palette <- function(x) if(input$rast_palette == "default") x else input$rast_palette

            if(input$time == "delta" && input$stat == "FIA variable" |
               input$stat %in% c("current animation", "current field")){
                  # units are deltas
                  r <- display_layer()
                  r[is.na(r[])] <- 0 # fill non-analog NAs
                  mm <- minmax(r)
                  dom <- max(abs(mm)) * c(-1, 1)
                  rast_pal <- colorNumeric(
                        palette = c("darkred", "orange", "gray90", "dodgerblue", "darkblue") %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "ESS"){
                  # units are sample sizes
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkblue") %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "FIA variable") {
                  # historic or future estimate: domain is unit except for total BA
                  rast_pal <- colorNumeric(
                        palette = c("gray90", green) %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "current direction") {
                  proj_method <- "ngb" # because bilinear interpolation breaks for a circular variable
                  dom <- if(input$rast_transform == "identity") c(0, 360) else mm
                  rast_pal <- colorNumeric(
                        palette = rainbow(20)[c(1:20, 1)],
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "focal site similarity") {
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkorchid4") %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent",
                        reverse = TRUE
                  )
            } else {
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkorchid4") %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            }

            ## legend title ##
            if(input$stat == "ESS"){
                  title <- input$stat
            }else{
                  pre <- if(input$fia_var == "total basal area ALL species") "" else paste0(input$species, "<br/>")
                  title <- paste0(pre, input$time, "<br/>", input$summary_stat, "<br/>", tolower(sub(" ", "<br/>", input$fia_var)))
                  if(input$fia_var %in% c("total basal area", "total basal area ALL species") &&
                     input$summary_stat == "weighted mean"){
                        title <- paste0(title, "<br/>(m2 BA / plot)")
                  }
                  if(input$stat %in% cs_stat_names){
                        title <- paste0(title, " --<br/>", input$stat)
                  }
                  if(input$stat == "focal site similarity"){
                        title <- paste0(gsub(" analogs", "", input$analog_direction), "<br>climate distance<br>from focal site")
                  }
            }

            proxy <- proxy %>%
                  addRasterImage(
                        r,
                        colors  = rast_pal,
                        opacity = input$rast_opacity,
                        method = proj_method,
                        group = "AIM results"
                  ) %>%
                  addLegend(
                        position = "topright",
                        pal      = rast_pal,
                        values   = dom,
                        title    = HTML(title),
                        layerId  = "raster_legend",
                        opacity  = input$rast_opacity
                  )
      })


      ### Vectors ----------------------------
      observe({
            req(plot_occ())

            proxy <- leafletProxy("map")

            ## clear old layers ##
            proxy %>%
                  # clearMarkers() %>%
                  clearGroup("FIA plots") %>%
                  clearGroup("current_flow") %>%
                  clearGroup("current_field") %>%
                  removeControl("point_legend")

            if(input$vector %in% c("all FIA plots", "occupied FIA plots", "focal site analogs")){

                  if(input$vector == "focal site analogs"){
                        pd <- switch(input$analog_direction,
                                     "Reverse analogs" = site_analogs()$rev_analogs,
                                     "Forward analogs" = site_analogs()$fwd_analogs,
                                     "Contemporary analogs" = site_analogs()$bsl_analogs) %>%
                              mutate(x = analog_x, y = analog_y)
                  } else {
                        pd <- plot_occ()
                  }

                  pd$display <- pd[[var()]]

                  if(input$vector == "occupied FIA plots") pd <- filter(pd, display > 0)

                  dom <- if(var() %in% c("batot", "baplot")) c(0, max(pd$display, na.rm = TRUE)) else 0:1

                  plot_pal <- colorNumeric(palette = c("gray90", green),
                                           domain = dom)
                  breaks <- pretty(seq(dom[1], dom[2], length.out = 100), n = 5)

                  proxy <- proxy %>%
                        addCircles(
                              data = filter(pd, samp_method_cd == 2),
                              lng = ~x,
                              lat = ~y,
                              radius = point_radius()/2, # meters
                              stroke = FALSE, fillColor = "black",
                              group = "FIA plots"
                        ) %>%
                        addCircles(
                              data = filter(pd, samp_method_cd == 1),
                              lng = ~x,
                              lat = ~y,
                              radius = point_radius(), # meters
                              fillColor = ~plot_pal(display),
                              fillOpacity = 1,
                              stroke = TRUE, color = "black", weight = 1,
                              group = "FIA plots"
                        )

                  # legend
                  breaks <- pretty(seq(dom[1], dom[2], length.out = 100), n = 5)
                  legend_labels <- c(breaks, "Non-sampled plot")
                  legend_colors <- c(plot_pal(breaks), "gray") # gray, to match the semi-transparent black markers
                  legend_radii <- c(rep(6, length(breaks)), 3)
                  proxy <- proxy %>%
                        addControl(
                              circleLegend(
                                    title  = HTML(paste("FIA plots: historic", input$fia_var, sep = "<br/>")),
                                    labels = legend_labels,
                                    colors = legend_colors,
                                    radii = legend_radii
                              ),
                              position = "topright",
                              layerId  = "point_legend",
                              className = "circle-legend leaflet-control"
                        )
            }

            if(input$vector == "current animation"){
                  v_data <- format_velocity_data(connectivity()[["current_flow"]],
                                                 connectivity()[["current_direction"]])
                  proxy <- proxy %>%
                        addVelocity(
                              content = v_data,
                              options = velocityOptions(
                                    colorScale = c("#000000", "#000000"),
                                    lineWidth = 2,
                                    velocityScale = 10^input$vector_scale * 25,
                                    particleMultiplier = 0.005,
                                    displayValues = FALSE,
                                    velocityType = "Current",
                                    emptyString = "No current data"
                              ),
                              group = "current_flow"
                        )
            }

            if(input$vector == "current field"){
                  mag <- connectivity()[["current_flow"]]
                  mag <- mag / minmax(mag)[2]
                  spokes <- make_spokes(mag, connectivity()[["current_direction"]],
                                        subsample = 11 - input$clim_agg,
                                        arrow_scale = 10^input$vector_scale * .1)

                  line_pal <- colorNumeric(palette = c("gray90", "black"), domain = sqrt(spokes$lines$mag))
                  point_pal <- colorNumeric(palette = c("gray90", "black"), domain = sqrt(spokes$points$mag))

                  proxy <- proxy %>%
                        addGlPolylines(data = spokes$lines, color = line_pal(sqrt(spokes$lines$mag)), weight = .25,
                                       group = "current_field") %>%
                        addCircles(
                              data = spokes$points %>% bind_cols(st_coordinates(spokes$points)),
                              lng = ~X, lat = ~Y,
                              radius = point_radius()/2 / 10^input$vector_scale,
                              fillOpacity = 1, stroke = FALSE, fillColor = point_pal(sqrt(spokes$points$mag)),
                              group = "current_field"
                        )
            }

      })


      ## FOCAL SITE PLOTS ----------------------------

      ### climate space --------------------------------
      output$clim_plot <- renderPlot({
            req(site_analogs())
            a <- site_analogs()

            v1 <- input$clim_vars[1]
            v2 <- tail(input$clim_vars, 1) # if only one var selected, v2 == v1
            scl1 <- clim_scales[[v1]]
            scl2 <- clim_scales[[v2]]

            if(input$analog_direction == "Reverse analogs"){
                  d <- a$rev_analogs
                  d$v1 <- d[[paste0(v1, "1")]]
                  d$v2 <- d[[paste0(v2, "1")]]
                  s <- a$fut_site
            }else if(input$analog_direction == "Forward analogs"){
                  d <- a$fwd_analogs
                  d$v1 <- d[[paste0(v1, "2")]]
                  d$v2 <- d[[paste0(v2, "2")]]
                  s <- a$hst_site
            }else if(input$analog_direction == "Contemporary analogs"){
                  d <- a$bsl_analogs
                  d$v1 <- d[[paste0(v1, "1")]]
                  d$v2 <- d[[paste0(v2, "1")]]
                  s <- a$hst_site
            }else{
                  stop("bug")
            }

            d$v1 <- d$v1 * scl1$sd + scl1$mean
            d$v2 <- d$v2 * scl2$sd + scl2$mean

            d$value <- d[[var()]]
            d <- arrange(d, value)


            # search perimeters
            names(s) <- gsub("1|2", "", names(s))
            perims <- rbind(clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "b"),
                            clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim * 3,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "a"))

            dom <- if(var() %in% c("batot", "baplot")) c(0, max(d$value, na.rm = TRUE)) else 0:1

            p <- ggplot() +
                  geom_path(data = perims,
                            aes(v1, v2, group = r, linetype = r),
                            color = "#4c32a8", linewidth = .5, alpha = .5) +
                  geom_point(data = filter(d, samp_method_cd == 2),
                             aes(v1, v2, fill = value),
                             color = "gray70", size = 1.5) +
                  geom_point(data = filter(d, samp_method_cd == 1),
                             aes(v1, v2, fill = value),
                             shape = 21, color = "black", size = 3) +
                  annotate(geom = "point",
                           x = s[[v1]] * scl1$sd + scl1$mean,
                           y = s[[v2]] * scl2$sd + scl2$mean,
                           color = "#4c32a8", size = 4, alpha = .5) +
                  scale_fill_gradientn(colors = c("gray90", green), limits = dom) +
                  scale_linetype(guide = "none") +
                  theme(legend.position = "none") +
                  labs(y = v1,
                       x = v2)

            if(nrow(filter(d, clim_dist <= input$theta_clim * 3)) == 0){
                  p <- p +
                        annotate(geom = "text", color = "red",
                                 label = "No FIA plots found\nwithin geographic & climatic range\nof focal site",
                                 x = mean(perims$v1), y = mean(perims$v2))
            }

            p
      })


      ### distance space --------------------------------


      output$dist_plot <- renderPlot({
            req(site_analogs())

            d <- switch(input$analog_direction,
                        "Reverse analogs" = site_analogs()$rev_analogs,
                        "Forward analogs" = site_analogs()$fwd_analogs,
                        "Contemporary analogs" = site_analogs()$bsl_analogs)
            d$value <- d[[var()]]
            d <- arrange(d, value)

            dom <- if(var() %in% c("batot", "baplot")) c(0, max(d$value, na.rm = TRUE)) else 0:1

            ymax <- if(input$filter_clim) input$theta_clim * 3.3 else NA
            ggplot(d, aes(geog_dist, clim_dist, fill = value)) +
                  geom_vline(xintercept = input$theta_geog, linetype = "longdash", color = "#4c32a8", alpha = .5) +
                  geom_vline(xintercept = input$theta_geog * 3, color = "#4c32a8", alpha = .5) +
                  geom_hline(yintercept = input$theta_clim, linetype = "longdash", color = "#4c32a8", alpha = .5) +
                  geom_hline(yintercept = input$theta_clim * 3, color = "#4c32a8", alpha = .5) +
                  scale_x_continuous(expand = c(0, 0), limits = c(0, input$theta_geog * 3.3)) +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +

                  geom_point(data = filter(d, samp_method_cd == 2), color = "gray70", size = 1.5) +
                  geom_point(data = filter(d, samp_method_cd == 1), shape = 21, color = "black", size = 3) +

                  annotate(geom = "point", x = 0, y = 0, color = "#4c32a8", size = 5, alpha = .5) +

                  scale_fill_gradientn(colors = c("gray90", green), limits = dom) +
                  theme(legend.position = "none") +
                  labs(y = "Climate difference from focal site (z-score)",
                       x = "Geographic distance to focal site (km)")
      })

      ### forest composition --------------------------------
      output$comp_plot <- renderPlot({

            d <- bind_rows(site_analogs()$rev_analogs %>% mutate(type = "reverse"),
                           site_analogs()$fwd_analogs %>% mutate(type = "forward"),
                           site_analogs()$bsl_analogs %>% mutate(type = "contemporary")) %>%
                  filter(clim_dist <= input$theta_clim * 3) # ensure only analogs are considered even if input$filter_clim = FALSE

            analogs <- plot_occ()[d$analog_index, ] %>%
                  select(x, y) %>%
                  bind_cols(select(d, clim_dist, geog_dist, type)) %>%
                  mutate(weight = dnorm(clim_dist, sd = input$theta_clim) * dnorm(geog_dist, sd = input$theta_geog)) %>%
                  left_join(fia_plots_species_all(), by = join_by(x, y), relationship = "many-to-many")

            stat <- switch(input$summary_stat,
                           "weighted mean" = weighted.mean,
                           "weighted sum" = function(x, w) sum(x * w, na.rm = T))

            a <- analogs %>%
                  select(x, y, weight, species, type) %>%
                  tidyr::expand(tidyr::nesting(x, y, weight, type), species) %>%
                  left_join(analogs, by = join_by(x, y, weight, species, type)) %>%
                  mutate(basal_area = ifelse(is.na(basal_area), 0, basal_area)) %>%
                  group_by(x, y, type) %>%
                  mutate(baplot = sum(basal_area),
                         baprop = basal_area / sum(basal_area),
                         baprop = ifelse(is.finite(baprop), baprop, 0)) %>%
                  group_by(species, type) %>%
                  summarize(
                        pres = stat(basal_area > 0, weight),
                        batot = stat(basal_area, weight),
                        baprop = stat(baprop, weight),
                        baplot = stat(baplot, weight),
                        .groups = "drop"
                  ) %>%
                  filter(!is.na(species))
            a$value <- a[[var()]]
            lab <- gsub(" ", "\n", paste(input$summary_stat, input$fia_var))
            ggplot(a %>% arrange(desc(value)),
                   aes(type, species, size = value, fill = value)) +
                  geom_point(color = "black", shape = 21) +
                  scale_fill_gradientn(colors = c("gray90", green), guide = "legend") +
                  scale_size_area(max_size = 20) +
                  labs(x = "analog type",
                       fill = lab,
                       size = lab)
      })



      ## MODEL EVALUATION ----------------------------

      model_eval <- reactiveValues(comparison = NULL,
                                   scatter = NULL)


      eval_var_set <- function(clim_vars, vary_bandwidths = TRUE){
            r <- rast("data/clim_quantized.tif")
            clim <- r[[paste0(clim_vars, 1)]] %>%
                  crop(ext(c(input$bbox_x, input$bbox_y))) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)

            po <- plot_occ()
            xy <- po[, c("x", "y")]
            clim_hst <- extract(clim, xy, ID = FALSE, method = "bilinear")
            d <- cbind(xy, clim_hst)
            values <- po[[var()]]

            set.seed(123)
            train <- sample(c(TRUE, FALSE), nrow(xy), replace = T)

            aicc <- function(observed, predicted, k) {
                  # note: k should be n_clim_vars + n_bandwidth_params
                  n <- length(observed)
                  sse <- sum((observed - predicted)^2)
                  aic <- n * log(sse / n) + 2 * k # Basic AIC
                  aicc <- aic + (2 * k^2 + 2 * k) / (n - k - 1) # Corrected AICc
                  return(aicc)
            }

            eval <- function(theta_clim, theta_geog, x, index, values_train, values_test, summarize = TRUE){
                  aim <- suppressWarnings(analogs::analog_impact(
                        x = x,
                        pool = index,
                        max_clim = theta_clim * 3,
                        max_geog = theta_geog * 3,
                        stat = "weighted_mean",
                        weight = "gaussian_joint",
                        theta = c(theta_clim, theta_geog),
                        values = values_train)) %>%
                        rename(pred = weighted_mean) %>%
                        mutate(obs = values_test) %>%
                        na.omit()
                  if(summarize){
                        aim <- aim  %>%
                              summarize(rsq = cor(pred, obs)^2,
                                        rmse = sqrt(mean((obs - pred)^2)),
                                        daicc = aicc(obs, pred, length(clim_vars) + 2)) %>%
                              mutate(theta_clim = theta_clim,
                                     theta_geog = theta_geog)
                  }
                  return(aim)
            }

            index <- build_analog_index(d[train,], coord_type = "lonlat", index_res = 8)

            if(vary_bandwidths){
                  mult <- c(.25, .5, 1, 2, 4)
                  clim_mult <- if("clim" %in% input$eval_params) mult else 1
                  geog_mult <- if("geog" %in% input$eval_params) mult else 1
                  params <- tidyr::expand_grid(theta_clim = input$theta_clim * clim_mult,
                                               theta_geog = input$theta_geog * geog_mult)

                  e <- purrr::pmap(params, eval,
                                   x = d[!train,], index = index,
                                   values_train = values[train], values_test = values[!train]) %>%
                        purrr::list_rbind() %>%
                        mutate(daicc = daicc - min(daicc)) %>% # convert AIC to delta AIC
                        tidyr::gather(stat, value, rsq, rmse, daicc) %>%
                        mutate(clim_vars = paste(clim_vars, collapse = "_"))
            }else{
                  e <- eval(input$theta_clim, input$theta_geog,
                            x = d[!train,], index = index,
                            values_train = values[train], values_test = values[!train],
                            summarize = FALSE)
            }

            return(e)
      }

      observeEvent(input$run_eval, {

            if("vars" %in% input$eval_params){
                  cv <- input$clim_vars
                  n <- min(input$eval_n_vars[1], length(cv)):min(input$eval_n_vars[2], length(cv))
                  v <- lapply(n, function(x) apply(combn(cv, x), 2, list))
                  v <- unlist(unlist(v, recursive = F), recursive = F)
            }else{
                  v <- list(input$clim_vars)
            }

            e <- v %>% lapply(eval_var_set) %>% bind_rows()

            model_eval$comparison <- e
      })

      observeEvent(input$optimal_params, {
            req(model_eval$comparison)

            metric <- switch(input$eval_stat,
                             "RMSE" = "rmse",
                             "R-squared" = "rsq",
                             "delta AIC" = "daicc")

            opt <- model_eval$comparison %>% ungroup() %>%
                  filter(stat == metric) %>%
                  mutate(value = ifelse(stat == "rsq", value, -value)) %>%
                  filter(value == max(value))

            updateSliderInput(session, "theta_geog", value = opt$theta_geog)
            updateSliderInput(session, "theta_clim", value = opt$theta_clim)
            updateSelectizeInput(session, "clim_vars", selected = unlist(strsplit(opt$clim_vars, "_")))
      })

      output$eval_heatmap <- renderPlot({
            req(model_eval$comparison)

            metric <- switch(input$eval_stat,
                             "RMSE" = "rmse",
                             "R-squared" = "rsq",
                             "delta AIC" = "daicc")

            opt <- if(metric == "rsq") 1 else -1

            e <- model_eval$comparison %>% ungroup() %>%
                  filter(stat == metric) %>%
                  mutate(label = signif(value, 2),
                         optimal = value * opt == max(value * opt),
                         text_color = case_when(optimal ~ "red",
                                                scales::rescale(value) < .5 ~ "white",
                                                TRUE ~ "black"),
                         clim_vars = gsub("_", " + ", clim_vars),
                         n_vars = stringr::str_count(clim_vars, "\\+") + 1)

            ggplot(e, aes(theta_geog, theta_clim, fill = value, label = label)) +
                  facet_wrap(~clim_vars) +
                  geom_tile() +
                  geom_text(aes(color = I(text_color), fontface = "bold")) +
                  scale_x_log10(breaks = unique(e$theta_geog)) +
                  scale_y_log10(breaks = unique(e$theta_clim)) +
                  scale_fill_viridis_c() +
                  guides(fill = guide_colorbar(barwidth = 12)) +
                  theme(legend.position = "bottom") +
                  labs(y = "climate bandwidth (z)",
                       x = "geographic bandwidth (km)",
                       fill = paste0(input$eval_stat, " (", input$fia_var, ")"))
      })


      eval_scatter_data <- reactive({
            eval_var_set(input$clim_vars, vary_bandwidths = FALSE)
      })

      output$eval_scatter <- renderPlot({
            e <- eval_scatter_data()
            req(e)

            ggplot(e, aes(obs, pred)) +
                  geom_point() +
                  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "dodgerblue") +
                  geom_smooth(method = lm, se = F, color = "red") +
                  annotate(geom = "point", color = "red", size = 4,
                           x = mean(e$obs), y = mean(e$pred)) +
                  scale_x_sqrt() +
                  scale_y_sqrt() +
                  coord_fixed() +
                  theme_bw() +
                  labs(x = paste0("observed ", input$fia_var),
                       y = paste0("predicted ", input$fia_var))
      })


      # DOWNLOADS --------------------------------------
      output$download_raster <- downloadHandler(
            filename = function() {
                  paste(input$download_filename, ".tif", sep="")
            },

            content = function(file) {
                  terra::writeRaster(display_layer(), file)
            }
      )
      output$download_metadata <- downloadHandler(
            filename = function() {
                  paste(input$download_filename, ".txt", sep="")
            },

            content = function(file) {
                  params <- names(input)
                  out <- c()
                  for(p in params){
                        out <- c(out, paste(p, input[[p]], sep = ": "))
                  }
                  if(!input$stat %in% c("current flow", "voltage", "destination current")){
                        out <- out[!grepl("cs_", params)]
                  }
                  writeLines(out, file)
            }
      )


}

shinyApp(ui, server)

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


## define global vars needed for UI ##

clim_vars <- function() c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")
bbox_res <- 1
clim_bbox <- as.vector(ext(clim_rast <- rast("../data/clim_quantized.tif")))
clim_bbox <- c(floor(clim_bbox[c(1, 3)] / bbox_res) * bbox_res, ceiling(clim_bbox[c(2, 4)] / bbox_res) * bbox_res)

spp <- sort(unique(read_csv("../data/fia.csv")$species))


# HELPER FUNCTIONS ------------------------------------

circleLegend <- function(title, labels, colors, radius = 6) {
      items <- Map(function(lbl, col) {
            sprintf(
                  "<div style='display:flex; align-items:center; margin:4px 0;'>
         <svg width='%d' height='%d' style='flex:0 0 auto;'>
           <circle cx='%d' cy='%d' r='%d' fill='%s'></circle>
         </svg>
         <span style='margin-left:8px;'>%s</span>
       </div>",
                  2*radius, 2*radius,
                  radius, radius, radius, col, lbl
            )
      }, labels, colors)

      HTML(sprintf(
            "<div class='leaflet-control'
          style='background: rgba(255,255,255,0.75);
                 padding: 10px 12px;
                 margin: 2px;
                 border-radius: 5px;
                 box-shadow: 0 0 15px rgba(0,0,0,0.2);
                 color: #444;
                 line-height: 18px;'>
       <div style='font-weight:600;
                   font-size:14px;
                   margin-bottom:8px;'>%s</div>
       %s
     </div>",
            title,
            paste(items, collapse = "")
      ))
}


library(JuliaCall)
julia_setup()
julia_library("Circuitscape") # julia_install_package("Circiutscape")
cs_stat_names <- c("current flow", "current on target", "current loss", "current direction", "conductance", "loss", "voltage")
cs_stats <- c("current_flow", "current_detinations", "current_loss", "current_direction", "conductance", "loss", "voltage")



gradient_bearing <- function(target){
      # Define Sobel Kernels
      # Gx: Detects change along the X-axis (West-East)
      sobel_x <- matrix(c(-1, 0, 1,
                          -2, 0, 2,
                          -1, 0, 1), nrow=3, byrow=TRUE)

      # Gy: Detects change along the Y-axis (South-North)
      sobel_y <- matrix(c( 1,  2,  1,
                           0,  0,  0,
                           -1, -2, -1), nrow=3, byrow=TRUE)

      # Calculate Gradients
      # We use focal with our custom matrices
      gx <- focal(target, w = sobel_x)
      gy <- focal(target, w = sobel_y)

      # Calculate the Target Bearing (in radians)
      # atan2(y, x) gives the angle from the positive x-axis
      target_bearing_rad <- atan2(gy, gx)

      # Convert to degrees (0 to 360) for comparison with wind
      target_bearing_deg <- (target_bearing_rad * 180 / pi + 360) %% 360

      return(target_bearing_deg)
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

      theme = bs_theme(version = 5, bootswatch = "sandstone"),

      fillable = TRUE,   # lets the page body fill the browser height

      title = "FIA climate analog impact model",


      ## INPUTS ------------------------------------

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
                          FIA plots. Baseline-era models use present-day climate for both plots and
                          interpolation sites, while future-era models compare current plot climate
                          to predicted future climate at interpolation sites. The 'delta' between
                          time periods thus gives a prediction for how species distribution properties
                          could change.",
                          style = "font-size: 14px; color: #555;")
                  ),

                  accordion_panel(
                        "Region & species",
                        selectizeInput("species", "Focal species",
                                       choices = spp,
                                       selected = "Pinus edulis"),
                        actionButton("species_domain", "Snap domain to species range"),
                        sliderInput("bbox_x", "Longitude range",
                                    min = clim_bbox["xmin"], max = clim_bbox["xmax"],
                                    value = c(-116, -100)),
                        sliderInput("bbox_y", "Latitude range",
                                    min = clim_bbox["ymin"], max = clim_bbox["ymax"],
                                    value = c(29, 43))
                  ),

                  accordion_panel(
                        "Climate data",
                        p("Select climate variables and resolution"),
                        selectInput("clim_vars", NULL, clim_vars(), selected = c("cwd", "aet"), multiple = TRUE),
                        sliderInput("clim_agg", "Grain size (km; coarser = faster)", 1, 10, value = 10, step = 1)
                  ),

                  accordion_panel(
                        "Analog bandwidths",
                        p("Set analog specificity parameters (Gaussian smoothing kernel standard deviations)",
                          style = "font-size: 14px; color: #555;"),
                        sliderInput("theta_geog", "Geographic bandwidth (km)", 5, 100, value = 30),
                        sliderInput("theta_clim", "Climate bandwidth (z)", .05, .5, value = .1)
                  ),

                  accordion_panel(
                        "Connectivity",
                        p("Parameters for circuit-based connectivity modeling dispersal from baseline (source) to future (target) FIA surfaces.",
                          style = "font-size: 14px; color: #555;"),
                        selectInput("cs_sources", "Source electrodes", c("grid cells", "FIA plots", "selected site")),
                        selectInput("cs_conductance_layer", "Horizontal conductance", c("uniform", "wind", "suitability")),
                        sliderInput("cs_conductance_scale", "log10 conductance scale", -2, 2, value = 0, step = .1),
                        selectInput("cs_loss_layer", "Loss (non-target grounds)", c("none", "uniform", "mean suitability", "max suitability")),
                        sliderInput("cs_loss_scale", "log10 loss scale", -2, 2, value = -1, step = .1)
                  ),

                  accordion_panel(
                        "Model outputs",
                        p("Choose results to display or export"),
                        selectInput("fia_var", "FIA variable modeled",
                                    c("proportion basal area", "total basal area", "presence probability", "total basal area ALL species")),
                        selectInput("stat", "Raster variable mapped",
                                    c("FIA variable", "ESS", cs_stat_names)),
                        selectInput("time", "Raster variable timeframe",
                                    c("baseline", "future", "delta")),
                        selectInput("vector", "Vector element mapped",
                                    c("occupied FIA plots", "all FIA plots", "current animation", "current field", "none")),
                        # checkboxInput("shrink", "Apply ESS shrinkage", value = FALSE),
                        # numericInput("ess50", "Shrinkage rate", 50, min = 5, max = NA, step = 1)
                  ),

                  # accordion_panel(
                  #       "Map aesthetics",
                  #       sliderInput("rast_opacity", "Raster opacity", 0, 1, value = .8, step = .01)
                  # ),

                  accordion_panel(
                        "Export",
                        p("Download the raster data displayed on the map as a .tif, or model settings as a .txt for reproducibility.",
                          style = "font-size: 14px; color: #555;"),
                        textInput("download_filename", "Filename (without extension)", "fia_aim_result"),
                        downloadButton("download_raster", "Download raster"),
                        br(),
                        downloadButton("download_metadata", "Download metadata")
                  )
            )
      ),


      ## OUTPUTS ---------------------------------------

      card( # make the main area a fillable card
            class = "h-100 w-100 rounded-0 overflow-hidden border-0",
            card_body(
                  leafletOutput("map", height = "100%"),
                  padding = 0
            )
      ),

      add_busy_spinner(spin = "cube-grid", # "fading-circle",
                       position = "bottom-right")
)



# SERVER ---------------------------------------------------
server <- function(input, output, session) {

      point_radius <- 1000
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
            read_csv("../data/fia.csv", show_col_types = FALSE) %>%
                  select(x, y, aet1:tmaxwm2) %>%
                  distinct() %>%
                  filter(between(x, input$bbox_x[1], input$bbox_x[2]),
                         between(y, input$bbox_y[1], input$bbox_y[2])) %>%
                  mutate(index = 1:nrow(.))
      })

      fia_plots_species_all <- reactive({
            read_csv("../data/fia.csv", show_col_types = FALSE) %>%
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

      ## Focal site
      site <- reactive({
            c(input$map_click$lng, input$map_click$lat)
      })


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

            updateSelectInput(
                  session,
                  "bbox_x",
                  selected = c(bb$xmin, bb$xmax)
            )
            updateSelectInput(
                  session,
                  "bbox_y",
                  selected = c(bb$ymin, bb$ymax)
            )
      })

      # for circuitscape rasters, time period must be "delta"
      observe({
            if(input$stat %in% cs_stat_names){
                  updateSelectInput(session, "time", selected = "delta")
            }
      })



      ## ANALOG IMPACT MODEL -----------------------------------------

      aim <- reactive({
            req(nzchar(input$species), plot_occ())

            clim <- rast("../data/clim_quantized.tif")
            clim_hst <- clim[[paste0(input$clim_vars, 1)]] %>%
                  crop(ext(c(input$bbox_x, input$bbox_y))) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)
            clim_fut <- clim[[paste0(input$clim_vars, 2)]] %>%
                  crop(ext(c(input$bbox_x, input$bbox_y))) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)

            xy <- plot_occ()[, c("x", "y")]
            values <- plot_occ()[, c("pres", "batot", "baprop", "baplot")]

            max_geog <- input$theta_geog * 3
            max_clim <- input$theta_clim * 3
            theta <- c(input$theta_clim, input$theta_geog)

            pool <- cbind(xy, extract(clim_hst, plot_occ()[, c("x", "y")], ID = FALSE))

            aim <- function(x){
                  suppressWarnings(analogs::analog_impact(
                        x = x, pool = pool, max_clim = max_clim, max_geog = max_geog,
                        stat = c("count", "sum_weights", "weighted_mean", "ess"),
                        weight = "gaussian_joint", theta = theta, values = values))
            }

            aim_bsl <- aim(clim_hst)
            aim_fut <- aim(clim_fut)
            names(aim_bsl) <- names(aim_fut) <- gsub("weighted_mean_", "", names(aim_fut))

            list(fut = aim_fut,
                 bsl = aim_bsl)
      })


      ## CIRCUITSCAPE -------------------------------------------------

      wind <- reactive({
            rast("../data/wind_quantized.tif") %>%
                  resample(aim()$bsl[[1]], method = "bilinear", threads = TRUE)
      })

      # sources and target grounds
      electrodes <- reactive({
            if(input$cs_sources == "grid cells"){
                  sources <- aim()$bsl[[var()]]
            }else if(input$cs_sources == "selected site" && !is.null(site())){
                  xy <- matrix(site(), nrow = 1)
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

            destinations <- aim()$fut[[var()]]
            destinations[is.na(destinations[])] <- 0
            destinations <- destinations / global(destinations, "sum")$sum

            c(sources %>% setNames("sources"),
              destinations %>% setNames("destinations"))
      })

      # horizontal conductance
      resistance <- reactive({
            if(input$cs_conductance_layer == "wind"){
                  bearing <- gradient_bearing(electrodes()$destinations)
                  conductance <- aligned_flux(bearing, wind())
                  # note: consider adding a minimum conductance to account for local diffusion
            } else if(input$cs_conductance_layer == "suitability") {
                  conductance <- sqrt(electrodes()$sources * electrodes()$destinations)
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

      # vertical current loss
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

      connectivity <- reactive({

            grounds <- (electrodes()$destinations + leakage()) %>% setNames("grounds")

            # write rasters to disk as ASCII grids
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
            current <- rast("circuitscape/result_curmap.asc") %>% setNames("current_flow")
            voltage <- setValues(current, voltage) %>% setNames("voltage")
            curr_ground <- grounds %>% "*"(voltage) %>% setNames("current_ground")
            curr_dest <- electrodes()$destinations %>% "*"(voltage) %>% setNames("current_detinations")
            curr_loss <- leakage() %>% "*"(voltage) %>% setNames("current_loss")

            curr_dir <- gradient_bearing(-voltage) %>% setNames("current_direction")

            # Clean temp files
            # file.remove(list.files("circuitscape", full.names = T))

            c(current + 0, # force minmax calc
              voltage, curr_dest, curr_loss, resistance()$conductance, leakage(), curr_dir) %>%
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

                  # add opacity slider
                  addControl(
                        sliderInput("rast_opacity", "Raster opacity", 0, 1, value = .8, step = .01, width = "150px"),
                        position = "bottomright"
                  )
      })

      ### Focal site -----------------------------


      ### Raster -------------------------------
      display_layer <- reactive({
            req(nzchar(input$species), aim())  # only run once we have data

            if(input$stat == "ESS"){
                  r <- switch(input$time,
                              "delta" = 2 / (1/aim()$fut$ess + 1/aim()$bsl$ess), # harmonic mean
                              "baseline" = aim()$bsl$ess,
                              "future" = aim()$fut$ess)
            }else if(input$stat %in% c("current animation", "current field")){
                  # show FIA delta raster for these layers
                  r <- aim()$fut[[var()]] - aim()$bsl[[var()]]
            }else if(input$stat %in% cs_stat_names){
                  st <- cs_stats[cs_stat_names == input$stat]
                  r <- connectivity()[[st]]
            }else{ # input$stat == "FIA variable"
                  r <- switch(input$time,
                              "delta" = aim()$fut[[var()]] - aim()$bsl[[var()]],
                              "baseline" = aim()$bsl[[var()]],
                              "future" = aim()$fut[[var()]])
            }

            # if(input$shrink && input$stat == "FIA variable"){
            #       # works for delta, and for estimate if we assume that no nearby samples means no forest
            #       ess <- switch(input$time,
            #                     delta = 2 / (1/aim()$fut$ess + 1/aim()$bsl$ess), # harmonic mean
            #                     baseline = aim()$bsl$ess,
            #                     future = aim()$fut$ess)
            #       k <- input$ess50 # ESS at which the effect is shrunk by 50% (default was 100)
            #       shrinkage <- ess / (ess + k)
            #       #prior_mean <- 0 # shrink toward this value
            #       r <- r * shrinkage #+ shrinkage * prior_mean
            # }

            r
      })

      observe({
            req(display_layer())

            proxy <- leafletProxy("map")

            ## clear old layers ##
            proxy %>%
                  clearImages() %>%
                  removeControl("raster_legend")

            ## add raster ##

            ## palettes ##
            mm <- minmax(display_layer())
            if(input$time == "delta" && input$stat == "FIA variable" |
               input$stat %in% c("current animation", "current field")){
                  # units are deltas
                  dom <- max(abs(mm)) * c(-1, 1)
                  rast_pal <- colorNumeric(
                        palette = c("darkred", "orange", "gray90", "dodgerblue", "darkblue"),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "ESS"){
                  # units are sample sizes
                  dom <- c(0, mm[2])
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkblue"),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "FIA variable") {
                  # baseline or future estimate: domain is unit except for total BA
                  # dom <- if(var() %in% c("batot", "baplot")) c(0, mm[2]) else 0:1
                  dom <- c(0, mm[2])
                  rast_pal <- colorNumeric(
                        palette = c("gray90", green),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "current direction") {
                  dom <- c(0, 360)
                  rast_pal <- colorNumeric(
                        palette = rainbow(10),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else {
                  # other connectivity variable
                  dom <- c(0, mm[2])
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkorchid4"),
                        domain = dom,
                        na.color = "transparent"
                  )
            }

            ## legend title ##
            if(input$stat == "ESS"){
                  title <- input$stat
            }else{
                  title <- paste0(input$species, ":<br/>", input$time, " ",
                                  tolower(sub(" ", "<br/>", input$fia_var)))
                  if(input$fia_var %in% c("total basal area", "total basal area ALL species")){
                        title <- paste0(title, "<br/>(m2 BA / plot)")
                  }
                  if(input$stat %in% cs_stat_names){
                        title <- paste0(title, " --<br/>", input$stat)
                  }
            }

            proxy <- proxy %>%
                  addRasterImage(
                        display_layer(),
                        colors  = rast_pal,
                        opacity = input$rast_opacity,
                        # project = FALSE,
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
                  clearMarkers() %>%
                  clearGroup("FIA plots") %>%
                  clearGroup("current_flow") %>%
                  clearGroup("current_field") %>%
                  removeControl("point_legend")

            if(input$vector %in% c("all FIA plots", "occupied FIA plots")){

                  pd <- plot_occ()
                  pd$display <- pd[[var()]]
                  if(input$vector == "occupied FIA plots") pd <- filter(pd, display > 0)

                  dom <- if(var() %in% c("batot", "baplot")) c(0, max(pd$display)) else 0:1
                  plot_pal <- colorNumeric(palette = c("gray90", green),
                                           domain = dom)
                  breaks <- pretty(seq(dom[1], dom[2], length.out = 100), n = 5)

                  proxy <- proxy %>%
                        addCircles(
                              data = filter(pd, samp_method_cd == 2),
                              lng = ~x,
                              lat = ~y,
                              radius = point_radius/3, # meters
                              stroke = FALSE, fillColor = "black",
                              group = "FIA plots"
                        ) %>%
                        addCircles(
                              data = filter(pd, samp_method_cd == 1),
                              lng = ~x,
                              lat = ~y,
                              radius = point_radius, # meters
                              fillColor = ~plot_pal(display),
                              fillOpacity = 1,
                              stroke = TRUE, color = "black", weight = 1,
                              group = "FIA plots"
                        )

                  proxy <- proxy %>%
                        addControl(
                              circleLegend(
                                    title  = HTML(paste(input$fia_var, "(FIA plots, baseline)", sep = "<br/>")),
                                    labels = breaks,
                                    colors = plot_pal(breaks),
                                    radius = 6
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
                                    velocityScale = 0.25,
                                    particleMultiplier = 0.001
                              ),
                              group = "current_flow"
                        )
            }

            if(input$vector == "current field"){
                  mag <- connectivity()[["current_flow"]]
                  mag <- mag / minmax(mag)[2]
                  spokes <- make_spokes(mag, connectivity()[["current_direction"]],
                                        subsample = 11 - input$clim_agg,
                                        arrow_scale = .1)

                  pl <- colorNumeric(palette = c("gray90", "black"), domain = sqrt(spokes$lines$mag))

                  proxy <- proxy %>%
                        addGlPolylines(data = spokes$lines, color = pl(sqrt(spokes$lines$mag)), weight = .25,
                                       group = "current_field") %>%
                        addCircles(
                              data = spokes$points %>% bind_cols(st_coordinates(spokes$points)),
                              lng = ~X, lat = ~Y,
                              radius = point_radius/3, fillOpacity = 1,
                              stroke = FALSE, fillColor = pl(sqrt(spokes$points$mag)),
                              group = "current_field"
                        )
            }


            # ## add layer controls ##
            # proxy <- proxy %>%
            #       addLayersControl(
            #             overlayGroups = c("AIM results", "FIA plots"),
            #             options = layersControlOptions(collapsed = FALSE)
            #       )
      })


      output$plot1 <- renderPlot({
            ggplot(mtcars, aes(mpg, cyl)) + geom_point()
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

library(shiny)
library(leaflet)
library(ggplot2)
library(terra)
library(analogs)
library(dplyr)
library(readr)
library(bslib)
library(shinybusy)


clim_vars <- function() c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")
bbox_res <- 1
clim_bbox <- as.vector(ext(clim_rast <- rast("../data/clim_quantized.tif")))
clim_bbox <- c(floor(clim_bbox[c(1, 3)] / bbox_res) * bbox_res, ceiling(clim_bbox[c(2, 4)] / bbox_res) * bbox_res)

spp <- sort(unique(read_csv("../data/fia.csv")$species))

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
cs_stat_names <- c("current flow", "current on target", "current loss", "resistance", "loss", "voltage")
cs_stats <- c("current_flow", "current_detinations", "current_loss", "resistance", "loss", "voltage")


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
                                       # selected = "Pinus edulis"),
                                       selected = "Pinus contorta"),
                        selectInput("fia_var", "FIA variable",
                                    c("proportion basal area", "total basal area", "presence probability",
                                      "total basal area ALL species")),

                        # p("Specify geographic bounding box"),
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
                        "Model outputs",
                        p("Choose results to display or export"),
                        selectInput("time", "Timeframe", c("baseline", "future", "delta")),
                        selectInput("stat", "Statistic", c("estimate", "ESS", cs_stat_names)),
                        checkboxInput("shrink", "Apply ESS shrinkage", value = FALSE),
                        numericInput("ess50", "Shrinkage rate", 50, min = 5, max = NA, step = 1)
                  ),

                  # accordion_panel(
                  #       "Connectivity",
                  #       p("This module uses Circuitscape to model connectivity between baseline and future surfaces.",
                  #         style = "font-size: 14px; color: #555;"),
                  #       sliderInput("cs_conductance", "log10 conductance (1/resitance)", -2, 2, value = 0),
                  #       sliderInput("cs_loss", "log10 conductance leakage", -2, 2, value = -1)#,
                  #       # actionButton("cs_click", "Run Circuitscape")
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


# Server ---------------------------------------------------

server <- function(input, output, session) {

      point_radius <- 1000
      green <- "darkgreen"


      ## Data setup -----------------------

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


      ## Snap domain to species range -----------------------

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





      ## Base map -----------------------------------------------
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
                        sliderInput("rast_opacity", "Raster opacity",
                                    0, 1, value = .8, step = .01, width = "150px"),
                        position = "bottomright"
                  )
      })


      ## AIM computation -----------------------------------------
      aim <- reactive({
            req(nzchar(input$species), plot_occ())
# if(input$bbox_x[1] < -120) browser()
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

            aim_bsl <- suppressWarnings(analogs::analog_impact(
                  x = clim_hst, pool = pool, max_clim = max_clim, max_geog = max_geog,
                  weight = "gaussian_joint", theta = theta, values = values))

            aim_fut <- suppressWarnings(analogs::analog_impact(
                  x = clim_fut, pool = pool, max_clim = max_clim, max_geog = max_geog,
                  weight = "gaussian_joint", theta = theta, values = values))

            # aim_bsl[is.na(aim_bsl)] <- 0
            # aim_fut[is.na(aim_fut)] <- 0

            names(aim_bsl) <- names(aim_fut) <- gsub("weighted_mean_", "", names(aim_fut))

            list(fut = aim_fut,
                 bsl = aim_bsl)
      })

      ## Circuitscape -------------------------------------------------
      connectivity <- reactive({

            # data setup
            var <- switch(input$fia_var,
                          "proportion basal area" = "baprop",
                          "total basal area" = "batot",
                          "presence probability" = "pres",
                          "total basal area ALL species" = "baplot")

            sources <- aim()$bsl[[var]] %>% setNames("sources")
            land_mask <- sources
            destinations <- aim()$fut[[var]] %>% setNames("destinations")
            sources[is.na(sources[])] <- 0
            destinations[is.na(destinations[])] <- 0

            v <- 1 # v <- 10^input$cs_conductance
            conductance <- setValues(sources, v) %>% setNames("conductance")
            resistance <- (1/conductance) %>% setNames("resistance")

            # survival proportional to the (geometric) mean of baseline and future suitabilities
            # loss <- setValues(sources, 10^input$cs_loss) %>% setNames("loss")
            loss <- (1 - sqrt(sources * destinations)) %>% setNames("loss")

            grounds <- (destinations + loss) %>% setNames("grounds")

            # write rasters to disk as ASCII grids
            writeRaster(resistance, "circuitscape/resistance.asc", overwrite = TRUE)
            writeRaster(sources, "circuitscape/sources.asc", overwrite = TRUE)
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
            curr_dest <- destinations %>% "*"(voltage) %>% setNames("current_detinations")
            curr_loss <- loss %>% "*"(voltage) %>% setNames("current_loss")

            c(current + 0, # force minmax calc
              voltage, curr_dest, curr_loss, resistance, loss) %>%
                  mask(land_mask)
      })


      ## Raster for display/download --------------------------------------------
      display_layer <- reactive({
            req(nzchar(input$species), aim())  # only run once we have data

            var <- switch(input$fia_var,
                          "proportion basal area" = "baprop",
                          "total basal area" = "batot",
                          "presence probability" = "pres",
                          "total basal area ALL species" = "baplot")

            if(input$stat == "ESS"){
                  r <- switch(input$time,
                              "delta" = 2 / (1/aim()$fut$ess + 1/aim()$bsl$ess), # harmonic mean
                              "baseline" = aim()$bsl$ess,
                              "future" = aim()$fut$ess)
            }else if(input$stat %in% cs_stat_names){
                  # browser()
                  st <- cs_stats[cs_stat_names == input$stat]
                  r <- connectivity()[[st]]
            }else{ # input$stat == "estimate"
                  r <- switch(input$time,
                              "delta" = aim()$fut[[var]] - aim()$bsl[[var]],
                              "baseline" = aim()$bsl[[var]],
                              "future" = aim()$fut[[var]])
            }

            if(input$shrink && input$stat == "estimate"){
                  # works for delta, and for estimate if we assume that no nearby samples means no forest
                  ess <- switch(input$time,
                                delta = 2 / (1/aim()$fut$ess + 1/aim()$bsl$ess), # harmonic mean
                                baseline = aim()$bsl$ess,
                                future = aim()$fut$ess)
                  k <- input$ess50 # ESS at which the effect is shrunk by 50% (default was 100)
                  shrinkage <- ess / (ess + k)
                  #prior_mean <- 0 # shrink toward this value
                  r <- r * shrinkage #+ shrinkage * prior_mean
            }

            r
      })


      ## Add data to map -------------------------------
      observe({
            req(plot_occ(), display_layer())  # only run once we have data

            proxy <- leafletProxy("map")

            ## clear old layers ##
            proxy %>%
                  clearMarkers() %>%
                  clearImages() %>%
                  removeControl("point_legend") %>%
                  removeControl("raster_legend")


            ## add raster ##

            var <- switch(input$fia_var,
                          "proportion basal area" = "baprop",
                          "total basal area" = "batot",
                          "presence probability" = "pres",
                          "total basal area ALL species" = "baplot")
            if(input$species == "Quercus agrifolia") browser()
            ## palettes ##
            mm <- minmax(display_layer())
            if(input$time == "delta" && input$stat == "estimate"){
                  # units are deltas
                  dom <- max(abs(mm)) * c(-1, 1)
                  rast_pal <- colorNumeric(
                        palette = c("darkred", "orange", "gray90", "dodgerblue", "darkblue"),
                        domain  = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "ESS"){
                  # units are sample sizes
                  dom <- c(0, mm[2])
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkblue"),
                        domain  = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "estimate") {
                  # baseline or future estimate: domain is unit except for total BA
                  dom <- if(var %in% c("batot", "baplot")) c(0, mm[2]) else 0:1
                  rast_pal <- colorNumeric(
                        palette = c("gray90", green),
                        domain  = dom,
                        na.color = "transparent"
                  )
            } else {
                  # connectivity variable
                  # dom <- c(0, mm[2])
                  dom <- c(mm[1], mm[2])
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkorchid4"),
                        domain  = dom,
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


            ## add FIA points ##
            pd <- plot_occ()
            pd$display <- pd[[var]]

            dom <- if(var %in% c("batot", "baplot")) c(0, max(pd$display)) else 0:1
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

            ## add layer controls ##
            proxy <- proxy %>%
                  addLayersControl(
                        overlayGroups = c("AIM results", "FIA plots"),
                        options = layersControlOptions(collapsed = FALSE)
                  )
      })


      ## Download handler --------------------------------------
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

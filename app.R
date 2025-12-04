library(shiny)
library(leaflet)
library(ggplot2)
library(terra)
library(analogs)
library(dplyr)
library(readr)


# helper functions ---------------------

pareto_optimal <- function(
            x # a two-column matrix of clim_dist and geog_dist values for a focal site
){
      x <- as.matrix(x)
      idx <- 1:nrow(x)
      x <- cbind(x, idx)

      frontier_idx <- vector()
      f <- x

      while(nrow(f) > 0){
            i <- which.min(f[,2])
            frontier_idx <- c(frontier_idx, f[i,3])
            f <- f[f[,1] < f[i,1], , drop = F]
      }

      idx %in% frontier_idx
}

on_pareto_front <- function(x) {
      x <- as.matrix(x)
      stopifnot(ncol(x) == 2)

      # Columns: clim_dist, geog_dist
      # Order by geog_dist (col 2), then clim_dist (col 1)
      o <- order(x[, 2], x[, 1])
      xs <- x[o, , drop = FALSE]

      best_clim <- Inf
      on_front_sorted <- logical(nrow(xs))

      for (k in seq_len(nrow(xs))) {
            clim_k <- xs[k, 1]
            if (clim_k <= best_clim) {
                  best_clim <- clim_k
                  on_front_sorted[k] <- TRUE
            }
      }

      # Map back to original order
      on_front <- logical(nrow(x))
      on_front[o] <- on_front_sorted
      on_front
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

clim_vars <- function() c("tmean", "ppt", "aet", "cwd")

# set_up_data <- function(){
#       if(!dir.exists("data")) dir.create("data")
#       if(!file.exists("data/clim.tif")) download.file("", "data/clim.tif")
#       if(!file.exists("data/fia.csv")) download.file("", "data/fia.csv")
# }



# UI ---------------------------------------------------

ui <- fluidPage(
      tags$style(HTML("
    #map .leaflet-tile {
      filter: grayscale(100%) brightness(100%) opacity(35%);
      -webkit-filter: grayscale(100%) brightness(100%) opacity(35%);
    }
  ")),


      titlePanel("FIA climate analogs"),

      fluidRow(
            column(width = 2,
                   # selectInput("type", "Analog type",
                   #             c("forward (outbound)", "reverse (inbound)", "baseline")),
                   selectInput("mode", "Query type",
                               c("Nearest geographic neighbors", "Best climate matches", "Pareto front"))),

            column(width = 2,
                   selectInput("v1", "Climate variable", clim_vars(), "cwd"),
                   selectInput("v2", "Second climate variable", clim_vars(), "aet")),

            column(width = 2, sliderInput("max_clim", "Max climate difference (sigma)", .1, 2, value = .5)),
            column(width = 2, sliderInput("max_geog", "Max distance (km)", 50, 1000, value = 300)),
            column(width = 2, sliderInput("max_k", "Max analog matches", 1, 100, value = 50)),

            column(width = 2, selectInput("plot_weight", "Weight",
                                          c("inverse geographic distance", "inverse climate distance", "uniform")))

      ),

      fluidRow(
            column(
                  width = 4,
                  h4("Geography (click map to select focal site)"),
                  leafletOutput("map", height = "600px")
            ),
            column(
                  width = 3,
                  h4("Climate"),
                  plotOutput("analog_plot", height = "600px")
            ),
            column(
                  width = 3,
                  h4("Analog similarity"),
                  plotOutput("dist_plot", height = "600px")
            ),
            column(
                  width = 2,
                  h4("Forest composition"),
                  plotOutput("forest_plot", height = "600px")
            )
      )
)



# Server ---------------------------------------------------

server <- function(input, output, session) {

      point_radius <- 1500

      type_pal <- data.frame(type = c("baseline", "forward", "reverse"),
                             color = c("gold", "blue", "red"))

      green <- scales::muted("darkgreen", l = 15)



      ## Reactive containers for data -----------------------------------
      data <- reactiveValues(
            climate_raster = NULL,  # base / current climate raster
            fia_points = NULL,  # base / current FIA plot data
            clim_scales = NULL,
            fia_idx_hst = NULL,
            fia_idx_fut = NULL,
            focal_site = NULL, # data for clicked location
            fia_analogs = NULL, # results of FIA analog query for clicked location
            fia_candidates = NULL
      )


      ## Load data -----------------------------------------
      observe({
            # set_up_data()
            data$climate_raster <- rast("data/clim.tif")
            data$clim_scales <- readRDS("data/clim_scales.rds")
            data$fia_points <- read_csv("data/fia.csv", show_col_types = FALSE)
      })


      ## Base map -----------------------------------------------
      output$map <- renderLeaflet({
            req(data$fia_points)

            leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
                  addProviderTiles(providers$OpenTopoMap,
                                   options = providerTileOptions(className = "grayscale")) %>%

                  # zoom to FIA extent
                  fitBounds(lng1 = min(data$fia_points$x),
                            lat1 = min(data$fia_points$y),
                            lng2 = max(data$fia_points$x),
                            lat2 = max(data$fia_points$y))
      })


      ## Draw FIA points on the map -------------------------------
      observe({
            req(data$fia_points)  # only run once we have data

            proxy <- leafletProxy("map")

            # Clear old layers
            proxy %>%
                  clearMarkers()


            # add FIA points
            points <- data$fia_points %>% select(x, y) %>% distinct()
            if (!is.null(data$fia_points)) {
                  proxy <- proxy %>%
                        addCircles(
                              data = points,
                              lng = ~x,
                              lat = ~y,
                              radius = point_radius, # meters
                              fillColor = green,
                              fillOpacity = 0.4,
                              stroke = FALSE
                        )

                  proxy <- proxy %>%
                        addLegend(
                              position = "topright",
                              pal = colorFactor(green, "FIA plots"),
                              values = "FIA plots"
                        )
            }
      })


      vars <- reactive({c(input$v1, input$v2)})

      ## Build search indices for focal vars ----------------------
      observe({
            data$fia_idx_hst <- build_analog_index(
                  data$fia_points[, c("x", "y", paste0(input$v1, 1), paste0(input$v2, 1))])
            data$fia_idx_fut <- build_analog_index(
                  data$fia_points[, c("x", "y", paste0(input$v1, 2), paste0(input$v2, 2))])
      })

      ## Analog query on click or slider change ---------------------------------
      observe({
            req(input$map_click)
            lat <- input$map_click$lat
            lon <- input$map_click$lng
            # lon <- -116
            # lat <- 45

            max_clim <- input$max_clim
            max_geog <- input$max_geog
            mode <- switch(input$mode,
                           "Nearest geographic neighbors" = "knn_geog",
                           "Best climate matches" = "knn_clim",
                           "Pareto front" = "all",
                           stop("invalid mode"))
            max_k <- if(mode != "all") input$max_k else NULL

            # focal site climate
            focal_clim <- terra::extract(data$climate_raster, matrix(c(lon, lat), nrow = 1))
            focal_hst <- cbind(x = lon, y = lat, focal_clim[, paste0(vars(), 1)])
            focal_fut <- cbind(x = lon, y = lat, focal_clim[, paste0(vars(), 2)])

            # store clicked location including climate
            data$focal_site <- data.frame(
                  lon = lon,
                  lat = lat,
                  value = 1  # placeholder column
            ) %>% bind_cols(focal_clim)

            # FIA analogs
            fia_query <- function(x, pool){
                  a <- analogs::analog_search(x, pool, mode = mode, k = max_k,
                                              max_clim = max_clim, max_geog = max_geog)
                  if(input$mode == "Pareto front"){
                        # return the pareto-optimal subset of "all" analogs passing filters
                        on_front <- a %>% select(clim_dist, geog_dist) %>% on_pareto_front()
                        return(a[on_front,])
                  }else{
                        return(a)
                  }
            }
            fia_analogs_bsl <- fia_query(focal_hst, data$fia_idx_hst)
            fia_analogs_fwd <- fia_query(focal_hst, data$fia_idx_fut)
            fia_analogs_rev <- fia_query(focal_fut, data$fia_idx_hst)

            fia_analogs <- bind_rows(
                  data$fia_points[fia_analogs_bsl$analog_index, ] %>% mutate(type = "baseline"),
                  data$fia_points[fia_analogs_fwd$analog_index, ] %>% mutate(type = "forward"),
                  data$fia_points[fia_analogs_rev$analog_index, ] %>% mutate(type = "reverse")
            )

            dists <- bind_rows(fia_analogs_bsl[,c("clim_dist", "geog_dist")],
                               fia_analogs_fwd[,c("clim_dist", "geog_dist")],
                               fia_analogs_rev[,c("clim_dist", "geog_dist")])
            fia_analogs <- bind_cols(fia_analogs, dists)


            a <- fia_analogs %>%
                  # rename(t = t1, p = p1) %>%
                  mutate(v1 = fia_analogs[[paste0(input$v1, 1)]],
                         v2 = fia_analogs[[paste0(input$v2, 1)]]) %>%
                  select(-(tmean1:cwd2)) %>%
                  mutate(era = "baseline")
            a <- fia_analogs %>%
                  # rename(t = t2, p = p2) %>% select(-t1, -p1) %>%
                  mutate(v1 = fia_analogs[[paste0(input$v1, 2)]],
                         v2 = fia_analogs[[paste0(input$v2, 2)]]) %>%
                  select(-(tmean1:cwd2)) %>%
                  mutate(era = "future") %>%
                  bind_rows(a)
            a <- a %>%
                  filter(type == "forward" & era == "future" |
                               type == "reverse" & era == "baseline" |
                               type == "baseline" & era == "baseline") %>%
                  mutate(type = factor(type, levels = c("reverse", "baseline", "forward")))

            # store FIA analog results
            data$fia_analogs <- a
      })


      ## Candidate FIA plots within max_geog radius -----------------
      observe({
            req(data$focal_site)
            set.seed(123)
            cand <- data$fia_points %>%
                  # select(x, y, t1, p1, t2, p2) %>%
                  select(x, y, tmean1:cwd2) %>%
                  distinct() %>%
                  mutate(geog_dist = geosphere::distHaversine(
                        cbind(x, y),
                        cbind(data$focal_site$lon, data$focal_site$lat)),
                        geog_dist = geog_dist / 1000
                  ) %>%
                  filter(geog_dist <= input$max_geog) %>%
                  sample_n(min(nrow(.), 10000))

            cand <- cand %>%
                  mutate(clim_dist = sqrt((
                        cand[[paste0(input$v1, 2)]] - data$focal_site[[paste0(input$v1, 1)]])^2 +
                              (cand[[paste0(input$v2, 2)]] - data$focal_site[[paste0(input$v2, 1)]])^2))

            data$fia_candidates <- cand
      })


      ## Map update based on query results ------------------------
      observe({
            req(data$focal_site)
            req(data$fia_analogs)

            pal <- colorFactor(type_pal$color, domain = type_pal$type)

            proxy <- leafletProxy("map")

            # Clear any previous “result” markers first
            proxy <- proxy %>%
                  clearGroup("focal_site") %>%
                  clearGroup("max_geog") %>%
                  clearGroup("fia_analogs") %>%
                  clearGroup("analog_links") %>%
                  removeControl("analog_legend")

            if(is.finite(input$max_geog)){
                  poly <- geo_circle(data$focal_site[,c("lon", "lat")], input$max_geog)
                  proxy <- proxy %>%
                        addPolygons(lng <- poly[,1], lat <- poly[,2],
                                    color = "black", weight = 2, fill = FALSE,
                                    group = "max_geog")
            }

            # Add line segments
            segments <- data$fia_analogs %>%
                  filter(type != "baseline") %>%
                  mutate(lon0 = data$focal_site$lon[1],
                         lat0 = data$focal_site$lat[1])
            for (i in seq_len(nrow(segments))) {
                  proxy <- proxy %>% addPolylines(
                        lng = c(segments$lon0[i], segments$x[i]),
                        lat = c(segments$lat0[i], segments$y[i]),
                        color = pal(segments$type[i]),
                        weight = .5,
                        opacity = 1,
                        group = "analog_links"
                  )
            }

            # Add markers
            proxy <- proxy %>%
                  addCircleMarkers(
                        data = data$focal_site,
                        lng = ~lon,
                        lat = ~lat,
                        radius = 10,
                        color = "black",
                        fillOpacity = 1,
                        stroke = TRUE,
                        weight = 2,
                        group = "focal_site"
                  ) %>%
                  addCircleMarkers(
                        data = data$fia_analogs,
                        lng = ~x,
                        lat = ~y,
                        radius = 6,
                        color = "black",
                        fillColor = ~pal(type),
                        fillOpacity = 1,
                        stroke = TRUE,
                        weight = 2,
                        group = "fia_analogs"
                  )

            # Add legend
            proxy <- proxy %>%
                  addLegend(
                        position = "topright",
                        pal = pal,
                        values = data$fia_analogs$type,
                        title = "Analog type",
                        opacity  = 1,
                        layerId = "analog_legend"
                  )

            proxy
      })


      ## plot: climate space --------------------------------
      output$analog_plot <- renderPlot({
            req(data$focal_site, data$fia_analogs)

            v11 <- paste0(input$v1, 1)
            v12 <- paste0(input$v1, 2)
            v21 <- paste0(input$v2, 1)
            v22 <- paste0(input$v2, 2)

            mean <- c(data$clim_scales[[input$v1]]$mean,
                      data$clim_scales[[input$v2]]$mean)
            sd <- c(data$clim_scales[[input$v1]]$sd,
                    data$clim_scales[[input$v2]]$sd)

            a <- data$fia_analogs %>%
                  select(type, era, v1, v2) %>%
                  distinct() %>%
                  mutate(v1 = v1 * sd[1] + mean[1],
                         v2 = v2 * sd[2] + mean[2])


            p <- data$focal_site %>%
                  select(v1 = all_of(v11),
                         v2 = all_of(v21)) %>%
                  mutate(era = "baseline")
            p <- data$focal_site %>%
                  select(v1 = all_of(v12),
                         v2 = all_of(v22)) %>%
                  mutate(era = "future") %>%
                  bind_rows(p)
            p <- p %>%
                  mutate(v1 = v1 * sd[1] + mean[1],
                         v2 = v2 * sd[2] + mean[2])

            # search perimeters
            perims <- bind_rows(
                  clim_circle(c(data$focal_site[, v11], data$focal_site[, v21]),
                              input$max_clim, mean, sd) %>% mutate(era = "baseline"),
                  clim_circle(c(data$focal_site[, v12], data$focal_site[, v22]),
                              input$max_clim, mean, sd) %>% mutate(era = "future"))

            # background FIA plots
            bg <- data$fia_candidates
            bg <- bg %>%
                  mutate(v11 = bg[[v11]] * sd[1] + mean[1],
                         v21 = bg[[v21]] * sd[2] + mean[2],
                         v12 = bg[[v12]] * sd[1] + mean[1],
                         v22 = bg[[v22]] * sd[2] + mean[2])

            # build ggplot
            q <- .25
            ggplot() +
                  geom_density_2d_filled(data = bg %>% mutate(era = "baseline"),
                                         aes(v11, v21, fill = NULL, linetype = era),
                                         breaks = c(q, 1), contour_var = "ndensity",
                                         fill = green, color = green, alpha = .3, linewidth = .2) +
                  geom_density_2d_filled(data = bg %>% mutate(era = "future"),
                                         aes(v12, v22, fill = NULL, linetype = era),
                                         breaks = c(q, 1), contour_var = "ndensity",
                                         fill = green, color = green, alpha = .3, linewidth = .2) +
                  geom_polygon(data = perims,
                               aes(v1, v2, group = era, linetype = era),
                               color = "gray40", linewidth = .5, fill = NA, alpha = .5) +
                  geom_segment(data = filter(a, type == "forward"),
                               aes(xend = v1, yend = v2),
                               x = p$v1[p$era == "baseline"],
                               y = p$v2[p$era == "baseline"],
                               color = "blue") +
                  geom_segment(data = filter(a, type == "reverse"),
                               aes(x = v1, y = v2),
                               xend = p$v1[p$era == "future"],
                               yend = p$v2[p$era == "future"],
                               color = "red") +
                  geom_point(data = a,
                             aes(v1, v2, shape = era, fill = type),
                             size = 3) +
                  geom_point(data = p,
                             aes(v1, v2, shape = era),
                             fill = "black", color = "white", size = 5) +
                  scale_shape_manual(values = c(21, 24)) +
                  scale_fill_manual(values = type_pal$color,
                                    breaks = type_pal$type, guide = "none") +
                  theme_bw() +
                  theme(legend.position = "bottom") +
                  labs(shape = "Time period", linetype = "Time period",
                       x = input$v1,
                       y = input$v2)
      })


      ## plot: distances --------------------------------
      output$dist_plot <- renderPlot({
            req(data$fia_analogs)

            p <- ggplot(data$fia_analogs,
                        aes(geog_dist, clim_dist, fill = type, color = type)) +
                  geom_vline(xintercept = input$max_geog) +
                  geom_hline(yintercept = input$max_clim) +
                  geom_point(data = data$fia_candidates,
                             color = green, fill = NA,
                             size = 2, alpha = .5, shape = 16)

            if(input$mode == "Pareto front") p <- p + geom_line(color = "black", linewidth = 1) + geom_line(linewidth = .5)

            p + geom_point(shape = 21, color = "black", size = 3) +
                  scale_fill_manual(values = type_pal$color, breaks = type_pal$type, guide = "none") +
                  scale_color_manual(values = type_pal$color, breaks = type_pal$type, guide = "none") +
                  scale_x_continuous(expand = c(0, 0), limits = c(0, max(data$fia_analogs$geog_dist * 1.25))) +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, max(data$fia_analogs$clim_dist * 1.25))) +
                  scale_alpha_continuous(range = c(0, 1)) +
                  theme_bw() +
                  labs(y = "Climate difference from focal site (sigma)",
                       x = "Geographic distane to focal site (km)")
      })


      ## plot: forest attributes --------------------------------
      output$forest_plot <- renderPlot({
            req(data$fia_analogs)

            a <- data$fia_analogs

            a$weight <- switch(input$plot_weight,
                               "inverse geographic distance" = 1/a$clim_dist,
                               "inverse climate distance" = 1/a$geog_dist,
                               "uniform" = 1)

            a <- a %>%
                  mutate(species = paste(genus, species)) %>%
                  group_by(type, species) %>%
                  summarize(basal_area = weighted.mean(pi*(dia/2)^2, weight) / 1550) %>% # sq in to sq m
                  group_by(type) %>%
                  mutate(total = sum(basal_area)) %>%
                  arrange(desc(total))

            ggplot(a, aes(type, species, fill = type, size = basal_area)) +
                  geom_point(color = "black", shape = 21) +
                  scale_fill_manual(values = type_pal$color, breaks = type_pal$type, guide = "none") +
                  scale_size(range = c(0, 20),
                             limits = c(0, NA),
                             breaks = signif(max(a$basal_area) * c(.01, .1, 1), 1)) +
                  theme_bw() +
                  theme(legend.position = "bottom") +
                  labs(x = "Analog type",
                       y = NULL,
                       size = "Mean basal area,\nbaseline era (m^2)")
      })
}

shinyApp(ui, server)

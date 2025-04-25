####################
# Title: 4_figures.R
# Date Created: 4/17/2025
# Author: Chen Chen, Lara Schwarz and Tim B. Frankland
# Purpose: Plot figures for the long-term wildfire smoke on mortality KPSC project
####################

library(data.table)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

indir1 <- "" ## work directory

## set the stage and read in data for figure 2, 4, S3, and S4
#######
metrics <- c("Mean daily wildfire PM2.5",
             "Mean daily wildfire PM2.5 during peak week",
             "Number of days with wildfire PM2.5 > 0",
             "Number of weeks with wildfire PM2.5 > 5",
             "Number of smoke waves",
             "Mean daily non-WF PM2.5")
exposure_labels <- c(expression("Mean daily wildfire PM"[2.5]), 
                     expression("Mean daily wildfire PM"[2.5]~"during peak week"),
                     expression("Number of days with wildfire PM"[2.5]~">0"),
                     expression("Number of weeks with wildfire PM"[2.5]~">5"),
                     expression("Number of smoke waves"),
                     expression("Mean daily non-WF"[2.5]))

## Effect estimates in table 2 were gathered from output files and transformed so that
## they represent odds ratio per 5th to 95th percentile change in exposure.
out <- fread(file.path(indir1, "results", "Table S2.csv"))
out$exposure <- factor(out$exposure, levels = metrics)

colrs <- brewer.pal(n = 6, name = "Set2")
dodge <- position_dodge(width=0.6)
#######

## Figure 2
#######
nn <- 5
bar <- out[out$exposure!=metrics[6] & out$type=="all", ]
loc_con <- 1:2 ## used in two panels
loc_day <- 3:5 ## used in two panels
or90_ylim <- c(min(bar[1:nn, .(or90, or90_ul, or90_ll)]) * 0.99, 
               max(bar[1:nn, .(or90, or90_ul, or90_ll)]) * 1.01)
p1 <- ggplot(bar[bar$metric_type=="concentration",], 
             aes(y=or90, x=exposure, ymax=or90_ul, ymin=or90_ll, col=exposure)) + 
  geom_hline(yintercept = 1, col="darkgrey", linewidth=1.2) +
  geom_point(position=dodge, size=3, shape = 16) + geom_errorbar(position=dodge, width=0.2, linewidth=1.2) + 
  scale_color_manual(values = colrs[loc_con], labels = exposure_labels[loc_con]
  ) +
  labs(title="Concentration metrics", x="", 
       y = "", 
       col = "Exposure metric")+ 
  scale_y_continuous(trans = "log", limits = or90_ylim) +
  theme_bw()+
  theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.position = "none", legend.text.align = 0
  )
p2 <- ggplot(bar[bar$metric_type=="exceedance",], 
             aes(y=or90, x=exposure, ymax=or90_ul, ymin=or90_ll, col=exposure)) + 
  geom_hline(yintercept = 1, col="darkgrey", linewidth=1.2) +
  geom_point(position=dodge, size=3, shape = 17) + geom_errorbar(position=dodge, width=0.2, linewidth=1.2) + 
  # scale_y_continuous(trans = "log", breaks = c(1, 2, 5, 10, 20, 30, 50)) + 
  scale_color_manual(values = colrs[loc_day], labels = exposure_labels[loc_day]
  ) +
  labs(title="Exceedance metrics", x="", 
       y = "", 
       col = "Exposure metrics")+ 
  scale_y_continuous(trans = "log", limits = or90_ylim) +
  theme_bw()+
  theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.position = "none", legend.text.align = 0
  )
legend_f <- get_plot_component(
  ggplot(bar[bar$type=="all", ], 
         aes(y=or90, x=exposure, ymax=or90_ul, ymin=or90_ll, col=exposure, shape = exposure)) + 
    geom_hline(yintercept = 1, col="darkgrey", linewidth=1.2) +
    geom_point(position=dodge, size=3) + geom_errorbar(position=dodge, width=0.2, linewidth=1.2) + 
    scale_color_manual(values = colrs[c(loc_con, loc_day)], labels = exposure_labels[c(loc_con, loc_day)]
    ) +
    scale_shape_manual(values = c(16, 16, 17, 17, 17), labels = exposure_labels[c(loc_con, loc_day)]) +
    guides(color=guide_legend(nrow=3,byrow=T), shape = guide_legend(nrow=3,byrow=T)) +
    labs(title="", x="", 
         y = expression(""), 
         col = "Exposure metrics", 
         shape = "Exposure metrics")+ 
    theme_bw()+
    theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "bottom", legend.text.align = 0
    )
  , "guide-box", return_all = TRUE)[[3]]
ylab <- ggplot()+
  geom_text(aes(x=0,y=0),label="Odds ratio per 5th to 95th percentile\nincrease in exposure metric",
            size=6,angle=90)+
  theme_void()

tiff(file.path(indir1, "figures", "publication", "Figure2.tiff"), compression ="zip",
# png(file.path(indir1, "figures", "publication", "Figure2.png"),
     width=9,height=6,units="in", res = 600, bg="white",
     family="sans")
p <- plot_grid(p1, p2, labels=NULL, ncol = 2)
p.lab <- plot_grid(ylab, p, nrow=1, rel_widths = c(0.07,1))
plot_grid(p.lab, legend_f, ncol=1, rel_heights = c(1, 0.2))
dev.off()
#######

## Figure 4
#######
nn <- 5
locs <- 1:5 ## used in one panel
bar <- out[out$exposure!=metrics[6], ]
bar <- bar[!bar$subgroup %in% c("other-known", "other-unknown"), ]

## clean up labels
bar$label<- factor(bar$subgroup,
                   levels = c("all", "age74below", "age75up", 
                              "asian", "black", "hispanic", "white", "other",
                              "female", "male", "pov15below", "pov15up"
                   ),
                   labels = c("All", "Age < 75", "Age >or= 75", 
                              "Asian", "Black", "Hispanic", "White", "Other",
                              "Female", "Male", "Pov < 15%", "Pov >or= 15%"
                   ))
bar$label <- as.character(bar$label)
bar$cq_p <- as.numeric(bar$cq_p)
temp <- bar
temp$star <- ifelse(temp$cq_p<0.1, "*", "")
temp$star[is.na(temp$star)] <- ""

bar$type <- factor(bar$type, levels = c("all", "sex", "age", "race", "poverty"),
                   labels = c("All", "Sex", "Age", "Race", "Poverty"))
or90_ylim <- c(min(bar[1:nn, .(or90, or90_ul, or90_ll)]) * 0.99, max(bar[1:nn, .(or90, or90_ul, or90_ll)]) * 1.01)

tiff(file.path(indir1, "figures", "publication", "Figure4.tiff"), compression ="zip",
# png(file.path(indir1, "figures", "publication", "Figure4.png"),
     width=9,height=6,units="in", res = 600, bg="white",
     family="sans")
ggplot(temp[temp$subgroup!="other", ], 
       aes(y=or90, x=label, ymax=or90_ul, ymin=or90_ll, col=exposure, label=star, shape=exposure)) + 
  geom_hline(yintercept = 1, col="darkgrey", linewidth=1.2) +
  geom_point(position=dodge, size=3) + 
  geom_errorbar(position=dodge, width=0.2, linewidth=1.2) + 
  geom_text(position=dodge, hjust=0, vjust=-1, show.legend = FALSE, size=4) +
  facet_grid(.~type, scales = "free", space = "free_x") +
  scale_color_manual(values = colrs[locs], labels = exposure_labels[locs]
  ) +
  scale_shape_manual(values = c(16, 16, 17, 17, 17), labels = exposure_labels[locs]) +
  # scale_y_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log") +
  labs(title="", x="Subgroup", 
       y = "Odds ratio per 5th to 95th percentile\nincrease in exposure metric", 
       col = "Exposure metrics",
       shape = "Exposure metrics") + 
  guides(color=guide_legend(nrow=3,byrow=T),  shape = guide_legend(nrow=3,byrow=T)) +
  theme_bw()+
  theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "inside", 
        legend.position = "bottom", legend.text.align = 0
  )
dev.off()
#######

## Figure S3
#######
bar <- out[out$exposure == metrics[6], ]
bar <- bar[!bar$subgroup %in% c("other-known", "other-unknown"), ]
## clean up labels
bar$label<- factor(bar$subgroup,
                   levels = c("all", "age74below", "age75up", 
                              "asian", "black", "hispanic", "white", "other",
                              "female", "male", "pov15below", "pov15up"
                   ),
                   labels = c("All", "Age < 75", "Age >or= 75", 
                              "Asian", "Black", "Hispanic", "White", "Other",
                              "Female", "Male", "Pov < 15%", "Pov >or= 15%"
                   ))
bar$label <- as.character(bar$label)
bar$cq_p <- as.numeric(bar$cq_p)
bar$label <- ifelse(bar$cq_p<0.1, paste0(bar$label, "*"), bar$label) ## this does not work in colored figure
bar$label[is.na(bar$label)] <- "All"

bar$type <- factor(bar$type, levels = c("all", "sex", "age", "race", "poverty"),
                   labels = c("All", "Sex", "Age", "Race", "Poverty"))

png(file.path(indir1, "figures", "publication", "FigureS3.png"),
    width=6,height=4,units="in", res = 600, bg="white",
    family="sans")
ggplot(bar, aes(y=or90, x=label, ymax=or90_ul, ymin=or90_ll)) +
  geom_hline(yintercept = foo$or90[foo$type=="all"], col="lightgrey", linewidth=1, linetype = 2) +
  geom_hline(yintercept = 1, col="darkgrey", linewidth=1.2) +
  geom_point(position=dodge, size=3) + geom_errorbar(position=dodge, width=0.2, linewidth=1.2) +
  facet_grid(.~type, scales = "free", space = "free_x") +
  scale_y_continuous(trans = "log") +
  labs(title=exposure_labels[6], x="Subgroup",
       y = "Odds ratio per 5th to 95th percentile\nincrease in exposure metric") +
  theme_bw()+
  theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "inside"
  )
dev.off()
#######

## Figure S4
#######
foo.unknown <- out[out$type %in% c("all", "race"), ]
foo.unknown$subgroup <- factor(foo.unknown$subgroup, 
                               levels = c("all", "asian", "black", "hispanic", "white", 
                                          "other", "other-unknown", "other-known"),
                               labels = c("All", "Asian", "Black", "Hispanic", "White", 
                                          "Other", "Other-unknown", "Other-known"))

nn <- 5
for (i in 1:nn) {
  
  exposure_ <- metrics[i]
  foo <- foo.unknown[foo.unknown$exposure == exposure_, ]
  
  assign(paste0("por90_", i), 
         ggplot(foo, aes(y=or90, x=subgroup, ymax=or90_ul, ymin=or90_ll)) + 
           geom_hline(yintercept = foo$or90[foo$subgroup=="all"], col="lightgrey", linewidth=1, linetype = 2) +
           geom_hline(yintercept = 1, col="darkgrey", linewidth=1.2) +
           geom_point(position=dodge, size=3) + geom_errorbar(position=dodge, width=0.2, linewidth=1.2) + 
           scale_y_continuous(trans = "log1p") +
           labs(title=exposure_labels[i], x="", 
                y = "") + 
           theme_bw()+
           theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5),
                 axis.text.x = element_text(angle = 45, hjust = 1),
                 strip.background = element_blank(),
                 strip.placement = "inside"
           ))
}

xlab <- ggplot()+
  geom_text(aes(x=0,y=0),label="Subgroups",size=7)+
  theme_void()
ylab_or90 <- ggplot()+
  geom_text(aes(x=0,y=0),label="Odds ratio per 5th to 95th percentile increase in exposure metric",size=7,angle=90)+
  theme_void()

png(file.path(indir1, "figures", "publication", "figureS4.png"),
    width=6,height=15,units="in", res = 600, bg="white",
    family="sans")
p <- plot_grid(por90_1, por90_2, por90_3, por90_4, por90_5, labels="AUTO", nrow = nn)
p.lab <- plot_grid(ylab_or90, p, nrow=1, rel_widths = c(0.05,1))
plot_grid(p.lab, xlab, ncol=1, rel_heights = c(1, 0.03))
dev.off()

#######

#######
## figures below require reading in of other sets of data
#######

## Figure 1
#######
# Load and prepare wildfire smoke data
wf_data <- readRDS(file.path(indir1, "data", "r_kpsc_wfpm_average_rolling_seasonal_8724.RData"))

# Load CA census tracts (2010)
CA_census_tracts <- tracts(state = "CA", year = 2010, class = "sf")
kaiser_socal_counties <- c("06037", "06059", "06071", "06065", "06111", 
                           "06073", "06025", "06083", "06029", "06079")

# Filter to relevant counties
CA_county_outlines <- CA_census_tracts %>%
  filter(substr(GEOID10, 1, 5) %in% kaiser_socal_counties)

# Dissolve county boundaries into one
county_outlines_sf <- CA_county_outlines %>%
  st_union() %>%
  st_cast("POLYGON")

# Filter for years 2009–2019
wf_data_filtered <- wf_data %>%
  mutate(year = as.numeric(substr(date, 1, 4))) %>%
  filter(year >= 2009 & year < 2020)

# Define variables of interest
variables <- c("rolling_wf_average", "rolling_peak", "rolling_non_zero", 
               "rolling_freq", "rolling_smoke_waves", "rolling_nonwf_average")

# Average over time per tract
wf_data_avg <- wf_data_filtered %>%
  group_by(geoid) %>%
  summarise(across(all_of(variables), ~ mean(.x, na.rm = TRUE), .names = "{.col}"))

# Merge with spatial data
CA_census_tracts$geoid <- as.numeric(CA_census_tracts$GEOID10)
KPSC_wildfire_smoke <- merge(CA_census_tracts, wf_data_avg, by = "geoid")

# Create union of catchment area for possible plotting
catchment_area <- st_union(KPSC_wildfire_smoke)

# Create thematic maps

map_vars <- c("rolling_wf_average", "rolling_peak", "rolling_non_zero", 
              "rolling_freq", "rolling_smoke_waves")
legend_titles <- c(expression("Mean daily wildfire PM"[2.5]), 
                   expression("mean peak week"), 
                   expression("non-zero days"),
                   expression("weeks over five"), 
                   expression("smoke waves"))
colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")

maps <- lapply(1:5, function(i) {
  ggplot(KPSC_wildfire_smoke) +
    geom_sf(aes_string(fill = map_vars[i]), color = NA) +
    geom_sf(data = county_outlines_sf, fill = NA, color = "black", size = 0.8) +
    scale_fill_gradient(low = "white", high = colors[i]) +
    theme_minimal() +
    labs(fill = legend_titles[i]) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.margin = unit(rep(0.5, 4), "cm"),
      legend.position = "right",
      legend.box.spacing = unit(1, "cm")
    )
})

# Combine maps into a single plot
combined_maps <- wrap_plots(maps, ncol = 2, nrow = 3)

# Correlation heatmap
KPSC_corr_data <- wf_data_avg %>%
  select(all_of(map_vars))

corr_matrix <- cor(KPSC_corr_data, use = "complete.obs")
rownames(corr_matrix) <- legend_titles
colnames(corr_matrix) <- legend_titles

# Prepare for ggplot
corr_data <- melt(corr_matrix)
corr_data$color <- rep(colors, each = length(colors))


# Create parseable labels as strings
axis_labels <- c(
  "Mean~daily~wildfire~PM[2.5]", 
  "\"mean peak week\"",
  "\"non-zero days\"",
  "\"weeks over five\"",
  "\"smoke waves\""
)

# Make sure order matches levels in Var1 and Var2
corr_data$Var1 <- factor(corr_data$Var1, levels = legend_titles)
corr_data$Var2 <- factor(corr_data$Var2, levels = legend_titles)

# Plot
p1 <- ggplot(corr_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 5, color = "black") +
  scale_fill_gradient2(high = "#8e8e8e", low = "#FF6F61", mid = "white", midpoint = 0) +
  scale_x_discrete(labels = parse(text = axis_labels)) +
  scale_y_discrete(labels = parse(text = axis_labels)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, color = colors, face = "bold"),
    axis.text.y = element_text(size = 15, color = colors, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank()
  )

# Save output to PNG
png(file.path(indir1, "figure", "publication", "combined_plot_v2.png"), width = 14, height = 10, units = "in", res = 300)
wrap_plots(maps, ncol = 2, nrow = 3) + p1
dev.off()
#######

## Figure 2
#######
colrs <- brewer.pal(n = 5, name = "Set2")
colrs.ribbon <-  rgb(t(col2rgb(colrs)), alpha=100, maxColorValue = 255)

wf_pm_pred <- readRDS(file.path(indir1, "results", "wf_pm_pred_mort_dlmn.rds")) ### UPDATE ME
wf_pm_pred$col <- "wf_pm_pred"
p1 <- ggplot(wf_pm_pred, aes(x=mean_wf_pm, y=or, col=col))+  
  geom_ribbon(aes(ymin=or.ll, ymax=or.ul),fill=colrs.ribbon[1]) +  
  geom_line() + 
  geom_hline(yintercept = 1, col="black") +
  theme_minimal(base_size=12) + ### UPDATE ME
  scale_color_manual(values = colrs[1]) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        legend.position = "none", 
        panel.grid.minor=element_blank(), panel.background=element_blank()) +
  scale_x_continuous(expression(paste("Mean wildfire P",M[2][.][5])),expand=c(0,0)) + 
  scale_y_continuous("Odds ratio relative to zero exposure", trans="log"
                     , breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
  ) +
  geom_rug(sides="b")  

md_pw_pred <- readRDS(file.path(indir1, "results", "md_pw_pred_mort_dlmn.rds")) ### UPDATE ME
md_pw_pred$col <- "md_pw_pred"
p2 <- ggplot(md_pw_pred, aes(x=mean_daily_peak_week, y=or, col=col))+  
  geom_ribbon(aes(ymin=or.ll, ymax=or.ul),fill=colrs.ribbon[2]) +  
  geom_line() + 
  geom_hline(yintercept = 1, col="black") +
  theme_minimal(base_size=12) + ### UPDATE ME
  scale_color_manual(values = colrs[2]) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        legend.position = "none", 
        panel.grid.minor=element_blank(), panel.background=element_blank()) +
  scale_x_continuous("Mean peak week",expand=c(0,0)) + 
  scale_y_continuous("Odds ratio relative to zero exposure", trans="log"
                     , breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
  ) +
  geom_rug(sides="b")  

nzd_pred <- readRDS(file.path(indir1, "results", "nzd_pred_mort_dlmn.rds")) ### UPDATE ME
nzd_pred$col <- "non_zero_days"
p3 <- ggplot(nzd_pred, aes(x=non_zero_days, y=or, col=col))+  
  geom_ribbon(aes(ymin=or.ll, ymax=or.ul),fill=colrs.ribbon[3]) +  
  geom_line() + 
  geom_hline(yintercept = 1, col="black") +
  theme_minimal(base_size=12) + ### UPDATE ME
  scale_color_manual(values = colrs[3]) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        legend.position = "none", 
        panel.grid.minor=element_blank(), panel.background=element_blank()) +
  scale_x_continuous("Non-zero days",expand=c(0,0)) + 
  scale_y_continuous("Odds ratio relative to zero exposure", trans="log"
                     , breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
  ) +
  geom_rug(sides="b")  

wks_gt5_pred <- readRDS(file.path(indir1, "results", "wks_gt5_pred_mort_dlmn.rds")) ### UPDATE ME
wks_gt5_pred$col <- "weeks_gt_5"
p4 <- ggplot(wks_gt5_pred, aes(x=weeks_gt_5, y=or, col=col))+  
  geom_ribbon(aes(ymin=or.ll, ymax=or.ul),fill=colrs.ribbon[4]) +  
  geom_line() + 
  geom_hline(yintercept = 1, col="black") +
  theme_minimal(base_size=12) + ### UPDATE ME
  scale_color_manual(values = colrs[4]) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        legend.position = "none", 
        panel.grid.minor=element_blank(), panel.background=element_blank()) +
  scale_x_continuous("Weeks over five",expand=c(0,0)) + 
  scale_y_continuous("Odds ratio relative to zero exposure", trans="log"
                     , breaks = c(0.6, 1.0, 1.2, 1.4, 1.6)
  ) +
  geom_rug(sides="b") 

smk_wvs_pred <- readRDS(file.path(indir1, "results", "smk_wvs_pred_mort_dlmn.rds")) ### UPDATE ME
smk_wvs_pred$col <- "smoke_waves"
p5 <- ggplot(smk_wvs_pred, aes(x=smoke_waves, y=or, col=col))+  
  geom_ribbon(aes(ymin=or.ll, ymax=or.ul),fill=colrs.ribbon[5]) +  
  geom_line() + 
  geom_hline(yintercept = 1, col="black") +
  theme_minimal(base_size=12) + ### UPDATE ME
  scale_color_manual(values = colrs[5]) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        legend.position = "none", 
        panel.grid.minor=element_blank(), panel.background=element_blank()) +
  scale_x_continuous("Smoke waves",expand=c(0,0)) + 
  scale_y_continuous("Odds ratio relative to zero exposure", trans="log"
                     , breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
  ) +
  geom_rug(sides="b") 

png(file.path(indir1, "figures", "publication", "figure3.png"),  
    width=8.3, height=11.7,units="in", res = 600, bg="white",
    family="sans")
plot_grid(p1, p2, p3, p4, p5, labels="AUTO", nrow = 3, byrow = T)
dev.off()
#######

## Figure S2
#######
## single figure for non-wildfire pm
non_wf_pm_pred <- readRDS(file.path(indir1, "results", "non_wf_pm_pred_mort_dlmn.rds"))
spline_plot<-non_wf_pm_pred %>%                                                                                                                                      
  ggplot(aes(x=mean_non_wf_pm,y=rr)) +                                                                                                                              geom_ribbon(aes(ymin=ll,ymax=ul),fill="grey90") +  
  geom_line() + 
  theme_minimal(base_size=14) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), panel.background=element_blank()) +
  scale_x_continuous(expression(paste("Mean Non-Wildfire P",M[2][.][5])),expand=c(0,0)) +
  scale_y_continuous("Odds of Mortality", trans="log"
                     , breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)) +
  geom_rug(sides="b")  

png(file.path(indir1, "figures", "publication", "figureS2.png"),
    width=6, height=4,units="in", res = 600, bg="white",
    family="sans")
print(spline_plot)
dev.off()
#######
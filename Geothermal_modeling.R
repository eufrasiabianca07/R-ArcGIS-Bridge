###Geothermal potential modeling and correlation analysis###
###Esri Indonesia Webinar Session###


# Call all the libraries needed
library(arcgisbinding)
library(sp)
library(spdep)
library (reshape2)
library (ggplot2)
library (ggmap)
library (reshape)

# Check ArcGIS licensing
arc.check_product()

# Open, select and enrich the data
enrich_df <- arc.open(path = "C:\\Users\\ebadiatmiko\\Documents\\ArcGIS\\Projects\\Geothermal Potential\\Geothermal Potential.gdb\\Geothermal_potential_enriched")
enrich_select_df <- arc.select(object = enrich_df, fields = c('OBJECTID', 'SUM_VALUE', 'hotspring_dist', 'caldera_dist', 'volcanic_dist', 'alteration_dist', 'resistivity', 'faults_dist'))
enrich_spdf <- arc.data2sp(enrich_select_df)
head(enrich_spdf@data)

# Change field names
col_names <- c("OBJECTID", "potential_index", "hotspring_dist", "caldera_dist", "volcanic_dist", "alteration_dist", "resistivity", "faults_dist")
colnames(enrich_spdf@data) <- col_names
head(enrich_spdf@data)

# Calculate Empirical Bayesian rate
n <- enrich_spdf@data$potential_index
x <- enrich_spdf@data$hotspring_dist
EB <- EBest (x,n)
p <- EB$raw
b <- attr(EB, "parameters")$b
a <- attr(EB, "parameters")$a
v <- a + (b/x)
v[v < 0] <- b/x
z <- (p - b)/sqrt(v)
enrich_spdf@data$EB_Rate <- z
arcgis_df <- arc.sp2data(enrich_spdf)

# Save output to spatial format
arc.write("C:\\Users\\ebadiatmiko\\Documents\\ArcGIS\\Projects\\Geothermal Potential\\Geothermal Potential.gdb\\Geothermal_EBrate", arcgis_df, shape_info = arc.shapeinfo(enrich_df))

# Call the EB rate calculation output
rate_df <- arc.open("C:\\Users\\ebadiatmiko\\Documents\\ArcGIS\\Projects\\Geothermal Potential\\Geothermal Potential.gdb\\Geothermal_EBrate")
rate_select_df <- arc.select(rate_df, fields = c("OBJECTID", "potential_index", "hotspring_dist", "caldera_dist", "volcanic_dist", "alteration_dist", "resistivity", "faults_dist", "EB_Rate"))
rate_spdf <- arc.data2sp(rate_select_df)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}
#
reorder_cormat <- function(cormat) {
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat) / 2)
  hc <- hclust(dd)
  cormat <- cormat [hc$order, hc$order]
}

# Create correlation matrix
corr_sub <- rate_spdf@data [ c ("potential_index", "resistivity", "faults_dist", "volcanic_dist", "alteration_dist", "caldera_dist", "EB_Rate")]
cormax <- round (cor(corr_sub), 2)
upper_tri <- get_upper_tri (cormax)
melted_cormax <- melt (upper_tri, na.rm = TRUE)
cormax <- reorder_cormat (cormax)
upper_tri <- get_upper_tri (cormax)
melted_cormax <- melt (upper_tri, na.rm = TRUE)
head(melted_cormax)
ggheatmap <- ggplot (melted_cormax, aes (X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2 (low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name = "Pearson\nCorrelation") +
  theme_minimal() + # minimal theme
  theme (axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed()
print (ggheatmap)
ggheatmap +
  geom_text (aes (X2, X1, label = value), color = "black", size = 4) +
  theme (
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c (1, 0),
   # legend.position = c (0.6, 0.7),
    legend.direction = "horizontal") +
 guides (fill = guide_colorbar (barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
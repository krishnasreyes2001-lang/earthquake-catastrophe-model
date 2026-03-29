#install.packages(c("sf", "rnaturalearth", "rnaturalearthdata"))
#install.packages("geosphere")
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)

options(scipen = 999)
######load all 4 data sets#####
EQ_76_96 = read.csv("\\query 1976-96.csv")
EQ_97_06 = read.csv("\\query 1997-06.csv")
EQ_07_16 = read.csv("\query 2007-16.csv")
EQ_17_26 = read.csv("\\query 2017-26.csv")

EQ = rbind(EQ_17_26, EQ_07_16, EQ_97_06, EQ_76_96)

colnames(EQ)

#####Data insights#####
EQ$time =  as.POSIXct(EQ$time, tz = "UTC")
EQ = EQ %>% mutate(Year = as.integer(format(time,"%Y")))
EQ = EQ %>% mutate("DepthClass" = ifelse(depth<70, "Shallow (<70 km)", ifelse(depth<300, "Intermediate (70–300 km)", "Deep (>300 km)")))
EQ_sf = st_as_sf(EQ, coords = c("longitude", "latitude"), crs = 4326)
world = ne_countries(scale = "medium", returnclass = "sf")

thresholds_mag = seq(2.5, floor(max(EQ$mag)*2)/2,by=0.5)
magnitude_thresholds_counts = 
  EQ %>% 
  crossing(threshold = thresholds_mag) %>%
  group_by(Year, threshold) %>% 
  summarise(Thres = sum(mag >= threshold, na.rm = TRUE), .groups = "drop") %>% 
  mutate(colname = paste0("MAG Greater/Equal ", threshold)) %>%
  select(-threshold) %>%
  pivot_wider(names_from = colname, values_from = Thres, values_fill = 0) %>%
  arrange(Year)

thresholds_depth = seq(floor(min(EQ$depth)*2)/2, floor(max(EQ$depth)*2)/2,by=40)
depth_thresholds_counts = 
  EQ %>% 
  crossing(threshold = thresholds_depth) %>%
  group_by(Year, threshold) %>% 
  summarise(Thres = sum(depth >= threshold, na.rm = TRUE), .groups = "drop") %>%
  mutate(colname = paste0("Depth Greater/Equal ", threshold)) %>%
  select(-threshold) %>%
  pivot_wider(names_from = colname, values_from = Thres, values_fill = 0) %>%
  arrange(Year)

year_counts = 
  EQ %>% 
  group_by(Year) %>% 
  summarise(Total_EQs = n()) %>%
  ungroup() %>%
  arrange(Year)

depth_counts = EQ %>% 
  group_by(Year, DepthClass) %>% 
  summarise(Counts = n(), .groups = "drop") %>% 
  pivot_wider(names_from = DepthClass, values_from = Counts, values_fill = 0) %>% 
  arrange(Year)

ggplot(EQ, aes(x = Year)) +
  geom_bar(color = "white", fill = "black") +
  labs(title = "Count of earthquakes", xlab = "Year", ylab = "Count") +
  theme_minimal()

ggplot(EQ, aes(mag)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "white") + 
  labs(title = "Histogram of EQ magnitude", x = "Magnitude", y = "Frequency") + 
  theme_minimal()

ggplot(EQ, aes(mag)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "white") + 
  scale_y_log10() +
  labs(title = "Histogram of EQ magnitude Log Transformed", x = "Magnitude", y = "Frequency") + 
  theme_minimal()

ggplot(EQ, aes(depth)) +
  geom_histogram(binwidth = 10, fill = "black", color = "white") + 
  labs(title = "Histogram of EQ depth", x = "Depth", y = "Frequency") + 
  theme_minimal()

bbox_eq = st_bbox(EQ_sf)
sf_use_s2(FALSE)
world_valid = st_make_valid(world)
world_crop = st_crop(world_valid, bbox_eq)
ggplot() +
  geom_sf(data = world_crop, fill = "gray95", color = "gray60", linewidth = 0.2) +
  geom_point(data = EQ, aes(x = longitude, y = latitude, color = mag), alpha = 0.6, size = 0.7) + 
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Earthquake Geographical Distribution", x = "Longitude", y = "Latitude") +
  coord_sf(expand = FALSE) +
  theme_minimal()
  
ggplot(EQ, aes(x = depth, y = mag)) +
  geom_point(alpha = 0.4) + 
  coord_cartesian(xlim = c(0,300)) +
  labs(title = "Depth VS Magnitude Scatter Plot", x = "Depth", y = "Magnitude") + 
  theme_classic()

cor(EQ$depth, EQ$mag)

nrow(EQ)
sum(is.na(EQ$mag))
table(EQ$magType)


####Below what magnitude is my data unreliable?####
mag_bins = seq(floor(min(EQ$mag)), ceiling(max(EQ$mag)), by = 0.1)

mag_cum = EQ %>% mutate(mag_bin = floor(mag*10)/10) %>% count(mag_bin) %>% arrange(desc(mag_bin)) %>% mutate(cum_count = cumsum(n))

ggplot(mag_cum, aes(x = mag_bin, y = cum_count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  labs(
    title = "Cumulative magnitude frequency distribution",
    x = "Magnitute",
    y = "Cummulative count (log scale)"
  ) +
  theme_minimal()

Mc <- 4.2 #lets remove magnitude below 4.2 because many low magnitude earthquakes are missing from observations.

EQ_model = EQ %>% filter(mag >= Mc)

yearly_model_count = EQ_model %>% count(Year)

lambda_hat = mean(yearly_model_count$n)
lambda_hat

#estimating the a and b values of Gutenberg–Richter model
thresholds = seq(Mc, floor(max(EQ_model$mag)*2)/2, by = 0.5)

mag_counts = data.frame(
  thresholds = thresholds,
  count = sapply(thresholds, function(x) {
    sum(EQ_model$mag >= x)
  })
)

mag_counts$log_counts = log10(mag_counts$count)

model = lm(log_counts ~ thresholds, data = mag_counts)

a_value = coef(model)[1]
b_value = -coef(model)[2]

####Simulation of earthquakes####
n_years = 1000
set.seed(123)
n_events_per_yesr = rpois(n_years, lambda_hat)
head(n_events_per_yesr)

#Simulate magnitude for each year and each earthquake from that year

simulated_events = data.frame()

for(i in 1:n_years) {
  n_events = n_events_per_yesr[i]
  if(n_events == 0) next
  mags = sample(EQ_model$mag, n_events, replace = TRUE)
  
  #sample location and depth from the actual data
  
  sample_rows = EQ_model %>% sample_n(n_events, replace = TRUE)
  
  temp = data.frame(
    Year = i,
    mag = mags,
    longitude = sample_rows$longitude,
    latitude = sample_rows$latitude,
    depth = sample_rows$depth
  )
  simulated_events = rbind(simulated_events, temp)
}

nrow(simulated_events)

ggplot(simulated_events, aes(mag)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "white") + 
  scale_y_log10() +
  labs(
    title = "Histogram of simulated Magnitudes (Y axis is scaled)",
    x = "Magnitude",
    y = "Frequency"
  )

ggplot() +
  geom_point(data = simulated_events, aes(longitude, latitude), alpha = 0.2)

####Modeling PGA in Chennai, India####
#Approximate latitude and longitude of Chennai
site_lon <- 80.2707
site_lat <- 13.0827

#epicenter distances from Chennai in KM 
simulated_events$distance_km = distHaversine(
  matrix(c(simulated_events$longitude, simulated_events$latitude), ncol = 2),
  c(site_lon, site_lat)
)/1000

simulated_events$HC_distance = sqrt(simulated_events$distance_km^2 + simulated_events$depth^2)

#Simulating PGA
compute_pga = function(M, X) {
  log10_A = -1.072 + 0.3903 * M - 1.21 * log10(X + exp(0.5873 * M))
  A = 10^(log10_A)
  return(A)
}

simulated_events$PGA = compute_pga(simulated_events$mag, simulated_events$HC_distance)

summary(simulated_events$PGA)

ggplot(simulated_events, aes(PGA)) +
  geom_histogram(bins = 50, fill = "black", color = "white") +
  scale_x_log10() +
  labs(
    title = "Histogram of PGA",
    x = "PGA",
    y = "Frequency"
  ) +
  theme_minimal()

#Maximum PGA for the year
annual_max_pga = simulated_events %>% group_by(Year) %>% summarise(MAX_pga = max(PGA,na.rm = TRUE)) %>% ungroup()

summary(annual_max_pga$MAX_pga)

#Generate an hazard curve
pga_thresholds = seq(0, max(annual_max_pga$MAX_pga), length.out = 100)

hazard_curve = data.frame(
  PGA = pga_thresholds,
  Exceedence_Probability = sapply(pga_thresholds, function(x) {
    mean(annual_max_pga$MAX_pga >= x)
  })
)

ggplot(hazard_curve, aes(x = PGA, y = Exceedence_Probability)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Hazard Curve for Chennai",
    x = "PGA",
    y = "Annual Exceedence Probability"
  )

ggplot(hazard_curve, aes(x = PGA, y = Exceedence_Probability)) +
  geom_line() +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Hazard Curve for Chennai (Y Axis Scaled)",
    x = "PGA",
    y = "Annual Exceedence Probability"
  )

hazard_curve = hazard_curve %>% mutate(Return_Period = 1/Exceedence_Probability)

#Due to lack of publicly available insured exposure data, a data set was constructed based on reasonable assumptions regarding building density, value distribution, and construction type.

#####Generated Exposure data#####
set.seed(123)

exposure = data.frame(
  asset_id = 1:1000,
  longitude = runif(1000, 79.8, 80.5),
  latitude = runif(1000, 12.8, 13.3),
  value = runif(1000, 10^5, 10^6),
  building_type = sample(
    c("Masonry", "RC"),
    1000,
    replace = TRUE,
    prob = c(0.65, 0.35)
  )
)

summary(exposure)
table(exposure$building_type)

#####Vulnerability Model#####
#We’ll use a simple deterministic function returning damage ratio
vulnerability_function = function(pga, building_type) {
  if(building_type == "Masonry") {
    if(pga < 0.005) return(0.00)
    if(pga < 0.010) return(0.02)
    if(pga < 0.020) return(0.08)
    if(pga < 0.030) return(0.15)
    if(pga < 0.050) return(0.30)
    return(0.60)
  }
  if(building_type == "RC") {
    if(pga < 0.005) return(0.00)
    if(pga < 0.010) return(0.01)
    if(pga < 0.020) return(0.04)
    if(pga < 0.030) return(0.08)
    if(pga < 0.050) return(0.18)
    return(0.35)
  }
}

vulnerability_function(0.015, "Masonry")
vulnerability_function(0.015, "RC")

####Annual Losses####
#We assume that the max pga for the year is experienced uniformly across all assets in Chennai

annual_loss = annual_max_pga %>% rowwise() %>% mutate(Portfolio_Loss = sum(sapply(1:nrow(exposure), function(i){
  dr = vulnerability_function(MAX_pga, exposure$building_type[i])
  exposure$value[i]*dr
}))) %>% ungroup()

summary(annual_loss$Portfolio_Loss)

ggplot(annual_loss, aes(Portfolio_Loss)) +
  geom_histogram(bins = 40, fill="black", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Annual Portfolio Loss",
    x = "Annual Loss",
    y = "Frequency"
  )

loss_thresholds = seq(0, max(annual_loss$Portfolio_Loss), length.out = 100)
loss_curve = data.frame(
  Loss = loss_thresholds,
  Exceed_prob = sapply(loss_thresholds, function(x) {
    mean(annual_loss$Portfolio_Loss >= x)
  })
)

ggplot(loss_curve, aes(Loss, Exceed_prob)) + 
  geom_point() +
  geom_line() + 
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Loss Exceedance Probabilities",
    x = "Loss",
    y = "Probability that Loss exceeds X"
  )



AAL = mean(annual_loss$Portfolio_Loss)
format(AAL, big.mark = ",", scientific = FALSE)
format(quantile(annual_loss$Portfolio_Loss, c(0.9, 0.95, 0.99, 0.995)), big.mark = ",", scientific = FALSE)

































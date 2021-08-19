library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

dives <- read.csv("../data/KBB_Ventana_dives.csv") %>%
  mutate(datetime = mdy_hm(index_recorded_timestamp),
         date = date(datetime))
head(dives)

counts <- dives %>%
  filter(activity=="transect") %>%
  mutate(depth_meters = round(depth_meters, -2)) %>%
  group_by(date, depth_meters, concept) %>%
  summarize(n=n())


relevant <- c("Actinopterygii", "Euphausiacea", "ghost tail", "Cyclothone", 
  "Eusergestes similis", "Amphipoda", "Mysida", "Pasiphaea", "Teuthoidea",
  "Leuroglossus stilbius", "Paralepididae", "Bathylagidae", "Sternoptychidae",
  "Pleuronectiformes 2", "Octopoteuthis deletron", "Cranchiidae", "Sebastes",
  "Melanostigma pammelas", "Histioteuthis heteropsis", "Phronima sedentaria",
  "Cystisoma", "Sergestoidea", "Sternoptyx diaphana", "Paraphronima",
  "Lagenorhynchus obliquidens", "Oplophoridae", "Parvilux ingens",
  "Euchirella bitumida", "Scina", "Idiacanthus antrostomus", "Stomiidae",
  "Tactostoma macropus", "Lanceolidae", "Gnathophausia ingens", "Gonatus",
  "Chiroteuthis calyx", "Sergestidae", "Munneurycope murrayi", "Idiacanthus",
  "Myctophidae", "Merluccius productus")

prey <- c("Actinpterygii", "Bathylagidae", "Chiroteuthis calyx", "Cranchiidae",
          "Cyclothone", "Gonatus", "Histioteuthis heteropsis", "Idiacanthus", 
          "Idiacanthus antrostomus", "Lanceolidae", "Leuroglossus stilbius",
          "Melanostigma pammelas", "Merluccius productus", "Myctophidae", 
          "Octopoteuthis deletron", "Parvilux ingens", "Pleuronectiformes 2",
          "Sebastes", "Sternoptychidae", "Sternoptyx diaphana", "Stomiidae",
          "Teuthoidea")

groups <- data.frame(
  concept = prey,
  Group = c("Small mesopelagics", "Small mesopelagics", "Squid", "Squid",
            "Small mesopelagics", "Squid", "Squid", "Small mesopelagics", 
            "Small mesopelagics", "Small mesopelagics", "Small mesopelagics",
            "Small mesopelagics", "Hake", "Small mesopelagics",
            "Squid", "Small mesopelagics", "Small mesopelagics",
            "Rockfish", "Small mesopelagics", "Small mesopelagics", "Small mesopelagics",
            "Squid"))

dives <- left_join(dives, groups)

dives %>%
  filter(activity=="transect", concept %in% prey) %>%
  mutate(depth_meters = round(depth_meters, -1)) %>%
  ggplot(aes(x=depth_meters, fill=Group)) + 
    geom_bar(width=20) +
    coord_flip() + 
    scale_x_reverse("Depth (m)", limits=c(800, 0)) + 
    ylab("Count") +
    theme_minimal() +
    facet_wrap(~date, nrow=1)
ggsave("../graphics/rov_counts.png", h=9/3, w=16/2)


layerbase <- read.csv("../data/layerbase.csv")
ggplot(layerbase, aes(x=layer_base, y=..density..)) + 
  geom_histogram() +
  coord_flip() +
  scale_x_reverse("Layer depth (m)", limits=c(800, 0)) + 
  ylab("Probability") +
  theme_minimal()
ggsave("../graphics/layerbase_histogram.png", h=9/3, w=9/5)



profiles <- dives %>%
  filter(activity=="descend", light_transmission > 0) %>%
  mutate(depth = round(depth_meters, -1)) %>%
  group_by(date, depth) %>%
  summarize(temperature = mean(temperature_celsius),
            salinity = mean(salinity),
            oxygen = mean(oxygen_ml_per_l),
            transmission = mean(light_transmission)) %>%
  arrange(depth) %>%
  ungroup()

profiles %>%
  pivot_longer(c(-date, -depth)) %>%
  ggplot(aes(y=depth, x=value, color=as.factor(date))) +
    geom_path() +
    scale_y_reverse() +
    facet_wrap(~name, scales="free_x")

profiles %>%
  filter(date == "2019-07-23") %>%
  ggplot(aes(y=depth, x=oxygen)) + 
  geom_path() + scale_y_reverse()
profiles %>%
  filter(date == "2019-07-23", depth %in% c(150, 250))

# 
k = 0.08
light = 2000 * exp(-k*0:800)

png("../graphics/demo_ctd_profile.png", w=16/4, h=7, units="in", res=150)
  par(mar=c(7, 4, 4, 2))
  plot(-depth ~ temperature, filter(profiles, date=="2019-02-07"), ty='l',
       xlab="", ylab="Depth (m)")
  mtext("Temperature (C)", side=1, line=2)
  par(new=T)
  plot(-(0:800) ~ light, filter(profiles, date=="2019-02-07"), ty='l', col="red",
       xaxt="n", xlab="", yaxt="n", ylab="")
  axis(3, col='red', col.axis="red",las=1)
  mtext(expression(Light~(mu*mol~m^-2~s^-1)), side=3, line=2, col="red")
  par(new=T)
  plot(-depth ~ oxygen, filter(profiles, date=="2019-02-07"), ty='l', col="blue",
       xaxt="n", xlab="", yaxt="n", ylab="")
  axis(1, col="blue", col.axis="blue", line=3.5)  
  mtext("Oxygen (mL/L)", side=1, line=5.5, col="blue")
dev.off()

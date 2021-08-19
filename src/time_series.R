library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(viridis)
library(scales)

TZ.OFFSET <- hours(8)
DATA.DIR <- "../data/"
GRAPHICS.DIR <- "../graphics/"

local.tz <- "America/Los_Angeles"

clicks <- arrow::read_feather(paste0(DATA.DIR, "clicks_10min.feather")) %>%
  # mutate(datetime = datetime_round) %>% 
  select(datetime_round, n)

clicks %>%
  filter(datetime_round > ymd_hm("20190210 1700"), 
         datetime_round < ymd_hm("20190211 1200")) %>%
  ggplot(aes(x=datetime_round, y=n)) + geom_line()

# Center of mass, excluding backscatter from epipelagic schools.
# Not using the "layerbase" variable from this data frame, because these
# metrics were calculated on just the upper 500 m (to better capture the
# change in CM due to migration, by excluding the deep non-migrating layer), so
# this layerbase is artificially squashed
center_of_mass <- read.csv("../data/center_of_mass.csv", stringsAsFactors = F) %>%
  mutate(datetime_round = round_date(ymd_hms(datetime), "10 minutes")) %>%
  group_by(datetime_round) %>%
  summarize(center_of_mass = mean(cm),
            layertop = mean(layertop))

layerbase <- read.csv("../data/layerbase.csv", stringsAsFactors = F) %>%
  mutate(datetime_round = ymd_hms(datetime_round) - hours(8))

dlayerbase <- layerbase %>%
  select(datetime_round, layer_base, center_of_mass) %>%
  arrange(datetime_round) %>%
  mutate(dlayer_base = c(0, diff(layer_base)),
         dcenter_of_mass = c(0, diff(center_of_mass))) %>%
  ungroup()

clicks.diff <- clicks  %>%
  arrange(datetime_round) %>%
  mutate(dn = c(0, diff(n))) %>%
  ungroup() %>%
  left_join(dlayerbase, by="datetime_round") %>%
  filter(!is.na(layer_base))

t1 <- "20190210 1200"
t2 <- "20190211 1200"

p1 <- clicks.diff %>%
  filter(datetime_round > ymd_hm(t1),
         datetime_round < ymd_hm(t2)) %>%
  ggplot(aes(x=datetime_round, y=n)) + geom_line() + geom_point()

p2 <- clicks.diff %>%
  filter(datetime_round > ymd_hm(t1),
         datetime_round < ymd_hm(t2)) %>%
  ggplot(aes(x=datetime_round, y=layer_base)) + geom_line() + geom_point() + 
   scale_y_reverse()
gridExtra::grid.arrange(p1, p2, nrow=2)



dclicks <- clicks.diff$dn
dlayer <- clicks.diff$dlayer_base
click.ccf <- ccf(ts(dclicks, deltat=10), ts(dlayer, deltat=10), 
                 na.action=na.pass, lag.max=6)


#######################
# Clicking bout-by-bout
#######################
label_bouts <- function(bouts) {
  bouts.rle <- rle(bouts)
  bout.labels <- as.numeric(bouts.rle$values)
  bout.labels[bouts.rle$values] <- 1:sum(bouts.rle$values)
  return(inverse.rle(list(lengths=bouts.rle$lengths, values=bout.labels)))
}

dilate <- function(x, left=1, right=1) {
  n <- length(x)
  xnew <- rep(FALSE, n)
  for (i in 1:right) {
    xnew[i] <- max(x[1:i])
  }
  for (i in right:(n-left)) {
    xnew[i] <- max(x[(i-right):(i+left)])
  }
  for (i in (n-left):n) {
    xnew[i] <- max(x[(i-right):n])
  }
  return(xnew)
}

right.lag <- 4

clicks1 <- clicks %>%
  arrange(datetime_round) %>%
  mutate(bout_thresh = quantile(n, 0.95),
         bout = dilate(n > bout_thresh, left=1, right=right.lag) == 1,
         # bout = n > bout_thresh,
         bout_label = label_bouts(bout)) %>%
  ungroup()


clicks1 <- clicks1 %>%
  filter(bout) %>%
  left_join(layerbase) %>%
  group_by(bout_label) %>%
  mutate(bout_time = as.numeric(datetime_round - min(datetime_round)),
         rel_depth = layer_base - first(layer_base)) %>%
  arrange(bout_time) %>%
  filter(!is.na(year), bout_time <= 250*60) %>%
  ungroup() %>%
  drop_na()


ggplot(clicks1, aes(x=bout_time/60, y=rel_depth, group=bout_label)) +
  geom_line(alpha=0.25) + 
  geom_point(aes(size=n), alpha=0.25) 


# Average change in depth over course of bout
clicks1 %>%
  group_by(bout_time) %>%
  summarize(se = sd(rel_depth, na.rm=T) / n(),
            q25 = quantile(-rel_depth, 0.25, na.rm=T),
            q75 = quantile(-rel_depth, 0.75, na.rm=T),
            q40 = quantile(-rel_depth, 0.40, na.rm=T),
            q60 = quantile(-rel_depth, 0.60, na.rm=T),
            rel_depth = median(-rel_depth, na.rm=T),
            n=n()) %>%
  ungroup() %>%
  ggplot(aes(x=bout_time/60, y=rel_depth)) + 
  # geom_line(aes(x=bout_time/60, y=-rel_depth, group=bout_label), data=clicks, alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(ymin=q25, ymax=q75), alpha=0.2) +
  geom_ribbon(aes(ymin=q40, ymax=q60), alpha=0.2) +
  geom_line() + 
  geom_point() +
  ylab("Depth change (m)") + xlab("Time since start of bout (min)") +
  scale_color_discrete("Click type") +
  theme_minimal() + xlim(0, 100) + ylim(-25, 15)
ggsave("../graphics/depth_change_during_bout.png", h=9/2, w=16/2)  


bouts <- clicks1 %>%
  group_by(bout_label) %>%
  summarize(bout_start = min(datetime_round), 
            bout_end = max(datetime_round) + minutes(10),
            bout_start_utc = min(datetime_round) + hours(8),
            bout_end_utc = max(datetime_round) + minutes(10) + hours(8),
            nclicks = sum(n),
            layerdist = sum(abs(diff(rel_depth)), na.rm=T),
            maxchange = max(rel_depth, na.rm=T),
            init_dive_speed = nth(rel_depth, 2) / 20,
            distance = sum(abs(diff(rel_depth)))) %>%
  ungroup() %>%
  mutate(duration = as.numeric(bout_end - bout_start, units="mins"),
         date = date(bout_start),
         hour = hour(bout_start) + minute(bout_start)/60)

hist(bouts$duration)
mean(bouts$duration)
max(bouts$duration)

hist(bouts$maxchange)
mean(bouts$maxchange)
max(bouts$maxchange)
quantile(bouts$maxchange, c(60, 65, 75, 80, 85, 90, 95, 99)/100)

hist(bouts$init_dive_speed, 50)
mean(bouts$init_dive_speed)
median(bouts$init_dive_speed)
max(bouts$init_dive_speed)

hist(bouts$distance)
median(bouts$distance)
mean(bouts$distance)
max(bouts$distance)


daily.bouts <- bouts %>%
  mutate(date=date(bout_start)) %>%
  group_by(date) %>%
  summarize(nbouts=n()) %>%
  ungroup()
hist(daily.bouts$nbouts)
mean(daily.bouts$nbouts)

mean(daily.bouts$nbouts) * mean(bouts$distance)

bouts %>%
  mutate(daytime = (hour > 6) & (hour < 18)) %>%
  group_by(daytime) %>%
  summarize(nclicks = sum(nclicks)) %>%
  mutate(p = nclicks / sum(nclicks))



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


clicks <- arrow::read_feather(paste0(DATA.DIR, "click_counts_bio.feather")) %>%
  mutate(datetime = time_round, n = as.numeric(n)) %>%
  select(datetime, n)
tz(clicks$datetime) <- local.tz # Have to do this, otherwise comparisons w/ ymd_hm don't work


layerbase <- read.csv("../data/layerbase.csv") %>%
  arrange(datetime_round) %>%
  mutate(datetime = ymd_hms(datetime_round)) %>%
  select(datetime, layerbase=layer_base)


read.echo.df <- function(filename) {
  echo.df <- arrow::read_feather(filename) %>%
    filter(depth > 0)
  echo.df$Sv_mean[echo.df$Sv_mean < -90] <- -90 
  return(echo.df)
}


plot.echogram <- function(echo.df, plot.clicks=FALSE, depth.limits=NULL, time.limits=NULL, 
                          col.limits=c(-90, -60), date_breaks="1 hour", proportion=0.25) {
  
  if (is.null(time.limits)) {
    tmin <- min(echo.df$datetime_round)
    tmax <- max(echo.df$datetime_round)
  } else {
    tmin <- time.limits[1]
    tmax <- time.limits[2]
  }
  
  if (is.null(depth.limits)) {
    dmin <- min(echo.df$depth)
    dmax <- max(echo.df$depth)
    depth.limits <- c(dmax, dmin)
  } else {
    dmax <- max(depth.limits)
    dmin <- min(depth.limits)
  }
  drange <- dmax - dmin
  
  
  p <- ggplot() +
    geom_raster(aes(x=datetime_round, y=depth, fill=Sv_mean), echo.df) + 
    scale_fill_viridis("Sv (dB)", option="magma", limits=col.limits, oob=squish) +
    scale_x_datetime("Time", limits=c(tmin, tmax), date_breaks=date_breaks, date_labels="%H:%M", 
                     timezone=local.tz, expand=c(0, 0)) +
    ggtitle(date(first(echo.df$datetime))) +
    theme_bw()
  
  if (plot.clicks) {  
    clicks.rescaled <- clicks %>%
      filter(datetime > tmin, datetime < tmax)
    nmax = max(clicks.rescaled$n)
    clicks.rescaled <- mutate(clicks.rescaled, depth = (-drange * proportion * n / nmax) + dmax - 1)
    
    scaler <- function(x) (x - dmax) * -nmax / drange / proportion
    
    # maxtick <- 10^ceiling(log10(nmax))
    breaks <- axisTicks(c(0, nmax), log=FALSE, nint=3) #c(0, maxtick/2, maxtick)
    labels <- scientific(breaks)
    return(p + geom_line(aes(x=datetime, y=depth), data=clicks.rescaled, 
                         col="white") +
             scale_y_reverse("Depth (m)", expand=c(0, 0), limits=depth.limits,
                             sec.axis = sec_axis(scaler, name="", breaks=breaks, labels=labels)))
  } else {
    return(p + scale_y_reverse("Depth (m)", limits=depth.limits, expand=c(0, 0)))
  }
}

avoidance.files <- list.files(DATA.DIR, pattern="avoidance-*", full.names = T)
avoidance.files <- avoidance.files[startsWith(avoidance.files, paste0(DATA.DIR, "avoidance-"))]
avoidance.files <- avoidance.files[endsWith(avoidance.files, ".feather")]
auv.avoidance.files <- list.files(DATA.DIR, pattern="auv-avoidance-", full.names = T)
auv.avoidance.files <- auv.avoidance.files[endsWith(auv.avoidance.files, ".feather")]
school.avoidance.files <- list.files(DATA.DIR, pattern="school-avoidance-", full.names = T)
school.avoidance.files <- school.avoidance.files[endsWith(school.avoidance.files, ".feather")]

for (f in school.avoidance.files) {
  print(f)
  echo.df <- read.echo.df(paste0(DATA.DIR, f))
  p <- plot.echogram(echo.df, plot.clicks = F)
  figname <- paste0(GRAPHICS.DIR, gsub(".feather", ".png", basename(f)))
  ggsave(figname, p, w=4, h=3)
}


for (f in avoidance.files) {
  figname <- paste0(GRAPHICS.DIR, gsub(".feather", ".png", basename(f)))
  print(figname)
  echo.df <- read.echo.df(paste0(DATA.DIR, f))
  p <- plot.echogram(echo.df, plot.clicks = T)
  ggsave(figname, p, w=4.5, h=3)
}


# Figure 1A
f <- avoidance.files[3]
echo.df <- read.echo.df(paste0(DATA.DIR, f))
p <- plot.echogram(echo.df, plot.clicks = T, depth.limits=c(700, 300),
                   time.limits=c(ymd_hm("2019-03-06 22:00", tz=local.tz), ymd_hm("2019-03-07 01:00", tz=local.tz)))
figname <- paste0(GRAPHICS.DIR, "Figure-1A.png")
ggsave(figname, p, w=4.5, h=3)

# Figure 1B
f <- avoidance.files[6]
echo.df <- read.echo.df(paste0(DATA.DIR, f))
p <- plot.echogram(echo.df, plot.clicks = T, depth.limits=c(200, 10), date_breaks="20 min",
                   time.limits=c(ymd_hm("2019-05-03 20:30", tz=local.tz), ymd_hm("2019-05-03 21:30", tz=local.tz)))
figname <- paste0(GRAPHICS.DIR, "Figure-1B.png")
ggsave(figname, p, w=4.5, h=3)

# Figure 1C
f <- avoidance.files[2]
echo.df <- read.echo.df(paste0(DATA.DIR, f))
p <- plot.echogram(echo.df, plot.clicks = T, date_breaks="30 min",
                   time.limits=c(ymd_hm("2019-02-11 02:00", tz=local.tz), ymd_hm("2019-02-11 03:15", tz=local.tz)))
figname <- paste0(GRAPHICS.DIR, "Figure-1C.png")
ggsave(figname, p, w=4.5, h=3)

# Figure 1D
f <- auv.avoidance.files[1]
echo.df <- read.echo.df(paste0(DATA.DIR, f))
p <- plot.echogram(echo.df, plot.clicks = F, depth.limits = c(275, 75), date_breaks="20 min",
                   time.limits=c(ymd_hm("2019-10-07 14:31", tz=local.tz), ymd_hm("2019-10-07 15:30", tz=local.tz)))
figname <- paste0(GRAPHICS.DIR, "Figure-1D.png")
ggsave(figname, p, w=4., h=3)

# Figure 1E (wider format)
f <- paste0(DATA.DIR, "deimos-2019-02-11-med-res.feather")
figname <- paste0(GRAPHICS.DIR, "Figure-1E.png")
echo.df <- read.echo.df(paste0(DATA.DIR, f)) %>%
  mutate(datetime_round = floor_date(datetime_round, "1 minutes"))
p <- plot.echogram(echo.df, plot.clicks = T, date_breaks="2 hour") +
  geom_line(aes(x=datetime, y=layerbase), data=layerbase, color="white", linetype=2)
figname <- paste0(GRAPHICS.DIR, "Figure-1E.png")
ggsave(figname, p, w=8, h=3)


echo.df <- read.echo.df(paste0(DATA.DIR, "avoidance-2019-03-07-0600.feather"))
# plot.echogram(echo.df)
plot.echogram(echo.df, T)


#######################
# Migration transitions
#######################

# each set of dates is outer start, inner start, inner end, outer end
dt_bounds <- list(
  c(ymd_hm("2019-04-01 00:00"), ymd_hm("2019-04-06 00:00"), 
    ymd_hm("2019-04-09 00:00"), ymd_hm("2019-04-14 00:00")),
  c(ymd_hm("2019-08-16 00:00"), ymd_hm("2019-08-21 00:00"), 
    ymd_hm("2019-08-24 00:00"), ymd_hm("2019-08-29 00:00"))
)

cm.hourly <- read.csv("../data/center_of_mass.csv") %>%
  mutate(datetime = ymd_hms(datetime))

for (bounds in dt_bounds) {
  s2 <- bounds[1]
  s1 <- bounds[2]
  e1 <- bounds[3]
  e2 <- bounds[4]
  yearmonth <- format(s1, "%Y-%m")
  print(yearmonth)
  
  print("Loading data...")
  echo.mres <- read.echo.df(paste0(DATA.DIR, "deimos-", yearmonth, "-med-res.feather")) %>%
    mutate(datetime_round = floor_date(datetime_round, "2 minutes") - period(8, "hour"))

  echo <- read.echo.df(paste0(DATA.DIR, "deimos-", yearmonth, ".feather")) %>%
    mutate(datetime_round = floor_date(datetime_round, "10 minutes") - period(8, "hour"))
  
  
  print("Plotting...")
  p1 <- ggplot() +
    geom_raster(aes(x=datetime_round, y=depth, fill=Sv_mean), 
                filter(echo.mres, datetime > s1, datetime_round < e1)) + 
    geom_line(aes(x=datetime, y=layertop_smooth), color="white",
              data=filter(cm.hourly, datetime > s1, datetime < e1)) +
    scale_fill_viridis("Sv (dB)", option="magma", limits=c(-85, -55), oob=squish) +
    scale_x_datetime("Time", date_breaks="4 hour", date_labels="%H:%M", limits=c(s1, e1), expand=c(0, 0)) +
    scale_y_reverse("Depth (m)", expand=c(0, 0), limits=c(400, 10)) +
    theme_bw() + 
    theme(panel.background = element_rect(fill="black"), panel.grid=element_line(color="black"))
  
  p2 <- ggplot() +
    geom_raster(aes(x=datetime_round, y=depth, fill=Sv_mean), 
                filter(echo, datetime > s2, datetime_round < e2)) + 
    geom_line(aes(x=datetime, y=layertop_smooth), color="white",
              data=filter(cm.hourly, datetime > s2, datetime < e2)) +
    scale_fill_viridis("Sv (dB)", option="magma", limits=c(-85, -55), oob=squish) +
    # geom_vline(xintercept=c(s1, e1), color="white", linetype=2) +
    geom_rect(aes(xmin=s1, xmax=e1, ymin=0, ymax=400), color="white", linetype=2, fill=NA) +
    scale_x_datetime("Date", date_breaks="1 day", date_labels="%m-%d", limits=c(s2, e2), expand=c(0, 0)) +
    scale_y_reverse("Depth (m)", expand=c(0, 0)) +
    theme_bw() + 
    theme(panel.background = element_rect(fill="black"), panel.grid=element_line(color="black"))
  
  p <- gridExtra::grid.arrange(p1, p2, nrow=2, heights=c(1.5, 1))
  ggsave(paste0(GRAPHICS.DIR, "dvm-transition-", yearmonth, ".png"), p, width=9, height=4.5)
}


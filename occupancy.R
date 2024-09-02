# Code to fit an occupancy model to describe the distribution of wild species in
# the city of Toronto, in response to environmental and anthropogenic variables

library(spOccupancy)
library(tidyverse)

# Import and arrange data
occ_data_raw <- read_csv("../../Data/main_database_2207.csv")
env_data_raw <- read_csv("../../Data/station_covars.csv")

occ_data <- filter(occ_data_raw, !common_name %in% c("empty", "horse", "small rodent", "not listed", 
                                                     "bird", "mink", "groundhog", "chipmunk", "squirrel")) %>% 
  mutate(DateTime = lubridate::mdy_hms(DateTime))
ggplot(occ_data, aes(year_week,common_name))+geom_count(aes(size = after_stat(prop), group = common_name))

# Fit occupancy model
# check which cameras to include
siteswcovs <- env_data_raw[complete.cases(env_data_raw),"site_name"]
# Create occupancy model matrices
y_df <- filter(occ_data, year_week<=53) %>% 
  inner_join(siteswcovs) %>% 
  count(season, site_name, common_name) %>% 
  complete(nesting(season, site_name), common_name, fill = list(n=0)) %>% 
  complete(season, site_name, common_name) %>% 
  arrange(season,  site_name, common_name)
nsp <- length(unique(y_df$common_name))
nst <- length(unique(y_df$site_name))
nocc <- max(y_df$season)
nsp*nst*nocc
mody <- array(data = as.numeric(y_df$n>0), dim = c(nsp, nst, nocc), 
              dimnames = list(species = unique(y_df$common_name), sites = unique(y_df$site_name), occ = 1:nocc))
covs_df <- distinct(y_df,site_name) %>% left_join(env_data_raw) %>% 
  select(site_name, x, y, Tree_PA_500, DEM_500, NDVI_500, forest_dist) %>% 
  rename(trees = Tree_PA_500, elev = DEM_500, ndvi = NDVI_500, fdist = forest_dist)
temp <- matrix(rep(c(-4.4, -4.4, 0.3, 5.8, 12.3, 18.2, 21.9, 21.5, 18.2, 11.4, 5, -0.6), each=dim(mody)[2]), nrow = dim(mody)[2])
effort <- filter(y_df, n>0) %>% distinct(season,site_name) %>% mutate(a=1) %>% 
  pivot_wider(names_from = season, values_from = a) %>% 
  arrange(site_name) %>% 
  column_to_rownames("site_name") %>% 
  # mutate(across(where(is.numeric), \(x) as.numeric(x>0))) %>% 
  as.matrix()
detcov <- list(temp = temp*effort)
moddata <- list(y = mody, 
                occ.covs = covs_df[,c("trees","elev", "ndvi", "fdist")] ,
                det.covs = detcov,
                coords = covs_df[,c("x", "y")])
# fit model
m1 <- lfMsPGOcc(occ.formula = ~scale(trees)+scale(elev)+scale(ndvi)+scale(fdist),
                det.formula = ~temp,
                data = moddata, 
                n.factors = 2,n.burn = 3000, n.thin = 2, n.chains = 3, n.samples = 5000)
summary(m1)

# plot estimates
betas <- m1$beta.samples
tibble(coef = colnames(betas),
           mean = colMeans(betas),
           lq = apply(betas,2,quantile, prob = 0.025),
           uq = apply(betas,2,quantile, prob = 0.975)
) %>% separate_wider_delim(coef, "-", names = c("parameter", "species")) %>% 
  filter(parameter!="(Intercept)") %>% 
  mutate(parameter = gsub("scale(|)", "", parameter)) %>% 
  ggplot(aes(y = species, x = mean))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_linerange(aes(xmin = lq, xmax = uq))+
  geom_point()+
  facet_wrap(~parameter)

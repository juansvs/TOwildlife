# Code to fit an occupancy model to describe the distribution of wild species in
# the city of Toronto, in response to environmental and anthropogenic variables

library(spOccupancy)
library(tidyverse)

# Import and arrange data --------
occ_data_raw <- read_csv("../../Data/main_database_2207.csv")
env_data_raw <- read_csv("../../Data/station_covars.csv")

occ_data <- filter(occ_data_raw, !common_name %in% c("empty", "horse", "small rodent", "not listed", 
                                                     "bird", "mink", "groundhog", "chipmunk", "squirrel")) %>% 
  mutate(DateTime = lubridate::mdy_hms(DateTime))
ggplot(occ_data, aes(year_week,common_name))+geom_count(aes(size = after_stat(prop), group = common_name))

# Fit occupancy model
# check which cameras to include
siteswcovs <- env_data_raw[complete.cases(env_data_raw),"site_name"]
# Create model matrices -----------
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
# fit model --------
m1 <- lfMsPGOcc(occ.formula = ~scale(trees)+scale(elev)+scale(ndvi)+scale(fdist),
                det.formula = ~temp,
                data = moddata, 
                n.factors = 2,n.burn = 3000, n.thin = 2, n.chains = 3, n.samples = 5000)
summary(m1)

# Null model
mnull <- lfMsPGOcc(occ.formula = ~1,
                   det.formula = ~1,
                   data = moddata, 
                   n.factors = 2,n.burn = 3000, n.thin = 2, n.chains = 3, n.samples = 5000)
summary(mnull)
waicOcc(mnull)
waicOcc(m1)

m2 <- lfMsPGOcc(occ.formula = ~scale(elev)+scale(ndvi)+scale(fdist),
                det.formula = ~temp,
                data = moddata, 
                n.factors = 2,n.burn = 3000, n.thin = 2, n.chains = 3, n.samples = 5000)

m3 <- lfMsPGOcc(occ.formula = ~scale(ndvi)+scale(fdist),
                det.formula = ~temp,
                data = moddata, 
                n.factors = 2,n.burn = 3000, n.thin = 2, n.chains = 3, n.samples = 5000)
m3b <- lfMsPGOcc(occ.formula = ~scale(ndvi)+scale(fdist),
                det.formula = ~1,
                data = moddata, 
                n.factors = 2,n.burn = 3000, n.thin = 2, n.chains = 3, n.samples = 5000)
m4 <- lfMsPGOcc(occ.formula = ~scale(trees)+scale(fdist),
                det.formula = ~temp,
                data = moddata, 
                n.factors = 2,n.burn = 3000, n.thin = 2, n.chains = 3, n.samples = 5000)
waicOcc(m4)
# m3 is better based on waic = 2385.99
summary(m3)

# plot estimates --------
betas <- m3$beta.samples
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



#' If we focus just on canids, ie foxes and coyotes, we could use another 
#' model that explicitly estimates the joint probabilities of use, rather 
#' than the correlation between occupancy among species. The Rota et al. 
#' (2016) model allows to do this, and to estimate how different covariates 
#' affect the marginal and joint probabilities. This model is implemented
#' in the unmarked package

##### Rota model ####
library(unmarked)

canidy <- lapply(c("coyote", "fox"), \(x) mody[dimnames(mody)[[1]]==x,,])
umdata <- unmarkedFrameOccuMulti(y = canidy, #A list of S (species) elements, each one a M(sites)xJ(occasions) matrix
                                 siteCovs = covs_df,  # a dataframe of site covariates
                                 obsCovs = detcov# named list of dataframe of covariates that vary within sites
)

rt_m1 <- occuMulti(detformulas = c('~temp', '~temp'), # one per sp
          stateformulas = c('~ndvi+fdist','~ndvi+fdist','~1'), # columns correspond to 10, 01, and 11
          data =umdata)
rt_m1 
plot(rt_m1)

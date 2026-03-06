## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
library(skater)

## -----------------------------------------------------------------------------
famfile <- system.file("extdata", "3gens.fam", package="skater", mustWork=TRUE)
fam <- read_fam(famfile)
fam

## -----------------------------------------------------------------------------
peds <- fam2ped(fam)
peds

## -----------------------------------------------------------------------------
peds %>% 
  dplyr::filter(fid=="testped1") %>% 
  tidyr::unnest(cols=data)

## -----------------------------------------------------------------------------
peds$ped[[1]]

## ----plotped, fig.width=8, fig.height=8, fig.align="center"-------------------
plot(peds$ped[[1]], mar=c(1,4,1,4))

## ----eval=FALSE---------------------------------------------------------------
# plot_pedigree(peds$ped, file="3gens.ped.pdf")

## -----------------------------------------------------------------------------
ped2kinpair(peds$ped[[1]])

## -----------------------------------------------------------------------------
kinpairs <- 
  peds %>% 
  dplyr::mutate(pairs=purrr::map(ped, ped2kinpair)) %>% 
  dplyr::select(fid, pairs) %>% 
  tidyr::unnest(cols=pairs)
kinpairs

## -----------------------------------------------------------------------------
dibble()

## -----------------------------------------------------------------------------
dibble(max_degree = 5)

## -----------------------------------------------------------------------------
kin2degree(.25, max_degree=3)

## -----------------------------------------------------------------------------
kin2degree(.0312, max_degree=3)

## -----------------------------------------------------------------------------
kin2degree(.0312, max_degree=5)

## -----------------------------------------------------------------------------
# Get two pairs from each type of relationship we have in kinpairs:
kinpairs_subset <- 
  kinpairs %>% 
  dplyr::group_by(k) %>% 
  dplyr::slice(1:2)
kinpairs_subset

# Infer degree out to third degree relatives:
kinpairs_subset %>% 
  dplyr::mutate(degree=kin2degree(k, max_degree=3))

## -----------------------------------------------------------------------------
# Function to randomly flip levels of a factor (at 20%, by default)
randomflip <- function(x, p=.2) ifelse(runif(length(x))<p, sample(unique(x)), x)

# Infer degree (truth/target) using kin2degree, then randomly flip 20% of them
set.seed(42)
kinpairs_inferred <- kinpairs %>% 
  dplyr::mutate(degree_truth=kin2degree(k, max_degree=3)) %>% 
  dplyr::mutate(degree_truth=as.character(degree_truth)) %>%
  dplyr::mutate(degree_truth=tidyr::replace_na(degree_truth, "unrelated")) %>% 
  dplyr::mutate(degree_inferred=randomflip(degree_truth))
kinpairs_inferred

## -----------------------------------------------------------------------------
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth)

## -----------------------------------------------------------------------------
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth) %>% 
  purrr::pluck("Table")

## -----------------------------------------------------------------------------
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth, 
                 longer = TRUE) %>% 
  purrr::pluck("Other") %>% 
  tidyr::spread(Class, Value) %>% 
  dplyr::relocate(Average, .after=dplyr::last_col()) %>% 
  dplyr::mutate_if(rlang::is_double, signif, 2) %>% 
  knitr::kable()

## -----------------------------------------------------------------------------
hapibd_fp <- system.file("extdata", "GBR.sim.ibd.gz", package="skater", mustWork=TRUE)
hapibd_seg <- read_ibd(hapibd_fp, source = "hapibd")
hapibd_seg

## -----------------------------------------------------------------------------
gmapfile <- system.file("extdata", "sexspec-avg-min.plink.map", package="skater", mustWork=TRUE)
gmap <- read_map(gmapfile)
gmap

## -----------------------------------------------------------------------------
ibd_dat <- ibd2kin(.ibd_data=hapibd_seg, .map=gmap)
ibd_dat

## -----------------------------------------------------------------------------
pedsim_fp <- system.file("extdata", "GBR.sim.seg.gz", package="skater", mustWork=TRUE)
pedsim_seg <- read_ibd(pedsim_fp, source = "pedsim")
pedsim_seg

## -----------------------------------------------------------------------------
ibd1_dat <- ibd2kin(.ibd_data=pedsim_seg$IBD1, .map=gmap, type="IBD1")
ibd2_dat <- ibd2kin(.ibd_data=pedsim_seg$IBD2, .map=gmap, type="IBD2")
dplyr::bind_rows(ibd1_dat,ibd2_dat) %>%
  dplyr::group_by(id1,id2) %>%
  dplyr::summarise(kinship = sum(kinship), .groups = "drop")


library(rbi, warn.conflicts=FALSE)
library(rbi.helpers)
library(magrittr, warn.conflicts=FALSE)
suppressPackageStartupMessages(library(tidyverse))
library(pander)
library(lubridate, warn.conflicts=FALSE)

save_file <- function(libbi_data, file_path){
    out <- tryCatch(
      {
        message("Saving rds file.")
        save_libbi(libbi_data, name=file_path)
      },
      error=function(cond){
        message("Save rds failed with error:")
        message(cond)
      },
      finally={
        invisible()
      }
  )
    return(out)
}

# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
print(commandArgs(trailingOnly=FALSE))
# 
# end_date should be as.Date(args[1], format='%m-%d')
end_date <- args[1]
state_abbr <- args[2]
iter <- args[3]
epi_model <- args[5]

if (state_abbr == "TN"){
    population <- 6833174
    state <- 'Tennessee'
} else {
    population <- 19453561
    state <- 'New York'
}

print(state_abbr)
print(end_date)
path <- '/home/cades/volume/covid/'

in_file = paste(path, end_date, '/', state, '-covid.txt', sep="")
print(paste("Reading data from ", in_file, sep=""))

v <- read.csv(in_file, sep="\t", header=FALSE) %>%
  rowSums()

y <- data.frame(value=v) %>%
  mutate(time = seq(7, by=7, length.out=n())) %>%
  dplyr::select(time, value)

min_particles <- 5000

input_lst <- list(N=population)  
end_time <- max(y$time)
obs_lst <- list(y=y %>% dplyr::filter(time<=end_time))

if (args[4] == 'vary'){
    model_file <- paste(epi_model, "_vary_p.bi", sep='')
} else {
    model_file <- paste(epi_model, "_const_p.bi", sep='')
}

path <- ### your path here ###
model_path <- paste(path, model_file, sep="")
print(paste("Reading model from ", model_path, sep=""))
bi_model <- libbi(model_path)

posterior <- sample(bi_model, end_time=end_time, input=input_lst,
                    obs=obs_lst, nsamples=5000, nparticles=min_particles,
                    proposal='prior', nthreads=32, thin=1) %>%
  adapt_proposal(min=0.125, max=0.4) %>%
  sample(nsamples=50000)

print("Saving Data files.")
save(sample_obs, posterior,  file=paste(path, state_abbr, '_', epi_model, '_', end_date, '-', iter, '.RData', sep=""))
save_file(posterior, name=paste(path, state_abbr, '_', epi_model, '_', end_date, '.RDS', sep=""))

#####  make predictions  ######

print("Making predictions...")

pred_end <- end_time + (7 * 4)
pred_start <- end_time - (4 * 7)
pred_bi <- predict(posterior, start_time=pred_start, end_time=pred_end,
                   output_every=7, target='prediction', nsamples=25000,
                   with=c("transform-obs-to-state"))

print("Saving Data files.")
save(pred_bi, file=paste(path, state_abbr, '_', epi_model, '_preds-onemonth_', epi_model, '_', '.RData', sep=""))
save_file(pred_bi, name=paste(path, state_abbr, '_', epi_model, '_preds-onemonth_', epi_model, '_', end_date, '.RDS', sep=""))

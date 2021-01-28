
# setup netlogo
load_netlogo <- function(x,nl.path,nl.jarname,model.path) {
  # start netlogo
  library(RNetLogo)
  NLStart(nl.path, nl.jarname=nl.jarname,gui=F)
  NLLoadModel(paste(nl.path,model.path,sep=""))
}

# setup netlogo and place all run results in dataframe
sim_netlogo <- function(parameters.list) {
  # libraries
  library(stringr)
  library(dplyr)
  library(reshape2)
  
  population_from_distribution <- ifelse(parameters.list[1]==1,TRUE,FALSE)
  use_capture_time <- ifelse(parameters.list[7]==1,TRUE,FALSE)
  
  # select parameter values for the run
  NLCommand(paste0("set population-from-distribution ",population_from_distribution))
  NLCommand(paste0("set population-mean ",parameters.list[2]))
  NLCommand(paste0("set population ",parameters.list[3]))
  NLCommand(paste0("set capture-line-length ",parameters.list[4]))
  NLCommand(paste0("set num-to-capture ",parameters.list[5]))
  NLCommand(paste0("set capture-time-line ",parameters.list[6]))
  NLCommand(paste0("set use-capture-time ",use_capture_time))
  NLCommand(paste0("set total-captures ",parameters.list[8]))
  NLCommand(paste0("set time-between-captures ",parameters.list[9]))
  
  # set default "flocking" settings
  NLCommand("set vision 5.0")
  NLCommand("set minimum-separation 1.0")
  NLCommand("set max-align-turn 5.0")
  NLCommand("set max-cohere-turn 3.0")
  NLCommand("set max-separate-turn 1.5")
  
  # setup new run of the model  
  NLCommand("setup")
  
  # run model while number of current captures is less than total captures selected by user  
  NLDoCommandWhile("num-captures < total-captures + 1", "go")
  # collect results by agent
  data1 <- NLGetAgentSet(c("ticks","who","xcor","ycor","(word capture-tags)"),"turtles", as.data.frame=FALSE)
    
  return(data1)
}

nlToDataFrame <- function(netlogoOut,parameters.df){
  nRuns = nrow(netlogoOut)
  nAgents = sapply(netlogoOut[,1], length)
  
  out = data.frame(runID = rep(1:nRuns, nAgents), population_size = rep(c(nAgents), nAgents), do.call(data.frame, lapply(netlogoOut[1:5], unlist)))
  rownames(out) = NULL
  
  # add parameter information
  tmp <- parameters.df[rep(seq_len(nrow(parameters.df)), nAgents), ]
  out <- cbind(out,tmp)
  
  # rename columns
  out <- rename(out,last_tick=ticks, capture_tags=X.word.capture.tags.)
  
  # convert capture tags to better format
  # first remove brackets
  out$capture_tags<-gsub('\\[', '', out$capture_tags)
  out$capture_tags<-gsub('\\]', '', out$capture_tags)
  # split capture information into multiple columns
  #tmp<-str_split_fixed(out$capture_tags, " ", 2)
  tmp<-str_split_fixed(out$capture_tags, " ", n=Inf)
  # add new columns to data
  out<-cbind(out,tmp)
  # remove original capture columns
  out<-select(out,-c(capture_tags))
  # convert to long format
  out<-melt(out,id.vars=c("runID","population_size","last_tick","who","xcor","ycor",colnames(parameters.df)))
  #out$variable<-as.numeric(as.character(out$variable))
  out$value<-as.numeric(as.character(out$value))
  #out<-rename(out,capture_num=variable,captured=value)
  out<-rename(out,capture_num=value)
  out$capture_num[is.na(out$capture_num)]<-0
  #out<-filter(out,capture_num>0)
  out$captured<-ifelse(out$capture_num>0,1,0)
  
  return(out)
}

# function to quit netlogo
quit_netlogo <- function(x) {
  NLQuit()
}

# function to parallelize netlogo
run_netlogo <- function(parameters.df,nl.path,nl.jarname,model.path,runs,num_cores) {
  
  library(parallel)
  
  # detect the number of cores available
  processors <- num_cores
  # create cluster
  cl <- makeCluster(processors)
  #clusterExport(cl, list("nl.path","nl.jarname","model.path"), envir=environment())
  #tryCatch()
  
  # set variables for the start up process
  # load NetLogo in each processor/core
  invisible(parLapply(cl, 1:processors, load_netlogo, nl.path=nl.path, nl.jarname=nl.jarname, model.path=model.path))
  
  # create a row for every run to perform
  parameters.df <- parameters.df[rep(seq_len(nrow(parameters.df)), runs), ]
  
  # convert parameter information to list 
  parameters.list <- as.list(as.data.frame(t(parameters.df)))
  
  # run a simulation for each density value by calling parallel sapply
  result.par <- parSapply(cl, parameters.list, FUN=sim_netlogo)
  results <- data.frame(t(result.par))
  # add parameter information
  results <- cbind(results,parameters.df)
  
  # convert to easier to use format
  results <- nlToDataFrame(results,parameters.df)
  
  # quit NetLogo in each processor/core
  invisible(parLapply(cl, 1:processors, quit_netlogo))
  # stop cluster
  stopCluster(cl)
  
  return(results)
}

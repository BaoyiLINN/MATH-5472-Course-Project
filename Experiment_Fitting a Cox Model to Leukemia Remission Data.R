## Experiment: Fitting a Cox Model to Leukemia Remission Data

gehan <- read.table("https://data.princeton.edu/wws509/datasets/gehan.raw", 
                    header =  FALSE, sep ="", 
                    col.names = c( "group", "weeks", "remission"))

cox.model <- coxph(Surv(time = weeks, event = remission) ~ as.factor(group), 
                   data = gehan, method = "breslow")
summary(cox.model)
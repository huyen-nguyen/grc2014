#reads names and ages of time intervals
time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=1&limit=all")
time_int<-subset(time_int,time_int$level==5)

#finds maximum age of oldest user-specified interval
max_interval_ma<-subset(time_int$early_age,time_int$interval_name=="Bashkirian")

med_interval_ma<-subset(time_int$early_age,time_int$interval_name=="Changhsingian")

#finds minimum age of youngest user-specified interval
min_interval_ma<-subset(time_int$late_age,time_int$interval_name=="Callovian")

#reads occurrences based on specified taxa and age range
occurrences1<-read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=Brachiopoda,Bivalvia,Gastropoda,Porifera,Cnidaria,Echinodermata,Bryozoa,Trilobita,Ostracoda&min_ma=",med_interval_ma,"&max_ma=",max_interval_ma,"&show=time,phylo,geo,ident&limit=all",sep=""))
occurrences2<-read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=Brachiopoda,Bivalvia,Gastropoda,Porifera,Cnidaria,Echinodermata,Bryozoa,Trilobita,Ostracoda&min_ma=",min_interval_ma,"&max_ma=",med_interval_ma,"&show=time,phylo,geo,ident&limit=all",sep=""))
occurrences<-rbind(occurrences1,occurrences2)

#finds only those collections resolved to a single time interval
resolved_occs<-subset(occurrences,occurrences$cx_int_no %in% subset(time_int$interval_no,time_int$level==5))
  
#creates variable in data frame with stage number
resolved_occs$time_int<-resolved_occs$cx_int_no


#list of marine environments
carbonate_env<-c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
siliciclastic_env<-c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
marine_env<-c(carbonate_env,siliciclastic_env)

resolved_occs<-subset(resolved_occs,resolved_occs$environment %in% marine_env)

#FILE PREPARATION AND DATA CLEANING

  #file preparation and data cleaning for genus level analysis
  #deletes occurrences not resolved to at least genus level using classified and unclassified species
  cleaned_occs<-subset(resolved_occs,resolved_occs$matched_rank<=5)
  
  #deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
  cleaned_occs<-subset(cleaned_occs,cleaned_occs$genus_reso=="" | cleaned_occs$genus_reso=="n. gen.")
  
  #extracts genus name from matched_name string
  cleaned_occs$matched_name<-gsub(" .*","",cleaned_occs$matched_name)

  cleaned_occs$matched_name<-as.factor(cleaned_occs$matched_name)

#reads activity data
activity_level<-read.csv("https://github.com/mclapham/grc2014/raw/master/activity_level.csv")

#finds row in activity data frame corresponding to each occurrence
activity_row<-cbind(match(cleaned_occs$class,activity_level$taxon),match(cleaned_occs$order,activity_level$taxon),match(cleaned_occs$family,activity_level$taxon))
activity_row<-apply(activity_row,1,function(x) min(x,na.rm=T))

#adds activity quotient to occurrences
cleaned_occs$activity<-activity_level$activity_quotient[activity_row]

#activity score for each genus
genus_activity<-sapply(split(cleaned_occs$activity,cleaned_occs$matched_name),mean)

#counts number of occurrences per genus in each time interval
genus_matrix<-sapply(split(cleaned_occs,cleaned_occs$time_int),function(x) table(x$matched_name))

#orders genus matrix so that time intervals are in chronologic order
genus_matrix<-genus_matrix[,order(max(match(colnames(genus_matrix),time_int$interval_no))-match(colnames(genus_matrix),time_int$interval_no)+1)]

#LOGISTIC REGRESSION FOR SURVIVAL AS FUNCTION OF ACTIVITY
log_odds_ratio<-numeric(0)
log_odds_error<-numeric(0)
log_odds_p<-numeric(0)

#logistic regression based on modified three-timer extinction (three bin moving window)
for (i in 1:(ncol(genus_matrix)-2)) {
  moving_window<-genus_matrix[,seq(i,i+2)]
  moving_window<-subset(moving_window,moving_window[,1]>0 & moving_window[,2]>0) #examines cohort of taxa present in both time t-1 and time t
  extinct<-apply(moving_window,1,function(x) ifelse(x[3]>0,1,0)) #modified version of three-timer extinction; two-timers=extinct, three-timers=survive
  extinct_select<-data.frame(extinct,activity=genus_activity[match(names(extinct),names(genus_activity))])
  activity_glm<-glm(extinct~activity,data=extinct_select)
  log_odds_ratio[i]<-summary(activity_glm)$coefficients[2]
  log_odds_error[i]<-summary(activity_glm)$coefficients[4]
  log_odds_p[i]<-summary(activity_glm)$coefficients[8]
}

#finds top and midpoint ages and interval names
bin_top<-time_int$late_age[match(colnames(genus_matrix),time_int$interval_no)]
bin_midpt<-rowMeans(cbind(time_int$late_age[match(colnames(genus_matrix),time_int$interval_no)],time_int$early_age[match(colnames(genus_matrix),time_int$interval_no)]))
int_names<-time_int$interval_name[match(colnames(genus_matrix),time_int$interval_no)]

#compiles results into table with interval names and ages
activity_results<-data.frame(int=int_names[2:(length(int_names)-1)],top=bin_top[2:(length(bin_top)-1)],log_odds_ratio,log_odds_error,log_odds_p)

#prepares info for plot (axis limits and colors)
ylim_max<-max(activity_results$log_odds_ratio+1.96*activity_results$log_odds_error)
ylim_min<-min(activity_results$log_odds_ratio-1.96*activity_results$log_odds_error)
activity_results$color<-ifelse(activity_results$log_odds_p<0.05,ifelse(activity_results$log_odds_ratio>0,"red","blue"),"gray")

pdf("selectivity.pdf",width=24,height=16,pointsize=36)

plot(activity_results$top,activity_results$log_odds_ratio,ylim=c(ylim_min*1.2,ylim_max),xlim=rev(range(activity_results$top)),xlab="Age (Ma)",ylab="Log odds ratio",col=activity_results$color,pch=16,bty="n")
abline(h=0,lty=3)
segments(activity_results$top,activity_results$log_odds_ratio-1.96*activity_results$log_odds_error,activity_results$top,activity_results$log_odds_ratio+1.96*activity_results$log_odds_error,col=activity_results$color)

#adds geological timescale to plot
rect(bin_top[1:(length(bin_top)-2)],ylim_min*1.2,bin_top[2:(length(bin_top)-1)],1.2*ylim_min+0.05,col=paste(time_int$color[match(int_names,time_int$interval_name)])[2:length(bin_top)])
text(bin_midpt[2:(length(bin_midpt)-1)],ylim_min*1.2-0.007,strtrim(activity_results$int,-diff(bin_top[1:(length(bin_top)-1)])),pos=3,cex=0.45)

dev.off()



#adds activity quotient to occurrences
cleaned_occs$calcified<-activity_level$calcareous[activity_row]
cleaned_occs$mineralogy<-activity_level$calcite[activity_row]

#calcified for each genus
genus_calcified<-sapply(split(cleaned_occs$calcified,cleaned_occs$matched_name),mean)

#mineralogy for each genus
genus_mineralogy<-sapply(split(cleaned_occs$mineralogy,cleaned_occs$matched_name),mean)

#MULTIPLE LOGISTIC REGRESSION FOR SURVIVAL AS FUNCTION OF ACTIVITY AND SHELL MINERALOGY
calcified_odds_ratio<-numeric(0)
calcified_odds_error<-numeric(0)
calcified_odds_p<-numeric(0)

mineralogy_odds_ratio<-numeric(0)
mineralogy_odds_error<-numeric(0)
mineralogy_odds_p<-numeric(0)

#logistic regression for calcified and mineralogy based on modified three-timer extinction (three bin moving window)
for (i in 1:(ncol(genus_matrix)-2)) {
  moving_window<-genus_matrix[,seq(i,i+2)]
  moving_window<-subset(moving_window,moving_window[,1]>0 & moving_window[,2]>0) #examines cohort of taxa present in both time t-1 and time t
  extinct<-apply(moving_window,1,function(x) ifelse(x[3]>0,1,0)) #modified version of three-timer extinction; two-timers=extinct, three-timers=survive
  extinct_select<-data.frame(extinct,activity=genus_activity[match(names(extinct),names(genus_activity))],calcified=genus_calcified[match(names(extinct),names(genus_calcified))],mineralogy=genus_mineralogy[match(names(extinct),names(genus_mineralogy))])
  mineralogy_glm<-glm(extinct~activity+calcified+mineralogy,data=extinct_select)
  calcified_odds_ratio[i]<-summary(mineralogy_glm)$coefficients[3]
  calcified_odds_error[i]<-summary(mineralogy_glm)$coefficients[7]
  calcified_odds_p[i]<-summary(mineralogy_glm)$coefficients[15]
  mineralogy_odds_ratio[i]<-summary(mineralogy_glm)$coefficients[4]
  mineralogy_odds_error[i]<-summary(mineralogy_glm)$coefficients[8]
  mineralogy_odds_p[i]<-summary(mineralogy_glm)$coefficients[16]
}

#finds top and midpoint ages and interval names
bin_top<-time_int$late_age[match(colnames(genus_matrix),time_int$interval_no)]
bin_midpt<-rowMeans(cbind(time_int$late_age[match(colnames(genus_matrix),time_int$interval_no)],time_int$early_age[match(colnames(genus_matrix),time_int$interval_no)]))
int_names<-time_int$interval_name[match(colnames(genus_matrix),time_int$interval_no)]

#compiles results into table with interval names and ages
calcified_results<-data.frame(int=int_names[2:(length(int_names)-1)],top=bin_top[2:(length(bin_top)-1)],calcified_odds_ratio,calcified_odds_error,calcified_odds_p)
mineralogy_results<-data.frame(int=int_names[2:(length(int_names)-1)],top=bin_top[2:(length(bin_top)-1)],mineralogy_odds_ratio,mineralogy_odds_error,mineralogy_odds_p)


pdf("mineralogy.pdf",width=24,height=24,pointsize=36)
par(mfrow=c(2,1))
par(mar=c(3,4,1,4))

#removes rows where too few contrasts for logistic regression
calcified_results<-subset(calcified_results,is.na(calcified_results$calcified_odds_p)==0)

ylim_max<-max(calcified_results$calcified_odds_ratio+1.96*calcified_results$calcified_odds_error)
ylim_min<-min(calcified_results$calcified_odds_ratio-1.96*calcified_results$calcified_odds_error)
calcified_results$color<-ifelse(calcified_results$calcified_odds_p<0.05,ifelse(calcified_results$calcified_odds_ratio>0,"red","blue"),"gray")

plot(calcified_results$top,calcified_results$calcified_odds_ratio,ylim=c(ylim_min*1.2,ylim_max),xlim=rev(range(calcified_results$top)),xlab="Age (Ma)",ylab="Log odds ratio",col=calcified_results$color,pch=16,bty="n")
abline(h=0,lty=3)
segments(calcified_results$top,calcified_results$calcified_odds_ratio-1.96*calcified_results$calcified_odds_error,calcified_results$top,calcified_results$calcified_odds_ratio+1.96*calcified_results$calcified_odds_error,col=calcified_results$color)

#adds geological timescale
rect(bin_top[1:(length(bin_top)-2)],ylim_min*1.2,bin_top[2:(length(bin_top)-1)],1.2*ylim_min+0.1,col=paste(time_int$color[match(int_names,time_int$interval_name)])[2:length(bin_top)])
text(bin_midpt[2:(length(bin_midpt)-1)],ylim_min*1.2-0.08,strtrim(int_names[2:(length(int_names)-1)],-diff(bin_top[1:(length(bin_top)-1)])),pos=3,cex=0.45)

#removes rows where too few contrasts for logistic regression
mineralogy_results<-subset(mineralogy_results,is.na(mineralogy_results$mineralogy_odds_p)==0)

ylim_max<-max(mineralogy_results$mineralogy_odds_ratio+1.96*mineralogy_results$mineralogy_odds_error)
ylim_min<-min(mineralogy_results$mineralogy_odds_ratio-1.96*mineralogy_results$mineralogy_odds_error)
mineralogy_results$color<-ifelse(mineralogy_results$mineralogy_odds_p<0.05,ifelse(mineralogy_results$mineralogy_odds_ratio>0,"red","blue"),"gray")

plot(mineralogy_results$top,mineralogy_results$mineralogy_odds_ratio,ylim=c(ylim_min*1.2,ylim_max),xlim=rev(range(mineralogy_results$top)),xlab="Age (Ma)",ylab="Log odds ratio",col=mineralogy_results$color,pch=16,bty="n")
abline(h=0,lty=3)
segments(mineralogy_results$top,mineralogy_results$mineralogy_odds_ratio-1.96*mineralogy_results$mineralogy_odds_error,mineralogy_results$top,mineralogy_results$mineralogy_odds_ratio+1.96*mineralogy_results$mineralogy_odds_error,col=mineralogy_results$color)

#adds geological timescale
rect(bin_top[1:(length(bin_top)-2)],ylim_min*1.2,bin_top[2:(length(bin_top)-1)],1.2*ylim_min+0.07,col=paste(time_int$color[match(int_names,time_int$interval_name)])[2:length(bin_top)])
text(bin_midpt[2:(length(bin_midpt)-1)],ylim_min*1.2-0.06,strtrim(int_names[2:(length(int_names)-1)],-diff(bin_top[1:(length(bin_top)-1)])),pos=3,cex=0.45)

dev.off()

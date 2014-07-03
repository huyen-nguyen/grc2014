#reads extinction data from file produced at fossilworks.org
ext_data<-read.csv("https://github.com/mclapham/grc2014/raw/master/subsamp_ext.csv")

#reads time intervals from PBDB data services
time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=1&min_ma=160&max_ma=325")
time_int<-subset(time_int,time_int$level==5)

stage_dur<-diff(ext_data$Base_ma) #stage duration

#plot extinction rates over time
pdf("extinction.pdf",width=24,height=16,pointsize=36)

plot(ext_data$Base_ma[1:(length(ext_data$Base_ma)-1)],ext_data$extinction_3T[2:length(ext_data$extinction_3T)],type="o",pch=16,cex=1.25,col="red",lwd=3,xlim=c(305,165),ylim=c(-0.05,max(ext_data$extinction_3T,na.rm=T)),xlab="Age (Ma)",ylab="Extinction rate",bty="n")

#adds rectangles for geological timescale
rect(ext_data$Base_ma[1:(length(ext_data$Base_ma)-1)],-0.1,ext_data$Base_ma[2:length(ext_data$Base_ma)],0,col=paste(time_int$color)[2:length(ext_data$Base_ma)])
text(ext_data$Midpoint_ma[2:(length(ext_data$Midpoint_ma)-2)],-0.125,strtrim(time_int$interval_name[2:(length(time_int$interval_name)-2)],diff(ext_data$Base_ma)),pos=3,cex=0.48)

dev.off()

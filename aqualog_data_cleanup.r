options(scipen=999)



convert_aqualog<-function(){
library(tcltk)
#main_folder<-"/Users/prime/gdrive/uf_aqualog/uf_aqualog/20190408/20190408_Group1"
choose.dir <- function() {
	system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
			intern = FALSE, ignore.stderr = TRUE)
	p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
	return(ifelse(length(p), p, NA))
}

#choose folder where corrected data will be stored.
cat("\nchoose destination folder")
destination_folder<-choose.dir()

x<-"yes"
while(x=="yes"){
#choose folder where original aqualog data are stored.
cat("\nchoose folder where original aqualog data are stored.")
main_folder<-choose.dir()

abs_files<-list.files(main_folder,pattern="*ABS.dat")
#abs_file<-read.table("/Users/prime/gdrive/uf_aqualog/uf_aqualog/20190408/20190408_Group1/abs/B1S15S15ABS.csv")

#gsub(".*/(.*?)/$", "\\1", main_folder)
#dir.create(paste(destination_folder,gsub(".*/(.*?)/$", "\\1", main_folder),sep=""))
#dir.create(paste(destination_folder,gsub(".*/(.*?)/$", "\\1", main_folder),"/abs",sep=""))
#dir.create(paste(destination_folder,gsub(".*/(.*?)/$", "\\1", main_folder),"/eem",sep=""))
dir.create(paste(destination_folder,"/abs",sep=""))
dir.create(paste(destination_folder,"/eem",sep=""))
dir.create(paste(destination_folder,"/abs/",gsub(".*/(.*?)/$", "\\1", main_folder),sep=""))
dir.create(paste(destination_folder,"/eem/",gsub(".*/(.*?)/$", "\\1", main_folder),sep=""))
abs_cor_folder<-paste(destination_folder,"/abs/",gsub(".*/(.*?)/$", "\\1", main_folder),sep="")
eem_cor_folder<-paste(destination_folder,"/eem/",gsub(".*/(.*?)/$", "\\1", main_folder),sep="")


for(i in 1:length(abs_files)){
	abs_temp<-read.table(paste(main_folder,"/",abs_files[i],sep=""),sep="\t")
	abs_temp<-subset(abs_temp,!is.na(V2))
	write.table(abs_temp,paste(abs_cor_folder,"/",gsub("ABS.dat",".csv",abs_files[i],fixed=TRUE),sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}


sem_files<-list.files(main_folder,pattern="*SEM.dat")

for(i in 1:length(sem_files)){
	abs_temp<-read.table(paste(main_folder,"/",sem_files[i],sep=""),sep="\t")
	abs_col<-ncol(abs_temp)
	abs_temp2<-abs_temp[,2:abs_col]
	abs_temp2<-abs_temp2[,order(abs_temp2[1,])]
	abs_temp3<-cbind.data.frame(abs_temp[,1],abs_temp2)
	write.table(abs_temp3,paste(eem_cor_folder,"/",gsub("SEM.dat",".csv",sem_files[i],fixed=TRUE),sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}


bem_files<-list.files(main_folder,pattern="*BEM.dat")
bem_temp<-read.table(paste(main_folder,"/",bem_files[1],sep=""),sep="\t")
bem_col<-ncol(bem_temp)
bem_temp2<-bem_temp[,2:bem_col]
bem_temp2<-bem_temp2[,order(bem_temp2[1,])]
bem_temp3<-cbind.data.frame(bem_temp[,1],bem_temp2)
write.table(bem_temp3,paste(eem_cor_folder,"/",gsub("BEM.dat","_blank.csv",bem_files[i],fixed=TRUE),sep=""),row.names=FALSE,col.names=FALSE,sep=",")

cat("\nCorrect Another Folder?")
x<-readLines(n=1)
print(x)
}
}

convert_aqualog()
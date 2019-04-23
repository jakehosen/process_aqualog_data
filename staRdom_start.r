#check this tutorial for reference: https://cran.r-project.org/web/packages/staRdom/vignettes/PARAFAC_analysis_of_EEM.html

#Run these once!
#install.packages("devtools") # Run this only, if devtools is not installed already.
#devtools::install_github("jakehosen/stardom_aqualog")

choose.dir <- function() {
	system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
			intern = FALSE, ignore.stderr = TRUE)
	p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
	return(ifelse(length(p), p, NA))
}

library(staRdom)

  
#correction files
#excorfile <- system.file("extdata/CorrectionFiles/xc06se06n.csv",package="staRdom")
#excorfile<-"/Users/prime/gdrive/uf_aqualog/uf_aqualog/xcorrect_uf.csv"
cat("identify xcorrect file")
excorfile<-file.choose()
Excor <- data.table::fread(excorfile)
Excor<-Excor[order(Excor$V1),]

#emcorfile <- system.file("extdata/CorrectionFiles/mcorrs_4nm.csv",package="staRdom")
#emcorfile<-"/Users/prime/gdrive/uf_aqualog/uf_aqualog/mcorrect_uf.csv"
cat("identify xcorrect file")
emcorfile<-file.choose()
Emcor <- data.table::fread(emcorfile)
Emcor<-Emcor[order(Emcor$V1),]




data_folder<-choose.dir()
#data_folder<-"/Users/prime/gdrive/uf_aqualog/uf_aqualog/20190408/20190408_Group1/"
cores <- detectCores(logical = FALSE)
#folder <- system.file("extdata/EEMs/", package = "staRdom") # folder containing example EEMs
eem_folder<-paste(data_folder,"/eem",sep="") 
eem_list <- eem_read_csv(eem_folder) # in case you use your own data, just replace folder by a path. e.g. "C:/folder/another folder" and use eem_read if you do not have plain csv tables to import

#eem_list <- eem_read(eem_folder)

#eem_overview_plot(eem_list, spp=8)

#absorbance_path = system.file("extdata/absorbance", package = "staRdom") # load example data, set a path without using system.file to use your own data e.g. "C:/folder/another folder"


absorbance_path<-paste(data_folder,"/abs",sep="")
absorbance <- absorbance_read(absorbance_path) # load csv or txt tables in folder
absorbance <- abs_blcor(absorbance,wlrange = c(680,700))

#metatable <- system.file("extdata/metatable_dreem.csv",package = "staRdom") # path to example data, can be replaced by a path to your own data
meta <- read.table(paste(data_folder,"/metatable.csv",sep=""), header = TRUE, sep = ",", dec = ".", row.names = 1) # load data
dilution = "meta" # e.g. 1 for undiluted samples
dil_sample_name_correction = FALSE

#eem_listtemplate(eem_list, absorbance) %>%
#  write.csv(file="metatable.csv", row.names = FALSE)
  
  



# adjust range of EEMs to cover correction vectors
#eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1]))
#eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1][1],780))
eem_list <- eem_range(eem_list,ex = c(250,600), em = range(Emcor[,1][1],600))

eem_list <- eem_spectral_cor(eem_list,Excor,Emcor)

eem_list2 <- eem_extend2largest(eem_list, interpolation = 1, extend = FALSE, cores = cores)

# blank subtraction
eem_list3 <- eem_remove_blank(eem_list2)


  
#eem_overview_plot(eem_list2, spp=8)

eem_list4 <- eem_ife_correction(eem_list3,absorbance, cuvl = 5)




eem_list4b <- eem_raman_normalisation2(eem_list4, blank = "blank")
#eem_list <- eem_raman_normalisation(eem_list)

eem_list5 <- eem_extract(eem_list4b, c("nano", "miliq", "milliq", "mq", "blank","B1S13SPB04No28","B1S14SPB04No33","B1S14SPB04No29"),ignore_case = TRUE)

absorbance <- dplyr::select(absorbance, -dplyr::matches("nano|miliq|milliq|mq|blank|B1S13SPB04No28|B1S14SPB04No33|B1S14SPB04No29", ignore.case = TRUE))


remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(20,15,30,30)

eem_list_rem <- eem_rem_scat(eem_list5, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)
eem_overview_plot(eem_list_rem, spp=3)


eem_list_rem_int <- eem_interp(eem_list_rem, cores = cores, type = 1, extend = FALSE)
eem_overview_plot(eem_list_rem_int, spp=3)

dil_data <- meta["dilution"]

eem_list_dil <- eem_dilution(eem_list_rem_int,dil_data)

write.csv(cbind.data.frame(eem_biological_index(eem_list_dil),eem_coble_peaks(eem_list_dil),eem_fluorescence_index(eem_list_dil),eem_fluorescence_index(eem_list_dil)),paste(data_folder,"basic_fl_stats.csv",sep=""),row.names=FALSE)

# minimum and maximum of numbers of components
dim_min <- 3
dim_max <- 7

nstart <- 20 # number of similar models from which best is chosen
maxit = 5000 # maximum number of iterations in PARAFAC analysis
ctol <- 10^-5 # tolerance in PARAFAC analysis

# calculating PARAFAC models, one for each number of components
#pf1 <- eem_parafac(eem_list_rem_int, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("uncons", "uncons", "uncons"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)

exclude1<-list(sample=c("B1S15S15","B1S1S1","B1S2S2","B1S3S3","B1S4S4","B1S5S5"))
eem_list_dil2<-eem_exclude(eem_list_dil, exclude1)

# same model but using non-negative constraints
pf1n <- eem_parafac(eem_list_dil2, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)

# rescale B and C modes to a maximum fluorescence of 1 for each component
pf1 <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
pf1n <- lapply(pf1n, eempf_rescaleBC, newscale = "Fmax")



#compare models
eempf_compare(pf1n)
eempf_corplot(pf1n[[4]], progress = FALSE)
# plot leverage (nice plot)

# calculate leverage
cpl <- eempf_leverage(pf1n[[4]])
eempf_leverage_plot(cpl,qlabel=0.1)

exclude <- eempf_leverage_ident(cpl,qlabel=0.1)

eem_list_dil2_ex1 <- eem_exclude(eem_list_dil2, exclude)

pf2n <- eem_parafac(eem_list_dil2_ex1, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
eempf_compare(pf2n)

eempf_compare(pf1n)





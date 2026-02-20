#################
### Libraries ###
#################

library(ncdf4)
library(rworldmap)
library(RColorBrewer)
library(raster)
library(npreg)
library(fields)
library(RColorBrewer)






#################
### Datasets #####
#################

aegtoday<-raster("TCurMean30Sum_97ae.tif")
albtoday<-raster("TCurMean30Sum_97al.tif")
aeg26<-raster("rcp26_2050_HE_sum_97ae.tif")
aeg85<-raster("rcp85_2050_HE_sum_97ae.tif")
alb26<-raster("rcp26_2050_HE_sum_97al.tif")
alb85<-raster("rcp85_2050_HE_sum_97al.tif")


### Calculate max of albopictus and aegypti in each location
s <- stack(aegtoday, albtoday)
maxAedesMap <- calc(s, fun=function(x) max(x,na.rm=T))
s <- stack(aeg26, alb26)
maxAedesMap26 <- calc(s, fun=function(x) max(x,na.rm=T))
s <- stack(aeg85, alb85)
maxAedesMap85 <- calc(s, fun=function(x) max(x,na.rm=T))

combMax<-list()
combMax[[1]]<-maxAedesMap
combMax[[2]]<-maxAedesMap26
combMax[[3]]<-maxAedesMap85

rm(aegtoday)
rm(albtoday)
rm(aeg26)
rm(aeg85)
rm(alb26)
rm(alb85)

### Bring in population distribution
pop2050ssp2<-nc_open('population_ssp2soc_2.5min_annual_2006-2100.nc4')
fillvalue <- ncatt_get(pop2050ssp2, "number_of_people", "_FillValue")
lon <- ncvar_get(pop2050ssp2, "lon")
lat <- ncvar_get(pop2050ssp2, "lat", verbose = F)
t <- ncvar_get(pop2050ssp2, "time")
year2050<-which(floor(1661+t)==2050)
year2025<-which(floor(1661+t)==2025)

v3      <- pop2050ssp2$var[[1]]
varsize <- v3$varsize
ndims   <- v3$ndims

start <- rep(1,ndims)   # begin with start=(1,1,1,...,1)
start[ndims] <-  year2050
count <- varsize        # begin w/count=(nx,ny,nz,...,nt), reads entire var
count[ndims] <- 1
pop.array <- ncvar_get(pop2050ssp2, "number_of_people",start=start,count=count)
pop.array[pop.array == fillvalue$value] <- NA
mapssp2 <- raster(t(pop.array), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

start[ndims] <-  year2025
pop.array <- ncvar_get(pop2050ssp2, "number_of_people",start=start,count=count)
pop.array[pop.array == fillvalue$value] <- NA
mapssp2_2025 <- raster(t(pop.array), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

nc_close(pop2050ssp2)

pop2050ssp5<-nc_open('population_ssp5soc_2.5min_annual_2006-2100.nc4')
fillvalue <- ncatt_get(pop2050ssp5, "number_of_people", "_FillValue")
lon <- ncvar_get(pop2050ssp5, "lon")
lat <- ncvar_get(pop2050ssp5, "lat", verbose = F)
t <- ncvar_get(pop2050ssp5, "time")
year2050<-which(floor(1661+t)==2050)

v3      <- pop2050ssp5$var[[1]]
varsize <- v3$varsize
ndims   <- v3$ndims

start <- rep(1,ndims)   # begin with start=(1,1,1,...,1)
start[ndims] <-  year2050
count <- varsize        # begin w/count=(nx,ny,nz,...,nt), reads entire var
count[ndims] <- 1
pop.array <- ncvar_get(pop2050ssp5, "number_of_people",start=start,count=count)
pop.array[pop.array == fillvalue$value] <- NA
mapssp5 <- raster(t(pop.array), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

nc_close(pop2050ssp5)

popDat<-list()
popDat[[1]]<-mapssp2_2025
popDat[[2]]<-mapssp2
popDat[[3]]<-mapssp5

### Analysis of countries
world <- getMap(resolution = "low")
minPop<-200000
uniqueCountries<-as.character(world@data$NAME_SORT[which(world@data$POP_EST>=minPop)])
a<-c("Cyprus, Northern","Gaza Strip","Hong Kong S.A.R.","Macau S.A.R","West Bank","Somaliland","Kosovo","Taiwan","Vanuatu")
uniqueCountries<-uniqueCountries[-which(uniqueCountries%in%a)] ### Removes places where we don't have demography estimates or baseline epi estimates
EstCountries<-world[match(uniqueCountries,as.character(world@data$NAME_SORT)),]


### Demography of countries
demogDat<-read.csv("wcde_data_NAMESMATCHED.csv")
originalEstimates<-read.csv("epiStatusByCountry.csv")[,-1]
originalEstimates<-originalEstimates[match(uniqueCountries,originalEstimates[,1]),]

### Age categories and estimates of disease by age
maxAge<-100
AgesMax<-c(1,9,19,29,39,49,59,69,maxAge)
AgesMin<-c(0,1,10,20,30,40,50,60,70)
ageIndex<-c(1,9,10,10,10,10,10,10,maxAge-70)

probSymptom<-c(0.104,0.0471,0.0493,0.0574,0.0614,0.0617,0.0696,0.0995,0.143)
probArthralgia<-probSymptom/2
deaths<-c(50,2,3,1,7,11,27,39,158)
infecs<-c(44800,399000,403000,395000,357000,254000,192000,126000,90200)
IFR<-deaths/infecs
IFR<-c(0.112,0.000501,0.000745,0.000253,0.00196,0.00433,0.0140,0.031,0.175)/100

probSymptom2<-rep(probSymptom,ageIndex)
probArthralgia2<-rep(probArthralgia,ageIndex)
IFR2<-rep(IFR,ageIndex)

DiseaseMat<-cbind(probSymptom2,probArthralgia2,IFR2)





##################################
###### ANALYSIS ##################
##################################


### Put age data into right bins
agebins<-unique(demogDat$Age)
agebins<-agebins[-which(agebins=="All")]
agebinWidth<-rep(5,length(agebins))

popByageByCountryBase<-array(NaN,c(length(uniqueCountries),length(agebins),3))
for (i in 1:length(agebins)){
	for(j in 1:length(uniqueCountries)){
		a<-which(demogDat$Age==agebins[i]&demogDat$Area==uniqueCountries[j]&demogDat$Scenario=="SSP2"&demogDat$Year=="2020")
		popByageByCountryBase[j,i,1]<-demogDat$Population[a]

		a<-which(demogDat$Age==agebins[i]&demogDat$Area==uniqueCountries[j]&demogDat$Scenario=="SSP2"&demogDat$Year=="2050")
		popByageByCountryBase[j,i,2]<-demogDat$Population[a]

		a<-which(demogDat$Age==agebins[i]&demogDat$Area==uniqueCountries[j]&demogDat$Scenario=="SSP5"&demogDat$Year=="2050")
		popByageByCountryBase[j,i,3]<-demogDat$Population[a]
	}
}

popByageByCountry<-array(NaN,c(length(uniqueCountries),maxAge,3))
for (i in 1:length(uniqueCountries)){
	for (j in 1:3){
		tmp<-as.numeric(popByageByCountryBase[i,,j])
		popByage<-rep(tmp/5,each=5)
		popByageByCountry[i,1:(maxAge-1),j]<-popByage[1:(maxAge-1)]
		popByageByCountry[i,maxAge,j]<-(sum(popByage)-sum(popByage[1:(maxAge-1)]))
	}
}	

popByageByCountry<-popByageByCountry*1000
totByCountry<-apply(popByageByCountry,c(1,3),"sum")

### Population structure for sensitivity analyses
tmp1<-sweep(popByageByCountry[,,2],1,rowSums(popByageByCountry[,,2]),"/")
tmp2<-sweep(popByageByCountry[,,3],1,rowSums(popByageByCountry[,,3]),"/")
popByageByCountry2050_tot2025<-popByageByCountry
popByageByCountry2050_tot2025[,,2]<-sweep(tmp1,1,rowSums(popByageByCountry[,,1]),"*")
popByageByCountry2050_tot2025[,,3]<-sweep(tmp2,1,rowSums(popByageByCountry[,,1]),"*")
### Total is same as today but structure different

tmp1<-sweep(popByageByCountry[,,1],1,rowSums(popByageByCountry[,,1]),"/")
popByageByCountry_tot2050<-popByageByCountry
popByageByCountry_tot2050[,,2]<-sweep(tmp1,1,rowSums(popByageByCountry[,,2]),"*")
popByageByCountry_tot2050[,,3]<-sweep(tmp1,1,rowSums(popByageByCountry[,,3]),"*")
### Total pop is same as 2050 but structure like today



### Calculate max Aedes in each country weighted by population distribution
maxAegAlb<-rep(NaN,length(uniqueCountries),1)
for (i in 1:length(uniqueCountries)){
	map1<-world[which(as.character(world@data$NAME_SORT)%in%uniqueCountries[i]),]
	tmp<-extract(combMax[[1]],map1, cellnumbers=TRUE)
	tmp[[1]][,2][which(is.na(tmp[[1]][,2]))]<-0 ### Places with no suitability have NA - replace with 0
	tmp[[1]][,2][which(tmp[[1]][,2]==-Inf)]<-0
	coords<-xyFromCell(combMax[[1]], cell=tmp[[1]][,1])

	AedesMaxBase<-extract(combMax[[1]],coords)	
	popsize<-extract(popDat[[1]],coords)
	maxAegAlb[i]<-weighted.mean(tmp[[1]][,2],popsize,na.rm=TRUE)
	print(i)
}


### Look at relationship with ever outbreaks
originalEstimates$maxAegAlb<-maxAegAlb
maxmosq<-seq(2,12,2)
minmosq<-maxmosq-2
midmosq<-(minmosq+maxmosq)/2
probAegAlbwithCI<-matrix(NaN,length(maxmosq),3)
for(i in 1:length(maxmosq)){
	a<-which(originalEstimates$maxAegAlb>=minmosq[i]&originalEstimates$maxAegAlb<maxmosq[i])
	if(length(a)<4)next
	b<-prop.test(sum(originalEstimates$EverCases[a]=="TRUE"),length(a))
	probAegAlbwithCI[i,1]<-b$estimate
	probAegAlbwithCI[i,2:3]<-b$conf.int
}


### Regression of ever had outbreaks with population level average 
x<-originalEstimates$maxAegAlb
y<-originalEstimates$EverCases
xs<-seq(0,12,0.1)
smod <- sm(y ~ x)
pred<-predict(smod,newdata=data.frame(x=xs))


### Predict effective population size in each country
popSizeCountry<-monthsMosq<-effectivePopSize<-propAtRisk<-matrix(0,length(uniqueCountries),3)
for (i in 1:length(uniqueCountries)){
	map1<-world[which(as.character(world@data$NAME_SORT)%in%uniqueCountries[i]),]

	popBase<-extract(popDat[[1]],map1, cellnumbers=TRUE)
	pop26<-extract(popDat[[2]],map1, cellnumbers=TRUE)
	pop85<-extract(popDat[[3]],map1, cellnumbers=TRUE)

	popSizeCountry[i,1]<-sum(popBase[[1]][,2],na.rm=T)
	popSizeCountry[i,2]<-sum(pop26[[1]][,2],na.rm=T)
	popSizeCountry[i,3]<-sum(pop85[[1]][,2],na.rm=T)

	coords<-xyFromCell(popDat[[1]], cell=popBase[[1]][,1])

	AedesMaxBase<-extract(combMax[[1]],coords)
	AedesMax26<-extract(combMax[[2]],coords)
	AedesMax85<-extract(combMax[[3]],coords)
	AedesMaxBase[which(is.na(AedesMaxBase))]<-0
	AedesMax26[which(is.na(AedesMax26))]<-0
	AedesMax85[which(is.na(AedesMax85))]<-0

	exc<-which(AedesMax26==0&AedesMaxBase>10) ## Clear discrepancies where new versions of map have gaps. Only relevant in a few places.
	if(length(exc)>0){
		AedesMaxBase<-AedesMaxBase[-exc]
		AedesMax26<-AedesMax26[-exc]
		AedesMax85<-AedesMax85[-exc]
		popBase[[1]]<-popBase[[1]][-exc,]
		pop26[[1]]<-pop26[[1]][-exc,]
		pop85[[1]]<-pop85[[1]][-exc,]
	}

	if(sum(AedesMaxBase,na.rm=T)>0){
		pred1<-predict(smod,newdata=data.frame(x=AedesMaxBase))
		monthsMosq[i,1]<-weighted.mean(AedesMaxBase,popBase[[1]][,2],na.rm=TRUE)
		effectivePopSize[i,1]<-sum(popBase[[1]][,2]*pred1,na.rm=TRUE)
		propAtRisk[i,1]<-effectivePopSize[i,1]/sum(popBase[[1]][,2],na.rm=T)
	}
	if(sum(AedesMax26,na.rm=T)>0){
		pred1<-predict(smod,newdata=data.frame(x=AedesMax26))
		monthsMosq[i,2]<-weighted.mean(AedesMax26,pop26[[1]][,2],na.rm=TRUE)
		effectivePopSize[i,2]<-sum(pop26[[1]][,2]*pred1,na.rm=TRUE)
		propAtRisk[i,2]<-effectivePopSize[i,2]/sum(pop26[[1]][,2],na.rm=T)
	}
	if(sum(AedesMax85,na.rm=T)>0){
		pred1<-predict(smod,newdata=data.frame(x=AedesMax85))
		monthsMosq[i,3]<-weighted.mean(AedesMax85,pop85[[1]][,2],na.rm=TRUE)
		effectivePopSize[i,3]<-sum(pop85[[1]][,2]*pred1,na.rm=TRUE)
		propAtRisk[i,3]<-effectivePopSize[i,3]/sum(pop85[[1]][,2],na.rm=T)
		
	}
	print(i)
}


#### Estimate case numbers function
caseNumbers<-function(foi=0.02,yearsCirculating=20,diseaseMat=DiseaseMat,maxAge=maxAge, foi_base=0.02){

	ageBase<-1:maxAge
	ageBase[which(ageBase>yearsCirculating)]<-yearsCirculating
	infections<-exp(-foi_base*ageBase)*foi
	disease<-sweep(diseaseMat,1,infections,"*")

	outMat<-cbind(infections,disease)
	return(outMat)
}


#### Estimate burden by country
outputMat_Today<-outputMat_26<-outputMat_85<-outputMat_26_pop2025<-outputMat_85_pop2025<-outputMat_26_popStructureSens<-outputMat_85_popStructureSens<-outputMat_26_popSizeSens<-outputMat_85_popSizeSens<-matrix(NaN,length(uniqueCountries),4)
infecMod_epidemic<-caseNumbers(foi=0.016,yearsCirculating=30,foi_base=0.016,maxAge=maxAge)
infecMod_endemic<-caseNumbers(foi=0.021,yearsCirculating=30, foi_base=0.021,maxAge=maxAge)
for (i in 1:length(uniqueCountries)){
	epiStatus<-originalEstimates$status[which(originalEstimates$Country==uniqueCountries[i])]

	infecModTMP<-if(epiStatus=="Endemic")(infecMod_endemic)else(infecMod_epidemic)
	
	outputMat_Today[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry[i,,1]*propAtRisk[i,1],"*"))
	outputMat_26[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry[i,,2]*propAtRisk[i,2],"*"))
	outputMat_85[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry[i,,3]*propAtRisk[i,3],"*"))

	outputMat_26_pop2025[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry[i,,1]*propAtRisk[i,2],"*"))
	outputMat_85_pop2025[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry[i,,1]*propAtRisk[i,3],"*"))

	outputMat_26_popStructureSens[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry2050_tot2025[i,,2]*propAtRisk[i,1],"*"))
	outputMat_85_popStructureSens[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry2050_tot2025[i,,3]*propAtRisk[i,1],"*"))

	outputMat_26_popSizeSens[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry_tot2050[i,,2]*propAtRisk[i,1],"*"))
	outputMat_85_popSizeSens[i,]<-colSums(sweep(infecModTMP,1,popByageByCountry_tot2050[i,,3]*propAtRisk[i,1],"*"))
}
diffByCountry_26<-outputMat_26-outputMat_Today
diffByCountry_85<-outputMat_85-outputMat_Today


#Continental average burden
nameslist<-c("Europe","North America","South America","Asia","Africa","Australia")
outputMat_TodayCont<-outputMat_26Cont<-outputMat_85Cont<-outputMat_26_pop2025Cont<-outputMat_85_pop2025Cont<-outputMat_26_popStructureSensCont<-outputMat_85_popStructureSensCont<-outputMat_26_popSizeSensCont<-outputMat_85_popSizeSensCont<-matrix(NaN,6,4)
for(ii in 1:6){
	outputMat_TodayCont[ii,]<-colSums(outputMat_Today[which(EstCountries@data$REGION==nameslist[ii]),])
	outputMat_26Cont[ii,]<-colSums(outputMat_26[which(EstCountries@data$REGION==nameslist[ii]),])
	outputMat_85Cont[ii,]<-colSums(outputMat_85[which(EstCountries@data$REGION==nameslist[ii]),])

	outputMat_26_pop2025Cont[ii,]<-colSums(outputMat_26_pop2025[which(EstCountries@data$REGION==nameslist[ii]),])
	outputMat_85_pop2025Cont[ii,]<-colSums(outputMat_85_pop2025[which(EstCountries@data$REGION==nameslist[ii]),])
	outputMat_26_popStructureSensCont[ii,]<-colSums(outputMat_26_popStructureSens[which(EstCountries@data$REGION==nameslist[ii]),])
	outputMat_85_popStructureSensCont[ii,]<-colSums(outputMat_85_popStructureSens[which(EstCountries@data$REGION==nameslist[ii]),])
	outputMat_26_popSizeSensCont[ii,]<-colSums(outputMat_26_popSizeSens[which(EstCountries@data$REGION==nameslist[ii]),])
	outputMat_85_popSizeSensCont[ii,]<-colSums(outputMat_85_popSizeSens[which(EstCountries@data$REGION==nameslist[ii]),])
}
ratByCont_26<-outputMat_26Cont/outputMat_TodayCont
ratByCont_85<-outputMat_85Cont/outputMat_TodayCont

diffByCont_26<-outputMat_26Cont-outputMat_TodayCont
diffByCont_85<-outputMat_85Cont-outputMat_TodayCont

outputMat_TodayOverall<-colSums(outputMat_Today)
outputMat_26Overall<-colSums(outputMat_26)
outputMat_85Overall<-colSums(outputMat_85)
outputMat_26_pop2025Overall<-colSums(outputMat_26_pop2025)
outputMat_85_pop2025Overall<-colSums(outputMat_85_pop2025)
outputMat_26_popStructureSensOverall<-colSums(outputMat_26_popStructureSens)
outputMat_85_popStructureSensOverall<-colSums(outputMat_85_popStructureSens)
outputMat_26_popSizeSensOverall<-colSums(outputMat_26_popSizeSens)
outputMat_85_popSizeSensOverall<-colSums(outputMat_85_popSizeSens)

## Percentage change versus today
CC26<-(outputMat_26_pop2025Overall-outputMat_TodayOverall)/outputMat_TodayOverall
CC85<-(outputMat_85_pop2025Overall-outputMat_TodayOverall)/outputMat_TodayOverall
PopSize2050<-(outputMat_26_popSizeSensOverall-outputMat_TodayOverall)/outputMat_TodayOverall
PopStructure2050<-(outputMat_26_popStructureSensOverall-outputMat_TodayOverall)/outputMat_TodayOverall
combined26<-(outputMat_26Overall-outputMat_TodayOverall)/outputMat_TodayOverall
combined85<-(outputMat_85Overall-outputMat_TodayOverall)/outputMat_TodayOverall
plotdat<-100*cbind(CC26,CC85,PopSize2050,PopStructure2050,combined26,combined85)

## Same but at continental level
CC26Cont<-(outputMat_26_pop2025Cont-outputMat_TodayCont)/outputMat_TodayCont
CC85Cont<-(outputMat_85_pop2025Cont-outputMat_TodayCont)/outputMat_TodayCont
PopSize2050Cont<-(outputMat_26_popSizeSensCont-outputMat_TodayCont)/outputMat_TodayCont
PopStructure2050Cont<-(outputMat_26_popStructureSensCont-outputMat_TodayCont)/outputMat_TodayCont
combined26Cont<-(outputMat_26Cont-outputMat_TodayCont)/outputMat_TodayCont
combined85Cont<-(outputMat_85Cont-outputMat_TodayCont)/outputMat_TodayCont


#### Different FOI scenarios by continents
relFOI<-c(0.2,0.5,1,1.2,1.5, 2)
sensitivityDeaths85<-sensitivityCases85<-sensitivityDeaths26<-sensitivityCases26<-matrix(NaN,length(uniqueCountries),length(relFOI))
for(ii in 1:length(relFOI)){
	infecMod_epidemicTMP<-caseNumbers(foi=0.016*relFOI[ii],yearsCirculating=30,foi_base=0.016,maxAge=maxAge)
	infecMod_endemicTMP<-caseNumbers(foi=0.021*relFOI[ii],yearsCirculating=30, foi_base=0.021,maxAge=maxAge)
	for (i in 1:length(uniqueCountries)){
		epiStatus<-originalEstimates$status[which(originalEstimates$Country==uniqueCountries[i])]
		if(epiStatus=="Endemic"){
			tmp<-colSums(sweep(infecMod_endemicTMP,1,popByageByCountry[i,,2]*propAtRisk[i,2],"*"))
			sensitivityCases26[i,ii]<-tmp[2]
			sensitivityDeaths26[i,ii]<-tmp[4]

			tmp<-colSums(sweep(infecMod_endemicTMP,1,popByageByCountry[i,,3]*propAtRisk[i,3],"*"))
			sensitivityCases85[i,ii]<-tmp[2]
			sensitivityDeaths85[i,ii]<-tmp[4]
		}else{
			tmp<-colSums(sweep(infecMod_epidemicTMP,1,popByageByCountry[i,,2]*propAtRisk[i,2],"*"))
			sensitivityCases26[i,ii]<-tmp[2]
			sensitivityDeaths26[i,ii]<-tmp[4]

			tmp<-colSums(sweep(infecMod_epidemicTMP,1,popByageByCountry[i,,3]*propAtRisk[i,3],"*"))
			sensitivityCases85[i,ii]<-tmp[2]
			sensitivityDeaths85[i,ii]<-tmp[4]
		}
	}
	print(ii)
}

sensitivityDeaths85Cont<-sensitivityCases85Cont<-sensitivityDeaths26Cont<-sensitivityCases26Cont<-matrix(NaN,6,length(relFOI))
for(i in 1:length(relFOI)){
	for(ii in 1:6){
		sensitivityCases26Cont[ii,i]<-sum(sensitivityCases26[which(EstCountries@data$REGION==nameslist[ii]),i])
		sensitivityCases85Cont[ii,i]<-sum(sensitivityCases85[which(EstCountries@data$REGION==nameslist[ii]),i])
		sensitivityDeaths26Cont[ii,i]<-sum(sensitivityDeaths26[which(EstCountries@data$REGION==nameslist[ii]),i])
		sensitivityDeaths85Cont[ii,i]<-sum(sensitivityDeaths85[which(EstCountries@data$REGION==nameslist[ii]),i])
	}
}


##################
#### Vaccines ####
##################

foi_reduction_epidemic<-0.05
foi_reduction_endemic<-0.05
vaccAgeMin<-c(12, 12, 12, 12,1)
coverage<-c(0.5, 0.5, 0.5, 0.9,0.9)
protectionInfectionDirect<-0.4
protectionDisease<-c(0.7,0.98,0.7,0.7,0.98) ### Overall
protectionDiseaseNet<-(protectionDisease-protectionInfectionDirect)/(1-protectionInfectionDirect) #Reflects those infected who don't get sick

VaccinationDelay<-c(0.2,0.2,0.01,0.2,0) ##Proportion of infections that have already occurred during outbreak by the time you vaccinate (epidemic only)

infecMod_epidemic_V<-caseNumbers(foi=0.016*(1-foi_reduction_epidemic),yearsCirculating=30,foi_base=0.016,maxAge=maxAge)
infecMod_endemic_V<-caseNumbers(foi=0.021*(1-foi_reduction_endemic),yearsCirculating=30,foi_base=0.021,maxAge=maxAge)

#### Estimate by country
outputMat_TodayV<-outputMat_26V<-outputMat_85V<-array(NaN,c(length(uniqueCountries),4,length(vaccAgeMin)))
vaccinecoverage_today<-vaccinecoverage_26<-vaccinecoverage_85<-matrix(NaN,length(uniqueCountries),length(vaccAgeMin))
for(ii in 1:length(vaccAgeMin)){
	for (i in 1:length(uniqueCountries)){
		tmp<-popByageByCountry[i,,2]*propAtRisk[i,2]*coverage[ii]
		vaccinecoverage_26[i,ii]<-sum(tmp[vaccAgeMin[ii]:maxAge])
		tmp<-popByageByCountry[i,,3]*propAtRisk[i,3]*coverage[ii]
		vaccinecoverage_85[i,ii]<-sum(tmp[vaccAgeMin[ii]:maxAge])

		epiStatus<-originalEstimates$status[which(originalEstimates$Country==uniqueCountries[i])]

		covs<-rep(0,maxAge)
		if(epiStatus=="Endemic"){
			infecModTMP<-infecMod_endemic_V
			covs[vaccAgeMin[ii]:maxAge]<-coverage[ii]
		}else{
			infecModTMP<-infecMod_epidemic_V
			covs[vaccAgeMin[ii]:maxAge]<-coverage[ii]*(1-VaccinationDelay[ii])}

		tmp_vacc<-sweep(infecModTMP,1,popByageByCountry[i,,2]*propAtRisk[i,2]*covs*(1-protectionInfectionDirect),"*")
		tmp_unvacc<-sweep(infecModTMP,1,popByageByCountry[i,,2]*propAtRisk[i,2]*(1-covs),"*")
		tmp_vacc[vaccAgeMin[ii]:maxAge,2:4]<-tmp_vacc[vaccAgeMin[ii]:maxAge,2:4]*(1-protectionDiseaseNet[ii])
		tmp<-tmp_vacc+tmp_unvacc
		outputMat_26V[i,,ii]<-colSums(tmp)

		tmp_vacc<-sweep(infecModTMP,1,popByageByCountry[i,,3]*propAtRisk[i,2]*covs*(1-protectionInfectionDirect),"*")
		tmp_unvacc<-sweep(infecModTMP,1,popByageByCountry[i,,2]*propAtRisk[i,3]*(1-covs),"*")
		tmp_vacc[vaccAgeMin[ii]:maxAge,2:4]<-tmp_vacc[vaccAgeMin[ii]:maxAge,2:4]*(1-protectionDiseaseNet[ii])
		tmp<-tmp_vacc+tmp_unvacc
		outputMat_85V[i,,ii]<-colSums(tmp)
	}
}

outputMat_26OverallV<-apply(outputMat_26V,c(2,3),sum)
outputMat_85OverallV<-apply(outputMat_85V,c(2,3),sum)

V26<-sweep(sweep(outputMat_26OverallV,1,outputMat_26Overall,"-"),1,outputMat_26Overall,"/")
V85<-sweep(sweep(outputMat_85OverallV,1,outputMat_85Overall,"-"),1,outputMat_85Overall,"/")

TotalDose26<-colSums(vaccinecoverage_26)
TotalDose85<-colSums(vaccinecoverage_85)

V26PerDose<--100000*sweep(sweep(outputMat_26OverallV,1,outputMat_26Overall,"-"),2,TotalDose26,"/")
V85PerDose<--100000*sweep(sweep(outputMat_85OverallV,1,outputMat_85Overall,"-"),2,TotalDose85,"/")




#####################################
####### FIGURES #####################
#####################################

### PLOT FIGURE 1A
quartz(width=10,height=3.5)
par(mar=c(0,0,0,0))
plot(maxAedesMap,axes=F,zlim=c(0,12))
plot(world, add=T,lwd=0.1)


### PLOT FIGURE 1B-D
quartz(width=10,height=3.5)
par(mfrow=c(1,3))

par(mar=c(3,4,1,5))
plot(midmosq,probAegAlbwithCI[,1],pch=20,col="dark green",ylim=c(0,1),xlim=c(0,12),xlab="Months mosquito",ylab="Prob CHIKV",axes=F)
arrows(midmosq,probAegAlbwithCI[,2],midmosq,probAegAlbwithCI[,3],length=0,col="dark green")
lines(xs,pred)
axis(1)
axis(2,las=1)

plot(popSizeCountry[,1],popSizeCountry[,2],pch=20,col="dark blue",xlab="2025",ylab="2050",axes=F,log="xy")
points(popSizeCountry[,1],popSizeCountry[,3],pch=20,col="red2")
abline(a=0,b=1,lty=2)
axis(1)
axis(2,las=1)
legend("bottomright",legend=c("RCP2.6","RCP8.5"),pch=20,col=c("dark blue","red2"),bty="n",cex=1.5)

plot(monthsMosq[,1],monthsMosq[,2],pch=20,col="dark blue",xlab="2025",ylab="2050",axes=F,log="")
points(monthsMosq[,1],monthsMosq[,3],pch=20,col="red2")
abline(a=0,b=1,lty=2)
axis(1)
axis(2,las=1)
legend("bottomright",legend=c("RCP2.6","RCP8.5"),pch=20,col=c("dark blue","red2"),bty="n",cex=1.5)




### PLOT FIGURE 2
cols<-brewer.pal(9,"YlGn")
cols2<-colorRampPalette(cols,space="Lab")
Ncol<-50
cols<-(cols2(Ncol))
Ints<-seq(-0.00001,1.000001,length.out=Ncol+1)

cols2<-brewer.pal(9,"RdBu")
cols2<-colorRampPalette(cols2,space="Lab")
Ncol<-50
cols2<-rev(cols2(Ncol))
Ints2<-seq(-0.5,0.5,length.out=Ncol+1)

quartz(width=5,height=10)
par(mfrow=c(3,1))
par(mar=c(0.5,0.5,0.5,0.5))
plot(EstCountries,col=cols[findInterval(propAtRisk[,1],Ints)],lwd=0.05)
image.plot(1,1,zlim=c(0,1),col=cols,legend.only=TRUE)
range<-propAtRisk[,2]-propAtRisk[,1]
range[which(range<(-0.5))]<--0.5
plot(EstCountries,col=cols2[findInterval(range,Ints2)],lwd=0.05)
image.plot(1,1,zlim=c(-0.5,0.5),col=cols2,legend.only=TRUE)
range<-propAtRisk[,3]-propAtRisk[,1]
range[which(range<(-0.5))]<--0.5
plot(EstCountries,col=cols2[findInterval(range,Ints2)],lwd=0.05)
image.plot(1,1,zlim=c(-0.5,0.5),col=cols2,legend.only=TRUE)




### Figure 3A-B
cols<-brewer.pal(9,"PuOr")
cols2<-colorRampPalette(cols,space="Lab")
Ncol<-50
cols<-rev(cols2(Ncol))
Ints<-seq(-0.00001,1.000001,length.out=Ncol+1)

cols2<-brewer.pal(9,"RdBu")
cols2<-colorRampPalette(cols2,space="Lab")
Ncol<-50
cols2<-rev(cols2(Ncol))
Ints2<-exp(seq(-log(2),log(2),length.out=Ncol+1))
Ints3<-seq(-20001,20001,length.out=Ncol+1)
Ints4<-seq(-101,101,length.out=Ncol+1)

quartz(width=10,height=5)
par(mfrow=c(1,2))
par(mar=c(0.5,0.5,0.5,0.5))

range<-diffByCountry_26[,2]
range[which(range>20000)]<-20000
plot(EstCountries,col=cols[findInterval(range,Ints3)],lwd=0.05)

range<-diffByCountry_85[,2]
range[which(range>20000)]<-20000
plot(EstCountries,col=cols[findInterval(range,Ints3)],lwd=0.05)
image.plot(1,1,zlim=c(-20000,20000),col=cols,legend.only=TRUE)



### Figure 3C-F
ranks<-order(ratByCont_26[,2],decreasing=T)
nameslist[ranks]

quartz(width=10,height=4)
par(mfrow=c(1,4))
par(mar=c(3.5,5,0.5,0.5))

barplot(rbind(ratByCont_26[ranks,2],ratByCont_85[ranks,2]),col=c("dark blue","red2"),beside=T,axes=F,log="y",ylim=c(1,5))
legend("topleft",legend=c("RCP2.6", "RCP8.5"),fill=c("dark blue","red2"),bty="n",cex=1.5)
axis(2,las=1,at=c(1:5))
barplot(rbind(ratByCont_26[ranks,4],ratByCont_85[ranks,4]),col=c("dark blue","red2"),beside=T,axes=F,log="y",ylim=c(1,5))
axis(2,las=1,at=c(1:5))

barplot(rbind(outputMat_26Cont[ranks,2],outputMat_85Cont[ranks,2])/1e3,col=c("dark blue","red2"),beside=T,axes=F,log="",ylim=c(0,3000))
axis(2,las=1)
barplot(rbind(outputMat_26Cont[ranks,4],outputMat_85Cont[ranks,4])/1e3,col=c("dark blue","red2"),beside=T,axes=F,log="",ylim=c(0,10))
axis(2,las=1)




### Figure 4. Overall comparisons
plotdatCont_Inf<-100*rbind(CC26Cont[,1],CC85Cont[,1],PopSize2050Cont[,1],PopStructure2050Cont[,1],combined26Cont[,1],combined85Cont[,1])
plotdatCont_Case<-100*rbind(CC26Cont[,2],CC85Cont[,2],PopSize2050Cont[,2],PopStructure2050Cont[,2],combined26Cont[,2],combined85Cont[,2])
plotdatCont_Death<-100*rbind(CC26Cont[,4],CC85Cont[,4],PopSize2050Cont[,4],PopStructure2050Cont[,4],combined26Cont[,4],combined85Cont[,4])

plotdat<-100*cbind(CC26,CC85,PopSize2050,PopStructure2050,combined26,combined85)

cols<-brewer.pal(9,"YlGnBu")
cols2<-colorRampPalette(cols,space="Lab")
Ncol<-50
cols<-(cols2(Ncol))
Ints4<-seq(-200,200,length.out=Ncol+1)

zlim=c(-100,400)

quartz(width=8,height=3)
par(mfrow=c(1,3))
par(mar=c(3.5,1.2,3.5,1.2))
image.plot(1:6,1:6,plotdatCont_Inf[,ranks],col=cols,zlim=zlim,axes=F)
image.plot(1:6,1:6,plotdatCont_Case[,ranks],col=cols,zlim=zlim,axes=F)
image.plot(1:6,1:6,plotdatCont_Death[,ranks],col=cols,zlim=zlim,axes=F)


quartz(width=8,height=2)
par(mfrow=c(1,3))
par(mar=c(0,1.2,0.5,1.2))
barplot(plotdat[1,],axes=F,col=c("dark blue","red2","grey","grey","dark blue","red2"),density = c(20,20, NA,NA,NA,NA), angle = c(45,45,0,0,0,0),ylim=c(-0.1,1.5)*100)
axis(2,las=1)
barplot(plotdat[2,],axes=F,col=c("dark blue","red2","grey","grey","dark blue","red2"),density = c(20,20,NA,NA,NA,NA), angle = c(45,45,0,0,0,0),ylim=c(-0.1,1.5)*100)
barplot(plotdat[4,],axes=F,col=c("dark blue","red2","grey","grey","dark blue","red2"),density = c(20,20,NA,NA,NA,NA), angle = c(45,45,0,0,0,0),ylim=c(-0.1,1.5)*100)



###### PLOT FIGURE 5

quartz(width=8,height=6)
par(mfrow=c(2,3))
par(mar=c(0,3.5,0.5,0.5))

aa<-barplot(-V26[1,]*100,axes=F,col=c("grey"),ylim=c(0,1)*100)
points(aa[,1],-V85[1,]*100,pch=2,col="red2")
axis(2,las=1)
barplot(-V26[2,]*100,axes=F, col=c("grey"),ylim=c(0,1)*100)
points(aa[,1],-V85[2,]*100,pch=2,col="red2")
axis(2,las=1)
barplot(-V26[4,]*100,axes=F, col=c("grey"),ylim=c(0,1)*100)
points(aa[,1],-V85[4,]*100,pch=2,col="red2")
axis(2,las=1)

par(mar=c(3.5,3.5,0,0.5))
barplot(-V26PerDose[1,],axes=F,col=c("orange"),ylim=c(-600,0))
points(aa[,1],-V85PerDose[1,],pch=2,col="red2")
axis(2,las=1)
barplot(-V26PerDose[2,],axes=F, col=c("orange"),ylim=c(-100,0))
points(aa[,1],-V85PerDose[2,],pch=2,col="red2")
axis(2,las=1)
barplot(-V26PerDose[4,],axes=F, col=c("orange"),ylim=c(-0.4,0))
points(aa[,1],-V85PerDose[4,],pch=2,col="red2")
axis(2,las=1)



### Sensitivity analysis on FOI
pchs=c(1,2,3,15,16,17)

quartz(width=7,height=3.5)
par(mfrow=c(1,2))
par(mar=c(3,3,1,1))
plot(relFOI,sensitivityCases26Cont[1,],ylim=c(1,20000),type="n",log="xy",axes=F,xlim=c(0.2,2))
for (i in 1:6){
	# lines(relFOI,sensitivityCases85Cont[i,],col="red2",type="b",pch=pchs[i],lty=i+1)
	lines(relFOI,sensitivityCases26Cont[i,]/1e3,col="dark blue",type="b",pch=pchs[i],lty=1)
}
axis(1,at=c(0.2,0.5,1,2))
axis(2,las=1)
legend("topleft",legend=nameslist,pch=pchs,bty="n",col="dark blue")

plot(relFOI,sensitivityDeaths26Cont[1,],ylim=c(1,50000),type="n",log="xy",axes=F)
for (i in 1:6){
	lines(relFOI,sensitivityDeaths26Cont[i,],col="dark blue",type="b",pch=pchs[i],lty=1)
}
axis(1,at=c(0.2,0.5,1,2))
axis(2,las=1)









#######################
## Key estimates ######
#######################

## Pop at risk
colSums(effectivePopSize)
colSums(effectivePopSize)/colSums(effectivePopSize)[1]
colSums(effectivePopSize)/colSums(popSizeCountry)
colSums(popSizeCountry)-colSums(popSizeCountry)[3]
colSums(popSizeCountry)

outputMat_26Cont

## Number of countries with reduced mosquito suitability
sum(monthsMosq[,2]<monthsMosq[,1])
sum(monthsMosq[,3]<monthsMosq[,1])
1-sum(monthsMosq[,2]<monthsMosq[,1])/length(uniqueCountries)
1-sum(monthsMosq[,3]<monthsMosq[,1])/length(uniqueCountries)

weighted.mean(monthsMosq[,1],popSizeCountry[,1])
weighted.mean(monthsMosq[,2],popSizeCountry[,2])
weighted.mean(monthsMosq[,3],popSizeCountry[,3])
weighted.mean(monthsMosq[,2],popSizeCountry[,2])-weighted.mean(monthsMosq[,1],popSizeCountry[,1])
weighted.mean(monthsMosq[,3],popSizeCountry[,3])-weighted.mean(monthsMosq[,1],popSizeCountry[,1])



### Extract average changes by continent
tmp<-which(EstCountries@data$REGION=="Europe")
sum(effectivePopSize[tmp,2])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="North America")
sum(effectivePopSize[tmp,2])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="South America")
sum(effectivePopSize[tmp,2])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="Asia")
sum(effectivePopSize[tmp,2])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="Africa")
sum(effectivePopSize[tmp,2])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="Australia")
sum(effectivePopSize[tmp,2])/sum(effectivePopSize[tmp,1])


tmp<-which(EstCountries@data$REGION=="Europe")
sum(effectivePopSize[tmp,3])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="North America")
sum(effectivePopSize[tmp,3])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="South America")
sum(effectivePopSize[tmp,3])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="Asia")
sum(effectivePopSize[tmp,3])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="Africa")
sum(effectivePopSize[tmp,3])/sum(effectivePopSize[tmp,1])
tmp<-which(EstCountries@data$REGION=="Australia")
sum(effectivePopSize[tmp,2])/sum(effectivePopSize[tmp,1])


### Case burden
outputMat_26Overall/1e3
outputMat_85Overall/1e3
outputMat_26Overall/outputMat_TodayOverall
outputMat_85Overall/outputMat_TodayOverall

ratByCont_26[ranks,]-1

##Proportion by continent
sweep(outputMat_26Cont,2,colSums(outputMat_26Cont),"/")
sweep(outputMat_85Cont,2,colSums(outputMat_85Cont),"/")

### Force of infection sensitivity analysis
sensitivityCases26Cont
sensitivityDeaths26Cont
nameslist
relFOI


## Independent effects
CC26
CC26Cont
CC85
PopSize2050
PopSize2050Cont
PopStructure2050
PopStructure2050Cont

PopSize2050Cont

CC26Cont/combined26Cont
PopStructure2050Cont/combined26Cont
PopSize2050Cont/combined26Cont

tmp<-apply(cbind(CC26Cont[,1],PopSize2050Cont[,1],PopStructure2050Cont[,1]),1,which.max)
tmp[ranks]

tmp<-apply(cbind(CC26Cont[,2],PopSize2050Cont[,2],PopStructure2050Cont[,2]),1,which.max)
tmp[ranks]

tmp<-apply(cbind(CC26Cont[,4],PopSize2050Cont[,4],PopStructure2050Cont[,4]),1,which.max)
tmp[ranks]

tmp<-apply(cbind(CC26,PopSize2050,PopStructure2050),1,which.max)
tmp

## Vaccine impact
V26
TotalDose26

#Ancilleno Davis

#This r-code has been developed by Ancilleno Davis in partial completion of the 
#requirements for PhD in Ecology Evolution and Environmental Science at Miami University in Oxford, Ohio.

#The purpose of this Rcode is to import habitat classification data 
#generated in Google Earth Engine from LAndsat 8 OLI reflectance data,
#and determine correlations between eBird observations and habitat.
#for the island of Grand Bahama(GB), The Bahamas

#### Initial workspace parameters####
#set seed
set.seed(1981)

#Set my working directory and load necessary libraries
setwd("C:/Users/davisao2/Desktop/open source GIS/eBird r code")
#note this working directory should include:
#1: the geotiff of the raster data created using the Google Earth Engine habitat classification
#2: a spreadsheet of the eBird observations from the region of interest in a comma separated values text file.

#importing raster data
require(caret)#this has confusionMatrix for determining raster accuracy
require(cooccur)
require(data.table)#to allow rearrangment of columns by their headings
require(EcoSimR)
require("FactoMineR")#used for the Multiple Correspondence Analysis of the raster datasets
require("factoextra")#used for the Multiple Correspondence Analysis of the raster datasets
require(fmsb)
require(geosphere)
require(ggplot2)
require(ggpubr)
require(gstat)
require(lme4)
require(lmerTest)
require(plyr)
require(psych)
require(raster)
require(readr)
require(readxl)
require(reshape2)#contains dcast for reorganizing data into presence absence etc.
require(rgdal)#read in shapefiles
require(rgeos)
require(RStoolbox)
require(sp)

dev.off()#this resets the graphics window so that margins and other adjustments are returned to default

getmode <- function(v) { #creating this mode function allows me to get the most common value
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Assign the colors I would like to represent my habitat classes####

col7class=c(#include all 7 colors including water. This will be used to plot the raster image.
  "blue", #ocean 
  "dark green", #Pine Forest
  "brown", #Wetlands 
  "seashell2",#sand
  "gray", #Urban
  "Light green", #grass
  "yellow")#high water table plant communities
##Assign label names for the habitat classes####
legend=c("Water", 
         "Pine", 
         "Wetland", 
         "Sand", 
         "Urban", 
         "Grass", 
         "HWTC")
#import the Outline of Grand Bahama Island and reproject to WGS84####
GBOutline<-readOGR(".","Grand_Bahama_Outline")
print(proj4string(GBOutline)) #details of the projection for GBOutline
plot(GBOutline, 
     axes=T, #remove the lat long coordinates with F and add them with T
     border="black")

#Reproject GBOutline to WGS84
GBOutlineWGS84 <- 
  spTransform(GBOutline, #this command reprojects a spatial item  
              CRS("+proj=longlat +datum=WGS84"))#the new projection is WGS 84

plot(GBOutlineWGS84, 
     axes=TRUE, 
     border="black",
     las=1,
     #ylab="Latitude",
     # xlab="Longitude",
     main="Grand Bahama Contiguous Area WGS84")
###Import Rasters of habitat classification ####

#Random Forests classification map created using
#6 Landsat 8 OLI bands 2,3,4,5,6,7
#collected during 2017-01-01 through 2017-12-31

RF7Classes2017<-raster("RFclass6bands30m7classes2017.tif")
RF7Classes2017@crs #finds the coordinate reference system which is already WGS84
RF7Classes2017@extent #gives me the boundaries of the raster

plot(RF7Classes2017, 
     las=1, #rotates the y-axis labels horizontal
     ylab="Latitude",
     xlab="Longitude",
     main = "Grand Bahama Habitat Map")
#overlay the reprojected shapefile on the raster####
plot(GBOutlineWGS84, add=T)
#contiguous area habitat

#Clip area to the contiguous area of Grand Bahama Island####
GBcontiguoushab<-#this vector will hold our clipped raster
  mask(RF7Classes2017#This is the raster image to be clipped
       , GBOutlineWGS84)#This is the shape you are clipping to

plot(GBcontiguoushab); plot(GBOutlineWGS84, add=T) #Plot the trimmed raster and the outline
summary(GBcontiguoushab) #note all the clipped area is now NA
##Save the newly clipped raster to a file you can use in other software####
writeRaster(GBcontiguoushab,#our new rater image
            "GBcontiguoushab.tif"#the file name we want to save it to
            , format="GTiff") #the raster format we are using.
##display the clipped raster####
GBcontiguoushab.range<-cellStats(GBcontiguoushab, range)

plot(GBcontiguoushab,#plot this raster file
     col=col7class, #a list of colors make sure they are in the order of the classes you are plotting
     axis.args=list(
       at= #how do you want the axis labels spread out?
         seq(#sequential separations from the 
           GBcontiguoushab.range[1],#lowest end of the range
           GBcontiguoushab.range[2], #to the highest end of the range
           1),#with one unit of separation between each tick mark
       labels=legend, #at each tick mark, what do you want the label to say?
       cex.axis=1),
     axes=T ,
     main = "Grand Bahama Island, Bahamas habitat map", #title of the plot 
     ylab="Latitude", #left side 
     xlab="Longitude", #x axis label
     las=1 #rotate the y axis labels to horizontal
)
plot(GBOutlineWGS84, add=T)
legend('bottomright', legend , fill=col7class, border="black",
       col=col7class, bty='n', cex=1.5)
#summarize the number of pixels and the geographic area in each class####
GBhabitatdf<-as.data.frame(GBcontiguoushab)
habitatpixels<-count(GBhabitatdf, "RFclass6bands30m7classes2017")
habitatareakm2<- #a vector of the total area in each habitat type
  (habitatpixels$freq[1:7] #1:7 selects the columns with the classes 0-6 and ignores the NA values
  )*0.0009 #0.0009 is the number of km2 in a 30 x 30 m pixel (Landsat 8 imagery)
habitatareakm2
totalpixels<- #vector for the total pixels in the study area
  sum(habitatpixels$freq[1:7])
totalarea<-#vector of total area in study area
  sum(habitatareakm2)
percentarea<-habitatareakm2/totalarea #%of total area in each habitat type

habitat<-data.frame(legend,habitatpixels$freq[1:7],habitatareakm2,percentarea)

##Display barplots of the pixel area in each habitat type.####
BpGBHabitatPixels<-barplot(GBcontiguoushab, 
                           col=col7class,
                           axes=F,
                           main="Pixels per terrestrial habitat class on Grand Bahama 2017", 
                           xlab="Habitat type",
                           ylab="Number of Pixels",
                           horiz=FALSE,
                           ylim=c(0,400000),
                           names.arg=c("Water", 
                                       "Pine", 
                                       "Wetland", 
                                       "Sand", 
                                       "Urban", 
                                       "Grass", 
                                       "HWTC"),
                           las=1)
text(BpGBHabitatPixels, 
     habitatpixels$freq[1:7], 
     label=habitatpixels$freq[1:7], 
     pos = 3, 
     xpd = NA)
##Display barplots of the km2 area in each habitat type.####

BpGBHabitatkm2<-barplot(habitatareakm2, 
                        col=col7class,
                        axes=F,
                        main="Square Km per terrestrial habitat class on Grand Bahama 2017", 
                        xlab="Habitat type",
                        ylab="Area (sq.Km)",
                        horiz=FALSE,
                        ylim=c(0,400),
                        names.arg=c("Water", 
                                    "Pine", 
                                    "Wetland", 
                                    "Sand", 
                                    "Urban", 
                                    "Grass", 
                                    "HWTC"))
text(BpGBHabitatkm2, 
     habitatareakm2, 
     label=habitatareakm2, 
     pos = 3, 
     xpd = NA)



###Import confusion matrix generate in Google Earth Engine Code API and calculate Cohen's Kappa####
#Import the confusion matrix that was calculated by Google Earth Engine 
#your confusionmatrix should be in a csv with no header rows
#using validation data
confmatrix20177classes <-read_csv("7x7 class confusion matrix 2017.csv", 
                                  col_names = TRUE)
##View(confmatrix20177classes)
confmatrix20177classes<-as.data.frame(confmatrix20177classes)
rownames(confmatrix20177classes)<-
  colnames(confmatrix20177classes)<-
  c("Water",
    "Pine",
    "Wetland",
    "Sand",
    "Urban",
    "Grass",
    "HWTC")

#calculate the cohen's kappa coefficient to determine the agreement between 
#the classification and validation data
Kappa.test(confmatrix20177classes,y=NULL, conf.level=0.95)




# ###Import the data table including the eBird locations and visual classification of the site by Ancilleno Davis####
# VisualClassification <-read_csv("Visual Classification of eBird sites on Grand Bahama.csv", 
#                                   col_names = TRUE)
# VisualClassificationcoords<-data.frame(VisualClassification$LONGITUDE,VisualClassification$LATITUDE)
# ##Extract the habitat data from the raster and add it to the visual classification data frame
# VisualClassification$Habitat<-
#   extract(RF7Classes2017,
#           VisualClassificationcoords)
# VisualVsGEEHabitat<-data.frame(as.factor(VisualClassification$Habitat),as.factor(VisualClassification$`Visual verification`))
# colnames(VisualVsGEEHabitat)<-c("raster","visual")
# 
# summary(VisualVsGEEHabitat)
# 
# merge.test <- merge(GBeBird1988.2016, VisualClassification, by=c('LATITUDE', 'LONGITUDE'),all=TRUE)[,-1]
# #Calculate Cohen's Kappa on the generated raster and visual verification of habitat at eBird Localities
# cohen.kappa(VisualVsGEEHabitat)




###import all locations with Bahamian bird species####
####create data frame with all bird observation data for Grand Bahama####

BSeBirdNOV2017complete <- read_csv("ebd_BS_relNov-2017.csv", 
col_types = cols(
  `ALL SPECIES REPORTED` = col_factor(levels = c("0","1")), 
  APPROVED = col_factor(levels = c("0","1")), 
  `EFFORT AREA HA` = col_double(), 
  `HAS MEDIA` = col_factor(levels = c("0","1")), 
  `OBSERVATION COUNT` = col_integer(), 
  `OBSERVATION DATE` = col_date(format = "%m/%d/%Y"), 
  REVIEWED = col_factor(levels = c("0","1")), 
  `TIME OBSERVATIONS STARTED` = col_time(format = "%H:%M:%S")))

#subset to data around Grand Bahama Island
GBeBirdcomplete<-BSeBirdNOV2017complete[ #subset
  BSeBirdNOV2017complete$STATE_PROVINCE %in% 
    c('Freeport and West Grand Bahama', 'East Grand Bahama'),]
#Replace all NAs in the Observation Count column with 1
#The species was identified but there was no count of individuals at least one individual is assumed

GBeBirdcomplete$`OBSERVATION COUNT`[is.na(GBeBirdcomplete$`OBSERVATION COUNT`)]<-1

#add columns for the month, year and day of the week of the observations
GBeBirdcomplete$YEAR<-as.numeric(format(GBeBirdcomplete$`OBSERVATION DATE`, format = "%Y"))
GBeBirdcomplete$MONTH<-as.numeric(format(GBeBirdcomplete$`OBSERVATION DATE`, format = "%m"))
GBeBirdcomplete$day <- as.factor(weekdays(as.Date(GBeBirdcomplete$`OBSERVATION DATE`)))
GBeBirdcomplete$day  <- factor(GBeBirdcomplete$day, levels= c("Sunday", "Monday", 
                                       "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))

GBeBirdcomplete[order(GBeBirdcomplete$day), ]
View(GBeBirdcomplete)



# select date range ------------------------------------------------------------
GBeBirdDEC2016<-GBeBirdcomplete[ #subset
  (GBeBirdcomplete$`OBSERVATION DATE`> "1987-12-31" & 
     GBeBirdcomplete$`OBSERVATION DATE` < "2017-01-01"),] #Before this date2017

GBeBird1988.2016<-GBeBirdcomplete[ #subset
  (GBeBirdcomplete$`OBSERVATION DATE`> "1987-12-31" & #After this date 1987 BirdsCaribbean was founded in 1988
     GBeBirdcomplete$`OBSERVATION DATE` < "2017-01-01"),] #Before this date2017

#Create yearly subset 2016
GBeBird2016only<-GBeBirdcomplete[ #subset
  (GBeBirdcomplete$`OBSERVATION DATE`> "2015-12-31" & #After this date
     GBeBirdcomplete$`OBSERVATION DATE` < "2017-01-01"),] #Before this date2017



# select only surveys in which all species were recorded ---------------------------------------------
GBeBirdDEC2016<-GBeBirdDEC2016[ #subset
  GBeBirdDEC2016$`ALL SPECIES REPORTED`== "1",] #records with all species reported


###The GBeBirdDEC2016 now includes all ebird data up to December 31 2016 and since 1988#### 

##in which observers attempted to record all the species they were able to identify.####
View(GBeBirdDEC2016)
##remove unneeded columns####
GBeBirdDEC2016$`LAST EDITED DATE`<-
  GBeBirdDEC2016$`SUBSPECIES COMMON NAME`<-
  GBeBirdDEC2016$`SUBSPECIES SCIENTIFIC NAME`<-
  GBeBirdDEC2016$`BREEDING BIRD ATLAS CODE`<-
  GBeBirdDEC2016$`BREEDING BIRD ATLAS CATEGORY`<-
  GBeBirdDEC2016$`AGE/SEX`<-
  GBeBirdDEC2016$COUNTRY<-
  GBeBirdDEC2016$COUNTRY_CODE<-
  GBeBirdDEC2016$COUNTY<-
  GBeBirdDEC2016$SUBNATIONAL1_CODE<-
  GBeBirdDEC2016$SUBNATIONAL2_CODE<-
  GBeBirdDEC2016$`IBA CODE`<-
  GBeBirdDEC2016$`BCR CODE`<-
  GBeBirdDEC2016$`USFWS CODE`<-
  GBeBirdDEC2016$`ATLAS BLOCK`<-
  NULL

View(GBeBirdDEC2016)
####Compile BAH and DAS species lists####
#Compile list of Bahamian endemic and resident species of importance for comparison###
BAHspecieslist<-c("Bahama Mockingbird", 
                  "Bahama Swallow", 
                  "Bahama Warbler", 
                  "Bahama Woodstar", 
                  "Bahama Yellowthroat",
                  "Brown-headed Nuthatch",
                  "Common Ground-Dove", 
                  "Cuban Pewee",
                  "Hairy Woodpecker",
                  "La Sagra's Flycatcher",
                  "Loggerhead Kingbird",
                  "Northern Mockingbird", 
                  "Olive-capped Warbler",
                  "Pine Warbler",
                  "Red-legged Thrush",
                  "Western Spindalis")

#Compile list of Detroit Audubon society Species of concern
DASspecieslist<-c("American Kestrel",
                  "Barn Swallow",
                  "Belted Kingfisher",
                  "Blue-winged Teal",
                  "Canvasback",##Note there are no canvasback in the data set
                  "Common Tern",
                  "Gray Catbird",
                  "Great Blue Heron",
                  "Indigo Bunting",
                  "Killdeer",
                  "Kirtland's Warbler",##no Kirtland's Warblers were found on GB
                  "Lesser Scaup",
                  "Osprey",
                  "Palm Warbler",
                  "Piping Plover",
                  "Ring-necked Duck",
                  "Ruddy Turnstone",
                  "Sanderling",
                  "Sharp-shinned Hawk",
                  "Turkey Vulture"
)

#NB these lists are not comprehensive

###insert columns to classify all bird species as a Bahamian or DAS species ####
#based on being in BAHspecieslist or DASspecieslist 

GBeBirdDEC2016$BAHspecies<-GBeBirdDEC2016$`COMMON NAME`%in% BAHspecieslist
GBeBirdDEC2016$DASspecies<-GBeBirdDEC2016$`COMMON NAME`%in% DASspecieslist

####Create subset of points for bird species in the Bahamian species list only####
GBeBirdDEC2016BAH.sub<-
  subset(GBeBirdDEC2016, 
         GBeBirdDEC2016$BAHspecies=='TRUE')

#View(GBeBirdDEC2016BAH.sub)

#####Create subset of points for bird species in the Detroit Audubon Society species list only####
GBeBirdDEC2016DAS.sub<-
  subset(GBeBirdDEC2016, 
         GBeBirdDEC2016$DASspecies=='TRUE')

View(GBeBirdDEC2016DAS.sub)


#####Get habitat value for each location####
####create location dataframes for all observations and BAH or DAS species observations.####
GBeBirdDEC2016Localities<-
  data.frame(GBeBirdDEC2016$LOCALITY, 
             GBeBirdDEC2016$`LOCALITY ID`,
             GBeBirdDEC2016$LONGITUDE,
             GBeBirdDEC2016$LATITUDE)

GBeBirdDEC2016Localities.BAH<-
  data.frame(GBeBirdDEC2016BAH.sub$LOCALITY, 
             GBeBirdDEC2016BAH.sub$`LOCALITY ID`,
             GBeBirdDEC2016BAH.sub$LONGITUDE,
             GBeBirdDEC2016BAH.sub$LATITUDE)

GBeBirdDEC2016Localities.DAS<-
  data.frame(GBeBirdDEC2016DAS.sub$LOCALITY, 
             GBeBirdDEC2016DAS.sub$`LOCALITY ID`,
             GBeBirdDEC2016DAS.sub$LONGITUDE,
             GBeBirdDEC2016DAS.sub$LATITUDE)
#Extract location coordinates for All species and BAH or DAS species####  
GBeBirdDEC2016coords<-
  data.frame(GBeBirdDEC2016$LONGITUDE,
             GBeBirdDEC2016$LATITUDE)
GBeBirdDEC2016coords.BAH<-
  data.frame(GBeBirdDEC2016BAH.sub$LONGITUDE,
             GBeBirdDEC2016BAH.sub$LATITUDE)
GBeBirdDEC2016coords.DAS<-
  data.frame(GBeBirdDEC2016DAS.sub$LONGITUDE,
             GBeBirdDEC2016DAS.sub$LATITUDE)


#extract habitat from the raster at the points in the location dataframes####
#add habitat to a new column in the eBird data
GBeBirdDEC2016$Habitat<-
  extract(RF7Classes2017,
          GBeBirdDEC2016coords)
GBeBirdDEC2016BAH.sub$Habitat<-
  extract(RF7Classes2017,
          GBeBirdDEC2016coords.BAH)
GBeBirdDEC2016DAS.sub$Habitat<-
  extract(RF7Classes2017,
          GBeBirdDEC2016coords.DAS)

##Export data with habitat type to csv files####
write.csv(GBeBirdDEC2016, "GBeBirdDEC2016ALLhabitat.csv")#all the ebird data plus habitat data for All species
write.csv(GBeBirdDEC2016DAS.sub, "GBeBird1988toDEC2016DAShabitat.csv")#all ebird data plus habitat data for Detroit Audubon Species
write.csv(GBeBirdDEC2016BAH.sub, "GBeBird1988toDEC2016BAHhabitat.csv")#all ebird data plus habitat data for Bahamas Species



## Which Localities have no DAS species or Bah species reported?####
# #Have no BAH species reported
# unique(GBeBirdDEC2016$`LOCALITY ID`)
# NOBAHLocalities<-setdiff(unique(GBeBirdDEC2016$`LOCALITY ID`), unique(GBeBirdDEC2016BAH.sub$`LOCALITY ID`))
# #Have no DAS species reported
# setdiff(GBeBirdDEC2016$`LOCALITY ID`, GBeBirdDEC2016DAS.sub$`LOCALITY ID`)
# #Have BAH species but no DAS species
# setdiff(GBeBirdDEC2016BAH.sub$`LOCALITY ID`, GBeBirdDEC2016DAS.sub$`LOCALITY ID`)
# #Have DAS species but no BAH species
# setdiff( GBeBirdDEC2016DAS.sub$`LOCALITY ID`,GBeBirdDEC2016BAH.sub$`LOCALITY ID`)



####get the list of Observers that have seen Bahamian or DAS species or both, neither etc####
allObservers<-unique(GBeBirdDEC2016$`OBSERVER ID`)
BAHSpsObserversList<-unique(GBeBirdDEC2016BAH.sub$`OBSERVER ID`)
DASSpsObserversList<-unique(GBeBirdDEC2016DAS.sub$`OBSERVER ID`)
FocalSpsObserversList<-intersect(BAHSpsObserversList,DASSpsObserversList) #Have seen at least one DAS and one BAH sps

NoBAHObservers<-setdiff(allObservers, BAHSpsObserversList)
NoDASObservers<-setdiff(allObservers, DASSpsObserversList)
NoFocalObservers<-intersect((setdiff(allObservers, BAHSpsObserversList)), (setdiff(allObservers, DASSpsObserversList)) )
DASonlyObservers<-intersect(DASSpsObserversList, NoBAHObservers)
BAHonlyObservers<-intersect(BAHSpsObserversList, NoDASObservers)
####Calculate the number of Locations, Observers, Surveys and Species per habitat####
###Observers per habitat all ####
ObsAndHabitat<- data.frame(GBeBirdDEC2016$`OBSERVER ID`,GBeBirdDEC2016$Habitat)
ObsAndHabitat<-unique(ObsAndHabitat)
colnames(ObsAndHabitat)<-c( "Observer ID", "Habitat")
ObsAndHabitat$Habitat<-as.factor(ObsAndHabitat$Habitat)
ObsAndHabitat<-ObsAndHabitat[complete.cases(ObsAndHabitat),]
AllObservers<-unique(ObsAndHabitat$`Observer ID`) #= 346
ObsAndHabitat$Habitat<-as.factor(ObsAndHabitat$Habitat)
O<-count(ObsAndHabitat, Habitat)
habitat$Observers<-O$n


###Locations per habitat####
LocAndHabitat<- data.frame(
  GBeBirdDEC2016$LOCALITY,
  GBeBirdDEC2016$`LOCALITY ID`,
  GBeBirdDEC2016$Habitat)
LocAndHabitat<-unique(LocAndHabitat)
colnames(LocAndHabitat)<-c("Locality", "Locality ID", "Habitat")
LocAndHabitat$Habitat<-as.factor(LocAndHabitat$Habitat)
LocAndHabitat<-LocAndHabitat[complete.cases(LocAndHabitat),] #remove the locations with habitat NA (these were outside the bounds of the raster)
AllLocalities<-unique(LocAndHabitat$`Locality ID`) #=603
L<-count(LocAndHabitat, Habitat) #frequencies correspond to the number of localities per habitat
habitat$locations<-L$n
summary(LocAndHabitat$Habitat)
# 0   1   2   3   4   5   6 
# 31  78  49 120  92 111  18 

###Species per habitat####
SpsAndHabitat<- data.frame(
  GBeBirdDEC2016$`SCIENTIFIC NAME`,
  GBeBirdDEC2016$`COMMON NAME`,
  GBeBirdDEC2016$Habitat)
SpsAndHabitat<-unique(SpsAndHabitat)
colnames(SpsAndHabitat)<-c("SCI_NAME", "COM_NAME", "Habitat")

SpsAndHabitat<-SpsAndHabitat[complete.cases(SpsAndHabitat),] #remove the locations with habitat NA (these were outside the bounds of the raster)
AllSps<-unique(SpsAndHabitat$COM_NAME) #=208
S<-count(SpsAndHabitat, 'Habitat') #frequencies correspond to the number of species per habitat
habitat$Species<-S$n
SpsAndHabitat$Habitat<-as.factor(SpsAndHabitat$Habitat)
summary(SpsAndHabitat$Habitat)
# 0   1   2   3   4   5   6 
# 108 155 176 217 163 247 113

###Surveys per habitat####
SurAndHabitat<- data.frame(
  GBeBirdDEC2016$`SAMPLING EVENT IDENTIFIER`,
    GBeBirdDEC2016$Habitat)
SurAndHabitat<-unique(SurAndHabitat)
colnames(SurAndHabitat)<-c("Survey", "Habitat")
SurAndHabitat$Habitat<-as.factor(SurAndHabitat$Habitat)
SurAndHabitat<-SurAndHabitat[complete.cases(SurAndHabitat),] #remove the locations with habitat NA (these were outside the bounds of the raster)
AllSur<-unique(SurAndHabitat) #=208
Sur<-count(SurAndHabitat, "Habitat") #frequencies correspond to the number of localities per habitat
habitat$Surveys<-Sur$freq
summary(SurAndHabitat$Habitat)
# 0    1    2    3    4    5    6 
# 41  325  358 1303  472 1671   72 




#add the predicted number of DAS and BAH species in each habitat to the habitat table####
#There is no predicted number for Water, Grass or HWTC communities as these habitats are not terrestrial (water) or atypical

habitat$DASpr<-c(0,8,11,9,4,0,0) #Predicted # of DAS species per habitat

habitat$BAHpr<-c(0,14,4,2,7,0,0)#Predicted # of DAS species per habitat
habitat$BAHobs<-SBAH$freq
habitat$DASobs<-SDAS$freq
habitat

###Habitats per observer all####

HabPerObs<-count(ObsAndHabitat, "`Observer ID`")
colnames(HabPerObs)<-c("Observer ID", "Habitats visited")
summary(HabPerObs$`Habitats visited`)

###Observers per habitat DAS ####
ObsAndHabitatDAS<- data.frame(GBeBirdDEC2016DAS.sub$`OBSERVER ID`,GBeBirdDEC2016DAS.sub$Habitat)
ObsAndHabitatDAS<-unique(ObsAndHabitatDAS)
colnames(ObsAndHabitatDAS)<-c( "Observer ID", "Habitat")
ObsAndHabitatDAS$Habitat<-as.factor(ObsAndHabitatDAS$Habitat)
ObsAndHabitatDAS<-ObsAndHabitatDAS[complete.cases(ObsAndHabitatDAS),]
AllObserversDAS<-unique(ObsAndHabitatDAS$`Observer ID`) #= 346
ObsAndHabitatDAS$Habitat<-as.factor(ObsAndHabitatDAS$Habitat)
ODAS<-count(ObsAndHabitatDAS,"Habitat")
habitat$ObserversDAS<-ODAS$n
summary(ObsAndHabitatDAS$Habitat)


###Habitats per observer DAS####

HabPerObsDAS<-count(ObsAndHabitatDAS, "`Observer ID`")
colnames(HabPerObsDAS)<-c("Observer ID", "Habitats visited")
summary(HabPerObsDAS$`Habitats visited`)

###Species per habitat DAS####
SpsAndHabitatDAS<- data.frame(
  GBeBirdDEC2016DAS.sub$`SCIENTIFIC NAME`,
  GBeBirdDEC2016DAS.sub$`COMMON NAME`,
  GBeBirdDEC2016DAS.sub$Habitat)
SpsAndHabitatDAS<-unique(SpsAndHabitatDAS)
colnames(SpsAndHabitatDAS)<-c("SCI_NAME", "COM_NAME", "Habitat")
SpsAndHabitatDAS$Habitat<-as.factor(SpsAndHabitatDAS$Habitat)
SpsAndHabitatDAS<-SpsAndHabitatDAS[complete.cases(SpsAndHabitatDAS),] #remove the locations with habitat NA (these were outside the bounds of the raster)
DASSps<-unique(SpsAndHabitatDAS$COM_NAME) #=208
SDAS<-count(SpsAndHabitatDAS, "Habitat") #frequencies correspond to the number of localities per habitat
habitat$DASSpecies<-SDAS$n
summary(SpsAndHabitatDAS$Habitat)
###Surveys per habitatDAS####
SurAndHabitatDAS<- data.frame(
  GBeBirdDEC2016DAS.sub$`SAMPLING EVENT IDENTIFIER`,
  GBeBirdDEC2016DAS.sub$Habitat)
SurAndHabitatDAS<-unique(SurAndHabitatDAS)
colnames(SurAndHabitatDAS)<-c("Survey", "Habitat")
SurAndHabitatDAS$Habitat<-as.factor(SurAndHabitatDAS$Habitat)
SurAndHabitatDAS<-SurAndHabitatDAS[complete.cases(SurAndHabitatDAS),] #remove the locations with habitat NA (these were outside the bounds of the raster)
DASSur<-unique(SurAndHabitatDAS) #=208
SurDAS<-count(SurAndHabitatDAS, "Habitat") #frequencies correspond to the number of localities per habitat
habitat$DASSurveys<-SurDAS$n
summary(SurAndHabitatDAS$Habitat)

###Observers per habitat BAH ####
ObsAndHabitatBAH<- data.frame(GBeBirdDEC2016BAH.sub$`OBSERVER ID`,GBeBirdDEC2016BAH.sub$Habitat)
ObsAndHabitatBAH<-unique(ObsAndHabitatBAH)
colnames(ObsAndHabitatBAH)<-c( "Observer ID", "Habitat")
ObsAndHabitatBAH$Habitat<-as.factor(ObsAndHabitatBAH$Habitat)
ObsAndHabitatBAH<-ObsAndHabitatBAH[complete.cases(ObsAndHabitatBAH),]
AllObserversBAH<-unique(ObsAndHabitatBAH$`Observer ID`) #= 346
ObsAndHabitatBAH$Habitat<-as.factor(ObsAndHabitatBAH$Habitat)
OBAH<-count(ObsAndHabitatBAH, "Habitat")
habitat$ObserversBAH<-OBAH$n
summary(ObsAndHabitatBAH$Habitat)


###Habitats per observer BAH####

HabPerObsBAH<-count(ObsAndHabitatBAH, "`Observer ID`")
colnames(HabPerObsBAH)<-c("Observer ID", "Habitats visited")
summary(HabPerObsBAH)
length(HabPerObsBAH$`Observer ID`)

###Species per habitat BAH####
SpsAndHabitatBAH<- data.frame(
  GBeBirdDEC2016BAH.sub$`SCIENTIFIC NAME`,
  GBeBirdDEC2016BAH.sub$`COMMON NAME`,
  GBeBirdDEC2016BAH.sub$Habitat)
SpsAndHabitatBAH<-unique(SpsAndHabitatBAH)
colnames(SpsAndHabitatBAH)<-c("SCI_NAME", "COM_NAME", "Habitat")
SpsAndHabitatBAH$Habitat<-as.factor(SpsAndHabitatBAH$Habitat)
SpsAndHabitatBAH<-SpsAndHabitatBAH[complete.cases(SpsAndHabitatBAH),] #remove the locations with habitat NA (these were outside the bounds of the raster)
BAHSps<-unique(SpsAndHabitatBAH$COM_NAME)# 16
SBAH<-count(SpsAndHabitatBAH, "Habitat") #frequencies correspond to the number of localities per habitat
habitat$BAHSpecies<-SBAH$freq
summary(SpsAndHabitatBAH$Habitat)
# 0  1  2  3  4  5  6 
# 16 16 16 15 14 16 16 


###Surveys per habitat BAH####
SurAndHabitatBAH<- data.frame(
  GBeBirdDEC2016BAH.sub$`SAMPLING EVENT IDENTIFIER`,
  GBeBirdDEC2016BAH.sub$Habitat)
SurAndHabitatBAH<-unique(SurAndHabitatBAH)
colnames(SurAndHabitatBAH)<-c("Survey", "Habitat")
SurAndHabitatBAH$Habitat<-as.factor(SurAndHabitatBAH$Habitat)
SurAndHabitatBAH<-SurAndHabitatBAH[complete.cases(SurAndHabitatBAH),] #remove the locations with habitat NA (these were outside the bounds of the raster)
BAHSur<-unique(SurAndHabitatBAH) #=208
SurBAH<-count(SurAndHabitatBAH, "Habitat") #frequencies correspond to the number of localities per habitat
habitat$BAHSurveys<-SurBAH$freq
summary(SurAndHabitatBAH$Habitat)
# 0    1    2    3    4    5    6 
# 16  272  275  931  361 1391   64 

#                                                              ####
###GEPHI NETWORK DATA####
# #####
##create source target data for gephi networks of Observers and Habitats####
ObserverHabitatAll<-unique(
  data.frame(GBeBirdDEC2016$`OBSERVER ID`, 
             GBeBirdDEC2016$Habitat))
colnames(ObserverHabitatAll)<-c("Source", "Target")
ObserverHabitatAll<-na.omit(ObserverHabitatAll) #removes locations outside the study area

ObserverHabitatDAS<-unique(
  data.frame(GBeBirdDEC2016DAS.sub$`OBSERVER ID`, 
             GBeBirdDEC2016DAS.sub$Habitat))
colnames(ObserverHabitatDAS)<-c("Source", "Target")
ObserverHabitatDAS<-na.omit(ObserverHabitatDAS) #removes locations outside the study area

ObserverHabitatBAH<-unique(
  data.frame(GBeBirdDEC2016BAH.sub$`OBSERVER ID`, 
             GBeBirdDEC2016BAH.sub$Habitat))
colnames(ObserverHabitatBAH)<-c("Source", "Target")
ObserverHabitatBAH<-na.omit(ObserverHabitatBAH) #removes locations outside the study area


##create source target data for gephi networks of Species and Habitats####
SpeciesHabitatAll<-unique(
  data.frame(GBeBirdDEC2016$`COMMON NAME`, 
             GBeBirdDEC2016$Habitat))
colnames(SpeciesHabitatAll)<-c("Source", "Target")
SpeciesHabitatAll<-na.omit(SpeciesHabitatAll) #removes locations outside the study area

SpeciesHabitatDAS<-unique(
  data.frame(GBeBirdDEC2016DAS.sub$`COMMON NAME`, 
             GBeBirdDEC2016DAS.sub$Habitat))
colnames(SpeciesHabitatDAS)<-c("Source", "Target")
SpeciesHabitatDAS<-na.omit(SpeciesHabitatDAS) #removes locations outside the study area

SpeciesHabitatBAH<-unique(
  data.frame(GBeBirdDEC2016BAH.sub$`COMMON NAME`, 
             GBeBirdDEC2016BAH.sub$Habitat))
colnames(SpeciesHabitatBAH)<-c("Source", "Target")
SpeciesHabitatBAH<-na.omit(SpeciesHabitatBAH) #removes locations outside the study area


##Export networks to Gephi csv files####
##Observer Presence
write.csv(ObserverHabitatAll, "ObserverHabitatAll.csv")
write.csv(ObserverHabitatDAS, "ObserverHabitatDAS.csv")
write.csv(ObserverHabitatBAH, "ObserverHabitatBAH.csv")
write.csv(SpeciesHabitatAll, "SpeciesHabitatAll.csv")
write.csv(SpeciesHabitatDAS, "SpeciesHabitatDAS.csv")
write.csv(SpeciesHabitatBAH, "SpeciesHabitatBAH.csv")
SpeciesHabitatDAS

##Surveys per Habitat####
SurAndHabitat<- data.frame(GBeBirdDEC2016$`SAMPLING EVENT IDENTIFIER`,GBeBirdDEC2016$Habitat)
SurAndHabitat<-unique(SurAndHabitat)
colnames(SurAndHabitat)<-c( "Sampling Event Identifier", "Habitat")
SurAndHabitat$Habitat<-as.factor(SurAndHabitat$Habitat)
SurAndHabitat<-SurAndHabitat[complete.cases(SurAndHabitat),]
AllSurveys<-unique(SurAndHabitat$`Sampling Event Identifier`) #= 4915
SurAndHabitat$Habitat<-as.factor(SurAndHabitat$Habitat)
Sur<-count(SurAndHabitat,'Habitat')
habitat$Surveys<-Sur$freq
summary(SurAndHabitat$Habitat)
# 0    1    2    3    4    5    6 
# 41  330  359 1304  474 1680   72  

####Species per Habitat (all)####
SpsAndHabitat<-unique(
  data.frame(
    GBeBirdDEC2016$`SCIENTIFIC NAME`,
    GBeBirdDEC2016$`COMMON NAME`, 
    GBeBirdDEC2016$Habitat))
colnames(SpsAndHabitat)<-c( "Scientific Name","Common Name", "Habitat")
SpsAndHabitat$`Scientific Name`<-as.factor(SpsAndHabitat$`Scientific Name`)
SpsAndHabitat$`Common Name`<-as.factor(SpsAndHabitat$`Common Name`)
SpsAndHabitat$Habitat<-as.factor(SpsAndHabitat$Habitat)
SpsAndHabitat<-SpsAndHabitat[complete.cases(SpsAndHabitat),]
SpsAndHabitat
AllSpecies<-unique(data.frame(SpsAndHabitat$`Scientific Name`,SpsAndHabitat$`Common Name`)) #= 302
colnames(AllSpecies)<-c("Scientific Name", "Common Name")
Sps<-count(SpsAndHabitat,'Habitat')
habitat$Species<-Sps$freq
summary(SpsAndHabitat$Habitat)
# 0   1   2   3   4   5   6 
# 108 156 176 217 163 248 113  

####Species per Habitat BAH####
SpsAndHabitatBAH<-unique(
  data.frame(
    GBeBirdDEC2016BAH.sub$`SCIENTIFIC NAME`,
    GBeBirdDEC2016BAH.sub$`COMMON NAME`, 
    GBeBirdDEC2016BAH.sub$Habitat))
colnames(SpsAndHabitatBAH)<-c( "Scientific Name","Common Name", "Habitat")
SpsAndHabitatBAH$`Scientific Name`<-as.factor(SpsAndHabitatBAH$`Scientific Name`)
SpsAndHabitatBAH$`Common Name`<-as.factor(SpsAndHabitatBAH$`Common Name`)
SpsAndHabitatBAH$Habitat<-as.factor(SpsAndHabitatBAH$Habitat)
SpsAndHabitatBAH<-SpsAndHabitatBAH[complete.cases(SpsAndHabitatBAH),]
SpsAndHabitatBAH
BAHSpecies<-unique(data.frame(SpsAndHabitatBAH$`Scientific Name`,SpsAndHabitatBAH$`Common Name`)) #= 302
colnames(BAHSpecies)<-c("Scientific Name", "Common Name")
BAHSps<-count(SpsAndHabitatBAH,'Habitat')

HabPerBAHSps<-count(SpsAndHabitatBAH$`Scientific Name`)
colnames(HabPerBAHSps)<-c("Scientific Name", "Number of habitats")
SpsAndHabitatBAH
HabPerBAHSps

####Species per Habitat DAS####
SpsAndHabitatDAS<-unique(
  data.frame(
    GBeBirdDEC2016DAS.sub$`SCIENTIFIC NAME`,
    GBeBirdDEC2016DAS.sub$`COMMON NAME`, 
    GBeBirdDEC2016DAS.sub$Habitat))
colnames(SpsAndHabitatDAS)<-c( "Scientific Name","Common Name", "Habitat")
SpsAndHabitatDAS$`Scientific Name`<-as.factor(SpsAndHabitatDAS$`Scientific Name`)
SpsAndHabitatDAS$`Common Name`<-as.factor(SpsAndHabitatDAS$`Common Name`)
SpsAndHabitatDAS$Habitat<-as.factor(SpsAndHabitatDAS$Habitat)
SpsAndHabitatDAS<-SpsAndHabitatDAS[complete.cases(SpsAndHabitatDAS),]
SpsAndHabitatDAS
DASSpecies<-unique(data.frame(SpsAndHabitatDAS$`Scientific Name`,SpsAndHabitatDAS$`Common Name`)) #= 302
colnames(DASSpecies)<-c("Scientific Name", "Common Name")
DASSps<-count(SpsAndHabitatDAS,'Habitat')

HabPerDASSps<-count(SpsAndHabitatDAS$`Scientific Name`)
colnames(HabPerDASSps)<-c("Scientific Name", "Number of habitats")
SpsAndHabitatDAS
HabPerDASSps
length(DASspecieslist)
#now we have data to say how many habitats each observer and species has been recorded in
#we also have know how each habitat is represented in the data


##Check to ensure that you have at least one record for each of your species.
#### Create a table of the occurrence of each species and observers in each habitat or location####

##Species habitat matrices####
SpsHabMatrix<-dcast(GBeBirdDEC2016,
      `COMMON NAME` + #Only use the common name because EcoSimR and Cooccur need dataframes with just the single column and row headings.
       `SCIENTIFIC NAME`
      ~Habitat,
      value.var = "SAMPLING EVENT IDENTIFIER")
SpsHabMatrix<-na.omit(SpsHabMatrix)# remove lines with NA

rownames(SpsHabMatrix)<-SpsHabMatrix$`COMMON NAME`

#now that the rownames are applied you can remove the common names column and the NA's column
#SpsHabMatrix$`COMMON NAME`<-
  SpsHabMatrix$`SCIENTIFIC NAME` <-
  SpsHabMatrix$`NA`<-NULL 
  SpsHabMatrix$Surveys<-rowSums(SpsHabMatrix[,
                                             c('0', '1', '2','3','4', '5', '6')], na.rm=TRUE)
  
SpsHabMatrix
SpeciesHabPresence<-SpsHabMatrix
SpeciesHabPresence[SpeciesHabPresence > 0] <- 1
SpeciesHabPresence
#Add columns for percentages
###Initialize percentage columns for percent of surveys with the species in each habitat.####

SpsHabMatrix$Water<-
  SpsHabMatrix$Pine<-
  SpsHabMatrix$Wetland<-
  SpsHabMatrix$Sand<-
  SpsHabMatrix$Urban<-
  SpsHabMatrix$Grass<-
  SpsHabMatrix$HWTC<-0
options(digits = 2)
SpsHabMatrix$Water<-(SpsHabMatrix$`0`/SpsHabMatrix$Surveys)*100
SpsHabMatrix$Pine<-(SpsHabMatrix$`1`/SpsHabMatrix$Surveys)*100
SpsHabMatrix$Wetland<-(SpsHabMatrix$`2`/SpsHabMatrix$Surveys)*100
SpsHabMatrix$Sand<-(SpsHabMatrix$`3`/SpsHabMatrix$Surveys)*100
SpsHabMatrix$Urban<-(SpsHabMatrix$`4`/SpsHabMatrix$Surveys)*100
SpsHabMatrix$Grass<-(SpsHabMatrix$`5`/SpsHabMatrix$Surveys)*100
SpsHabMatrix$HWTC<-(SpsHabMatrix$`6`/SpsHabMatrix$Surveys)*100

SpsHabMatrix$totwater<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='0',])
SpsHabMatrix$totPine<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='1',])
SpsHabMatrix$totWetland<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='2',])
SpsHabMatrix$totSand<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='3',])
SpsHabMatrix$totUrban<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='4',])
SpsHabMatrix$totGrass<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='5',])
SpsHabMatrix$totHWTC<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='6',])
colnames(SpsHabMatrix)[2:8]<-c("Water.records", "Pine.records", "Wetland.records", "Sand.records", "Urban.records", "Grass.records", "HWTC.records")
SpsHabMatrix$percentallsurveys<-(SpsHabMatrix$Surveys)/4242
SpsHabMatrix$random<-as.numeric(0.5)
SpsHabMatrix$zero<-as.numeric(0)


#insertcolumns for total observers who have seen the species

View(SpsHabMatrix)

BAHSpsHabMatrix<-data.frame(SpsHabMatrix[SpsHabMatrix$`COMMON NAME` %in% BAHspecieslist,])

DASSpsHabMatrix<-data.frame(SpsHabMatrix[SpsHabMatrix$`COMMON NAME` %in% DASspecieslist,])
colnames(BAHSpsHabMatrix)<-colnames(DASSpsHabMatrix)<-colnames(SpsHabMatrix)
BAHprop<-Map(prop.test, x=BAHSpsHabMatrix$Pine.records,
                   n=BAHSpsHabMatrix$totPine,
                   p=BAHSpsHabMatrix$random, alternative = "two.sided", correct=TRUE)
BAHprop

Spsdetection<-data.frame(SpsHabMatrix$`COMMON NAME`)
colnames(Spsdetection)<-"Common Name"
Spsdetection$Water<-(SpsHabMatrix$`0`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='0',])))*100
Spsdetection$Pine<-(SpsHabMatrix$`1`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='1',])))*100
Spsdetection$Wetland<-(SpsHabMatrix$`2`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='2',])))*100
Spsdetection$Sand<-(SpsHabMatrix$`3`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='3',])))*100
Spsdetection$Urban<-(SpsHabMatrix$`4`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='4',])))*100
Spsdetection$Grass<-(SpsHabMatrix$`5`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='5',])))*100
Spsdetection$HWTC<-(SpsHabMatrix$`6`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='6',])))*100
BAHSpsdetection<-data.frame(Spsdetection[Spsdetection$`Common Name`%in% BAHspecieslist,])
DASSpsdetection<-data.frame(Spsdetection[Spsdetection$`Common Name`%in% DASspecieslist,])
write.csv(Spsdetection, "Grand Bahama Spsdetection by habitat type 1988 to 2016.csv")
write.csv(BAHSpsdetection, "Bahamian Sps detection by habitat type 1988 to 2016.csv")
write.csv(DASSpsdetection, "Detroit Audubon Society Sps detection by habitat type 1988 to 2016.csv")

##DAS Species Habitat Matrix

DASSpsHabMatrix<-dcast(GBeBirdDEC2016DAS.sub,
                       `COMMON NAME`# + Only use the common name because EcoSimR and Cooccur need dataframes with just the single column and row headings.
                       #`SCIENTIFIC NAME`
                       ~Habitat,
                       value.var = "SAMPLING EVENT IDENTIFIER")
rownames(DASSpsHabMatrix)<-DASSpsHabMatrix$`COMMON NAME`
DASSpsHabMatrix$`COMMON NAME` <-DASSpsHabMatrix$`NA`<-NULL
DASSpsHabMatrix

DASSpsHabMatrix$Surveys<-rowSums(DASSpsHabMatrix[,
                                           c('0', '1', '2','3','4', '5', '6')], na.rm=TRUE)
DASSpsHabMatrix$Water<-
  DASSpsHabMatrix$Pine<-
  DASSpsHabMatrix$Wetland<-
  DASSpsHabMatrix$Sand<-
  DASSpsHabMatrix$Urban<-
  DASSpsHabMatrix$Grass<-
  DASSpsHabMatrix$HWTC<-0
options(digits = 2)
DASSpsHabMatrix$Water<-(DASSpsHabMatrix$`0`/DASSpsHabMatrix$Surveys)*100
DASSpsHabMatrix$Pine<-(DASSpsHabMatrix$`1`/DASSpsHabMatrix$Surveys)*100
DASSpsHabMatrix$Wetland<-(DASSpsHabMatrix$`2`/DASSpsHabMatrix$Surveys)*100
DASSpsHabMatrix$Sand<-(DASSpsHabMatrix$`3`/DASSpsHabMatrix$Surveys)*100
DASSpsHabMatrix$Urban<-(DASSpsHabMatrix$`4`/DASSpsHabMatrix$Surveys)*100
DASSpsHabMatrix$Grass<-(DASSpsHabMatrix$`5`/DASSpsHabMatrix$Surveys)*100
DASSpsHabMatrix$HWTC<-(DASSpsHabMatrix$`6`/DASSpsHabMatrix$Surveys)*100

##detectability per habitat what percentage of surveys in the habitat showed the species
DASSpsHabMatrix$DetWater<-(DASSpsHabMatrix$`0`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='0',])))*100
DASSpsHabMatrix$DetPine<-(DASSpsHabMatrix$`1`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='1',])))*100
DASSpsHabMatrix$DetWetland<-(DASSpsHabMatrix$`2`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='2',])))*100
DASSpsHabMatrix$DetSand<-(DASSpsHabMatrix$`3`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='3',])))*100
DASSpsHabMatrix$DetUrban<-(DASSpsHabMatrix$`4`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='4',])))*100
DASSpsHabMatrix$DetGrass<-(DASSpsHabMatrix$`5`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='5',])))*100
DASSpsHabMatrix$DetHWTC<-(DASSpsHabMatrix$`6`/(nrow(SurAndHabitat[SurAndHabitat$Habitat=='6',])))*100

surveysperhabitat<-data.frame(summary(SurAndHabitat$Habitat))
DASSpeciesHabPresence<-DASSpsHabMatrix
DASSpeciesHabPresence[DASSpeciesHabPresence > 0] <- 1
DASSpeciesHabPresence

##Bahamas Species Habitat Matrix
BAHSpsHabMatrix<-dcast(GBeBirdDEC2016BAH.sub,
                       `COMMON NAME`# + Only use the common name because EcoSimR and Cooccur need dataframes with just the single column and row headings.
                       #`SCIENTIFIC NAME`
                       ~Habitat,
                       value.var = "SAMPLING EVENT IDENTIFIER")
rownames(BAHSpsHabMatrix)<-BAHSpsHabMatrix$`COMMON NAME`
BAHSpsHabMatrix$`COMMON NAME` <-BAHSpsHabMatrix$`NA`<-NULL
BAHSpsHabMatrix

BAHSpeciesHabPresence<-BAHSpsHabMatrix
BAHSpeciesHabPresence[BAHSpeciesHabPresence > 0] <- 1
BAHSpeciesHabPresence

#Observer Habitat Matrices####
ObsHabMatrix<-dcast(GBeBirdDEC2016,
                    `OBSERVER ID`# + Only use the common name because EcoSimR and Cooccur need dataframes with just the single column and row headings.
                    #`SCIENTIFIC NAME`
                    ~Habitat,
                    value.var = "SAMPLING EVENT IDENTIFIER")
rownames(ObsHabMatrix)<-ObsHabMatrix$`OBSERVER ID`
ObsHabMatrix$`OBSERVER ID` <-ObsHabMatrix$`NA`<-NULL
ObsHabPresence<-ObsHabMatrix
ObsHabPresence[ObsHabPresence > 0] <- 1
ObsHabPresence



##Species Location Matrices####
SpsLocMatrix<-dcast(GBeBirdDEC2016,
                    `COMMON NAME`#+
                      #`SCIENTIFIC NAME`
                    ~`LOCALITY ID`,
                    value.var = "SAMPLING EVENT IDENTIFIER")
rownames(SpsLocMatrix)<-SpsLocMatrix$`COMMON NAME`
SpsLocMatrix$`COMMON NAME`<-SpsLocMatrix$`NA`<-NULL
SpsLocMatrix[SpsLocMatrix >= 1] <- 1
SpsLocMatrix

##Observer Location Matrices####
ObsLocMatrix<-dcast(GBeBirdDEC2016,
                    `OBSERVER ID`#+
                    #`SCIENTIFIC NAME`
                    ~`LOCALITY ID`,
                    value.var = "SAMPLING EVENT IDENTIFIER")
rownames(ObsLocMatrix)<-ObsLocMatrix$`OBSERVER ID`
ObsLocMatrix$`OBSERVER ID`<-ObsLocMatrix$`NA`<-NULL
ObsLocMatrix[ObsLocMatrix >= 1] <- 1
ObsLocMatrix

##Observer Species Matrices####
ObsSpsMatrix<-dcast(GBeBirdDEC2016,
                    `OBSERVER ID`#+
                    #`SCIENTIFIC NAME`
                    ~`COMMON NAME`,
                    value.var = "SAMPLING EVENT IDENTIFIER")
rownames(ObsSpsMatrix)<-ObsSpsMatrix$`OBSERVER ID`
ObsSpsMatrix$`OBSERVER ID`<-ObsSpsMatrix$`NA`<-NULL
ObsSpsMatrix[ObsSpsMatrix >= 1] <- 1
ObsSpsMatrix

## Species Observer Matrices####
SpsObsMatrix<-dcast(GBeBirdDEC2016,
                    `COMMON NAME`#+
                    #`SCIENTIFIC NAME`
                    ~`OBSERVER ID`,
                    value.var = "SAMPLING EVENT IDENTIFIER")
rownames(SpsObsMatrix)<-SpsObsMatrix$`COMMON NAME`
SpsObsMatrix$`COMMON NAME`<-SpsObsMatrix$`NA`<-NULL
SpsObsMatrix[SpsObsMatrix >= 1] <- 1
SpsObsMatrix

####Run Co-occurrence analysis in EcoSimR::####
##All Species in all Habitats
ModelSpsHab<-cooc_null_model(speciesData=SpsHabMatrix, suppressProg=TRUE)
summary(ModelSpsHab)
plot(ModelSpsHab, type= "hist")
title(main = "Cooccurrence Null Model Sps and Habitat via EcoSimR")
plot(ModelSpsHab, type= "cooc")
title(main = "Cooccurrence Null Model Sps and Habitat via EcoSimR")
plot(ModelSpsHab, type= "burn_in")
title(main = "Cooccurrence Null Model Sps and Habitat via EcoSimR")

##All Observers in all Habitats
ModelObsHab<-cooc_null_model(speciesData=ObsHabMatrix, suppressProg=TRUE)
summary(ModelObsHab)
plot(ModelObsHab, type= "hist")
title(main = "Cooccurrence Null Model Observers and Habitat via EcoSimR")
plot(ModelObsHab, type= "cooc")
title(main = "Cooccurrence Null Model Observers and Habitat via EcoSimR")
plot(ModelObsHab, type= "burn_in")
title(main = "Cooccurrence Null Model Observers and Habitat via EcoSimR")


##All Species in all Locations##
ModelSpsLoc<-cooc_null_model(speciesData=SpsLocMatrix, suppressProg=TRUE)
summary(ModelSpsLoc)
plot(ModelSpsLoc, type= "hist")
title(main = "Cooccurrence Null Model Sps and Location via EcoSimR")
plot(ModelSpsLoc, type= "cooc")
title(main = "Cooccurrence Null Model Sps and Habitat via EcoSimR")
plot(ModelSpsLoc, type= "burn_in")
title(main = "Cooccurrence Null Model Sps and Habitat via EcoSimR")

##All Observers in all Locations##
ModelObsLoc<-cooc_null_model(speciesData=ObsLocMatrix, suppressProg=TRUE)
summary(ModelObsLoc)
plot(ModelObsLoc, type= "hist")
title(main = "Cooccurrence Null Model Observers and Location via EcoSimR")
plot(ModelObsLoc, type= "cooc")
title(main = "Cooccurrence Null Model Observers and Habitat via EcoSimR")
plot(ModelObsLoc, type= "burn_in")
title(main = "Cooccurrence Null Model Observers and Habitat via EcoSimR")

#THere are so many species and habitats that the plots are overwhelmed

####Run Cooccurrence analysis using the Cooccur package####
#Species Cooccurrence by habitat
SpsHabMatrix
CooccurSpsHAb<-cooccur(SpsHabMatrix, #this must be entered into a vector
        type= "spp_site",
        spp_names = TRUE,
        #prob="comb", #uses the combinatorics approach from Veech 2013 to calculate co-occurrence probabilities
        thresh = TRUE#,
        #only_effects = TRUE, #do not calculate probabilities to reduce calculation time
        #eff_standard = TRUE, #standardize the effect sizes
        #eff_matrix = TRUE #returns a distance matrix of the effect sizes
        )
CooccurSpsHAb

#Observer Cooccurrence by Habitat
ObsHabMatrix
CooccurObsHAb<-cooccur(ObsHabMatrix, #this must be entered into a vector
                       type= "spp_site",
                       spp_names = TRUE,
                       #prob="comb", #uses the combinatorics approach from Veech 2013 to calculate co-occurrence probabilities
                       thresh = TRUE#,
                       #only_effects = TRUE, #do not calculate probabilities to reduce calculation time
                       #eff_standard = TRUE, #standardize the effect sizes
                       #eff_matrix = TRUE #returns a distance matrix of the effect sizes
)
View(CooccurObsHAb$results)

#Species Cooccurrence by Location
CooccurSpsLoc<-cooccur(SpsLocMatrix, #this must be entered into a vector
                       type= "spp_site",
                       spp_names = FALSE,
                       #prob="comb", #uses the combinatorics approach from Veech 2013 to calculate co-occurrence probabilities
                       thresh = TRUE#,
                       #only_effects = TRUE, #do not calculate probabilities to reduce calculation time
                       #eff_standard = TRUE, #standardize the effect sizes
                       #eff_matrix = TRUE #returns a distance matrix of the effect sizes
)
CooccurSpsLoc

#Observer Cooccurrence by Location
CooccurObsLoc<-cooccur(ObsLocMatrix, #this must be entered into a vector
                       type= "spp_site",
                       spp_names = FALSE,
                       #prob="comb", #uses the combinatorics approach from Veech 2013 to calculate co-occurrence probabilities
                       thresh = TRUE#,
                       #only_effects = TRUE, #do not calculate probabilities to reduce calculation time
                       #eff_standard = TRUE, #standardize the effect sizes
                       #eff_matrix = TRUE #returns a distance matrix of the effect sizes
)
CooccurObsLoc

#Observer Cooccurrence by Species seen
CooccurObsSps<-cooccur(ObsSpsMatrix, #this must be entered into a vector
                       type= "spp_site",
                       spp_names = FALSE,
                       #prob="comb", #uses the combinatorics approach from Veech 2013 to calculate co-occurrence probabilities
                       thresh = TRUE#,
                       #only_effects = TRUE, #do not calculate probabilities to reduce calculation time
                       #eff_standard = TRUE, #standardize the effect sizes
                       #eff_matrix = TRUE #returns a distance matrix of the effect sizes
)
CooccurObsSps

#Species Cooccurrence by Observer
CooccurSpsObs<-cooccur(SpsObsMatrix, #this must be entered into a vector
                       type= "spp_site",
                       spp_names = FALSE,
                       #prob="comb", #uses the combinatorics approach from Veech 2013 to calculate co-occurrence probabilities
                       thresh = TRUE#,
                       #only_effects = TRUE, #do not calculate probabilities to reduce calculation time
                       #eff_standard = TRUE, #standardize the effect sizes
                       #eff_matrix = TRUE #returns a distance matrix of the effect sizes
)
CooccurSpsObs

####All Locations and the birds present by date, Observer, Survey and species#### 
GBeBirdDEC2016sps<-dcast(GBeBirdDEC2016, 
                         LOCALITY+
                           `LOCALITY ID`+
                           LONGITUDE+
                           LATITUDE+
                           `OBSERVATION DATE`+
                           `OBSERVER ID`+
                           `SAMPLING EVENT IDENTIFIER`#using Sampling event creates a presence absence type output
                         ~`COMMON NAME`, #each record has a column for the species common name
                         value.var = "OBSERVATION COUNT") #each column include


BahSPSwide<-dcast(GBeBirdDEC2016BAH.sub, 
                  LOCALITY+
                    `LOCALITY ID`+
                    LONGITUDE+
                    LATITUDE+
                    `OBSERVATION DATE`+
                    `OBSERVER ID`+
                    `SAMPLING EVENT IDENTIFIER`#using Sampling event creates a presence absence type output
                    ~`COMMON NAME`, #each record has a column for the species common name
                  value.var = "OBSERVATION COUNT") #each column include

DASSPSwide<-dcast(GBeBirdDEC2016DAS.sub, 
                  LOCALITY+
                    `LOCALITY ID`+
                    LONGITUDE+
                    LATITUDE+
                    `OBSERVATION DATE`+
                    `OBSERVER ID`+
                    `SAMPLING EVENT IDENTIFIER`
                  ~`COMMON NAME`, #each record has a column for the species common name
                  value.var = "OBSERVATION COUNT") #each Common name column includes the number of surveys in which the species was seen.
##Matrices for species in each habitat, location, survey, seen by observer, date and hour of day####
SPSHAB<-dcast(GBeBirdDEC2016,
        `COMMON NAME`
        ~Habitat,
        value.var = "OBSERVATION COUNT")

SPSLOC<-dcast(GBeBirdDEC2016,
              `COMMON NAME`
              ~`LOCALITY ID`,
              value.var = "OBSERVATION COUNT")
SPSLOC
SPSOBS<-dcast(GBeBirdDEC2016,
              `COMMON NAME`
              ~`OBSERVER ID`,
              value.var = "OBSERVATION COUNT")
SPSSUR<-dcast(GBeBirdDEC2016,
              `COMMON NAME`
              ~`SAMPLING EVENT IDENTIFIER`,
              value.var = "OBSERVATION COUNT")
SPSDATE<-dcast(GBeBirdDEC2016,
               `COMMON NAME`
               ~`OBSERVATION DATE`,
               value.var = "OBSERVATION COUNT")
SPSHOUR<-dcast(GBeBirdDEC2016,
               `COMMON NAME`
               ~`TIME OBSERVATIONS STARTED`,
               value.var = "OBSERVATION COUNT")

BahSPSLoc<-dcast(GBeBirdDEC2016BAH.sub, 
                    `LOCALITY ID`+ 
                   LOCALITY+
                   ~`COMMON NAME`, #each record has a column for the species common name
                     value.var = "SAMPLING EVENT IDENTIFIER") #each column includes the number of surveys in which the species was found


BahSPScounts<-dcast(GBeBirdDEC2016BAH.sub, 
                    LOCALITY+
                      `LOCALITY ID`+
                      LONGITUDE+
                      LATITUDE+
                      `OBSERVATION DATE`+
                      `OBSERVER ID`+
                      `SAMPLING EVENT IDENTIFIER`#using Sampling event creates a presence absence type output
                    ~`COMMON NAME`, #each record has a column for the species common name
                    fun.aggregate = sum, #x shows up as NA so we know the species was there but we do not know how many.
                    value.var = "OBSERVATION COUNT") #each column include

DASSPScounts<-dcast(GBeBirdDEC2016DAS.sub, 
                    LOCALITY+
                      LONGITUDE+
                      LATITUDE+
                      `OBSERVATION DATE`+
                      `OBSERVER ID`+
                      `SAMPLING EVENT IDENTIFIER`#using Sampling event creates a presence absence type output
                    ~`COMMON NAME`, #each record has a column for the species common name
                    fun.aggregate = sum, #x shows up as NA so we know the species was there but we do not know how many.
                    value.var = "OBSERVATION COUNT") #each Common name column includes the number of surveys in which the species was seen.
##The number of surveys at each location where the species was reported
BahSPSLoc<-dcast(GBeBirdDEC2016BAH.sub, 
                 LOCALITY+
                   `LOCALITY ID`
                 #using Sampling event creates a presence absence type output
                 ~`COMMON NAME`, #each record has a column for the species common name
                 value.var = "SAMPLING EVENT IDENTIFIER") #each column include

DASSPSLoc<-dcast(GBeBirdDEC2016DAS.sub, 
                 LOCALITY
                 #using Sampling event creates a presence absence type output
                 ~`COMMON NAME`, #each record has a column for the species common name
                 value.var = "SAMPLING EVENT IDENTIFIER") #each column include

##The number of surveys at each location where the observer participated
BahSPSObs<-dcast(GBeBirdDEC2016BAH.sub, 
                 LOCALITY
                 #using Sampling event creates an count of how many times the observer visited
                 ~`OBSERVER ID`, #each record has a column for the species common name
                 value.var = "SAMPLING EVENT IDENTIFIER") #each column include

DASSPSObs<-dcast(GBeBirdDEC2016DAS.sub, 
                 LOCALITY
                 #using Sampling event creates an count of how many times the observer visited
                 ~`OBSERVER ID`, #each record has a column for the species common name
                 value.var = "SAMPLING EVENT IDENTIFIER") #each column include
#extract the raster values at each location
#Bahamian species locations
BAHSPSxy<-data.frame(BahSPSwide$LONGITUDE,BahSPSwide$LATITUDE)
#DAS species locations
DASSPSxy<-data.frame(DASSPSwide$LONGITUDE,DASSPSwide$LATITUDE)
#All SPS locations

GBeBirdDEC2016xy<-data.frame(GBeBirdDEC2016sps$LONGITUDE,GBeBirdDEC2016sps$LATITUDE)

#make the species data into a spatial points data frame
coordinates(BahSPSwide) <- ~LONGITUDE+LATITUDE 
coordinates(DASSPSwide) <- ~LONGITUDE+LATITUDE


#Extract the values from each raster at the the ebird locations 

BAHSPSextract<-extract(RF7Classes2017,BAHSPSxy)
DASSPSextract<-extract(RF7Classes2017,DASSPSxy)
AllSPSextract<-extract(RF7Classes2017,GBeBirdDEC2016xy)
summary(GBeBirdDEC2016)

#Insert the habitat types into the PA and count data frames as a new column
#Make sure the values are converted to factors
BahSPSwide$habitat<-as.factor(BAHSPSextract)
BahSPScounts$habitat<-as.factor(BAHSPSextract)
DASSPSwide$habitat<-as.factor(DASSPSextract)
DASSPScounts$habitat<-as.factor(DASSPSextract)
GBeBirdDEC2016sps$habitat<-as.factor(AllSPSextract)
#Combine ebird observations, habitat classes and locations data into one data frame.
# Simple Bar Plot based on habitat by survey
AllLocHabitat <- table(GBeBirdDEC2016sps$habitat)
max(AllLocHabitat)
bpAll<-barplot(AllLocHabitat, main="Habitat type in all eBird Surveys on Grand Bahama", 
        xlab="Habitat type",horiz = FALSE,
        ylim=c(0,(1.1*max(AllLocHabitat))),
        names.arg=c("Water", 
                    "Pine", 
                    "Wetland", 
                    "Sand", 
                    "Urban", 
                    "Grass", 
                    "HWTC"))
text(bpAll, AllLocHabitat, label=AllLocHabitat, pos = 3, xpd = NA)

BAHLocHabitat <- table(BahSPSwide$habitat)
max(BAHLocHabitat)
bpBAH<-barplot(BAHLocHabitat, 
               main="Habitat type in eBird Surveys with Bahamian Species on Grand Bahama", 
               col=col7class,
               xlab="Habitat type",horiz = FALSE,
               ylim=c(0,(1.1*max(BAHLocHabitat))),
               names.arg=c("Water", 
                           "Pine", 
                           "Wetland", 
                           "Sand", 
                           "Urban", 
                           "Grass", 
                           "HWTC"))
text(bpBAH, BAHLocHabitat, label=BAHLocHabitat, pos = 3, xpd = NA)
dev.off()

DASLocHabitat <- table(DASSPSwide$habitat)
max(DASLocHabitat)
bpDAS<-barplot(DASLocHabitat, 
               col=col7class,
               main="Habitat type in eBird surveys with DAS species on Grand Bahama", 
               xlab="Habitat type",
               horiz = FALSE,
               ylim=c(0,(1.1*max(DASLocHabitat))),
               names.arg=c("Water", 
                           "Pine", 
                           "Wetland", 
                           "Sand", 
                           "Urban", 
                           "Grass", 
                           "HWTC"))
text(bpDAS, DASLocHabitat, label=DASLocHabitat, pos = 3, xpd = NA)

#combinedhabitat data
combinedHAB<-matrix(c(AllLocHabitat,
         BAHLocHabitat,
         DASLocHabitat),nr=3, byrow = TRUE)
#create a vector of the number of pixels in each habitat class

dev.off()
bpcombined<-barplot(combinedHAB, 
        main="Habitat type in eBird surveys on Grand Bahama", 
        xlab="Habitat type",
        ylab="Number of Surveys",
        ylim=c(0,(1.1*max(combinedHAB))),
        beside=T, 
        col=c("#a6cee3","#1f78b4","#b2df8a"),#these colors were chosen because they are a colorblind friendly palette 
        names.arg=c("Water", 
                    "Pine", 
                    "Wetland", 
                    "Sand", 
                    "Urban", 
                    "Grass", 
                    "HWTC"))
legend("topleft", c("All Species","BAH Species", "DAS Species"), pch=15, 
       col=c("#a6cee3","#1f78b4","#b2df8a"),#make sure these colors agree with those above
       
       bty="n")
text(bpcombined, combinedHAB, label=combinedHAB, pos = 3, xpd = NA)


####--------------------------------------------####
#Create tables of habitat and Locality only
LocPerHabAll<-dcast(GBeBirdDEC2016sps, 
                    LOCALITY+
                      LONGITUDE+
                      LATITUDE+
                      `OBSERVATION DATE`+
                      `OBSERVER ID`+
                      `SAMPLING EVENT IDENTIFIER`
                    ~`COMMON NAME`, #each record has a column for the species common name
                    value.var = "OBSERVATION COUNT")
LocPerHabDAS
LocPerHabBAH
LocPerHabAll





write.csv(BahSPSwide, file = "BAHSPSPAhabitat.csv")
write.csv(BahSPScounts, file = "BahSPScountshabitat.csv")
write.csv(DASSPSwide, file = "DASSPSPAhabitat.csv")
write.csv(DASSPScounts, file = "DASSPScountshabitat.csv")


#//////////////////////////////////////////////////////////////////////////////////////////####
#------STOP!! GO AND CHANGE THE PREVIOUS FILE NAME SO YOU DO NOT SAVE OVER IT AGAIN!!----####
#//////////////////////////////////////////////////////////////////////////////////////////####
#I added the month and year of the final entry "DEC2016"

#Generate a bar plot of the habitat types of the Bahamian species locations

summary(BahSPSwide$habitat)
summary(BahSPSwide)
  plot(BahSPSwide$habitat, 
       main="Bahamian Species Habitat",
       ylim=c(0,1600),
       ylab = "Count", col="steelblue", las = 2, labels=TRUE,
       xaxt = "n", 
       xlab='Habitat')
       axis(1,at=c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9), labels=c("Water", "Pine", "Wetland", "Sand", "Urban", "Grass", "HWTC"))
       
#Generate a bar plot of the habitat types of the DAS species locations
       
       summary(DASSPSwide$habitat)
       summary(DASSPSwide)
       plot(DASSPSwide$habitat, 
            main="DAS Species Habitat",
            ylim=c(0,1600),
            ylab = "Count", col="steelblue", las = 2,
            xaxt = "n", 
            xlab='Habitat')
       axis(1,at=c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9), labels=c("Water", "Pine", "Wetland", "Sand", "Urban", "Grass", "HWTC"))

#Generate a bar plot of the habitat types of all species locations       
              
       plot(GBeBirdDEC2016sps$habitat, 
            main="All Species Habitat",
            ylim=c(0,2000),
            ylab = "Count", col="steelblue", las = 2,
            xaxt = "n", 
            xlab='Habitat')
       axis(1,at=c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9), labels=c("Water", "Pine", "Wetland", "Sand", "Urban", "Grass", "HWTC"))
       
       
propAll<-AllLocHabitat/sum(AllLocHabitat)
proplandcover<-landcover/sum(landcover)
propBAH<-BAHLocHabitat/sum(BAHLocHabitat)
propDAS<-DASLocHabitat/sum(DASLocHabitat)

##reorganize data by habitat##
HabitatByLocSps<-dcast(GBeBirdDEC2016, 
                       Habitat+
                         LOCALITY+
                         LONGITUDE+
                         LATITUDE+
                         `OBSERVATION DATE`+
                         `OBSERVER ID`+
                         `SAMPLING EVENT IDENTIFIER`+
                         `DURATION MINUTES`
                       
                       ~`COMMON NAME`, #each record has a column for the species common name
                       value.var = "OBSERVATION COUNT")
HabitatByLocSps$Habitat<-as.factor(HabitatByLocSps$Habitat)

levels(HabitatByLocSps$Habitat)
ggboxplot(HabitatByLocSps, x = "Habitat", y = "`DURATION MINUTES`", 
          color = "Habitat", 
          palette = c("#40A4DF", 
                      "#04460E", 
                      "#B5AB40", 
                      "#00ff00",
                      "#FF0000",
                      "#52D061",
                      "#A76B29" ),
          order = c("0", "1", "2", "3", "4", "5", "6"),
          ylab = "Minutes", xlab = "Habitat")
kruskal.test(`DURATION MINUTES` ~ Habitat, data = HabitatByLocSps)
wilcox.test(`DURATION MINUTES`~Habitat,data=HabitatByLocSps)
pairwise.wilcox.test(HabitatByLocSps$`DURATION MINUTES`, HabitatByLocSps$Habitat, p.adjust.method = p.adjust.methods,
                     paired = FALSE)

minsperhab<-aggregate(HabitatByLocSps$`DURATION MINUTES`, 
          by=list(Category=HabitatByLocSps$Habitat), 
          FUN=sum, na.rm=TRUE)
minsperhab
locsperhab<-with(GBeBirdDEC2016, tapply(`LOCALITY ID`, Habitat, FUN = function(x) length(unique(x))))
survperhab<-with(HabitatByLocSps, tapply(`SAMPLING EVENT IDENTIFIER`, Habitat, FUN = function(x) length(unique(x))))
Obsperhab<-with(HabitatByLocSps, tapply(HabitatByLocSps$`OBSERVER ID`, Habitat, FUN = function(x) length(unique(x))))
spsperhab<-with(GBeBirdDEC2016, tapply(GBeBirdDEC2016$`COMMON NAME`, Habitat, FUN = function(x) length(unique(x))))
spsperhab
Obsperhab
surveysperhabitat<-aggregate(HabitatByLocSps$`DURATION MINUTES`, 
          by=list(Category=HabitatByLocSps$LOCALITY), 
          FUN=length, 
          na.rm=TRUE)
surveydata<-as.data.frame(
  unique(GBeBirdDEC2016
         [,c('Habitat',
             'LOCALITY ID',
             'OBSERVER ID',
             'SAMPLING EVENT IDENTIFIER',
             'DURATION MINUTES')]))
minsperhab<-aggregate(HabitatByLocSps$`DURATION MINUTES`, 
                      by=list(Category=HabitatByLocSps$Habitat), 
                      FUN=sum, na.rm=TRUE)

OBSperhab<-aggregate(GBeBirdDEC2016$`OBSERVER ID`, 
                      by=list(Category=GBeBirdDEC2016$Habitat, GBeBirdDEC2016$`LOCALITY ID`), 
                      FUN=length, na.rm=TRUE)

Locationdata<-dcast(GBeBirdDEC2016, 
                       Habitat+
                         LOCALITY+
                         LONGITUDE+
                         LATITUDE+
                         `OBSERVATION DATE`+
                         `OBSERVER ID`+
                         `SAMPLING EVENT IDENTIFIER`+
                         `DURATION MINUTES`
                       
                       ~`COMMON NAME`, #each record has a column for the species common name
                       value.var = "OBSERVATION COUNT")

plot(`DURATION MINUTES`~ Habitat, Locationdata)
#HabitatByObsSps
#HabitatBySurSps
#HabitatByObsSur
#HabitatByObsLoc
#HabitatBySurSps

###Locations per habitat
LocAndHabitat<- data.frame(GBeBirdDEC2016$LOCALITY,GBeBirdDEC2016$`LOCALITY ID`,GBeBirdDEC2016$Habitat)
LocAndHabitat<-unique(LocAndHabitat)
colnames(LocAndHabitat)<-c("Locality", "Locality ID", "Habitat")

LocAndHabitat$Habitat<-as.factor(LocAndHabitat$Habitat)
LocAndHabitat<-na.omit(LocAndHabitat)



#Is the reporting of the species different than expected in that habitat
prop.test(x, n, p = NULL,
          alternative = c("two.sided", "less", "greater"),
          conf.level = 0.95, correct = TRUE)
SpsHabMatrix$totwater<- nrow(SurAndHabitat[SurAndHabitat$Habitat=='0',])
SpsHabMatrix$percentallsurveys<-(SpsHabMatrix$Surveys)/4242
SpsHabMatrix$random<-0.5
SpsHabMatrix$zero<-0
prop.test(BAHSpsHabMatrix$WaterSurveys,
          BAHSpsHabMatrix$totwater,
          BAHSpsHabMatrix$random, alternative = "greater", correct=TRUE)
BAHSpsHabMatrix
DASSpsHabMatrix

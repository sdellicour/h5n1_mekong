## Using R to define most probable sampling polygons assigned to sequences involved in a continuous phylogeographic analysis performed with the program BEAST

This repo gathers input files and scripts related to our study entitled "Incorporating heterogeneous sampling probabilities in continuous phylogeographic inference — Application to H5N1 spread in the Mekong region" ([Dellicour *et al*. 2020](https://academic.oup.com/bioinformatics/article/36/7/2098/5650406#supplementary-data)). The present tutorial constitutes an example on how preparing a continuous phylogeographic inference using heterogeneous sampling priors associated with the tips.

The continuous phylogeographic method developed by Lemey *et al*. (2010) and implemented in the software BEAST (Drummond *et al*. 2012, Suchard *et al*. 2018) can be used to infer virus dispersal history in a continuous space. This method relies on the geographic coordinates associated with the sampling locations of virus sequences. Yet, frequently, only the more or less broad administrative area of origin of viral sequences is known (e.g. a broad administrative area "admin-1" or only the country of origin). The low precision of sampling locations can represent a notable limitation in the field of viral phylogeography. Indeed, when precise coordinates are unavailable for sampling locations, a common practice consists in either discarding the sequence or to take the centroid coordinates of the administrative area within which the sequence has been sampled. However, the latter option can be inappropriate because in many cases, taking the centroid coordinates can be irrelevant, e.g., regarding the very low or even unlikely probability that the sequence has been sampled around these coordinates. The present appendix describes the new method used to define smaller polygons within the overall administrative area, and to associate all these smaller "sub-polygons" with a distinct sampling probability. The generated collection of sub-polygons can then be specified in BEAST (Drummond *et al*. 2012, Suchard *et al*. 2018) to define the prior distribution of sampling coordinates associated with each sampled sequence. The idea behind this approach is to allow an appropriate use of virus sequences with low sampling precision in continuous phylogeographic reconstructions. The following approach thus proposes to include these sequences in the analysis but while using external information such as epidemiological records or host densities to specify a sampling probability assigned to some sub-polygons (admin-2 polygons) contained in the overall administrative area of origin (admin-1 polygon). We here describe how to use R to prepare such sub-polygons for the continuous phylogeographic analysis of the H5N1 clade 1 spread through the Mekong region. The data set is made of 214 RNA sequences (HA gene). While the admin-2 of origin is known for 109 of them, only the admin-1 polygon of origin if known for the remaining 105 sequences (see Figure 2 for the related sampling map of clade 1).

** Step 1: installing and loading required R libraries
```R
install.packages("ape"); library(ape)
install.packages("rgdal"); library(rgdal)
install.packages("rgeos"); library(rgeos)
install.packages("raster"); library(raster)
```

** Step 2: loading admin-1 and admin-2 borders (shapefiles)
```R
admin_1 = readOGR(dsn="./admin_shapefiles", layer="admin_1_Mekong")
admin_2 = readOGR(dsn="./admin_shapefiles", layer="admin_2_Mekong")
```
In addition, we also have to load a one column table gathering the admin-1 ID of each admin-2 polygons gathered in the "admin_2" SpatialPolygonsDataFrame:
```R
admin2_admin1ID = read.table("admin_shapefiles/Admin2_admin1ID.txt", header=T)
```

** Step 3: reading sampling data associated with each sequence from a fasta alignment
```R
fasta = read.dna("H5N1_clade1.fasta", format="fasta")
names = row.names(fasta)
sequenceIDs = matrix(nrow=length(names), ncol=1)
coordinates = matrix(nrow=length(names), ncol=2)
precisions = matrix(nrow=length(names), ncol=1)
years = matrix(nrow=length(names), ncol=1)
for (i in 1:length(names)) {
	sequenceIDs[i,1] = unlist(strsplit(names[i],"_"))[1]
	coordinates[i,] = unlist(strsplit(names[i],"_"))[4:5]
	p = unlist(strsplit(names[i],"_"))[6]
	precisions[i,1] = as.numeric(unlist(strsplit(p,"-admin"))[2])
	if (is.na(precisions[i,1])) precisions[j,1] = 1
	if (precisions[i,1] == 3) precisions[i,1] = 2
	years[i,1] = as.numeric(unlist(strsplit(names[i],"_"))[3])
}
```
The "precisions" values indicate the precision level of the different sampling locations, i.e. if we know the admin-2 or only the admin-1 administrative area of origin.

** Step 4: loading hosts (chicken and duck) density rasters
```R
chickens_density = raster("Mekong_rasters/Chickens_Mekong.asc")
ducks_density = raster("Mekong_rasters/Ducks_Mekong.asc")
```
By multiplying raster cell density values by their corresponding cell area, we can then transform density measures by estimates of absolute numbers of host individuals:
```R
chickens_number = chickens_density*raster::area(chickens_density)
ducks_number = ducks_density*raster::area(ducks_density)
hosts = chickens_number
hosts[] = hosts[]+ducks_number[]
```
Note that we here sum chicken and duck numbers to obtain an overall "host" number for each raster cell. Finally, we also log-transform resulting numbers in order to avoid providing an excessive importance to a few raster cells associated with high values.

** Step 5: loading H5N1 occurrence records  
These annual records are stored in different "csv" files and can be read by the following loop:
```R
records_all = c()
records_list = list()
for (year in 2003:2012) {
	records_list[[year-2002]] = read.csv(paste("H5N1_records/H5N1_",year,".csv",sep=""), header=T)[,cbind("LON","LAT")]
	records_all = rbind(records_all, records_list[[year-2002]])
}
```
See also Figure S1 for a graphical representation of these yearly records.

** Step 6: assigning a number of hosts to each admin-2 polygon
```R
hosts_counts = c()
records_counts_list = list()
for (i in 1:length(admin_2@polygons)) {
	hCounts = 0
	for (j in 1:length(admin_2@polygons[[i]]@Polygons)) {
		p = Polygon(coords)
		ps = Polygons(list(p),1)
		sps = SpatialPolygons(list(ps))
		r = mask(crop(hosts, coords), sps)
		hCounts = hCounts + sum(r[], na.rm=T)	
	}
	hosts_counts = c(hosts_counts, log(hCounts+1))	
}
```

** Step 7: assigning an annual number of H5N1 records to each admin-2 polygon
```R
for (h in 2:length(records_list)) {
	buffer = c()
	for (i in 1:length(admin_2@polygons)) {
		rCounts = 0; hCounts = 0
		for (j in 1:length(admin_2@polygons[[i]]@Polygons)) {
			coords = admin_2@polygons[[i]]@Polygons[[j]]@coords
			for (k in 1:dim(records_list[[h]])[1]) {
				if (point.in.polygon(records_list[[h]][k,1], records_list[[h]][k,2], coords[,1], coords[,2], mode.checked=FALSE) == 1) {
					rCounts = rCounts+1
				}
			}
		}
		buffer = c(buffer, rCounts)		
	}
	records_counts_list[[h]] = buffer	
}
```
In this step, we thus create a list of H5N1 record values for each year between 2003 and 2014.  

** Step 8: generating a collection of admin-2 polygons for each sampled sequence  
For sequences already assigned to an admin-2 polygon, the following script will simply generate a KML ("keyhole markup language") file specifying a single admin-2 polygon with an associated sampling probability equal to 1. While for sequences assigned to an admin-1 polygon, the script will generate a KML file specifying related admin-2 polygons each associated with a specific sampling probability (the sampling probabilities of the different admin-2 polygons contained in the admin-1 polygon summing to 1). If possible, the sampling probability assigned to each admin-2 polygon will be estimated according to the H5N1 outbreak records. In that case, the sampling probability assigned to an admin-2 polygon is simply the number of H5N1 outbreak records associated to this polygon and divided by the total number of H5N1 outbreak records in the admin-1 polygon. It is worth noting that to estimate these sampling probabilities we only consider H5N1 outbreak records for the sampling year of the considered sequence. When no annual H5N1 outbreak record is available for at least one of the admin-2 polygon, we can then use host incidence data (log-transformed host species numbers) to define the sampling probability to assign to these polygons. In the latter case, we follow the same logic as for H5N1 outbreak records, but this time defining the sampling probability as the ratio between log-transformed host species numbers in the admin-2 polygon and in the admin-1 polygon.
```R
admins = list()
admins[[1]] = admin_1
admins[[2]] = admin_2
for (i in 1:dim(coordinates)[1]) {
	admin = admins[[as.numeric(precisions[i])]]
	sampling_coord = cbind(as.numeric(coordinates[i,1]),as.numeric(coordinates[i,2]))
	sequenceID = sequenceIDs[i,1]
	index = 0
	for (j in 1:length(admin@polygons)) {
		for (k in 1:length(admin@polygons[[j]]@Polygons)) {
			if (point.in.polygon(sampling_coord[1], sampling_coord[2], admin@polygons[[j]]@Polygons[[k]]@coords[,1], admin@polygons[[j]]@Polygons[[k]]@coords[,2], mode.checked=FALSE) == 1) {
				index = j
			}
		}
	}
	if ((index == 0) & (precisions[i]==2)) {
		# print(paste("No initial administrative polygon admin-2 found for ",sequenceID,sep=""))
		counter = 0
		while (index == 0) {
			if (runif(1,0,1) < 0.5) {
				sampling_coord[1] = sampling_coord[1] - 0.001
			}	else	{
				sampling_coord[1] = sampling_coord[1] + 0.001
			}
			if (runif(1,0,1) < 0.5) {
				sampling_coord[2] = sampling_coord[2] - 0.001
			}	else	{
				sampling_coord[2] = sampling_coord[2] + 0.001
			}
			for (j in 1:length(admin@polygons)) {
				for (k in 1:length(admin@polygons[[j]]@Polygons)) {
					if (point.in.polygon(sampling_coord[1], sampling_coord[2], admin@polygons[[j]]@Polygons[[k]]@coords[,1], admin@polygons[[j]]@Polygons[[k]]@coords[,2], mode.checked=FALSE) == 1) {
						index = j
					}
				}
			}
		}
	}
	maxArea = 0; polIndex = 0
	for (j in 1:length(admin@polygons[[index]]@Polygons)) {
		if (maxArea < admin@polygons[[index]]@Polygons[[j]]@area) {
			maxArea = admin@polygons[[index]]@Polygons[[j]]@area
			polIndex = j
		}
	}
	pol1 = admin@polygons[[index]]@Polygons[[polIndex]]
	p = Polygon(pol1@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
	pol1 = sps; proj4string(pol1) = CRS("+init=epsg:4326")
	coords = admin@polygons[[index]]@Polygons[[polIndex]]@coords
	if (precisions[i] == 2) {
		sink(file=paste("H5N1_polygons/",sequenceID,".kml",sep=""))
		cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
		cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
		cat(paste("\t<polygon id=\"",paste(sequenceID,sep="_"),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
		cat("\t\t<coordinates>"); cat("\n")
		for (j in 1:dim(coords)[1]) {
			cat(paste("\t\t\t",coords[j,2],",",coords[j,1],",0",sep="")); cat("\n")
		}
		cat("\t\t</coordinates>"); cat("\n")
		cat("\t</polygon>"); cat("\n")
		cat("</kml>"); cat("\n")
		sink(NULL)
	}
	if (precisions[i] == 1) {
		indices = which(admin2_admin1ID == index)
		year_index = years[i,1]-2002
		records_selected = records_counts_list[[year_index]]
		oRecords = records_counts[indices]
		hRecords = hosts_counts[indices]
		records = c()
		if (sum(oRecords) > 0) {
			records = oRecords
		}	else	{
			records = hRecords
		}
		pols2_not_included = c()
		sink(file=paste("H5N1_polygons/",sequenceID,".kml",sep=""))
		cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
		cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
		firstPlot = TRUE
		for (j in 1:length(indices)) {
			fValue = records[j]/sum(records)						
			if (fValue > 0) {
				maxArea = 0; polIndex2 = 0
				for (k in 1:length(admin_2@polygons[[indices[j]]]@Polygons)) {
					if (maxArea < admin_2@polygons[[indices[j]]]@Polygons[[k]]@area) {
						maxArea = admin_2@polygons[[indices[j]]]@Polygons[[k]]@area
						polIndex2 = k
					}
				}
				coords = admin_2@polygons[[indices[j]]]@Polygons[[polIndex2]]@coords
				cat(paste("\t<polygon id=\"",paste(sequenceID,j,sep="_"),"\" samplingProbability=\"",fValue,"\">",sep=""))
				cat("\n")
				cat("\t\t<coordinates>"); cat("\n")
				for (l in 1:dim(coords)[1]) {
					cat(paste("\t\t\t",coords[l,2],",",coords[l,1],",0",sep="")); cat("\n")
				}
				cat("\t\t</coordinates>"); cat("\n")
				cat("\t</polygon>"); cat("\n")
				pol2 = admin_2@polygons[[indices[j]]]@Polygons[[polIndex2]]
				p = Polygon(pol2@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps)); 
				pol2 = sps; proj4string(pol2) = CRS("+init=epsg:4326")
				if (gIntersects(pol1, pol2) == TRUE) {
					pol3 = gIntersection(pol1, pol2)
					if (class(pol3) == "SpatialPolygons") {
						area2 = pol2@polygons[[1]]@area
						area3 = pol3@polygons[[1]]@area
						if (area3 < (0.5*area2)) {
							pols2_not_included = c(pols2_not_included, j)
						}
					}
				}
			}
		}
		cat("</kml>"); cat("\n")
		sink(NULL)
	}
}
```
In the end, a distinct KML file is thus generated per sampling sequence and saved in a folder here named "H5N1_polygons". In this case, the name of each KML file is simply the accession number of the corresponding sampling sequence (but could, of course, be any ID found in the sequence name).

** Step 9: generating a BEAST XML file based on a template file previously generated with BEAUti
This last step is optional as XML files can also be manually edited to add text parts linking the analysis to the different KML files generated above. As detailed in the resulting example file "H5N1_clade1.xml" available with this appendix, four new text blocks have to be added in the XML file to refer to these sampling polygons. For each sampled sequence: (i) a new "leafTraitParameter", (ii) a "flatGeoSpatialPrior" referencing the external KML file, and (iii) an "uniformGeoSpatialOperator" have to be added. In addition, "flatGeoSpatialPriors" have also to be listed in the operators section (see the content of the KML file for further details).
```R
xml = scan(file="BEAST_template.xml", what="", sep="\n", quiet=T)
directory = "H5N1_polygons"
sink(file="H5N1_clade1.xml")
for (i in 1:length(xml)) {
	cat(xml[i]); cat("\n")
	if (xml[i]=="\t</continuousDiffusionStatistic>") {
		cat("\n")
		for (j in 1:length(sequenceIDs)) {
			cat(paste("\t<leafTraitParameter id=\"",sequenceIDs[j],".trait\" taxon=\"",names[j],"\">",sep="")); cat("\n")
			cat(paste("\t\t<treeModel idref=\"treeModel\"/>",sep="")); cat("\n")
			cat(paste("\t\t<parameter idref=\"leaf.location\"/>",sep="")); cat("\n")
			cat(paste("\t</leafTraitParameter>",sep="")); cat("\n")
		}
		cat("\n")
		for (j in 1:length(sequenceIDs)) {
			cat(paste("\t<flatGeoSpatialPrior id=\"",sequenceIDs[j],"_polygons\" taxon=\"",names[j],"\" kmlFileName=\"",directory,"/",sequenceIDs[j],".kml\" inside=\"true\" union=\"true\" cache=\"true\">",sep="")); cat("\n")
			cat(paste("\t\t<data>",sep="")); cat("\n")
			cat(paste("\t\t\t<parameter idref=\"",sequenceIDs[j],".trait\"/>",sep="")); cat("\n")
			cat(paste("\t\t</data>",sep="")); cat("\n")
			cat(paste("\t</flatGeoSpatialPrior>",sep="")); cat("\n")
		}
		cat("\n")		
	}
	if (xml[i]=="\t\t</precisionGibbsOperator>") {
		cat("\n")
		for (j in 1:length(sequenceIDs)) {
			cat(paste("\t\t<uniformGeoSpatialOperator weight=\"0.01\">",sep="")); cat("\n")
			cat(paste("\t\t\t<parameter idref=\"",sequenceIDs[j],".trait\"/>",sep="")); cat("\n")
			cat(paste("\t\t\t<flatGeoSpatialPrior idref=\"",sequenceIDs[j],"_polygons\"/>",sep="")); cat("\n")
			cat(paste("\t\t</uniformGeoSpatialOperator>",sep="")); cat("\n")
		}
		cat("\n")
	}
	if (xml[i]=="\t\t\t\t<multivariateWishartPrior idref=\"location.precisionPrior\"/>") {
		cat("\n")
		cat("\t\t\t\t<geoDistributionCollection id=\"allGeoDistributions\">"); cat("\n")
		for (j in 1:length(sequenceIDs)) {
			cat(paste("\t\t\t\t<flatGeoSpatialPrior idref=\"",sequenceIDs[j],"_polygons\"/>",sep="")) 
			cat("\n")
		}
		cat("\t\t\t\t</geoDistributionCollection>"); cat("\n")
		cat("\n")
	}					
	if (xml[i]=="\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood\"/>") {
		if (xml[i-3]=="\t\t\t<strictClockBranchRates idref=\"branchRates\"/>") {
			cat("\n")
			for (j in 1:length(sequenceIDs)) {
				cat(paste("\t\t\t<leafTraitParameter idref=\"",sequenceIDs[j],".trait\"/>",sep="")); cat("\n")
			}
			cat("\n")
		}
	}	
}
sink(NULL)
```

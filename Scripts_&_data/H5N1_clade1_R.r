install.packages("ape"); library(ape)
install.packages("rgdal"); library(rgdal)
install.packages("rgeos"); library(rgeos)
install.packages("raster"); library(raster)

# Step 2: loading admin-1 and admin-2 borders (shapefiles)

admin_1 = shapefile("GADM_shapefiles/GADM_1_Mekong.shp")
admin_2 = shapefile("GADM_shapefiles/GADM_2_Mekong.shp")
admin2_admin1ID = read.table("GADM_shapefiles/Admin2_admin1ID.txt", header=T)

# Step 3: reading sampling data associated with each sequence from a fasta alignment

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

# Step 4: loading hosts (chicken and duck) density rasters

chickens_density = raster("Mekong_rasters/Chickens_Mekong.asc")
ducks_density = raster("Mekong_rasters/Ducks_Mekong.asc")

chickens_number = chickens_density*raster::area(chickens_density)
ducks_number = ducks_density*raster::area(ducks_density)
hosts = chickens_number
hosts[] = hosts[]+ducks_number[]

#Step 5: loading H5N1 occurrence records

records_all = c()
records_list = list()
for (year in 2003:2012) {
	records_list[[year-2002]] = read.csv(paste("H5N1_records/H5N1_",year,".csv",sep=""), header=T)[,cbind("LON","LAT")]
	records_all = rbind(records_all, records_list[[year-2002]])
}

# Step 6: assigning a number of hosts to each admin-2 polygon

hosts_counts = c()
for (i in 1:length(admin_2@polygons)) {
	hCounts = 0
	for (j in 1:length(admin_2@polygons[[i]]@Polygons)) {
		coords = admin_2@polygons[[i]]@Polygons[[j]]@coords
		p = Polygon(coords)
		ps = Polygons(list(p),1)
		sps = SpatialPolygons(list(ps))
		r = mask(crop(hosts, coords), sps)
		hCounts = hCounts + sum(r[], na.rm=T)	
	}
	hosts_counts = c(hosts_counts, log(hCounts+1))	
}

# Step 7: assigning the number of H5N1 records to each admin-2 polygon

records_counts = c()
for (i in 1:length(admin_2@polygons)) {
	rCounts = 0
	for (j in 1:length(admin_2@polygons[[i]]@Polygons)) {
		coords = admin_2@polygons[[i]]@Polygons[[j]]@coords
		for (k in 1:dim(records_all)[1]) {
			if (point.in.polygon(records_all[k,1], records_all[k,2], coords[,1], coords[,2], mode.checked=FALSE) == 1) {
				rCounts = rCounts+1
			}
		}
	}
	records_counts = c(records_counts, rCounts)		
}

# Step 8: generating a collection of admin-2 polygons for each sampled sequence

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

# Step 9: generating a BEAST XML file based on a template file previously generated with BEAUti

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


library(randomForest)
install.packages("FactoMineR")
library(FactoMineR)

rm(list = setdiff(ls(), lsf.str()))

#### Okay, now that we have done the k-means clustering let's loop through the results ###
### then we can do the random forest with this


setwd("/Users/carlschmidt/Desktop/xiaoke_2017_random_forest/oct_27th")


### go through the code to do the combinations, now ###
####Permutheus####
###Will need to change this to be specific to your local setup
totalMetab = read.table("/Users/carlschmidt/Downloads/xiaoke_metabolon_csv.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

###Will need to change this to be specific to your local setup
###Read in the liver data and subset
bigTable1 = read.table("/Users/carlschmidt/Downloads/All799_xiaoke_tsv.txt", sep = "\t", header = TRUE)

###Set interest genes back to the liver core enriched ones !
myTableOfInterest = read.table("/Users/carlschmidt/Downloads/liver_enriched_genes_xiaoke.txt", header = TRUE, sep = "\n", stringsAsFactors=FALSE)

#######MAY NEED TO CHANGE FOR YOUR SETUP DR SCHMIDT
###Note that you will need to make this directory on your computer to hold
#the output
###############

####For liver:
#we will need to get all of the liver libraries
#and then from that smaller table, subset from the RTE genes
#from the RTE genes
subsetLibsOfInterest <- c(1321,1322,1323,1324,1325,1326,1327,1328)
subsetLibsOfInterest2 <- c(1330,1331,1332,1333,1334,1335,1336,1337)

###temporarily join them both here !!!
subsetLibsOfInterest = c(subsetLibsOfInterest, subsetLibsOfInterest2)

###Declare the subtable as a frame
subsetClusterOfInterest = matrix(, nrow = length(bigTable1$Label), ncol = length(subsetLibsOfInterest))
subsetClusterOfInterestNames = as.vector(subsetLibsOfInterest)
colnames(subsetClusterOfInterest) = subsetClusterOfInterestNames
subsetClusterOfInterest = as.data.frame(subsetClusterOfInterest)

###Loop through the libs of interest
###adding to the table

###Join the names into a vector that will have all of the info
###for the names

##we need to isolate the set of samples representing the Liver
namesVec = vector()
for (i in 1:length(subsetLibsOfInterest))
{
  ###for the grep query
  s1 <- "^X"
  s2 <- subsetLibsOfInterest[i]
  
  toGrep = 0
  toGrep = paste(s1, s2, sep = "")
  print(grep(toGrep,names(bigTable1)))
  
  ##we are going to get the column of interest 
  ##and then we are going to subset and bind this
  ##to the table of interest
  toAddToSubset <- bigTable1[,grep(toGrep,names(bigTable1))]
  
  ##This is the table containing all of the libraries representing set "A"
  subsetClusterOfInterest[, which( colnames(subsetClusterOfInterest)== s2)] = toAddToSubset
  
  ###Put back the column with the names
  subsetClusterOfInterest$names = bigTable1$Label   
  
}  


myInterestGenes = vector()
myInterestGenes = myTableOfInterest[,1]
myInterestGenes = myInterestGenes[!myInterestGenes == "the"]
myInterestGenes = myInterestGenes[!myInterestGenes == "Libraries"] 

###Specify that these genes are enriched in liver
InterestGenesLiverEnriched = myInterestGenes

##now we've brought in the liver genes into a table !!
###Now, subset the liver and the metabs table !!!
#subset liver
liverRTEToJoin = subsetClusterOfInterest[subsetClusterOfInterest$names %in% InterestGenesLiverEnriched,]
liverRTEToJoin$names <- paste("Liver", liverRTEToJoin$names, sep = "_")
rownames(liverRTEToJoin) = liverRTEToJoin$names
liverRTEToJoin$names <- NULL

rownames(totalMetab) = totalMetab$Label
totalMetab$Label = NULL

colnames(totalMetab) = colnames(liverRTEToJoin) 

theCombinedGenesAndMetabs = rbind(totalMetab, liverRTEToJoin)  

#theCombinedGenes = liverRTEToJoin

###should get rid of all of the NA's before even doing the p-values (I think)
theCombinedGenesAndMetabs = theCombinedGenesAndMetabs[is.finite(rowSums(theCombinedGenesAndMetabs)) == TRUE,]

##in case that doesn't work, try this command as well##
theCombinedGenesAndMetabs = theCombinedGenesAndMetabs[rowSums(theCombinedGenesAndMetabs)!=0, ] 


#to use for the p-value
theCombinedGenesAndMetabsBeforeStandardizing = theCombinedGenesAndMetabs
theCombinedGenesAndMetabsBeforeStandardizing = theCombinedGenesAndMetabsBeforeStandardizing[rowSums(theCombinedGenesAndMetabsBeforeStandardizing)!=0, ] 


##save the combined genes and metabs table before standardizing as another variable
#we will want to save this for the pvalues.
theCombinedGenesAndMetabsBeforeStandardizing = theCombinedGenesAndMetabs

##now, do the scaling as agreed upon
for(i in 1:nrow(theCombinedGenesAndMetabs))
{
    #subtract the mean from each entry and then divide by the standard deviation
    theStandardizedVec = (as.numeric(theCombinedGenesAndMetabs[i,]) - mean(as.numeric(theCombinedGenesAndMetabs[i,])))/sd(as.numeric(theCombinedGenesAndMetabs[i,]))
    
    print(i)
    print("is i !!")
    
    print(theCombinedGenesAndMetabs[i,])
    print("is to be standardized !!!")
    
    print(theStandardizedVec)
    print("is the result of standardizing !!!")
    
    ##now, add the standardized vec back
    theCombinedGenesAndMetabs[i,] = theStandardizedVec 
}


theCombinedGenesAndMetabs$Compounds = NULL

##add the vector of the genes to the table, now
theCombinedGenesAndMetabs = theCombinedGenesAndMetabs[rowSums(theCombinedGenesAndMetabs)!= 0, ] 

###each entry in the combined genes and metabs will
#have an entry in the pValues table
#note that pValues table will have all our important statistical information
#compound, pvalue, gini and k-means cluster
pValues <- data.frame(Compound=character(length(rownames(theCombinedGenesAndMetabs))),pvalues=numeric(length(rownames(theCombinedGenesAndMetabs))),gini=numeric(length(rownames(theCombinedGenesAndMetabs))),cluster=numeric(length(rownames(theCombinedGenesAndMetabs))),ctrlmean=numeric(length(rownames(theCombinedGenesAndMetabs))),hsmean=numeric(length(rownames(theCombinedGenesAndMetabs))))
pValues$Compound = rownames(theCombinedGenesAndMetabs)

###should get rid of all of the NA's for doing the p-values (I think)
theCombinedGenesAndMetabsBeforeStandardizing = theCombinedGenesAndMetabsBeforeStandardizing[is.finite(rowSums(theCombinedGenesAndMetabsBeforeStandardizing)) == TRUE,]

for (i in 1:nrow(theCombinedGenesAndMetabsBeforeStandardizing))
{
    #get the gene name
    theGeneName = rownames(theCombinedGenesAndMetabsBeforeStandardizing)[i]
    
    #do the data look up
    theData = theCombinedGenesAndMetabsBeforeStandardizing[rownames(theCombinedGenesAndMetabsBeforeStandardizing) == theGeneName,]
    
    #now, do the t-test
    theResult = t.test(theData[1:5], theData[6:12], var.equal = FALSE)$p.value
    
    ###now, add the result to the table
    pValues[pValues$Compound == theGeneName,which(colnames(pValues) == "pvalues")] = theResult
    
    pValues[pValues$Compound == theGeneName,which(colnames(pValues) == "hsmean")] = mean(as.numeric(theData[1:5]))
    print(mean(as.numeric(theData[1:5])))
    print("is the data for hs")
    
    pValues[pValues$Compound == theGeneName,which(colnames(pValues) == "ctrlmean")] = mean(as.numeric(theData[6:12]))
    print(mean(as.numeric(theData[6:12])))
    print("is the data for ctrl")
    
    print(theResult)
    print("is the result !!!")
    
    print(theData)
    print("is the data !!!")
    print(theResult)
    print("is the result !!!")

}

###!!! Work with Xiaoke !!! ###
###https://rstudio-pubs-static.s3.amazonaws.com/33876_1d7794d9a86647ca90c4f182df93f0e8.html
###now, do the k-means clustering
set.seed(123456789) ## to fix the random starting clusters

##replace all the nan's and na's with zero that would have result while trying to do
##standardization with
theCombinedGenesAndMetabs = as.matrix(theCombinedGenesAndMetabs)
theCombinedGenesAndMetabs[is.nan(theCombinedGenesAndMetabs)] = 0
theCombinedGenesAndMetabs[is.na(theCombinedGenesAndMetabs)] = 0

##drop all of the rows that are simply zeroes
#execute the k means

theResults = kmeans(theCombinedGenesAndMetabs[,c("1322", "1323", "1325", "1326", "1328", "1330", "1331", "1332", "1333", "1334", "1336", "1337")], centers = 3, nstart = 10)

o=order(theResults$cluster)

#now, add the column to hold the names !!
#make sure it's back to a frame !!
theCombinedGenesAndMetabs = as.data.frame(theCombinedGenesAndMetabs)
theCombinedGenesAndMetabs$Compounds = rownames(theCombinedGenesAndMetabs)

###make the list of cluster assignments
#we will sift this table that contains all of the
#compounds and their cluster info by sending
#compounds associated with each cluster to the randomForest
#algorithm.  First, all compounds from cluster one

#then two and then three
tableToSift = data.frame(rownames(theCombinedGenesAndMetabs)[o],theResults$cluster[o])

    #Iterate through the random forest pipeline
    totalMetabII = theCombinedGenesAndMetabsBeforeStandardizing 

    i = 1

    toMakeRandomForest = tableToSift[tableToSift$theResults.cluster.o. == i,1]
    
    
    toMakeRandomForest = as.vector(toMakeRandomForest)
    
    
    theTarget = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
    
    ###now, we are going to make the random forest !!!
    #we will want to specify that we are working with the k means data and which cluster
    #we are working with
    theNameOfAnalysis = "kmeans"
    

    toModelII = theCombinedGenesAndMetabsBeforeStandardizing

    ##get A which will stay as raw data
    toSubset = toMakeRandomForest
  
    #get the vector of classifiers
    Target = theTarget
      

    cluster = i
            
    #include which iteration we are in
    theGoName = theNameOfAnalysis
    
    #this table is going to be of the form $Compounds, pValues = pvalue for the compound and then gini as $gini
    pValues = pValues
       
    #get the denominator
    ##get the experimental condition    
        
    subsetByInput = theCombinedGenesAndMetabs[rownames(theCombinedGenesAndMetabs) %in% toSubset,]

    predictor_names = rownames(subsetByInput)
    subsetByInput$Compounds = NULL

    theNameToPrint = paste(cluster, "_to_go_into_the_kMeans_raw_data_april1st.txt", sep = "_")

    
    #assuming that it is gene names row and libraries columns
    subsetByInput = t(subsetByInput)
    subsetByInput = as.data.frame(subsetByInput)
    subsetByInput["Target"] = Target

    ##note, also set the target variable now as zeroes and ones for whether it is heat stress or not:
    target = subsetByInput[,"Target"]
    target[target==0]="Heat Stress"
    target[target==1]="Control"
    target=as.factor(target)
 
    #predictor_data = t(subsetByInput)
    predictor_data = subsetByInput
    predictor_names = colnames(subsetByInput)
    colnames(predictor_data)=predictor_names

    #Now, we can execute the algorithm
    tmp = as.vector(table(target))
    num_classes = length(tmp)
    min_size = tmp[order(tmp,decreasing=FALSE)[1]]
    sampsizes = rep(min_size,num_classes)
        
    set.seed(1)
    rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)


    theSorted = head(sort(rf_output$importance[,colnames(rf_output$importance) == "MeanDecreaseGini"], decreasing = TRUE), n = 30) 
    theTableRawDataForSchmidt = theCombinedGenesAndMetabsBeforeStandardizing[rownames(theCombinedGenesAndMetabsBeforeStandardizing) %in% names(theSorted),]

            
    ###Get the table
    toModelII = theCombinedGenesAndMetabsBeforeStandardizing

    ##get A which will stay as raw data
    toSubset = toMakeRandomForest
  
    #get the vector of classifiers
    Target = theTarget
            
    #oh, also - get the cluster info as 8 for the function
    cluster = i
            
    #include which iteration we are in
    theGoName = theNameOfAnalysis
    
    #this table is going to be of the form $Compounds, pValues = pvalue for the compound and then gini as $gini
    pValues = pValues
       
    #get the denominator
    ##get the experimental condition    
        

    subsetByInput = theCombinedGenesAndMetabs[rownames(theCombinedGenesAndMetabs) %in% toSubset,]

    predictor_names = rownames(subsetByInput)
    subsetByInput$Compounds = NULL
    
    #pull out the top thirty compounds
    names(theSorted)
    
    res.pcaCors = PCA(cor(t(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],])) 
    
 
 
    res.pcaCors = PCA(cor(t(toReplaceWithNumbers)))

    #here, we store the names that we have replaced by numbers
    #in order to make the keys
    theVecOfActualNames = rownames(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],]) 
 
 
    theVecOfCompounds = vector()
    theVecOfCompounds = c(theVecOfCompounds, rownames(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))], ]))
 
 
    #write.table(dimdesc(res.pcaCors)[1], file = "dimDescPCA1.txt", sep = "\t")
    write.table(dimdesc(res.pcaCors)[1], file = "dimDescPCA_clus1_PC1.txt", sep = "\t")
    write.table(dimdesc(res.pcaCors)[2], file = "dimDescPCA_clus1_PC2.txt", sep = "\t")
    write.table(dimdesc(res.pcaCors)[3], file = "dimDescPCA_clus1_PC3.txt", sep = "\t")

    ### Also make with the names preserved in the PCA ### 
    res.pcaCorsWithNames = PCA(cor(t(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],])))
 
    
    write.table(dimdesc(res.pcaCorsWithNames)[1], file = "dimDescPCAWithNames_clus1_PC1.txt", sep = "\t")
    write.table(dimdesc(res.pcaCorsWithNames)[2], file = "dimDescPCAWithNames_clus1_PC2.txt", sep = "\t")
    write.table(dimdesc(res.pcaCorsWithNames)[3], file = "dimDescPCAWithNames_clus1_PC3.txt", sep = "\t")


    
    #try for i = 2
    i = 2    
    toMakeRandomForest = tableToSift[tableToSift$theResults.cluster.o. == i,1]
    

    toMakeRandomForest = as.vector(toMakeRandomForest)
    
    theTarget = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
    
    ###now, we are going to make the random forest !!!
    #we will want to specify that we are working with the k means data and which cluster
    #we are working with
    theNameOfAnalysis = "kmeans"
    
    print(i)
    print("is i")

    ### pay attention to these things which we are sending, above ###

    


    ###Get the table
    toModelII = theCombinedGenesAndMetabsBeforeStandardizing

    ##get A which will stay as raw data
    toSubset = toMakeRandomForest
  
    #get the vector of classifiers
    Target = theTarget

            
    cluster = i
            
    #include which iteration we are in
    theGoName = theNameOfAnalysis
    
    #this table is going to be of the form $Compounds, pValues = pvalue for the compound and then gini as $gini
    pValues = pValues
    subsetByInput = theCombinedGenesAndMetabs[rownames(theCombinedGenesAndMetabs) %in% toSubset,]

    predictor_names = rownames(subsetByInput)
    subsetByInput$Compounds = NULL

    theNameToPrint = paste(cluster, "_to_go_into_the_kMeans_raw_data_april1st.txt", sep = "_")

    #assuming that it is gene names row and libraries columns
    subsetByInput = t(subsetByInput)
    subsetByInput = as.data.frame(subsetByInput)
    subsetByInput["Target"] = Target

    ##note, also set the target variable now as zeroes and ones for whether it is heat stress or not:
    target = subsetByInput[,"Target"]
    target[target==0]="Heat Stress"
    target[target==1]="Control"
    target=as.factor(target)
 
    predictor_data = subsetByInput
    predictor_names = colnames(subsetByInput)
    colnames(predictor_data)=predictor_names

    #Now, we can execute the algorithm
    tmp = as.vector(table(target))
    num_classes = length(tmp)
    min_size = tmp[order(tmp,decreasing=FALSE)[1]]
    sampsizes = rep(min_size,num_classes)
        
    set.seed(1)
    rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)


    theSorted = head(sort(rf_output$importance[,colnames(rf_output$importance) == "MeanDecreaseGini"], decreasing = TRUE), n = 30)
    
    theVecOfCompounds = vector()
    theVecOfCompounds = c(theVecOfCompounds, rownames(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))], ]))    
    toReplaceWithNumbers = toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],]
 
    rownames(toReplaceWithNumbers) <- 1:nrow(toReplaceWithNumbers)
 
    res.pcaCors = PCA(cor(t(toReplaceWithNumbers)))

    #here, we store the names that we have replaced by numbers
    #in order to make the keys
    theVecOfActualNames = rownames(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],]) 
 
  
    write.table(dimdesc(res.pcaCors)[1], file = "dimDescPCA_clus2_PC1.txt", sep = "\t")
    write.table(dimdesc(res.pcaCors)[2], file = "dimDescPCA_clus2_PC2.txt", sep = "\t")
    write.table(dimdesc(res.pcaCors)[3], file = "dimDescPCA_clus2_PC3.txt", sep = "\t")
    
    ### Also make with the names preserved in the PCA ### 
    res.pcaCorsWithNames = PCA(cor(t(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],])))
 
    
    write.table(dimdesc(res.pcaCorsWithNames)[1], file = "dimDescPCAWithNames_clus2_PC1.txt", sep = "\t")
    write.table(dimdesc(res.pcaCorsWithNames)[2], file = "dimDescPCAWithNames_clus2_PC2.txt", sep = "\t")
    write.table(dimdesc(res.pcaCorsWithNames)[3], file = "dimDescPCAWithNames_clus2_PC3.txt", sep = "\t")
    
    
    
    #try for i = 3
    i = 3
    
    toMakeRandomForest = tableToSift[tableToSift$theResults.cluster.o. == i,1]
    

  
    toMakeRandomForest = as.vector(toMakeRandomForest)
    

    theTarget = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
    
    toModelII = theCombinedGenesAndMetabsBeforeStandardizing

    ##get A which will stay as raw data
    toSubset = toMakeRandomForest
  
    #get the vector of classifiers
    Target = theTarget

    #oh, also - get the cluster info as 8
    cluster = i
            
    #include which iteration we are in
    theGoName = theNameOfAnalysis
    
    #this table is going to be of the form $Compounds, pValues = pvalue for the compound and then gini as $gini
    pValues = pValues
    subsetByInput = theCombinedGenesAndMetabs[rownames(theCombinedGenesAndMetabs) %in% toSubset,]

    predictor_names = rownames(subsetByInput)
    subsetByInput$Compounds = NULL

    theNameToPrint = paste(cluster, "_to_go_into_the_kMeans_raw_data_april1st.txt", sep = "_")
    #write.table(subsetByInput, theNameToPrint, sep = "\t")
    
    #assuming that it is gene names row and libraries columns
    subsetByInput = t(subsetByInput)
    subsetByInput = as.data.frame(subsetByInput)
    subsetByInput["Target"] = Target

    ##note, also set the target variable now as zeroes and ones for whether it is heat stress or not:
    target = subsetByInput[,"Target"]
    target[target==0]="Heat Stress"
    target[target==1]="Control"
    target=as.factor(target)
 
    #predictor_data = t(subsetByInput)
    predictor_data = subsetByInput
    predictor_names = colnames(subsetByInput)
    colnames(predictor_data)=predictor_names

    #Now, we can execute the algorithm
    tmp = as.vector(table(target))
    num_classes = length(tmp)
    min_size = tmp[order(tmp,decreasing=FALSE)[1]]
    sampsizes = rep(min_size,num_classes)
        
    set.seed(1)
    rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)


    theSorted = head(sort(rf_output$importance[,colnames(rf_output$importance) == "MeanDecreaseGini"], decreasing = TRUE), n = 30)
    theVecOfCompounds = c(theVecOfCompounds, rownames(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))], ]))

    
    toReplaceWithNumbers = toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],]
 
    rownames(toReplaceWithNumbers) <- 1:nrow(toReplaceWithNumbers)
 
    res.pcaCors = PCA(cor(t(toReplaceWithNumbers)))


    #note that we take out the names for the readability
    write.table(dimdesc(res.pcaCors)[1], file = "dimDescPCA_clus3_PC1.txt", sep = "\t")
    write.table(dimdesc(res.pcaCors)[2], file = "dimDescPCA_clus3_PC2.txt", sep = "\t")
    write.table(dimdesc(res.pcaCors)[3], file = "dimDescPCA_clus3_PC3.txt", sep = "\t")

    #here, we store the names that we have replaced by numbers
    #in order to make the keys
    theVecOfActualNames = rownames(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],]) 
 
    ### Also make with the names preserved in the PCA ### 
    res.pcaCorsWithNames = PCA(cor(t(toModelII[rownames(toModelII) %in% names(theSorted)[-(which(names(theSorted) == "Target"))],])))
 
    
    write.table(dimdesc(res.pcaCorsWithNames)[1], file = "dimDescPCAWithNames_clus3_PC1.txt", sep = "\t")
    write.table(dimdesc(res.pcaCorsWithNames)[2], file = "dimDescPCAWithNames_clus3_PC2.txt", sep = "\t")
    write.table(dimdesc(res.pcaCorsWithNames)[3], file = "dimDescPCAWithNames_clus3_PC3.txt", sep = "\t")
 

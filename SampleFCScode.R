#This script shows how our FCS-BRM MI procedure works to generate m imputations of the incomplete dataset
#We use a simulated sample of the same form as our original dataset for this implementation
#The BRM used to impute includes geocoded, surname probabilities, family race indicators and certain interaction terms
library(pROC)
library(mix)
library(BayesLogit)

#Read in data
dataset<-read.csv("sampledata_git.csv")

geocoded<-c("WHITE_CY", "BLACK_CY", "AMERIND_CY", "HISPPOP_CY", "ASIANPACI_CY", "OTHER_CY")
surname<-c("prop_white", "prop_black", "prop_amerind", "prop_hisp", "prop_asianpa", "prop_other")
family<-c("white_fam", "black_fam", "amerind_fam", "hisp_fam", "asianpac_fam")
predictors<-c(geocoded, surname, family)

#input is dataset with race and set of predictors desired for BRM
#this function assumes geocoded and surname probabilities are missing as well
#dataset should have missing values
FCS_BRM<-function(data, predictors, interact=TRUE, burn_in=25, m=25){

    
#Selecting predictors of interest + Race
data<-data[ ,c("Race", predictors)]  
  
###############STEP: Transform geocoded and surname probabilities -> logit
a<-which(names(data)=="WHITE_CY") 
b<-which(names(data)=="prop_white") 

  
  #Logit transformation for probabilities
  for(j in 1:dim(data)[1]){
    for(i in (a+1):(a+5)){
      data[j, i]<-log((data[j, i])/data[j, a])
      data[j, i+6]<-log((data[j, i+6])/data[j, b])
    }
    
    #Fixing the NaN - occurs when both probabilities are 0
    for(k in 2:13){
      if(is.nan(data[j, k])){data[j, k]<-0.0015}
    }
    
  }
  
  #-Inf occurs for log(0/something) - category of interest is very low
  #Inf occurs for log(something/0) - reference category is very low
  data[data==-Inf]<-log(0.025/0.975)
  data[data==Inf]<-log(0.975/0.025)
  
  #Remove white reference categories for geocoded and surname probabilities
  data<-data[ , -c(a, b)]
  
###############STEP 0: Create starting dataset 
  
  #Characterizing missing patterns
  #A is missing race, geocode
  #B is missing race, surname
  #C is missing surname
  #D is missing geocode
  #E is missing race
  #F is complete
  miss_pat<-rep(NA, dim(data)[1])
  for(i in 1:dim(data)[1]){
    if(is.na(data$prop_black)[i]){miss_pat[i]<-"C"}
    if(is.na(data$BLACK_CY)[i]){miss_pat[i]<-"D"}
    if(is.na(data$Race)[i]){miss_pat[i]<-"E"}
    if(is.na(data$Race)[i]&is.na(data$BLACK_CY)[i]){miss_pat[i]<-"A"}
    if(is.na(data$Race)[i]&is.na(data$prop_black)[i]){miss_pat[i]<-"B"}
    if(complete.cases(data)[i]){miss_pat[i]<-"F"}
  }
  
  #Redefining geocoded, surname without white
  geocoded<-c("BLACK_CY", "AMERIND_CY", "HISPPOP_CY", "ASIANPACI_CY", "OTHER_CY")
  surname<-c("prop_black", "prop_amerind", "prop_hisp", "prop_asianpa", "prop_other")
  
  
  #For mix all the categorical variables must have categories 1, 2, .... and put in the first columns of the data
  data<-data[ ,c("Race", family, geocoded, surname)]
  #Form for indicators needed by Mix -> ex. binary indicator must be 1, 2 instead of 0, 1
  data[,family]<-data[,family]+1

  s<-prelim.mix(data, 6)
  thetahat<-em.mix(s)
  rngseed(12345678)
  start<-imp.mix(s, thetahat, data)
  data<-as.data.frame(start)
  data[,family]<-data[,family]-1

#### data is now the starting dataset

##################### The main loops begins here #####################
  
  #List to save datasets at the end
  dataset_save<-list()
  
  for(g in 1:(burn_in+m)){
    ####### STEP 1: Impute Race #######
    #set race equal to NA
    for(i in 1:dim(data)[1]){
      if(miss_pat[i]=="A"|miss_pat[i]=="B"|miss_pat[i]=="E"){data$Race[i]<-NA}
    }
    
    ################ Build BRM on data_train, use estimates to predict race for data_miss
    data_train<-data[!is.na(data$Race), ]
    data_miss<-data[is.na(data$Race), ]
    
    N<-nrow(data_train)
    P<-ncol(data_train)
    data_train$Race<-factor(data_train$Race)
    J<-nlevels(data_train$Race)
    if(interact==FALSE){X<-model.matrix(Race~., data=data_train)}
    # This is the modification that includes ALL two-way interactions
    if(interact==TRUE){X<-model.matrix(Race~.+BLACK_CY*black_fam+AMERIND_CY*amerind_fam+HISPPOP_CY*hisp_fam+ASIANPACI_CY*asianpac_fam+prop_black*black_fam+prop_amerind*amerind_fam+prop_hisp*hisp_fam+prop_asianpa*asianpac_fam+BLACK_CY*prop_black+AMERIND_CY*prop_amerind+HISPPOP_CY*prop_hisp+ASIANPACI_CY*prop_asianpa, data=data_train)
    }
    #y is now a matrix
    y.all<-model.matrix(~Race-1, data=data_train)
    #Excludes Race_Code1, which is white
    y<-y.all[, -1]
    
    out<-mlogit(y, X, samp=1, burn=10000, P.0=array(0.0001, dim=c(ncol(X), ncol(X), ncol(y))))
    beta_cat2<-as.vector(out$beta[ , , 1])
    beta_cat3<-as.vector(out$beta[ , , 2])
    beta_cat4<-as.vector(out$beta[ , , 3])
    beta_cat5<-as.vector(out$beta[ , , 4])
    
	#Cant't have outcome be all NAs to use model.matrix
	data_miss[,1]<-1
    if(interact==FALSE){X_miss<-model.matrix(Race~., data=data_miss)}
    # This is the modification that includes ALL two-way interactions
    if(interact==TRUE){X_miss<-model.matrix(Race~.+BLACK_CY*black_fam+AMERIND_CY*amerind_fam+HISPPOP_CY*hisp_fam+ASIANPACI_CY*asianpac_fam+prop_black*black_fam+prop_amerind*amerind_fam+prop_hisp*hisp_fam+prop_asianpa*asianpac_fam+BLACK_CY*prop_black+AMERIND_CY*prop_amerind+HISPPOP_CY*prop_hisp+ASIANPACI_CY*prop_asianpa, data=data_miss)
    }
    
    prediction<-c()
    #Going though every single person in this subset of dataset - predict race for each
    for(i in 1:dim(X_miss)[1]){
      cat2<-exp(sum(beta_cat2*X_miss[i, ]))
      cat3<-exp(sum(beta_cat3*X_miss[i, ]))
      cat4<-exp(sum(beta_cat4*X_miss[i, ]))
      cat5<-exp(sum(beta_cat5*X_miss[i, ]))
      cat1<-exp(sum(rep(0, dim(X_miss)[2])*X_miss[i, ]))
      #Computing p_ij's for each race category
      sumc<-cat1+cat2+cat3+cat4+cat5
      PosteriorProbs<-c(cat1, cat2, cat3, cat4, cat5)/sumc
      #generating a predicted race
      prediction[i]<-which(as.vector(rmultinom(1, 1, prob=PosteriorProbs))==1)
    }
    
    #Update dataset with new race predictions
    data$Race[is.na(data$Race)]<-prediction
    
    ####### STEP 2, 3: Impute Missing Geocoded, Surname Probabilities #######
    #Form needed for mix 
    data[,family]<-data[,family]+1
    
    for(i in 1:dim(data)[1]){
      if(miss_pat[i]=="A"|miss_pat[i]=="D"){data[i, geocoded]<-NA}
    }
    
    s<-prelim.mix(data, 6)
    thetahat<-em.mix(s)
    rngseed(123456+g)
    start<-imp.mix(s, thetahat, data)
    data<-as.data.frame(start)
    

   for(i in 1:dim(data)[1]){
      if(miss_pat[i]=="B"|miss_pat[i]=="C"){data[i, surname]<-NA}
    }
    s<-prelim.mix(data, 6)
    thetahat<-em.mix(s)
    rngseed(123456+g)
    start<-imp.mix(s, thetahat, data)
    data<-as.data.frame(start)
    
    
    #Returning to form needed for BRM
    data[,family]<-data[,family]-1
    
    #dataset is complete again, now we restart the process - set race equal to NA
    #Saving datasets onto a list after burn-in period has been completed
    if(g>burn_in){
      dataset_save[[g-burn_in]]<-data
    }
    print(g)
  }
   
 return(dataset_save) 
  }

imputed_datasets<-FCS_BRM(dataset, predictors, interact=TRUE)


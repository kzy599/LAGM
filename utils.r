library(AlphaSimR)
library(data.table)
library(optiSel)
library(lagm)
if (!require(jsonlite)) {
  install.packages("jsonlite")
  library(jsonlite)
}
source("ocs.R")
makeped = function(z){
  z[!(z == 0 | z == 1 | z == 2)] <- -9
  
  z[z == 2] <- 22
  z[z == 1] <- 12
  z[z == 0] <- 11

  dt <- as.data.table(z, keep.rownames = "ID")
  return(dt)
  }

   #make .map file required by plink 
  makemap = function(...){
  mapck = getSnpMap()
  mapck$id = mapck$chr
  mapck$chr = rownames(mapck)
  mapck$site = mapck$pos
  mapck$pos = rep(0,nrow(mapck))
  return(mapck)
  }

clnm = function(r,k){
  return(colnm = paste("r",r,"k",k,sep = "_"))
}

calmaf = function(pop){
  x_geno = pullSnpGeno(pop)
  maf = apply(x_geno,2,function(x){
    maf = sum(x)/(2*length(x))
    if(maf>0.5) maf = 1 - maf
    return(maf)
  })
  return(mean(maf))
}

makedt = function(phe_dt,r,k,effect){
  cn = clnm(r,k)
  trainID = phe_dt$ID[phe_dt[[cn]]=="train"]
  testID = phe_dt$ID[phe_dt[[cn]]!="train"]

  ID = c(trainID, testID)
  ef_dt= data.table(ID = ID, effect = effect)
  
  
  # Remove the first column (ID column) and return only the SNP information
  return(ef_dt)
}
read_geno <- function(dir_geno, phe_dt,r,k) {
  cn = clnm(r,k)
  trainID = phe_dt$ID[phe_dt[[cn]]=="train"]
  testID = phe_dt$ID[phe_dt[[cn]]!="train"]
  # Load the genotype data
  geno <- fread(dir_geno, sep = ",", header = TRUE)
  
  # Reorder rows based on trainID and testID
  geno_ordered <- geno[match(c(trainID, testID), geno[[1]]), ]
  
  # Remove the first column (ID column) and return only the SNP information
  return(geno_ordered[, -1, with = FALSE])
}


get_phe = function(phe_dt,trainID,r,k){
  cn = clnm(r,k)

  return(phe_dt$phe1[phe_dt[[cn]]=="train"])
}

makeGDE = function(x,type = "D",inv = FALSE){
  
  if(!is.matrix(x)) x <- as.matrix(x)
  
  if(type=="E"){

      M = x - 1
      
      E <- 0.5 * ((M %*% t(M)) * (M %*% t(M))) - 0.5 * ((M * M) %*% t(M * M))
      
      E <- E / (sum(diag(E)) / nrow(E))
      
      A = diag(1,nrow(E))
      
      E = E * 0.99 + A * 0.01

      
      if(inv) inverse <- solve(E) else inverse <- E

  }else if(type == "D"){

      P = apply(x,2,function(col){
      
        pi = sum(col)/(2*length(col))
      
        if(pi>0.5) pi = 1-pi
      
        return(pi)
      
      })

      W = apply(x,2,function(col){
          
          pi = sum(col)/(2*length(col))
          
          if(pi>0.5){
          
            pi = 1-pi
          
            AA_bool = which(col == 0)
          
            aa_bool = which(col == 2)
          
            col[AA_bool] = 2
          
            col[aa_bool] = 0
          
          }
          
          aa = -2 * (pi^2)
          
          Aa = 2 * pi * (1-pi)
          
          AA = -2 * ((1-pi)^2)
          
          aa_bool = which(col == 0)
          
          Aa_bool = which(col == 1)
          
          AA_bool = which(col == 2)
          
          col[AA_bool] = AA
          
          col[aa_bool] = aa
          
          col[Aa_bool] = Aa
          
          return(col)
      
      })

      D <- (W %*% t(W)) / sum((2 * P * (1 - P))^2)
      
      A = diag(1,nrow(D))
      
      D = D * 0.99 + A * 0.01
      
      if(inv) inverse <- solve(D) else inverse <- D

  }else{
      P = apply(x,2,function(col){
        
        pi = sum(col)/(2*length(col))
        
        if(pi>0.5) pi = 1-pi
        
        return(pi)
     
      })
      Z = apply(x,2,function(col){
        pi = sum(col)/(2*length(col))

        if(pi>0.5){
        
          pi = 1-pi
        
          AA_bool = which(col == 0)
        
          aa_bool = which(col == 2)
        
          col[AA_bool] = 2
        
          col[aa_bool] = 0
        
        }
        
        col = col - 2*pi
        
        return(col)
      })
    
      G = (Z %*% t(Z)) / sum((2 * P * (1 - P)))
      
      A = diag(1,nrow(G))

      G = G * 0.99 + A * 0.01
      
      if(inv) inverse <- solve(G) else inverse <- G

  }
  return(inverse)
}

calEBV = function(x,y,h2,type = "D"){
  train_len = length(y)

  whole_len = nrow(x)
  
  Z = matrix(0, nrow = train_len, ncol = whole_len)
  
  # Set the diagonal elements to 1
  diag_indices = 1:train_len
  
  Z[cbind(diag_indices, diag_indices)] = 1
  
  mat = makeGDE(x = x,type = type)
  
  if(type == "D"){
  
    lamb = (1 - h2*0.1)/(h2*0.1)
  
  }else if(type=="E"){
  
    lamb = (1 - h2*0.1)/(h2*0.1)
  
  }else{
  
    lamb = (1 - h2)/h2
  
  }
  #glamb = (1 - h2)/h2
  #dlamb = (1 - h2*0.1)/(h2*0.1)
  #elamb = (1 - h2*0.1)/(h2*0.1)

  y = y - mean(y)
  
  Z_t =  t(Z)
  
  effect  = solve(Z_t %*% Z + mat * lamb) %*% Z_t %*% y

  return(effect)
}



calac = function(phe_test,predictions_dt,model_name,a_dt=NULL){

  if(model_name %in% c("cnn_mt","fc_mt","mul_mt","cnn_soft","fc_soft","mul_soft","cnn_linear","fc_linear","mul_linear")){
  predictions = "Prediction_t"
}else if(model_name %in% c("cnn","fc","mul")){
  predictions = "Prediction"
}else if(model_name %in% c("gblup_mt","gblup")){
  predictions = "Gmat_b.GA"
}else{
  predictions = "Amat_b.PA"
}

if(!is.null(a_dt)){
  ebv1 = predictions_dt[[predictions]][match(phe_test$ID,predictions_dt$ID)] + a_dt$Gmat_b.GA[match(phe_test$ID,a_dt$ID)]
  model_name = paste(model_name,"_a",sep = "")
}else{
  ebv1 = predictions_dt[[predictions]][match(phe_test$ID,predictions_dt$ID)]
}

  GV1 = phe_test$gv1


  acc  = cor(GV1,ebv1)
  output = data.table(mname = model_name, ac = acc)
  return(output)

}



calWac = function(phe_test,predictions_dt,a_dt = NULL,model_name){

if(model_name %in% c("cnn_mt","fc_mt","mul_mt","cnn_soft","fc_soft","mul_soft","cnn_linear","fc_linear","mul_linear")){
  predictions = "Prediction_t"
}else if(model_name %in% c("cnn","fc","mul")){
  predictions = "Prediction"
}else{
  predictions = "Gmat_b.GA"
}

fac = data.table()
for( f in unique(phe_test$FamilyID)){
  phe_test_F = phe_test[FamilyID==f,]
  phe1_f = phe_test_F$phe1
  GV1_f = phe_test_F$gv1

  if(!is.null(a_dt)){

    ebv1_f = predictions_dt[[predictions]][match(phe_test_F$ID,predictions_dt$ID)] +a_dt$Gmat_b.GA[match(phe_test_F$ID,a_dt$ID)]

  }else{
    ebv1_f = predictions_dt[[predictions]][match(phe_test_F$ID,predictions_dt$ID)]
  }
  acg = cor(GV1_f,ebv1_f)
  acp = cor(phe1_f,ebv1_f)

  meanp = mean(phe1_f)
  meang = mean(GV1_f)
  meanebv = mean(ebv1_f)

  fac = rbind(fac, data.table(FamilyID = f, acg = acg, acp = acp,phe = meanp,gv = meang, ebv = meanebv ))
}

if(!is.null(a_dt)) model_name = paste(model_name,"_a",sep = "")

 output = data.table(mname = model_name,wac =mean(fac$acg) ,minac = min(fac$acg),maxac = max(fac$acg),sdac = sd(fac$acg),famacc = cor(fac$gv,fac$ebv))

 return(output)

}


calcvmean = function(dt){
  result <- dt[, lapply(.SD, mean), by = mname]
  return(result)
}


extractMeanSd = function(dt){
  result_mean <- dt[, lapply(.SD, mean), by = mname]

  result_sd <- dt[, lapply(.SD, sd), by = mname]


  bool = (colnames(result_sd)!="mname")

  colnames(result_sd)[bool] = paste(colnames(result_sd)[bool],"_sd",sep = "")

  result <- merge(result_mean, result_sd, by = "mname")
  
  return(result)
}


calinbnej = function(pop){

  nej <- function(snp_012_dt){
  calformula <- function(x){
    y <- sum(x!=1)/length(x)
    return(y)
  }
  inb <- apply(snp_012_dt, 1, calformula)
  return(inb)
}
  pop_select <- selectWithinFam(pop = pop,nInd = 20,use = "rand",simParam = SP)
  snp_012_dt = pullSnpGeno(pop_select)
  inbreeding_nej = mean(nej(snp_012_dt = snp_012_dt))
  return(inbreeding_nej)
}


calinblh= function(pop,snpfre){

  lh <- function(snp_012_dt,snpfre){
  calformula <- function(x){
    y <- (sum(x!=1) - sum(1 - 2 * snpfre * (1 - snpfre))) / (length(x) - sum(1 - 2 * snpfre * (1 - snpfre)))
    return(y)
  }
  inb <- apply(snp_012_dt, 1, calformula)
  return(inb)
}
  pop_select <- selectWithinFam(pop = pop,nInd = 20,use = "rand",simParam = SP)
  snp_012_dt = pullSnpGeno(pop_select)
  inbreeding_lh = mean(lh(snp_012_dt = snp_012_dt, snpfre = snpfre))
  return(inbreeding_lh)
}


calinbya2 <- function(snp_012_dt,snpfre){
  calformula <- function(x, snpfre){
    y <- mean((x^2 - ((1 + 2 * snpfre) * x) + 2 * (snpfre^2)) / (2 * snpfre * (1 - snpfre)), na.rm = TRUE)
    return(y)
  }
  inb <- apply(snp_012_dt, 1, calformula,snpfre = snpfre)
  return(inb)
}

calinb = function(pop,snpfre){
  pop_select <- selectWithinFam(pop = pop,nInd = 20,use = "rand",simParam = SP)
  snp_012_dt = pullSnpGeno(pop_select)
  inbreeding_ya2 = mean(calinbya2(snp_012_dt = snp_012_dt, snpfre = snpfre))
  return(inbreeding_ya2)
}

makejson = function(candidate,females,males,gen,gen_pre = NULL,He = NULL,A_matrix = NULL){

  cand_geno = pullSnpGeno(candidate)

  kinship_matrix <- makeGDE(x = cand_geno, type = "G", inv = FALSE)

  individual_ids <- candidate@id

  ebv_vector <- candidate@ebv[,1]

  # Female and Male IDs
  female_ids <- females@id

  male_ids <- males@id

  # Create a list containing all the data
  input_data <- list(
    individual_ids = individual_ids,
    ebv_vector = ebv_vector,
    kinship_matrix = kinship_matrix,
    female_ids = female_ids,
    male_ids = male_ids,
    geno_matrix = cand_geno,
    gen = gen,
    He = He,
    A_matrix = A_matrix,
    gen_pre = gen_pre
  )

  return(input_data)

}

calNeped = function(ped,keep_pop){
  ped_calparameters = rbind(ped,createPed(keep_pop))
  keep = keep_pop@id
  pedig <- prePed(ped_calparameters)
  pKin   <- pedIBD(pedig, keep.only = keep)
  Summary <- summary(pedig)
  id     <- keep
  x      <- Summary[Indiv %in% id]$equiGen
  N      <- length(x)
  n      <- (matrix(x, N, N, byrow = TRUE) + matrix(x, N, N, byrow = FALSE)) / 2
  deltaC <- 1 - (1 - pKin[id, id]) ^ (1 / n)
  Ne   <- 1 / (2 * mean(deltaC))
  return(Ne)
  
}

calInbped = function(ped,keep_pop){

  ped_calparameters = rbind(ped,createPed(keep_pop))

  keep = keep_pop@id
  
  Pedig <- prePed(ped_calparameters, keep=keep)

  Res   <- pedInbreeding(Pedig)
  
  inbreeding <- mean(Res$Inbr[Res$Indiv %in% keep])
  
  return(inbreeding)
  
}


createCandidate = function(pop){
  females <- selectWithinFam(pop,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  return(candidate)
}

makeEBV <- function(pop){
  needfam = unique(paste(pop@mother,pop@father,sep = "_"))

  pop_list = lapply(needfam,function(x){

          pop_v = pop[paste(pop@mother,pop@father,sep = "_") == x ]

          fgv = mean(pop_v@gv[,1])

          webv = pop_v@gv[,2] - mean(pop_v@gv[,2])

          Ebv = fgv + webv

          pop_v@ebv = matrix(Ebv,ncol = 1)

          pop_v@gv[,2] = Ebv - webv
        
        return(pop_v)
    })

    pop = mergePops(pop_list)

    return(pop)
}

createPed = function(pop){
  ped = data.table(id = pop@id, sire = pop@father,dam = pop@mother )
  colnames(ped) = c("id","sire","dam")
  return(ped)
}


calHe = function(pop){
  snp = pullSnpGeno(pop)
  He_loci = apply(snp,2,function(x){
    p <- sum(x, na.rm = TRUE)/(2*length(x))
    return(2*p*(1-p))
  })
  return(mean(He_loci))
}

calHo = function(pop){
  snp = pullSnpGeno(pop)
  Ho_loci = apply(snp,2,function(x){
    return(sum(x==1)/length(x))
  })
  return(mean(Ho_loci))
}

calAmatrix = function(ped,keep){
  
  Pedig <- prePed(ped, keep=keep)

  Amatrix <- makeA(ped)
  
  Amatrix <- as.matrix(Amatrix)

  tempBool_row = match(keep,rownames(Amatrix))
  tempBool_col = match(keep,colnames(Amatrix))
  
  AA<-Amatrix[tempBool_row,tempBool_col]
  
  return(AA)
}



run_hiblup = function(phename,trait_pos,addG="",domG="", out_prefix="G_hib"){

  #domG = ",Gmat_b.GD"
  #domA = ",Amat_b.PD"
  
  if(length(trait_pos)>1){

    trait_pos_str = paste(trait_pos, collapse = " ")
    model_for_trait = "--multi-trait"
    

  }else{
    trait_pos_str = as.character(trait_pos)
    model_for_trait = "--single-trait"
  }

      system(paste("hiblup",
        model_for_trait,
          "--pheno",
          phename,
            paste("--pheno-pos", trait_pos_str, sep = " ") ,
            paste("--xrm", paste(addG, domG, sep=","), sep=" "),
            "--vc-method AI",
            "--ai-maxit 30",
            "--threads 32",
            "--ignore-cove",
            paste("--out", out_prefix, sep = " "),sep = " "))

results_dt = fread(paste(out_prefix, ".rand", sep=""),sep = "\t")
  
dt_G = results_dt[, .(ID, Prediction = get(addG))]

if(domG != ""){
 dt_D = results_dt[, .(ID, Prediction = get(domG))]
}else{
  dt_D = results_dt[, .(ID, Prediction = NA)]
}
  
if(addG %flike%".GA"){
blup_models = list(
  list(dt=dt_G, model="gblup", comp="add"),
  list(dt=dt_D, model="gblup", comp="dom")
)
}else if(addG %flike%".PA"){
blup_models = list(
  list(dt=dt_G, model="pblup", comp="add"),
  list(dt=dt_D, model="pblup", comp="dom")
)
}
return(blup_models)
}


estimate_ebv = function(pop){
females_ref <- selectWithinFam(pop,nInd=nProgenyPerCross/4,use = "rand",sex = "F",trait = 1)
males_ref <- selectWithinFam(pop,nInd=nProgenyPerCross/4,use = "rand",sex = "M",trait = 1)
pop_ref = c(females_ref,males_ref)
genodt = pullSnpGeno(pop)
genodtput = makeped(z = genodt)
fwrite(genodtput,file="hib.ped",col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t",nThread = 8)
system("plink --file hib --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib")
system("hiblup --make-xrm --threads 32 --bfile hib --add --out Gmat_b")

phe = data.table(ID = pop@id,
                            phe1 = pop@pheno[,1],sire = pop@father,dam = pop@mother,
                            FamilyID = paste(pop@father,"_",pop@mother,sep = ""))
phe[!ID%in%pop_ref@id,phe1:=NA]
fwrite(phe,file = "phe.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
gblup_models = run_hiblup(phename = "phe.csv", trait_pos =which(colnames(phe)=="phe1"),addG="Gmat_b.GA",domG = "",out_prefix="G_hib")
ebv = copy(gblup_models[[1]]$dt)
pop@ebv = matrix(ebv$Prediction[match(pop@id,ebv$ID)],ncol=1)
return(pop[!pop@id %in% pop_ref@id])
}



estimate_ebv_pblup = function(pop,ped){

ped_temp = data.table(id = pop@id,sire = pop@father,dam = pop@mother)
colnames(ped_temp) = c("id","sire","dam")
ped_temp = rbind(ped,ped_temp)
fwrite(ped_temp,file = "ped_temp.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
system('hiblup --make-xrm --threads 32 --pedigree ped_temp.csv --add --out Amat_b')

phe = data.table(ID = pop@id,
                            phe1 = pop@pheno[,1],sire = pop@father,dam = pop@mother,
                            FamilyID = paste(pop@father,"_",pop@mother,sep = ""))
fwrite(phe,file = "phe_temp.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
pblup_models = run_hiblup(phename = "phe_temp.csv", trait_pos =which(colnames(phe)=="phe1"),addG="Amat_b.PA",domG = "",out_prefix="A_hib")
ebv = copy(pblup_models[[1]]$dt)
pop@ebv = matrix(ebv$Prediction[match(pop@id,ebv$ID)],ncol=1)
return(pop)
}



run_ocs_aqua = function(candidate,targetDegree,minInb,nCrosses,nProgenyPerCross,nDam,nSire,Fix_fmRatio=TRUE){
    if(Fix_fmRatio){
    pop <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        # maxMaleContribution = 4,
        # maxFemaleContribution = 2,
        targetDegree = targetDegree,
        use = "ebv",
        minInbreedingMating = minInb
  )
    }else{
    pop <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        # nFemalesMax = nDam,
        # nMalesMax = nSire,
        # equalizeFemaleContributions =TRUE,
        # equalizeMaleContributions = TRUE,
        maxMaleContribution = 4,
        maxFemaleContribution = 2,
        targetDegree = targetDegree,
        use = "ebv",
        minInbreedingMating = minInb
  )
    }

    return(pop)
}
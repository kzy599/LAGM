set.seed(42)
randomseed = sample(.Machine$integer.max,10)#.Machine$integer.max is 2147483647
for(r in 1:10){
rm(list = setdiff(ls(), c("r","randomseed")))
gc()
set.seed(randomseed[r])
nG = 5
dir = "/home/kangziyi/poster/BaseMapHaplo/"
mapFile = "mergedMap.csv"
haploFile = "mergedHaplo.csv"

source("utils.r")

num <- data.frame(chr = c("NC_088853.1",
                          "NC_088854.1",
                          "NC_088855.1",
                          "NC_088856.1",
                          "NC_088857.1",
                          "NC_088858.1",
                          "NC_088859.1",
                          "NC_088860.1",
                          "NC_088861.1",
                          "NC_088862.1"),
                  len = c(76070991,
                          61469542,
                          61039741,
                          57946171,
                          57274926,
                          56905015,
                          53672946,
                          51133819,
                          50364239,
                          37310742)
)

genMap = fread(paste(dir,mapFile,sep=""),sep = ",")

Haplo = fread(paste(dir,haploFile,sep=""),sep = ",",header = TRUE)

colnames(Haplo) = paste("Site_",colnames(Haplo),sep = "")

genMap[,site:=colnames(Haplo)]

map_list = list()
haplo_list = list()
qtl_pos = list()
snp_pos = list()

for(i in 1:nrow(num)){
    sites = genMap[chr == num$chr[i],site]

    Qtl_sites = genMap[chr == num$chr[i]&QTL==TRUE,site]

    Snp_sites = genMap[chr == num$chr[i]&QTL!=TRUE,site]

    map_list[[i]] = genMap[chr == num$chr[i],pos]
    map_list[[i]] = map_list[[i]]/num$len[i]

    tempIndex = match(sites,colnames(Haplo))

    haplo_list[[i]] = Haplo[,.SD,.SDcols = colnames(Haplo)[tempIndex] ]

    qtl_pos[[i]] = which(sites %in% Qtl_sites)

    snp_pos[[i]] = which(sites %in% Snp_sites)

}

founderPop <- newMapPop(genMap=map_list, haplotypes=haplo_list)
SP <- SimParam$new(founderPop)

SP$invalidQtl <- snp_pos
SP$invalidSnp <- qtl_pos

SP$addTraitA(nQtlPerChr = sapply(qtl_pos, length),
              mean=c(0,0),var = c(1,1),
              corA = matrix(c(1,0.3,0.3,1),nrow = 2))

SP$setVarE(h2 = c(0.3,0.3)) # 0.17 0.25 in the range of heritability for growth, meat yield, survival, etc
SP$addSnpChip(nSnpPerChr = sapply(snp_pos, length)) # all non-QTL SNPs saved from simulation
SP$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex

pop_founder = newPop(founderPop, simParam=SP)


nDam = 50
nSire = 25
nCrosses = 50
nProgenyPerCross = 100
nProgeny = nProgenyPerCross
pop <- selectCross(pop_founder,
                   nFemale = nDam,nMale = nSire,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)
ped_ne = data.table(id = c(unique(pop@mother,pop@father)),sire = NA,dam = NA)
colnames(ped_ne) = c("id","sire","dam")

# ped = data.table()
# phe = data.table()

for(g in c(1:nG)){
  
  # ped = rbind(ped,data.table(ID = pop@id,father = pop@father,mother = pop@mother))
  # phe = rbind(phe, data.table(ID = pop@id,env = rep(1,pop@nInd),
  #                             gv1 = pop@gv[,1],gv2 = pop@gv[,2],
  #                             phe1 = pop@pheno[,1],phe2 = pop@pheno[,2],
  #                             FamilyID = paste(pop@father,"_",pop@mother,sep = "")))
  if(g==1){
    pop <- makeEBV(pop)
    females <- selectWithinFam(pop,nInd=4,use = "ebv",sex = "F",trait = 1)

    males <- selectWithinFam(pop,nInd=2,use = "ebv",sex = "M",trait = 1)
  
    candidate = c(females,males)

    He = calHe(pop)
  }#else{
  #   pop_rl <- makeEBV(pop_rl)

  #   females <- selectWithinFam(pop_rl,nInd=4,use = "ebv",sex = "F",trait = 1)

  #   males <- selectWithinFam(pop_rl,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  #   candidate = c(females,males)
  # }
  # input_data = makejson(candidate = candidate,females = females,males = males,gen = ((nG+1)-g))
  # # Specify the output JSON file path
  # output_json_path <- "input_data.json"

  # # Write the JSON file
  # write_json(input_data, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

  # cat(paste("JSON input file has been generated at", output_json_path, "\n"))

  # system('bash -c "source activate tf-gpu && python optMating.py input_data.json"')

  # mating_plans <- fromJSON("breeding_pairs.json")[[1]]
  # mating_plans<- as.matrix(mating_plans)
  # #pop = selectInd(pop,nInd = 1000,use = "rand")

  # pop_rl <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)


if(g!=1){
  pop_rl2 <- makeEBV(pop_rl2)

  females <- selectWithinFam(pop_rl2,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  He = calHe(pop_rl2)

  ped_rl2 = rbind(ped_rl2,createPed(pop = candidate))

}else{

  ped_rl2 = rbind(ped_ne,createPed(pop = candidate))
}

  if(g == 1){
    input_data2 = makejson(candidate = candidate,females = females,males = males,gen = ((nG+1)-g),He = He)
    # input_data2 <- input_data

  }else{
    input_data2 <- makejson(candidate = candidate,females = females,males = males,gen = ((nG+1)-g),He = He)
  }
  # Specify the output JSON file path
  output_json_path <- "input_data2.json"

  # Write the JSON file
  write_json(input_data2, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

  cat(paste("JSON input file has been generated at", output_json_path, "\n"))

  system('bash -c "source activate tf-gpu && python optMatingP.py input_data2.json"')

  mating_plans <- fromJSON("breeding_pairs2.json")[[1]]
  mating_plans<- as.matrix(mating_plans)
  #pop = selectInd(pop,nInd = 1000,use = "rand")

  pop_rl2 <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)
  keep_rl2 <- selectWithinFam(pop = pop_rl2,nInd = 20,use = "rand",simParam = SP)


  if(g!=1){
    pop_tc <- makeEBV(pop_tc)
    candidate <- createCandidate(pop_tc)
    ped_tc = rbind(ped_tc,createPed(pop = candidate))
  }else{
    ped_tc = rbind(ped_ne,createPed(pop = candidate))
  }
  
  # pop_tc <- selectCross(candidate,
  #                    nFemale = nDam,nMale = nSire,
  #                    nCrosses = nCrosses,nProgeny = nProgenyPerCross,
  #                    use = "ebv",trait = 1,
  #                    simParam = SP)
  pop_tc <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 0,
        use = "ebv"
  )
  keep_tc <- selectWithinFam(pop = pop_tc,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    pop_ocs <- makeEBV(pop_ocs)
    candidate <- createCandidate(pop_ocs)
    ped_ocs = rbind(ped_ocs,createPed(pop = candidate))
  }else{
    ped_ocs = rbind(ped_ne,createPed(pop = candidate))
  }

  # candidate@ebv = matrix(candidate@gv[,2],ncol = 1)

  pop_ocs <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 45,
        use = "ebv"
  )
  keep_ocs <- selectWithinFam(pop = pop_ocs,nInd = 20,use = "rand",simParam = SP)

if(g!=1){
    pop_ran <- makeEBV(pop_ran)
    candidate <- createCandidate(pop_ran)
    ped_ran = rbind(ped_ran,createPed(pop = candidate))
  }else{
   ped_ran = rbind(ped_ne,createPed(pop = candidate))

  }
  pop_ran <- selectCross(candidate,
                     nFemale = nDam,nMale = nSire,
                     nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                     use = "rand",trait = 1,
                     simParam = SP)
  keep_ran <- selectWithinFam(pop = pop_ran,nInd = 20,use = "rand",simParam = SP)



  if(g!=1){
    pop_65 <- makeEBV(pop_65)
    candidate <- createCandidate(pop_65)
    ped_65 = rbind(ped_65,createPed(pop = candidate))
  }else{
    ped_65 = rbind(ped_ne,createPed(pop = candidate))
  }
  pop_65 <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 65,
        use = "ebv"
  )
  keep_65 <- selectWithinFam(pop = pop_65,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    pop_25 <- makeEBV(pop_25)
    candidate <- createCandidate(pop_25)
    ped_25 = rbind(ped_25,createPed(pop = candidate))
  }else{
    ped_25 = rbind(ped_ne,createPed(pop = candidate))
  }
  pop_25 <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 25,
        use = "ebv"
  )
  keep_25 = selectWithinFam(pop = pop_25,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    pop_90 <- makeEBV(pop_90)
    candidate <- createCandidate(pop_90)
    ped_90 = rbind(ped_90,createPed(pop = candidate))
  }else{
    ped_90 = rbind(ped_ne,createPed(pop = candidate))
  }
  # candidate@ebv[,1] = candidate@gv[,2]
  pop_90 <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 90,
        use = "ebv"
  )
  keep_90 = selectWithinFam(pop = pop_90,nInd = 20,use = "rand",simParam = SP)


  if(g!=1){
    pop_rate <- makeEBV(pop_rate)
    candidate <- createCandidate(pop_rate)
    ped_rate = rbind(ped_rate,createPed(pop = candidate))
  }else{
    ped_rate = rbind(ped_ne,createPed(pop = candidate))
  }
  # candidate@ebv[,1] = candidate@gv[,2]
  pop_rate <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        targetCoancestryRate = 0.01,
        use = "ebv"
  )
  keep_rate = selectWithinFam(pop = pop_rate,nInd = 20,use = "rand",simParam = SP)



  if(g!=1){
  pop_rl <- makeEBV(pop_rl)

  females <- selectWithinFam(pop_rl,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  He = calHe(pop_rl)

  ped_rl = rbind(ped_rl,createPed(pop = candidate))

}else{

  ped_rl = rbind(ped_ne,createPed(pop = candidate))
}

A_matrix = calAmatrix(ped = ped_rl,keep = candidate@id)

  if(g == 1){
    input_data = makejson(candidate = candidate,females = females,males = males,gen = ((nG+1)-g),He = He,A_matrix = A_matrix,gen_pre = g)
    # input_data2 <- input_data

  }else{
    input_data <- makejson(candidate = candidate,females = females,males = males,gen = ((nG+1)-g),He = He,A_matrix = A_matrix,gen_pre = g)
  }
  # Specify the output JSON file path
  output_json_path <- "input_data.json"

  # Write the JSON file
  write_json(input_data, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

  cat(paste("JSON input file has been generated at", output_json_path, "\n"))

  system('bash -c "source activate tf-gpu && python optMating.py input_data.json"')

  mating_plans <- fromJSON("breeding_pairs.json")[[1]]
  mating_plans<- as.matrix(mating_plans)
  #pop = selectInd(pop,nInd = 1000,use = "rand")

  pop_rl <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)
  keep_rl <- selectWithinFam(pop = pop_rl,nInd = 20,use = "rand",simParam = SP)

  if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}

pop_all = list(pop_tc,pop_ocs,pop_ran,pop_rl2,pop_65,pop_25,pop_90,pop_rate,pop_rl)
keep_all = list(keep_tc,keep_ocs,keep_ran,keep_rl2,keep_65,keep_25,keep_90,keep_rate,keep_rl)
ped_all = list(ped_tc,ped_ocs,ped_ran,ped_rl2,ped_65,ped_25,ped_90,ped_rate,ped_rl)
app = c("tc","ocs45","ran","rl2","ocs65","ocs25","ocs90","ocsrate","rl")
gv = c()
genetic = c()
genic = c()
inb = c()
inb_ped = c()
Ne_ped = c()
for(p in 1:length(pop_all)){
  gv = c(gv, mean(pop_all[[p]]@gv[,1]))
  genetic = c(genetic,varA(pop_all[[p]])[1,1])
  genic = c(genic,genicVarA(pop_all[[p]])[1])
  inb = c(inb, calinb(pop = pop_all[[p]],snpfre = snpfre_v))
  inb_ped = c(inb_ped,calInbped(ped = ped_all[[p]],keep_pop =keep_all[[p]]))
  Ne_ped = c(Ne_ped,calNeped(ped = ped_all[[p]],keep_pop =keep_all[[p]]))
}
  if(g == 1){
    output = data.table(app = app, gv = gv, genetic = genetic, genic = genic,inb = inb, Ne = Ne_ped,inb_ped = inb_ped,gen = rep(g,length(app)))
  }else{
    temp = data.table(app = app, gv = gv, genetic = genetic, genic = genic,inb = inb,Ne = Ne_ped,inb_ped = inb_ped, gen = rep(g,length(app)))

    output = rbind(output,temp)
  }

  # print(mean(pop_tc@gv[,2]))
  # print(mean(pop_ocs@gv[,2]))
  # print(mean(pop_rl@gv[,2]))
  # print(mean(pop_ran@gv[,2]))
}

setorder(output,gen)
fwrite(output,file = paste("output",r,".csv",sep = ""),sep = ",")
print(output)
}

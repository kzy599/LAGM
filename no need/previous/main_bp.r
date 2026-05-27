set.seed(42)
randomseed = sample(.Machine$integer.max,20)#.Machine$integer.max is 2147483647
for(r in 61:70){
rm(list = setdiff(ls(), c("r","randomseed")))
gc()
set.seed(randomseed[(r-60)])
nG = 20
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

  input_data2 <- makejson(candidate = candidate,females = females,males = males,gen = 7,He = He)

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
  pop_rl2_3 <- makeEBV(pop_rl2_3)

  females <- selectWithinFam(pop_rl2_3,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_3,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  He = calHe(pop_rl2_3)

  ped_rl2_3 = rbind(ped_rl2_3,createPed(pop = candidate))

}else{

  ped_rl2_3 = rbind(ped_ne,createPed(pop = candidate))
}

 input_data2_3 <- makejson(candidate = candidate,females = females,males = males,gen = 3,He = He)

  # Specify the output JSON file path
  output_json_path <- "input_data2_3.json"

  # Write the JSON file
  write_json(input_data2_3, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

  cat(paste("JSON input file has been generated at", output_json_path, "\n"))

  system('bash -c "source activate tf-gpu && python optMatingP.py input_data2_3.json"')

  mating_plans <- fromJSON("breeding_pairs2.json")[[1]]
  mating_plans<- as.matrix(mating_plans)
  #pop = selectInd(pop,nInd = 1000,use = "rand")

  pop_rl2_3 <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)
  
  keep_rl2_3 <- selectWithinFam(pop = pop_rl2_3,nInd = 20,use = "rand",simParam = SP)


if(g!=1){
  pop_rl2_5 <- makeEBV(pop_rl2_5)

  females <- selectWithinFam(pop_rl2_5,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_5,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  He = calHe(pop_rl2_5)

  ped_rl2_5 = rbind(ped_rl2_5,createPed(pop = candidate))

}else{

  ped_rl2_5 = rbind(ped_ne,createPed(pop = candidate))
}

  input_data2_5 <- makejson(candidate = candidate,females = females,males = males,gen = 5,He = He)

  # Specify the output JSON file path
  output_json_path <- "input_data2_5.json"

  # Write the JSON file
  write_json(input_data2_5, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

  cat(paste("JSON input file has been generated at", output_json_path, "\n"))

  system('bash -c "source activate tf-gpu && python optMatingP.py input_data2_5.json"')

  mating_plans <- fromJSON("breeding_pairs2.json")[[1]]
  mating_plans<- as.matrix(mating_plans)
  #pop = selectInd(pop,nInd = 1000,use = "rand")

  pop_rl2_5 <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)
  
  keep_rl2_5 <- selectWithinFam(pop = pop_rl2_5,nInd = 20,use = "rand",simParam = SP)



if(g!=1){
  pop_rl2_5m <- makeEBV(pop_rl2_5m)

  females <- selectWithinFam(pop_rl2_5m,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_5m,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  He = calHe(pop_rl2_5m)

  ped_rl2_5m = rbind(ped_rl2_5m,createPed(pop = candidate))

}else{

  ped_rl2_5m = rbind(ped_ne,createPed(pop = candidate))
}

  # input_data2_5m <- makejson(candidate = candidate,females = females,males = males,gen = calres(gen = g,Interval = 5),He = He)
  input_data2_5m <- makejson(candidate = candidate,females = females,males = males,gen = 1,He = He)
  # Specify the output JSON file path
  output_json_path <- "input_data2_5m.json"

  # Write the JSON file
  write_json(input_data2_5m, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

  cat(paste("JSON input file has been generated at", output_json_path, "\n"))

  system('bash -c "source activate tf-gpu && python optMatingP.py input_data2_5m.json"')

  mating_plans <- fromJSON("breeding_pairs2.json")[[1]]
  mating_plans<- as.matrix(mating_plans)
  #pop = selectInd(pop,nInd = 1000,use = "rand")

  pop_rl2_5m <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)
  
  keep_rl2_5m <- selectWithinFam(pop = pop_rl2_5m,nInd = 20,use = "rand",simParam = SP)





if(g!=1){
  pop_rl2_10m <- makeEBV(pop_rl2_10m)

  females <- selectWithinFam(pop_rl2_10m,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_10m,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  He = calHe(pop_rl2_10m)

  ped_rl2_10m = rbind(ped_rl2_10m,createPed(pop = candidate))

}else{

  ped_rl2_10m = rbind(ped_ne,createPed(pop = candidate))
}

  # input_data2_10m <- makejson(candidate = candidate,females = females,males = males,gen = calres(gen = g,Interval = 10),He = He)
  input_data2_10m <- makejson(candidate = candidate,females = females,males = males,gen = 2,He = He)
  # Specify the output JSON file path
  output_json_path <- "input_data2_10m.json"

  # Write the JSON file
  write_json(input_data2_10m, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

  cat(paste("JSON input file has been generated at", output_json_path, "\n"))

  system('bash -c "source activate tf-gpu && python optMatingP.py input_data2_10m.json"')

  mating_plans <- fromJSON("breeding_pairs2.json")[[1]]
  mating_plans<- as.matrix(mating_plans)
  #pop = selectInd(pop,nInd = 1000,use = "rand")

  pop_rl2_10m <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)
  
  keep_rl2_10m <- selectWithinFam(pop = pop_rl2_10m,nInd = 20,use = "rand",simParam = SP)


# if(g!=1){
#   pop_rl <- makeEBV(pop_rl)

#   females <- selectWithinFam(pop_rl,nInd=4,use = "ebv",sex = "F",trait = 1)

#   males <- selectWithinFam(pop_rl,nInd=2,use = "ebv",sex = "M",trait = 1)
  
#   candidate = c(females,males)

#   He = calHe(pop_rl)

#   ped_rl = rbind(ped_rl,createPed(pop = candidate))

# }else{

#   ped_rl = rbind(ped_ne,createPed(pop = candidate))
# }

# A_matrix = calAmatrix(ped = ped_rl,keep = candidate@id)

#   if(g == 1){
#     input_data = makejson(candidate = candidate,females = females,males = males,gen = ((nG+1)-g),He = He,A_matrix = A_matrix,gen_pre = g)
#     # input_data2 <- input_data

#   }else{
#     input_data <- makejson(candidate = candidate,females = females,males = males,gen = ((nG+1)-g),He = He,A_matrix = A_matrix,gen_pre = g)
#   }
#   # Specify the output JSON file path
#   output_json_path <- "input_data.json"

#   # Write the JSON file
#   write_json(input_data, path = output_json_path, pretty = TRUE, auto_unbox = TRUE)

#   cat(paste("JSON input file has been generated at", output_json_path, "\n"))

#   system('bash -c "source activate tf-gpu && python optMating.py input_data.json"')

#   mating_plans <- fromJSON("breeding_pairs.json")[[1]]
#   mating_plans<- as.matrix(mating_plans)
#   #pop = selectInd(pop,nInd = 1000,use = "rand")

#   pop_rl <- makeCross(pop = candidate, crossPlan = mating_plans, nProgeny = nProgenyPerCross, simParam = SP)
#   keep_rl <- selectWithinFam(pop = pop_rl,nInd = 20,use = "rand",simParam = SP)

if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}

pop_all = list(pop_rl2,pop_rl2_3,pop_rl2_5,pop_rl2_5m,pop_rl2_10m)
keep_all = list(keep_rl2,keep_rl2_3,keep_rl2_5,keep_rl2_5m,keep_rl2_10m)
ped_all = list(ped_rl2,ped_rl2_3,ped_rl2_5,ped_rl2_5m,ped_rl2_10m)
app = c("rl2","rl2_3","rl2_5","rl2_5m","rl2_10m")
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
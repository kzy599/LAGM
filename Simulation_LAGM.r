# remotes::install_local("/home/kangziyi/lagm_mating/lagmRcpp_original", force = TRUE)
set.seed(42)
randomseed = sample(.Machine$integer.max,20)#.Machine$integer.max is 2147483647
for(r in 1:20){
rm(list = setdiff(ls(), c("r","randomseed")))
gc()
set.seed(randomseed[r])
nG = 20
dir = "/home/kangziyi/poster/BaseMapHaplo/"

# feel free to contact me at:kangziyi1998@163.com
mapFile = "mergedMap.csv"
haploFile = "mergedHaplo.csv"

minInb = TRUE
Fix_fmRatio=FALSE

#default FALSE
rare_weight = FALSE #Experimental feature; disabled during testing/validation

if(Fix_fmRatio){
  f_min = 1L
  f_max = 1L
  m_min=2L
  m_max=2L
}else{
  f_min = 0L
  f_max = 2L
  m_min=0L
  m_max=4L
}

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

mapdt = makemap()
fwrite(mapdt, file = "hib.map",col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

nDam = 50
nSire = 25
nCrosses = 50
nProgenyPerCross = 100
nProgeny = nProgenyPerCross
pop <- selectCross(pop_founder,nFemale = nDam,nMale = nSire,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)
ped_ne = data.table(id = unique(c(pop@mother,pop@father)),sire = NA,dam = NA)
colnames(ped_ne) = c("id","sire","dam")


for(burn_in in c(1:5)){

  pop = estimate_ebv_pblup(pop = pop,ped=ped_ne)

  females <- selectWithinFam(pop,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)
  
  ped_ne = rbind(ped_ne,createPed(pop = candidate))

  pop = run_ocs_aqua(candidate = candidate,
  minInb=minInb,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 45,Fix_fmRatio=TRUE)

}


for(g in c(1:nG)){
  
  if(g==1){
    pop_candidate = estimate_ebv(pop = pop)

    females <- selectWithinFam(pop_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

    males <- selectWithinFam(pop_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
    candidate = c(females,males)

  }
if(g!=1){
  pop_rl2_candidate <- estimate_ebv(pop = pop_rl2)

  females <- selectWithinFam(pop_rl2_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)

  ped_rl2 = rbind(ped_rl2,createPed(pop = candidate))

}else{

  ped_rl2 = rbind(ped_ne,createPed(pop = candidate))
}

  pop_rl2 = lagm_mating(
  candidate = candidate,
  females = candidate[candidate@sex=="F"],
  males = candidate[candidate@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = ((nG+1)-g),
  female_min = rep(f_min, candidate[candidate@sex=="F"]@nInd),
  female_max = rep(f_max, candidate[candidate@sex=="F"]@nInd),
  male_min = rep(m_min, candidate[candidate@sex=="M"]@nInd),
  male_max = rep(m_max, candidate[candidate@sex=="M"]@nInd),
  diversity_mode = "genomic",
  base_diversity = 1,
  relationship_matrix = NULL,
  cooling_rate = 0.998, # 👉 配合高迭代次数，放缓降温
  stop_window = 10000,   # 10000次不进步则早停
  stop_eps = 1e-8,
  warmup_iter = 1000L,
  n_iter = 50000,
  n_pop = 300L,
  n_threads = 16L,
  n_progeny = nProgenyPerCross,
  rare_weight = rare_weight,
  sim_param = SP
)
  pop_rl2 = pop_rl2$offspring
  keep_rl2 <- selectWithinFam(pop = pop_rl2,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    pop_tc_candidate <- estimate_ebv(pop = pop_tc)
    candidate <- createCandidate(pop_tc_candidate)
    ped_tc = rbind(ped_tc,createPed(pop = candidate))
  }else{
    ped_tc = rbind(ped_ne,createPed(pop = candidate))
  }
  
  pop_tc <- run_ocs_aqua(candidate = candidate,
  minInb=minInb,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 0,Fix_fmRatio=Fix_fmRatio)

  keep_tc <- selectWithinFam(pop = pop_tc,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    pop_ocs_candidate <- estimate_ebv(pop = pop_ocs)
    candidate <- createCandidate(pop_ocs_candidate)
    ped_ocs = rbind(ped_ocs,createPed(pop = candidate))
  }else{
    ped_ocs = rbind(ped_ne,createPed(pop = candidate))
  }


  pop_ocs <- run_ocs_aqua(candidate = candidate,
  minInb=minInb,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 45,Fix_fmRatio=Fix_fmRatio)

  keep_ocs <- selectWithinFam(pop = pop_ocs,nInd = 20,use = "rand",simParam = SP)

if(g!=1){
    pop_ran_candidate <- estimate_ebv(pop = pop_ran)
    candidate <- createCandidate(pop_ran_candidate)
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
    pop_65_candidate <- estimate_ebv(pop = pop_65)
    candidate <- createCandidate(pop_65_candidate)
    ped_65 = rbind(ped_65,createPed(pop = candidate))
  }else{
    ped_65 = rbind(ped_ne,createPed(pop = candidate))
  }
  pop_65 <- run_ocs_aqua(candidate = candidate,
  minInb=minInb,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 65,Fix_fmRatio=Fix_fmRatio)
  keep_65 <- selectWithinFam(pop = pop_65,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    pop_25_candidate <- estimate_ebv(pop = pop_25)
    candidate <- createCandidate(pop_25_candidate)
    ped_25 = rbind(ped_25,createPed(pop = candidate))
  }else{
    ped_25 = rbind(ped_ne,createPed(pop = candidate))
  }
  pop_25 <- run_ocs_aqua(candidate = candidate,
  minInb=minInb,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 25,Fix_fmRatio=Fix_fmRatio)

  keep_25 = selectWithinFam(pop = pop_25,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    pop_90_candidate <- estimate_ebv(pop = pop_90)
    candidate <- createCandidate(pop_90_candidate)
    ped_90 = rbind(ped_90,createPed(pop = candidate))
  }else{
    ped_90 = rbind(ped_ne,createPed(pop = candidate))
  }
  pop_90 <- run_ocs_aqua(candidate = candidate,
  minInb=minInb,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 90,Fix_fmRatio=Fix_fmRatio)

  keep_90 = selectWithinFam(pop = pop_90,nInd = 20,use = "rand",simParam = SP)


  if(g!=1){
    pop_rate_candidate <- estimate_ebv(pop = pop_rate)
    candidate <- createCandidate(pop_rate_candidate)
    ped_rate = rbind(ped_rate,createPed(pop = candidate))
  }else{
    ped_rate = rbind(ped_ne,createPed(pop = candidate))
  }

  if(Fix_fmRatio){
  pop_rate <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        nFemalesMax = nDam,
        nMalesMax = nSire,
        equalizeFemaleContributions =TRUE,
        equalizeMaleContributions = TRUE,
        # maxMaleContribution = 4,
        # maxFemaleContribution = 2,
        targetCoancestryRate = 0.01,
        use = "ebv",
        minInbreedingMating = minInb
  )
  }else{
      pop_rate <- ocs(
        pop = candidate,
        nCrosses = nCrosses,
        nProgenyPerCross = nProgenyPerCross,
        # nFemalesMax = nDam,
        # nMalesMax = nSire,
        # equalizeFemaleContributions =TRUE,
        # equalizeMaleContributions = TRUE,
        maxMaleContribution = 4,
        maxFemaleContribution = 2,
        targetCoancestryRate = 0.01,
        use = "ebv",
        minInbreedingMating = minInb
  )
  }


  keep_rate = selectWithinFam(pop = pop_rate,nInd = 20,use = "rand",simParam = SP)



  if(g!=1){
  pop_rl_candidate <- estimate_ebv(pop = pop_rl)

  females <- selectWithinFam(pop_rl_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)


  ped_rl = rbind(ped_rl,createPed(pop = candidate))

  }else{

  ped_rl = rbind(ped_ne,createPed(pop = candidate))
}
  
G_matrix = makeGDE(x = pullSnpGeno(candidate), type="G",inv=FALSE)
pop_rl = lagm_mating(
  candidate = candidate,
  females = candidate[candidate@sex=="F"],
  males = candidate[candidate@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = ((nG+1)-g),
  female_min = rep(f_min, candidate[candidate@sex=="F"]@nInd),
  female_max = rep(f_max, candidate[candidate@sex=="F"]@nInd),
  male_min = rep(m_min, candidate[candidate@sex=="M"]@nInd),
  male_max = rep(m_max, candidate[candidate@sex=="M"]@nInd),
  diversity_mode = "relationship",
  base_diversity = 1,
  relationship_matrix = G_matrix,
  cooling_rate = 0.998, # 👉 配合高迭代次数，放缓降温
  stop_window = 10000,   # 10000次不进步则早停
  stop_eps = 1e-8,
  warmup_iter = 1000L,
  n_iter = 50000,
  n_pop = 300L,
  n_threads = 16L,
  n_progeny = nProgenyPerCross,
  rare_weight = rare_weight,
  sim_param = SP
)
  pop_rl = pop_rl$offspring
  keep_rl <- selectWithinFam(pop = pop_rl,nInd = 20,use = "rand",simParam = SP)

  if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}

if(g==1){
  app = c("tc","ocs45","ran","rl2","ocs65","ocs25","ocs90","ocsrate","rl")
  gv = rep(mean(pop@gv[,1]), length(app))
  genetic = rep(varA(pop)[1,1], length(app))
  genic = rep(genicVarA(pop)[1], length(app))
  inb = rep(calinb(pop = pop,snpfre = snpfre_v), length(app))
  keep_pop = selectWithinFam(pop,nInd = 20,use = "rand",simParam = SP)
  inb_ped = rep(calInbped(ped = ped_ne,keep_pop =keep_pop), length(app))
  Ne_ped = rep(calNeped(ped = ped_ne,keep_pop =keep_pop), length(app))
  He = rep(calHe(pop), length(app))
  Ho = rep(calHo(pop), length(app))
  output = data.table(app = app, gv = gv, genetic = genetic, genic = genic,inb = inb, Ne = Ne_ped,inb_ped = inb_ped,He = He,Ho = Ho,gen = rep(0,length(app)))
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
He = c()
Ho = c()

for(p in 1:length(pop_all)){
  gv = c(gv, mean(pop_all[[p]]@gv[,1]))
  genetic = c(genetic,varA(pop_all[[p]])[1,1])
  genic = c(genic,genicVarA(pop_all[[p]])[1])
  inb = c(inb, calinb(pop = pop_all[[p]],snpfre = snpfre_v))
  inb_ped = c(inb_ped,calInbped(ped = ped_all[[p]],keep_pop =keep_all[[p]]))
  Ne_ped = c(Ne_ped,calNeped(ped = ped_all[[p]],keep_pop =keep_all[[p]]))
  He = c(He,calHe(pop_all[[p]]))
  Ho = c(Ho,calHo(pop_all[[p]]))
}
  if(g == 1){

    temp = data.table(app = app, gv = gv, genetic = genetic, genic = genic,inb = inb,Ne = Ne_ped,inb_ped = inb_ped,He = He,Ho = Ho, gen = rep(g,length(app)))
    
    output = rbind(output,temp)


  }else{
    temp = data.table(app = app, gv = gv, genetic = genetic, genic = genic,inb = inb,Ne = Ne_ped,inb_ped = inb_ped,He = He,Ho = Ho, gen = rep(g,length(app)))

    output = rbind(output,temp)
  }


}

setorder(output,gen)
fwrite(output,file = paste("output",r,".csv",sep = ""),sep = ",")
print(output)
}

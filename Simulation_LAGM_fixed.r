set.seed(42)
randomseed = sample(.Machine$integer.max,20)#.Machine$integer.max is 2147483647
for(r in 1:20){
rm(list = setdiff(ls(), c("r","randomseed")))
gc()
set.seed(randomseed[r])
nG = 20
dir = "/home/kangziyi/poster/BaseMapHaplo/"

#feel free to contact me at: kangziyi1998@163.com
mapFile = "mergedMap.csv"
haploFile = "mergedHaplo.csv"

minInb = TRUE
Fix_fmRatio=FALSE
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

SP$setVarE(h2 = c(0.3,0.3))# don't worried about this, just the first trait was used in this simulation
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
pop <- selectCross(pop_founder,
                   nFemale = nDam,nMale = nSire,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)
ped_ne = data.table(id = c(unique(pop@mother,pop@father)),sire = NA,dam = NA)
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
  lookahead_generations = 7,
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
  sim_param = SP
)
  pop_rl2 = pop_rl2$offspring
  keep_rl2 <- selectWithinFam(pop = pop_rl2,nInd = 20,use = "rand",simParam = SP)


if(g!=1){
  pop_rl2_3_candidate = estimate_ebv(pop_rl2_3)

  females <- selectWithinFam(pop_rl2_3_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_3_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)


  ped_rl2_3 = rbind(ped_rl2_3,createPed(pop = candidate))

}else{

  ped_rl2_3 = rbind(ped_ne,createPed(pop = candidate))
}


  pop_rl2_3 = lagm_mating(
  candidate = candidate,
  females = candidate[candidate@sex=="F"],
  males = candidate[candidate@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = 3,
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
  sim_param = SP
)
  pop_rl2_3 = pop_rl2_3$offspring
  
  keep_rl2_3 <- selectWithinFam(pop = pop_rl2_3,nInd = 20,use = "rand",simParam = SP)


if(g!=1){
  pop_rl2_5_candidate = estimate_ebv(pop_rl2_5)

  females <- selectWithinFam(pop_rl2_5_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_5_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)


  ped_rl2_5 = rbind(ped_rl2_5,createPed(pop = candidate))

}else{

  ped_rl2_5 = rbind(ped_ne,createPed(pop = candidate))
}

  pop_rl2_5 = lagm_mating(
  candidate = candidate,
  females = candidate[candidate@sex=="F"],
  males = candidate[candidate@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = 5,
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
  sim_param = SP
)
  pop_rl2_5 = pop_rl2_5$offspring
  
  keep_rl2_5 <- selectWithinFam(pop = pop_rl2_5,nInd = 20,use = "rand",simParam = SP)



if(g!=1){
  pop_rl2_5m_candidate = estimate_ebv(pop_rl2_5m)

  females <- selectWithinFam(pop_rl2_5m_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_5m_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)


  ped_rl2_5m = rbind(ped_rl2_5m,createPed(pop = candidate))

}else{

  ped_rl2_5m = rbind(ped_ne,createPed(pop = candidate))
}

  pop_rl2_5m = lagm_mating(
  candidate = candidate,
  females = candidate[candidate@sex=="F"],
  males = candidate[candidate@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = 1,
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
  sim_param = SP
)
  pop_rl2_5m = pop_rl2_5m$offspring
  
  keep_rl2_5m <- selectWithinFam(pop = pop_rl2_5m,nInd = 20,use = "rand",simParam = SP)





if(g!=1){
  pop_rl2_10m_candidate = estimate_ebv(pop_rl2_10m)

  females <- selectWithinFam(pop_rl2_10m_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_10m_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate = c(females,males)


  ped_rl2_10m = rbind(ped_rl2_10m,createPed(pop = candidate))

}else{

  ped_rl2_10m = rbind(ped_ne,createPed(pop = candidate))
}


  pop_rl2_10m = lagm_mating(
  candidate = candidate,
  females = candidate[candidate@sex=="F"],
  males = candidate[candidate@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = 2,
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
  sim_param = SP
)
  pop_rl2_10m = pop_rl2_10m$offspring
  
  keep_rl2_10m <- selectWithinFam(pop = pop_rl2_10m,nInd = 20,use = "rand",simParam = SP)

if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}

if(g==1){
  app = c("rl2","rl2_3","rl2_5","rl2_5m","rl2_10m")
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
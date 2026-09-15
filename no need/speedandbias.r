# remotes::install_local("/home/kangziyi/lagm_mating/lagmRcpp_original", force = TRUE)
set.seed(42)
randomseed = sample(.Machine$integer.max,20)#.Machine$integer.max is 2147483647
for(r in 2:2){
rm(list = setdiff(ls(), c("r","randomseed")))
gc()
set.seed(randomseed[r])
nG = 20
dir = "/home/kangziyi/poster/BaseMapHaplo/"
mapFile = "mergedMap.csv"
haploFile = "mergedHaplo.csv"
minInb = FALSE
Fix_fmRatio=FALSE
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

calbias = function(pop){
y <- pop@pheno[, 1]
x <- pop@ebv[,1]
bias_slope <- cov(y, x, use = "complete.obs") / var(x, na.rm = TRUE)
return(bias_slope)
}

calbias_true = function(pop){
y <- pop@gv[, 1]
x <- pop@ebv[,1]
bias_slope <- cov(y, x, use = "complete.obs") / var(x, na.rm = TRUE)
return(bias_slope)
}

calacc = function(pop){
  return(cov(pop@ebv[,1],pop@gv[,1]))
}

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

  # pop@ebv = matrix(pop@pheno[,1],ncol = 1)

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
  
  # ped = rbind(ped,data.table(ID = pop@id,father = pop@father,mother = pop@mother))
  # phe = rbind(phe, data.table(ID = pop@id,env = rep(1,pop@nInd),
  #                             gv1 = pop@gv[,1],gv2 = pop@gv[,2],
  #                             phe1 = pop@pheno[,1],phe2 = pop@pheno[,2],
  #                             FamilyID = paste(pop@father,"_",pop@mother,sep = "")))
  if(g==1){
    # pop <- makeEBV(pop)
    pop_candidate = estimate_ebv(pop = pop)
    # pop_candidate = makingGVpop(pop = pop)

    females <- selectWithinFam(pop_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

    males <- selectWithinFam(pop_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
    candidate_300 = c(females,males)

    females <- selectWithinFam(pop_candidate,nInd=12,use = "ebv",sex = "F",trait = 1)

    males <- selectWithinFam(pop_candidate,nInd=6,use = "ebv",sex = "M",trait = 1)
  
    candidate_900 = c(females,males)

    females <- selectWithinFam(pop_candidate,nInd=16,use = "ebv",sex = "F",trait = 1)

    males <- selectWithinFam(pop_candidate,nInd=8,use = "ebv",sex = "M",trait = 1)
  
    candidate_1200 = c(females,males)
    # He = calHe(pop)

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
  # pop_rl2 <- makeEBV(pop_rl2)
  pop_rl2_candidate <- estimate_ebv(pop = pop_rl2)
  # pop_rl2_candidate = makingGVpop(pop = pop_rl2)

  females <- selectWithinFam(pop_rl2_candidate,nInd=4,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_candidate,nInd=2,use = "ebv",sex = "M",trait = 1)
  
  candidate_300 = c(females,males)


  females <- selectWithinFam(pop_rl2_candidate,nInd=12,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_candidate,nInd=6,use = "ebv",sex = "M",trait = 1)
  
  candidate_900 = c(females,males)



  females <- selectWithinFam(pop_rl2_candidate,nInd=16,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_rl2_candidate,nInd=8,use = "ebv",sex = "M",trait = 1)
  
  candidate_1200 = c(females,males)

  # He = calHe(pop_rl2)

  ped_rl2 = rbind(ped_rl2,createPed(pop = candidate_300))

}else{

  ped_rl2 = rbind(ped_ne,createPed(pop = candidate_300))
}
  start_time <- Sys.time()
  pop_rl2 = lagm_mating(
  candidate = candidate_1200,
  females = candidate_1200[candidate_1200@sex=="F"],
  males = candidate_1200[candidate_1200@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = (nG-g+1),
  female_min = rep(f_min, candidate_1200[candidate_1200@sex=="F"]@nInd),
  female_max = rep(f_max, candidate_1200[candidate_1200@sex=="F"]@nInd),
  male_min = rep(m_min, candidate_1200[candidate_1200@sex=="M"]@nInd),
  male_max = rep(m_max, candidate_1200[candidate_1200@sex=="M"]@nInd),
  diversity_mode = "genomic",
  base_diversity = 1,
  relationship_matrix = NULL,
  cooling_rate = 0.998, # 👉 配合高迭代次数，放缓降温
  stop_window = 1000,   # 1000次不进步则早停
  stop_eps = 1e-8,
  warmup_iter = 200L,
  n_iter = 10000,
  n_pop = 100L,
  n_threads = 16L,
  n_progeny = nProgenyPerCross,
  rare_weight = rare_weight,
  sim_param = SP
)
 end_time <- Sys.time()
 elapsed_time_1200_rl2 <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  pop_rl2 = lagm_mating(
  candidate = candidate_900,
  females = candidate_900[candidate_900@sex=="F"],
  males = candidate_900[candidate_900@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = (nG-g+1),
  female_min = rep(f_min, candidate_900[candidate_900@sex=="F"]@nInd),
  female_max = rep(f_max, candidate_900[candidate_900@sex=="F"]@nInd),
  male_min = rep(m_min, candidate_900[candidate_900@sex=="M"]@nInd),
  male_max = rep(m_max, candidate_900[candidate_900@sex=="M"]@nInd),
  diversity_mode = "genomic",
  base_diversity = 1,
  relationship_matrix = NULL,
  cooling_rate = 0.998, # 👉 配合高迭代次数，放缓降温
  stop_window = 1000,   # 1000次不进步则早停
  stop_eps = 1e-8,
  warmup_iter = 200L,
  n_iter = 10000,
  n_pop = 100L,
  n_threads = 16L,
  n_progeny = nProgenyPerCross,
  rare_weight = rare_weight,
  sim_param = SP
)
 end_time <- Sys.time()
 elapsed_time_900_rl2 <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  pop_rl2 = lagm_mating(
  candidate = candidate_300,
  females = candidate_300[candidate_300@sex=="F"],
  males = candidate_300[candidate_300@sex=="M"],
  n_crosses=nCrosses,
  lookahead_generations = (nG-g+1),
  female_min = rep(f_min, candidate_300[candidate_300@sex=="F"]@nInd),
  female_max = rep(f_max, candidate_300[candidate_300@sex=="F"]@nInd),
  male_min = rep(m_min, candidate_300[candidate_300@sex=="M"]@nInd),
  male_max = rep(m_max, candidate_300[candidate_300@sex=="M"]@nInd),
  diversity_mode = "genomic",
  base_diversity = 1,
  relationship_matrix = NULL,
  cooling_rate = 0.998, # 👉 配合高迭代次数，放缓降温
  stop_window = 1000,   # 1000次不进步则早停
  stop_eps = 1e-8,
  warmup_iter = 200L,
  n_iter = 10000,
  n_pop = 100L,
  n_threads = 16L,
  n_progeny = nProgenyPerCross,
  rare_weight = rare_weight,
  sim_param = SP
)
 end_time <- Sys.time()
 elapsed_time_300_rl2 <- difftime(end_time, start_time, units = "secs")

  pop_rl2 = pop_rl2$offspring
 
  keep_rl2 <- selectWithinFam(pop = pop_rl2,nInd = 20,use = "rand",simParam = SP)

  if(g!=1){
    # pop_25 <- makeEBV(pop_25)
    pop_25_candidate <- estimate_ebv(pop = pop_25)
    # pop_25_candidate = makingGVpop(pop= pop_25)

  candidate_300 <- createCandidate(pop_25_candidate)


  females <- selectWithinFam(pop_25_candidate,nInd=12,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_25_candidate,nInd=6,use = "ebv",sex = "M",trait = 1)
  
  candidate_900 = c(females,males)



  females <- selectWithinFam(pop_25_candidate,nInd=16,use = "ebv",sex = "F",trait = 1)

  males <- selectWithinFam(pop_25_candidate,nInd=8,use = "ebv",sex = "M",trait = 1)
  
  candidate_1200 = c(females,males)

    ped_25 = rbind(ped_25,createPed(pop = candidate_300))
  }else{
    ped_25 = rbind(ped_ne,createPed(pop = candidate_300))
  }
  start_time <- Sys.time()
  pop_25 <- run_ocs_aqua(candidate = candidate_1200,
  minInb=TRUE,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 45,Fix_fmRatio=Fix_fmRatio)
  end_time <- Sys.time()
  elapsed_time_1200_25 <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  pop_25 <- run_ocs_aqua(candidate = candidate_900,
  minInb=TRUE,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 45,Fix_fmRatio=Fix_fmRatio)
  end_time <- Sys.time()
  elapsed_time_900_25 <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  pop_25 <- run_ocs_aqua(candidate = candidate_300,
  minInb=TRUE,nCrosses = nCrosses,
  nProgenyPerCross=nProgenyPerCross,nDam=nDam,nSire=nSire,
  targetDegree = 45,Fix_fmRatio=Fix_fmRatio)
  end_time <- Sys.time()
  elapsed_time_300_25 <- difftime(end_time, start_time, units = "secs")

  keep_25 = selectWithinFam(pop = pop_25,nInd = 20,use = "rand",simParam = SP)

  if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}


app = c("LAGM","ocs45")
candidate_scale = c("n_300","n_900","n_1200")
app_sec = rep(app,each=3)
ncand = rep(candidate_scale,2)
time_secs = c(elapsed_time_300_rl2,elapsed_time_900_rl2,elapsed_time_1200_rl2,elapsed_time_300_25,elapsed_time_900_25,elapsed_time_1200_25)
gblup_bias = c()
gblup_biasP_true = c()
gblup_acc = c()

if(g != 1) pop_all = list(pop_rl2_candidate,pop_25_candidate)
for(p in c(1:2)){
  if(g==1) gblup_bias = c(gblup_bias, calbias(pop_candidate)) else gblup_bias = c(gblup_bias, calbias(pop_all[[p]]))
  if(g==1) gblup_acc = c(gblup_acc, calacc(pop_candidate)) else gblup_acc = c(gblup_acc, calacc(pop_all[[p]]))
  if(g==1) gblup_biasP_true = c(gblup_biasP_true, calbias_true(pop_candidate)) else gblup_biasP_true = c(gblup_biasP_true, calbias_true(pop_all[[p]]))

}

  if(g == 1){

    output = data.table(app = app, gen = rep(g,length(app)),gblup_bias = gblup_bias, gblup_acc = gblup_acc,gblup_bias_true = gblup_biasP_true)

    output_sec = data.table(app = app_sec, candidate_scale = ncand,gen = rep(g,length(app_sec)),time_secs = time_secs)

  }else{
    temp = data.table(app = app, gen = rep(g,length(app)),gblup_bias = gblup_bias, gblup_acc = gblup_acc,gblup_bias_true = gblup_biasP_true)

    temp_sec = data.table(app = app_sec, candidate_scale = ncand,gen = rep(g,length(app_sec)),time_secs = time_secs)

    output = rbind(output,temp)

    output_sec = rbind(output_sec,temp_sec)
  }

  # print(mean(pop_tc@gv[,2]))
  # print(mean(pop_ocs@gv[,2]))
  # print(mean(pop_rl@gv[,2]))
  # print(mean(pop_ran@gv[,2]))
}

setorder(output,gen)
setorder(output_sec,gen)
fwrite(output,file = paste("outputbias",r,".csv",sep = ""),sep = ",")
fwrite(output_sec,file = paste("outputsec",r,".csv",sep = ""),sep = ",")
print(output)
}

# model <- lm(gain_std ~ lost_diversity * app, data = dt_plot)

# # 查看整体方差分析表
# anova(model)

# library(emmeans)

# # 计算并比较各个策略的斜率 (trend)
# # var = "Diversity" 告诉软件我们要比较的是关于 Diversity 的斜率
# slope_comparison <- emtrends(model, pairwise ~ app, var = "lost_diversity")

# # 1. 查看每个策略的具体斜率数值和置信区间
# print(slope_comparison$emtrends)

# # 2. 查看两两比较的显著性结果（P值）
# print(slope_comparison$contrasts)

rm(list = ls())
library(data.table)
library(AlphaSimR)
#easy population
founderPop = quickHaplo(nInd=1000, nChr=10, segSites=2000, inbred=FALSE)
#set genetic parameters
SP = SimParam$new(founderPop)
SP$restrSegSites(minSnpFreq = 0.05,overlap = TRUE)
SP$addTraitA(100,mean=c(0,0),var = c(1,1),corA = matrix(c(1,0.8,0.8,1),nrow = 2))
SP$setVarE(h2=c(0.41,0.41))#遗传力和偏差
SP$addSnpChip(1250)#55kSNP芯片
SP$setSexes("yes_sys")#按照一个雌一个雄来分配个体的性别
pop_founder = newPop(founderPop, simParam=SP)
pop <- selectCross(pop_founder,
                   nFemale = 50,nMale = 25,
                   nCrosses = 50,nProgeny = 50,
                   use = "gv",
                   simParam = SP)

candidates = pop_founder[unique(c(pop@mother,pop@father))]

calrel = function(pop,candidates,pop_previous){
  parents_geno = pullSnpGeno(candidates)
  ped_map = data.table(id = pop@id,sire = pop@father,dam = pop@mother,fam = paste(pop@father,pop@mother,sep = "_"))
  rel_current = c()
  snp_012_dt = pullSnpGeno(pop_previous)
  snpfre <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
  for(f in unique(ped_map$fam)){
    sire = unique(ped_map[fam ==f,sire])
    dam = unique(ped_map[fam ==f,dam])
    sire_geno  = parents_geno[rownames(parents_geno)==sire,]
    dam_geno  = parents_geno[rownames(parents_geno)==dam,]

    sire_fre = sire_geno - 2*snpfre
    dam_fre = dam_geno - 2*snpfre

    rel = mean(sire_fre*dam_fre/(2*snpfre*(1-snpfre)))/2 #VR2

    # rel = (sum(sire_fre*dam_fre)/sum((2*snpfre*(1-snpfre))))/2#VR1
    
    rel_current = c(rel_current,rel)
  }
  rel_current = mean(rel_current)
  
  #using pop for prediction accuracy of He
  #using pop_previous checking the retention rate per generation
  #using pop_funder for the Ne calculation
  #rel_previous = calinb(pop = pop,snpfre = snpfre)

  return(rel_current)
}


calHe_predicted = function(pop,candidates,pop_previous){
  parents_geno = pullSnpGeno(candidates)
  ped_map = data.table(id = pop@id,sire = pop@father,dam = pop@mother,fam = paste(pop@father,pop@mother,sep = "_"))
  He_current = c()
  for(f in unique(ped_map$fam)){
    sire = unique(ped_map[fam ==f,sire])
    dam = unique(ped_map[fam ==f,dam])
    sire_geno  = parents_geno[rownames(parents_geno)==sire,]
    dam_geno  = parents_geno[rownames(parents_geno)==dam,]
    sire_fre = sire_geno/2
    dam_fre = dam_geno/2
    He_current = c(He_current,mean(sire_fre+dam_fre-2*sire_fre*dam_fre))
  }
  He_current = mean(He_current)
  
  #using pop for prediction accuracy of He
  #using pop_previous checking the retention rate per generation
  #using pop_funder for the Ne calculation
  He_previous = apply(pullSnpGeno(pop), 2, function(x){
    # p = sum(x)/(2*length(x))
    # q = 1-p
    # return(2*p*q)
    return(sum(x==1)/length(x))
  })
  He_previous = mean(He_previous)
  
  Ne = 1/(2*(1- (He_current/He_previous)^(1/g)))
  
  return(He_current)
}

calper = function(x){
  retentionrate = c()
  for(i in 1:(length(x)-1)){
   retentionrate = c(retentionrate,x[i+1]/x[i])
  }
  return(retentionrate)
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

caldf = function(pop1,pop2,snpfre){
  df = (1-calinb(pop = pop2,snpfre = snpfre_v))/(1-calinb(pop = pop1,snpfre = snpfre_v))
  return(df)
}
caldf_pre = function(rel_predicted){
  per = c()
  for(i in 1:(length(rel_predicted)-1)){
   per = c(per, (1-rel_predicted[i+1])/(1-rel_predicted[i]))
  }
  return(per)
}

if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}

He_predicted = c()
He = c()
deltaF = c()
rel_predicted = c()
for(g in c(1:5)){
  cand = selectWithinFam(pop,nInd = 6,use = "gv")
  pop_temp <- selectCross(cand,
                     nFemale = 50,nMale = 25,
                     nCrosses = 50,nProgeny = 50,
                     use = "gv",
                     simParam = SP)
  
  candidates = pop[unique(c(pop_temp@mother,pop_temp@father))]
  
  He_predicted = c(He_predicted,calHe_predicted(pop = pop_temp,candidates = candidates,pop_previous = pop))
  rel_predicted = c(rel_predicted,calrel(pop = pop_temp,candidates = candidates,pop_previous = pop))
  He_temp = apply(pullSnpGeno(pop_temp), 2, function(x){
    return(sum(x==1)/(length(x)))
  })

  deltaF = c(deltaF,caldf(pop1 = pop,pop2 = pop_temp,snpfre = snpfre_v))

  He = c(He,mean(He_temp))
  
  pop = pop_temp
}

1/(2*(1-calper(He)))

g = 0
mean(abs(c((He[1]*(calper(He)[1])^(4+g))/He[5+g],
(He[2]*(calper(He)[2])^(3+g))/He[5+g],
(He[3]*(calper(He)[3])^(2+g))/He[5+g],
(He[4]*(calper(He)[4])^(1+g))/He[5+g])-1))

mean(abs(c((He[1]*(calper(He_predicted)[1])^(4+g))/He[5+g],
(He[2]*(calper(He_predicted)[2])^(3+g))/He[5+g],
(He[3]*(calper(He_predicted)[3])^(2+g))/He[5+g],
(He[4]*(calper(He_predicted)[4])^(1+g))/He[5+g])-1))

1/(2*(1-(He[5]/He[1])^(1/5)))
1/(2*(1-(He[4]/He[1])^(1/4)))
1/(2*(1-(He[3]/He[1])^(1/3)))
1/(2*(1-(He[2]/He[1])^(1/2)))
mean(1/(2*(1-calper(He))))

deltaF

mean(c((He[1]*(deltaF[2])^(4+g))/He[5+g],
(He[2]*(deltaF[3])^(3+g))/He[5+g],
(He[3]*(deltaF[4])^(2+g))/He[5+g],
(He[4]*(deltaF[5])^(1+g))/He[5+g]))-1

mean(abs(c((He[1]*(caldf_pre(rel_predicted)[1])^(4+g))/He[5+g],
(He[2]*(caldf_pre(rel_predicted)[2])^(3+g))/He[5+g],
(He[3]*(caldf_pre(rel_predicted)[3])^(2+g))/He[5+g],
(He[4]*(caldf_pre(rel_predicted)[4])^(1+g))/He[5+g])-1))




He_predicted
He
ck = data.table(He_predicted,He)
cor(He_predicted,He)
Ne = c()
for(g in c(1:5)){
  decayrate= ck$He_predicted[g]
  Ne = c(Ne,1/(2*(1- (decayrate)^(1/g))))
}
ck
Ne

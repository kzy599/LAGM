rm(list = ls())
gc()
library(data.table)
library(dplyr)
library(multcompView)
getValue = function(rep,ap){
    value_c = c()
    for( r in rep){
        out = fread(paste("output",r,".csv",sep=""),sep = ",")
        value_c = c(value_c,out[app==ap&gen== 20,gv]/out[app==ap&gen== 20,inb_ped])
    }
    return(value_c)
}

ap = "rl2"

geno_old = c(11:20)
geno = c(21:30)
rel = c(31:40)
geno_new = c(41:50)

value_geno = getValue(rep = geno,ap=ap)

value_rel = getValue(rep = rel,ap=ap)

value_geno_new = getValue(rep = geno_new,ap=ap)

value_geno_old = getValue(rep = geno_old,ap=ap)

t.test(x = value_geno,y = value_rel)
t.test(x = value_geno_new,y = value_rel)
t.test(x = value_geno,y = value_geno_new)
t.test(x = value_geno,y = value_geno_old)

dt = rbind(data.table(type = "value_geno",value = value_geno),
data.table(type = "value_rel",value = value_rel),
data.table(type = "value_geno_new",value = value_geno_new),
data.table(type = "value_geno_old",value = value_geno_old))

g = "value"
aov_res <- aov(as.formula(paste(g, "~ type")), data = dt)
tukey_res <- TukeyHSD(aov_res,conf.level = 0.95)
letters <- multcompLetters4(aov_res, tukey_res)

tukey_res
letters


line_x = function(x){
  return(5*x)
}

line_y = function(y){
  return(9*y)
}

for(i in 1:10){
print(calpercent(x = line_y(i),y = line_x(i)))
}


calpercent = function(x,y){
  x = round(x,2)
  y = round(y,2)
  return(round(((x - y)/y)*100,1))
}

#bp_h
#ocs
ce5 = c(63.5,53.84,49.84,58.9,45.34)
ce10 = c(47.26,45.95,58.54,64.59,49.08)
ce15 = c(42.39,41.91,54.33,59.17,48.19)
ce20 = c(37.31,37.74,50.25,56.22,47.6)
calpercent(x = max(ce5),y = ce5)
calpercent(x = max(ce10),y = ce10)
calpercent(x = max(ce15),y = ce15)
calpercent(x = max(ce20),y = ce20)
calpercent(x = 73.97,y = ce5)
calpercent(x = 65.6,y = ce10)
calpercent(x = 57.16,y = ce15)
calpercent(x = 51.76,y = ce20)


calpercent(x = 6.22,y = c(7,5.36))
calpercent(x = 0.04,y = 0.07)
calpercent(x = 0.11,y = c(0.13,0.09))


calpercent(y = 38.27,x = ce5)
calpercent(y = 34.11,x = ce10)
calpercent(y = 29.76,x = ce15)
calpercent(y = 27.29,x = ce20)
calpercent(x = 73.97,y = 38.27)
calpercent(x = 65.6,y = 34.11)
calpercent(x = 57.16,y = 29.76)
calpercent(x = 51.76,y = 27.29)

calpercent(x = 0.05,y = 0.07)
calpercent(x = 0.07,y = 0.1)
calpercent(x = 0.09,y = 0.12)
calpercent(y = c(0.1,0.19,0.28,0.36), x = c(0.07,0.13,0.19,0.25))

#bp
calpercent = function(x,y){
  x = round(x,2)
  y = round(y,2)
  return(round(((x - y)/y)*100,1))
}

gg5 = c(2.2,2.29,2.58,2.81,3.01)
gg10 = c(4.02,4.34,4.76,5.07,5.42)
gg15 = c(5.61,5.93,6.45,6.89,7.39)
gg20 = c(6.83,7.21,7.79,8.29,8.81)

calpercent(x = max(gg20),y = gg20)

inb5  = c(0.03,0.04,0.03,0.03,0.03)
inb10 = c(0.07,0.09,0.08,0.07,0.07)
inb15 = c(0.15,0.13,0.12,0.11,0.10)
inb20 = c(0.2,0.17,0.16,0.14,0.14)
calpercent(x = max(inb20),y = inb20)



ce15 = c(50.22, 56.21)
ce20 = c(44.78, 50.93)

ce20 = c(44.78, 49.98)
calpercent(x = max(ce20),y = ce20)

calpercent(x = c(82.45,82),y = c(66.58,58.39))


calpercent(x = c(0.07,0.1,0.14),y = c(0.06,0.09,0.12))
calpercent(x = c(61.37,54.97,49.98),y = c(63.33,61,56.22))



#discussion 45 vs 65 in inbreeding and genetic gain acorss horizons
#extrem compared to explain rate of inbreeding across horizons
calpercent(x = c(2.73,5.07,6.97,8.41), y = c(2.02,3.83,5.47,6.8))
calpercent(x = c(0.04,0.08,0.13,0.17), y = c(0.04,0.06,0.09,0.12))

calpercent(x = c(3.58,6.38,8.39,9.8), y = c(1.19,2.41,3.46,4.43))
calpercent(x = c(0.1,0.19,0.28,0.36), y = c(0.03,0.05,0.07,0.09))


#comparaing LAGM with GOCS at G5 in bp
calpercent(x = c(76.45,82.45), y = c(66.58,46.29))

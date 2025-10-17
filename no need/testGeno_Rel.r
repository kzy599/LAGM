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


#bp
calpercent(x = 47.88,y =38.24)

gg5 = c(3.39,3.15,3.1,2.78,2.1)
gg10 = c(5.91,5.35,5.18,4.67,4.42)
gg15 = c(7.54,6.83,6.54,5.95,5.71)
gg20 = c(8.57,7.94,7.58,6.9,6.66)


calpercent(x = max(gg10),y = gg10)

inb5  = c(0.05,0.04,0.03,0.03)
inb10 = c(0.11,0.09,0.08,0.07,0.07)
inb15 = c(0.15,0.13,0.12,0.11,0.10)
inb20 = c(0.2,0.17,0.16,0.14,0.13)
calpercent(x = max(inb5),y = inb5)

ce5 = c(87.88, 69.2)
ce10 = c(65, 56.41)
ce15 = c(55.94, 49.38)
ce20 = c(49.46, 44.02)
calpercent(x = max(ce5),y = ce5)


calpercent(x = c(5.35,5.18), y = 5.92)
calpercent(x = c(6.83,6.54), y = 7.97)
calpercent(x = c(7.94,7.58), y = 9.52)
calpercent(x = 8.57, y = 9.52)


calpercent(x = c(4.67,4.42), y = 5.07)
calpercent(x = c(5.95,5.71), y = 6.97)
calpercent(x = c(6.9,6.66), y = 8.41)


calpercent(x = c(0.04,0.05), y = 0.07)
calpercent(x = c(0.08,0.11), y = 0.13)
calpercent(x = c(0.12,0.15), y = 0.19)
calpercent(x = c(0.16,0.2), y = 0.25)

calpercent(x = c(0.2), y = 0.17)


#2- and 3- generation compared to G0CS45 and 65
calpercent(x = c(3.15,3.1), y = 2.73)
calpercent(x = c(7.94,7.58), y = 8.41)

calpercent(x = c(0.04), y = 0.07)
calpercent(x = c(0.08,0.09), y = 0.13)
calpercent(x = c(0.12,0.13), y = 0.19)
calpercent(x = c(0.16,0.17), y = 0.25)


#5- and 7- generation compared to G0CS45 and 65
calpercent(x = c(0.13,0.14), y = 0.17)
calpercent(x = c(0.13,0.14), y = 0.12)


a=-1
b=4
c=5
x1 = (-b-sqrt(b^2-4*a*c))/(2*a)

library(data.table)
library(ggplot2)
curveA = function(x){
  return(-x^2+1)
}
minmax = function(a){
  return((a-min(a))/(max(a)-min(a)))
}
plotcurve = function(x,y){
  
df <- data.frame(x = x, y = y)

p =   ggplot() +
    geom_line(data = df, aes(x, y), color = "red", linewidth = 1) +  # Curve
    geom_point(data = df, aes(x, y), size = 3, color = "blue") +    # Points
    labs(title = "Transformed Coordinates Plot",
         x = "x (possibly transformed)",
         y = "y (possibly transformed)") +
    theme_minimal()

return(p)

}

x <- seq(0, 1, length.out = 100)  # 100 points for smooth curve
y <- curveA(x) + 10  # Apply curveA to these x-values
# x = rnorm(mean = 8,sd = 2, n = 100)
# y = curveA(x)
calpercent(x = y[10], y = y[20])
calpercent(x = (y-mean(y))[1], y= (y-mean(y))[20])
calpercent(x = minmax(y)[1], y= minmax(y)[20])
calpercent(x = y[10]*10, y= y[20]*10)

p1 = plotcurve(x = x,y = y)
p2 = plotcurve(x = minmax(x),y = minmax(y))
p1 = plotcurve(x = x,y = y-mean(y))

p3 = plotcurve(x = minmax(x^7),y = minmax(y))
p4 = plotcurve(x = x^7,y = y)
ggsave("Line1.pdf", p1 , width = 15, height = 6, dpi = 300)
ggsave("Line2.pdf", p2 , width = 15, height = 6, dpi = 300)
ggsave("Line3.pdf", p3 , width = 15, height = 6, dpi = 300)
ggsave("Line4.pdf", p4 , width = 15, height = 6, dpi = 300)


p5 = plotcurve(x = minmax(x),y = y)
ggsave("Line5.pdf", p5 , width = 15, height = 6, dpi = 300)

p6 = plotcurve(x = minmax(x^7),y = y)
ggsave("Line6.pdf", p6 , width = 15, height = 6, dpi = 300)

p7 = plotcurve(x = x^7,y = y-mean(y))
ggsave("Line7.pdf", p7 , width = 15, height = 6, dpi = 300)



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
calpercent(x = c(2.77,5.06,7,8.41), y = c(2.04,3.74,5.36,6.8))
calpercent(x = c(0.04,0.09,0.13,0.17), y = c(0.04,0.06,0.09,0.12))

calpercent(x = c(3.6,6.38,8.26,9.8), y = c(1.19,2.41,3.52,4.43))
calpercent(x = c(0.1,0.19,0.28,0.36), y = c(0.03,0.05,0.07,0.09))


#comparaing LAGM with GOCS at G5 in bp
calpercent(x = c(76.45,82.45), y = c(66.58,46.29))
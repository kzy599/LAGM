rm(list = ls())
library(data.table)
library(ggplot2)
library(ggsci)
library(dplyr)
library(multcompView)
calpercent = function(x,y){
  x = round(x,2)
  y = round(y,2)
  return(round(((x - y)/y)*100,1))
}
calpercent(x = 47.88,y =38.24)
caldeltaF_ped = function(output,opt){
  inb_dt = output[app == opt,.(inb_ped,gen)]
  setorder(inb_dt,gen)
  inb = inb_dt$inb_ped
  deltaF = c()
  for(f in 1:(length(inb)-1)){
    f_t = inb[f]
    f_tp1 = inb[f+1]
    deltaF = c(deltaF,(f_tp1-f_t)/(1-f_t))
  }
  return(round(mean(deltaF)*100,2))
}

caldeltaF = function(output,opt){
  inb_dt = output[app == opt,.(inb,gen)]
  setorder(inb_dt,gen)
  inb = inb_dt$inb
  deltaF = c()
  for(f in 1:(length(inb)-1)){
    f_t = inb[f]
    f_tp1 = inb[f+1]
    deltaF = c(deltaF,(f_tp1-f_t)/(1-f_t))
  }
  return(round(mean(deltaF)*100,2))
}

theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          #legend.title=element_blank(),
          #legend.position=c(0.5, 0.95),#图例在绘图区域的位置
          #legend.position="none",
          legend.position="right",
          #legend.direction = "horizontal",
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}

caleff = function(output,type,opt){
  out = output[app==opt,]
  
  genicVariance = out[[type]]
  
  genicall <- sqrt(out[[type]][1])
  
  genicper <- 1 - sqrt(out[[type]])/sqrt(genicall)
  
  gg0 = out$gg0
  
  ggst <- gg0/genicall
  
  dt_cv = data.table(meanGZ = ggst , sdGenicZ = genicper )
  
  fitmodel = lm(meanGZ ~ sdGenicZ, data=dt_cv )
  
  genicEff = coef(fitmodel)[2]
  
  return(genicEff)
}

calInbeff = function(output,opt){
  out = output[app==opt,]
  
  inbreeding = out$inb
  
  gg0_inb = out$gg0
  
  return(gg0_inb/inbreeding)
}

calInbeff_ped = function(output,opt){
  out = output[app==opt,]
  
  inbreeding = out$inb_ped
  
  gg0_inb = out$gg0
  
  return(gg0_inb/inbreeding)
}

calInbeff_ped_interval = function(output,opt,interval){
  out = output[app==opt,]

  gen_inb = as.numeric(out$gen)
  
  inbreeding = out$inb_ped
  
  gg0_inb = out$gg0

  # split into groups of length `interval`
  groups_inb  <- split(inbreeding, (gen_inb - 1) %/% interval)
  groups_gg0  <- split(gg0_inb,   (gen_inb - 1) %/% interval)

  # subtract the first element of each group
  groups_inb_rel <- lapply(groups_inb,  function(g) g - g[1])
  groups_gg0_rel <- lapply(groups_gg0,  function(g) g - g[1])

  # efficiency per group = ratio
  groups_eff <- mapply(
    function(gain, inb) gain / inb,
    groups_gg0_rel, groups_inb_rel,
    SIMPLIFY = FALSE
  )

  # unlist to restore generation-level vector
  efficiency <- unlist(groups_eff, use.names = FALSE)
  
  efficiency[is.na(efficiency)] = 0
  
  return(efficiency)
}

extractMeanSd = function(dt){
  result_mean <- dt[, lapply(.SD, mean), by = .(app,gen)]
  
  result_sd <- dt[, lapply(.SD, sd), by = .(app,gen)]
  
  
  bool = (!colnames(result_sd)%in%c("app","gen"))
  
  colnames(result_sd)[bool] = paste(colnames(result_sd)[bool],"_sd",sep = "")
  
  result <- merge(result_mean, result_sd, by = c("app","gen"))
  
  return(result)
}


#==========================================Begaining==========================================================================
main_dir = "/home/kangziyi/RLmating/"
# dir = c("gas20")
dir = c("gas20","Gen20")
alldir = paste(main_dir,dir,sep = "")
gen = c(5,10,15,20)
extractMeanSd = function(dt){

  result_mean <- dt[, lapply(.SD, mean), by = .(app,gen)]
  
  result_sd <- dt[, lapply(.SD, sd), by = .(app,gen)]
  
  
  bool = (!colnames(result_sd)%in%c("app","gen"))
  
  colnames(result_sd)[bool] = paste(colnames(result_sd)[bool],"_sd",sep = "")
  
  result <- merge(result_mean, result_sd, by = c("app","gen"))
  
  return(result)
}

plotDT = function(pv,pvse,dt){
output = data.table()
for(t in pv){
  locus = which(pv == t)
  sd_temp = pvse[locus]
  output_temp = data.table(app = dt$app,gen = dt$gen ,type = t,value =dt[[t]],valuese = dt[[sd_temp]])
  colnames(output_temp) = c("app","gen","type","value","valuese")
  output= rbind(output,output_temp)
}
return(output)
}

DT_test = function(rep_out,Vtype,gen,ck_dt){
output = data.table()
for(b in gen){
for(g in Vtype){
  aov_res <- aov(as.formula(paste(g, "~ app")), data = rep_out[gen == b,])
  tukey_res <- TukeyHSD(aov_res,conf.level = 0.95)
  letters <- multcompLetters4(aov_res, tukey_res)
  output = rbind(output,data.table(app = names(letters$app$Letters),Letter = letters$app$Letters,type = g,gen = b))
}
}
merged_dt <- merge(output, ck_dt, by = c("app", "gen","type"), all = FALSE)
return(merged_dt)
}
rename_dt = function(zt1){
zt1 = zt1[app!="rl",]
zt1[app=="rl2",app:="LAGM7"]
zt1[app=="rl2_3",app:="LAGM3"]
zt1[app=="rl2_5",app:="LAGM5"]
zt1[app=="rl2_5m",app:="LAGM1"]
zt1[app=="rl2_10m",app:="LAGM2"]
zt1[app=="tc",app:="TC"]
zt1[app=="ocs25",app:="OCS25"]
zt1[app=="ocs45",app:="OCS45"]
zt1[app=="ocs65",app:="OCS65"]
zt1[app=="ocs90",app:="OCS90"]
zt1[app=="ocsrate",app:="OCSrate"]
zt1[app=="ran",app:="Random"]
zt1[type=="gg0",type:="Gain"]
zt1[type=="Ne",type:="EffectiveSize"]
zt1[type=="inb_ped",type:="Inbreeding"]
zt1[type=="inbEff_ped",type:="Efficiency"]
zt1[type=="inbEff_ped_5",type:="Efficiency_5"]
zt1[type=="inbEff_ped_10",type:="Efficiency_10"]
return(zt1)
}
rep_out = data.table()
for(d in alldir){

setwd(d)

if(dir[which(alldir == d)] == "Gen20") nrep = c(1:20)
if(dir[which(alldir == d)] == "gas20") nrep = c(1:20) # c(21:30,41:50) c(40:60)
for(r in nrep){

output = fread(paste("output",r,".csv",sep = ""),sep = ",")

#command + shift + c
# output[,Index:=NULL]

colnames(output) = c("app","gg0","Va","genicVa","inb","Ne","inb_ped","gen")

if(d == "/home/kangziyi/RLmating/Gen20") output = output[!app%in%c("rl2","rl"),]

opt = unique(output$app)

for(o in opt){
vaeff = caleff(output = output, type = "Va",opt = o)[[1]]
geniceff = caleff(output = output, type = "genicVa",opt = o)[[1]]
deltaF = caldeltaF(output = output, opt = o)
deltaF_ped = caldeltaF_ped(output = output, opt = o)
inbeff = calInbeff(output = output, opt = o)
inbeff_ped = calInbeff_ped(output = output, opt = o)

inbeff_ped_5 = calInbeff_ped_interval(output = output, opt = o,interval = 5)
inbeff_ped_10 = calInbeff_ped_interval(output = output, opt = o,interval = 10)
output[app == o,inbEff_ped_5:=inbeff_ped_5]
output[app == o,inbEff_ped_10:=inbeff_ped_10]

output[app == o,vaEff:=vaeff]
output[app == o,genicEff:=geniceff]
output[app == o,inbrate:=deltaF]
output[app == o,inbrate_ped:=deltaF_ped]
output[app == o,inbEff:=inbeff]
output[app == o,inbEff_ped:=inbeff_ped]
}
rep_out = rbind(rep_out,output)

}

}

setwd("/home/kangziyi/RLmating/gas20")

pv = c("gg0","inb_ped","Ne","inbEff_ped","inbEff_ped_5","inbEff_ped_10")
pvse = paste(pv,"_sd",sep = "")

dt = extractMeanSd(rep_out)

ck_dt_line = plotDT(pv = c("gg0","inb_ped","Ne","inbEff_ped","inbEff_ped_5","inbEff_ped_10"),pvse = pvse,dt = dt)

ck_dt = DT_test(rep_out,Vtype = unique(ck_dt_line$type),gen = gen , ck_dt = ck_dt_line[gen%in%c(5,10,15,20),])

ck_dt = rename_dt(ck_dt)

ck_dt$app = factor(ck_dt$app,level = c("LAGM7","LAGM5","LAGM3","LAGM2","LAGM1","TC","OCSrate","OCS25","OCS45","OCS65","OCS90","Random"))
ck_dt$BP = paste(ck_dt$gen,"-generation",sep = "")
ck_dt$BP = factor(ck_dt$BP, level = c("5-generation","10-generation","15-generation","20-generation"))
                           
y_lab_name = "Inbreeding"#efficiency_interval meaningless
ck_dt$type = factor(ck_dt$type,level = c("Gain","Inbreeding","EffectiveSize","Efficiency","Efficiency_5","Efficiency_10"))
if(y_lab_name == "Gain"){
ck_dt[, letter_y_test := value + 3 * valuese, by = type]
ck_dt[, letter_y := value + 1.1 * valuese, by = type]
}else if(y_lab_name == "Inbreeding"){
ck_dt[, letter_y_test := value + 6 * valuese, by = type]
ck_dt[, letter_y := value + 1.1 * valuese, by = type]
ck_dt[app == "TC", letter_y_test := value + 3 * valuese]
}else{
ck_dt[, letter_y_test := value + 3 * valuese, by = type]
ck_dt[, letter_y := value + 1.1 * valuese, by = type]
# ck_dt[app =="GAS"&BP == "5-generations", letter_y_test := value + 2 * valuese, by = type]
}

base_cols <- ggsci::pal_lancet()(9)  # Lancet 原始 9 色
my_cols   <- grDevices::colorRampPalette(base_cols)(length(unique(ck_dt$app)))
names(my_cols) <- unique(ck_dt$app)
P <- ggplot(data = ck_dt[type == y_lab_name,], aes(x = app, y = value, group = app, fill = app)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠，position_dodge里的width 要与外面的width相等
  geom_text(aes(label = round(value,2),y = letter_y), position = position_dodge(width = 0.9),vjust = 0, size = 5) +  # 添加数值标签
  xlab("Mating strategies") +
  ylab("Value") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank(),
    #设置分面标题，要在facet_wrap之前用
     strip.text = element_text(
      size = 14,           # 字体大小
      face = "bold",       # 加粗（"plain", "italic", "bold.italic"）
      color = "black"      # 颜色
    )) +
  geom_errorbar(aes(ymin = value - valuese, ymax = value + valuese),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  # scale_fill_lancet()+
  scale_fill_manual(values = my_cols, drop = FALSE)+
  #添加position = position_dodge(width = 0.9)如果一个x分类内有多个组，要确保组不重叠，width 要与geom_bar里面的width相等
  geom_text(data = ck_dt[type == y_lab_name,], aes(label = Letter, y = letter_y_test), vjust = 0,size = 5) + facet_wrap(~BP, scales = "free_x")+
  ggtitle("Inbreeding across mating strategies in various breeding programs") +
  theme(
    plot.title = element_text(
      size = 20,          # 大字号
      face = "bold",      # 加粗
      hjust = 0.5,        # 居中
      color = "black"     # 颜色（可选）
    ),
    # 调整x轴标签字体
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    # 调整y轴标签字体
    axis.title.y = element_text(size = 16, face = "bold", color = "black")
  )

ggsave(paste("Figrue_bar_",y_lab_name,".pdf",sep = ""), P , width = 22, height = 10, dpi = 300)


ck_dt_line = rename_dt(ck_dt_line)
nG = 20
c("LAGM1","LAGM2","LAGM3","LAGM5","LAGM7","OCS45","OCS25")
P<- ggplot(data = ck_dt_line[type%in%c("Gain","Inbreeding")&app%in%c("LAGM2","LAGM7","OCS45","OCS65"),],aes(x=gen,y=value,group = app, color = app))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("Generation")+
  ylab("")+
  theme_zg()+theme(legend.title=element_blank())+
  #geom_ribbon(aes(ymin = value -valuese, ymax = value+valuese), alpha = 0.4)+labs(fill="Mat")+scale_fill_aaas()
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=value-valuese,
                    ymax=value+valuese),
                width=0.05,alpha = 0.5)+labs(color="Mat")+scale_color_lancet()+scale_x_continuous(limits = c(1,nG),breaks = seq(1,nG,4)) + facet_wrap(~type, scales = "free")+
ggtitle("Genetic gain and Inbreeding across 20 generations") +
  theme(
    plot.title = element_text(
      size = 20,          # 大字号
      face = "bold",      # 加粗
      hjust = 0.5,        # 居中
      color = "black"     # 颜色（可选）
    ),
    # 调整x轴标签字体
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    # 调整y轴标签字体
    axis.title.y = element_text(size = 16, face = "bold", color = "black")
  )
ggsave("Figrue_gain_inbreeding_line_27.pdf", P , width = 15, height = 6, dpi = 300)


y_lab_name = "Gain"
dt = ck_dt[type==y_lab_name,.(app,BP,value,valuese)]
# dt[, col_id := paste(type, dataset, sep = "_")]
dt[, value_display := sprintf("%.2f ± %.2f", value, valuese)]

wide_dt <- dcast(dt, app ~ BP, value.var = "value_display")

# ordered_cols = c("MODE")
# for(i in c("ST","MT_continuous","MT_binary")){
#    ordered_cols = c(ordered_cols,paste(i,c("Linear-A","Low-ADE","Med-ADE","Hig-ADE"),sep = '_'))
# }
ordered_cols = c("app","5-generations","10-generations","15-generations","20-generations")

wide_dt <- wide_dt[, ..ordered_cols]
colnames(wide_dt)[1] = "Mating_strategies" 
# Adjust header map accordingly
# header_map <- data.frame(
#   col_keys = ordered_cols,
#   Type = c("Model", rep(c("ST", "MT_continuous", "MT_binary"), each = 4)),
#   Datasets = c(" ", rep(c("Linear-A","Low-ADE","Med-ADE","Hig-ADE"), 3)),
#   stringsAsFactors = FALSE
# )

ft <- flextable(wide_dt)
# ft <- set_header_df(ft, mapping = header_map, key = "col_keys")
# ft <- merge_h(ft, part = "header")  # merge group headers
# ft <- merge_v(ft, j = "MODE", part = "body")  # optional: merge MODE column
ft <- theme_booktabs(ft)
ft <- autofit(ft)
ft <- fontsize(ft, size = 9, part = "all")        # smaller font
ft <- width(ft, width = 0.5)   

save_as_html(ft, path = "strategies.html")

library(officer)

doc <- read_docx() |> 
  body_add_par(paste(y_lab_name,"accuracy across models and simulated datasets",sep = " "), style = "heading 1") |> 
  body_add_flextable(ft)

print(doc, target = "strategies.docx")





setorder(ck_dt,BP)
y_lab_name = "Gain"
x = ck_dt[app == "OCS25"&type == y_lab_name,value]
y = ck_dt[app == "OCS90"&type == y_lab_name,value]

calpercent(x,y)

bp = "20-generation"
x = ck_dt[app == "GAS"&type == y_lab_name&BP == bp,value]
y = ck_dt[app == "OCS90"&type == y_lab_name&BP == bp,value]

calpercent(x,y)

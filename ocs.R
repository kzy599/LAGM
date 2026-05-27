G_A_ibd <- function(pop,pedmore){

    AnimalID <- pop@id
    
    ped_need = data.table(id = pop@id, sir = pop@father, dam = pop@mother)

    colnames(ped_need) = c("id","sir","dam")

    pedfinal <-rbind(pedmore,ped_need)
    pedfinal <-visPedigree:: tidyped(pedfinal,AnimalID)
    pedfinal<- pedfinal[,.(Ind,Sire,Dam)]

    Amatrix <- makeA(pedfinal)
    
    Amatrix <- as.matrix(Amatrix)

    tempBool_row = match(AnimalID,rownames(Amatrix))
    tempBool_col = match(AnimalID,colnames(Amatrix))
    
    AA<-Amatrix[tempBool_row,tempBool_col]
    
    AA <- cbind(pop@id,AA)
    
    colnames(AA) <- NULL
    
    rownames(AA) <- NULL
    
    return(AA)
}
get_os = function(){
  sysinf = Sys.info()
  if(!is.null(sysinf)){
    os = sysinf[['sysname']]
  }else{ # Sys.info not set up
    os = .Platform$OS.type
    if(grepl("^darwin", R.version$os)){
      os = "darwin"
    }else if(grepl("linux-gnu", R.version$os)){
      os = "linux"
    }
  }
  os = tolower(os)
  stopifnot(os=="windows" | os=="darwin" | os=="linux")
  return(os)
}

runProgram = function(progName, specFileName, runDir = "runDir") {
  #set up files
  dir.create(runDir,recursive = FALSE,showWarnings = FALSE)
  spec_file_path = file.path(runDir, specFileName)
  file.copy(specFileName, spec_file_path, overwrite = TRUE, recursive=FALSE)
  os = get_os()
  runCommand = paste(progName," ",specFileName,sep="")
  old_dir = getwd()
  setwd(runDir)
  system(runCommand)
  setwd(old_dir)
}

runAlphaImpute = function(pedigree, genotypes) {
  runDir = "runDir"
  dir.create(runDir,recursive = FALSE,showWarnings = FALSE)
  ped_path = file.path(runDir,"pedigree.txt")
  write.table(pedigree,file=ped_path, row.names=F, col.names=F, quote=F)
  write.table(cbind(pedigree[,1], genotypes),file=paste(runDir,"/","genotypes.txt",sep=""), row.names=F, col.names=F, quote=F)
  
  spec <-ImputeParam$new()
  spec$write_out_spec(runDir)
  runProgram("AlphaImpute", "AlphaImputeSpec.txt",runDir)
  return(readInAI(path=runDir))
}

fixGenotypes = function(genotypes,pop) {
  ids = genotypes[,1]
  newLoci = sapply(pop@id, function(id) {
    which(ids == id)
  })
  genotypes = genotypes[newLoci,-1]
  return(genotypes)
}

readInAI = function(path=".") {
  genotypes <- as.matrix(fread(paste(path,"/Results/ImputeGenotypes.txt",sep="")))
  phase <- as.matrix(fread(paste(path,"/Results/ImputePhase.txt",sep="")))
  gdosages <- as.matrix(fread(paste(path,"/Results/ImputeGenotypeProbabilities.txt",sep="")))
  return(list(genotypes = genotypes, phase = phase, gdosages=gdosages))
}

#' @title Optimal Contribution Selection
#'
#' @description Perform Optimal Contribution Selection via AlphaMate
#' .
#' @param pop population
#' @param nCrosses number of matings/crosses
#' @param nFemalesMax maximum number of females
#' @param nMalesMax maximum number of   males
#' @param minFemaleContribution minimum number of matings/crosses per female
#' @param maxFemaleContribution maximum number of matings/crosses per female
#' @param minMaleContribution minimum number of matings/crosses per male
#' @param maxMaleContribution maximum number of matings/crosses per male
#' @param targetDegree targeted trigonometric degrees between genetic gain and group coancestry
#' @param targetCoancestryRate targeted rate of group coancestry
#' @param nProgenyPerCross number of progeny per mating/cross
#' @param use character specifiying, which type of criterion to use, either "pheno", "ebv", or "gv"
#'
#' @export
ocs <- function(pop, nCrosses, usePed = FALSE, ped_parents = NULL,
                nFemalesMax = NULL, minFemaleContribution = NULL, maxFemaleContribution = NULL, equalizeFemaleContributions = NULL,
                nMalesMax   = NULL, minMaleContribution   = NULL, maxMaleContribution   = NULL, equalizeMaleContributions   = NULL,
                targetDegree = NULL, targetCoancestryRate = NULL,
                nProgenyPerCross, use,
                MateAllocation       = TRUE,
                minInbreedingMating  = FALSE) {   # ← 新增

  # ---- Prepare data ----
  runDir <- "runDir"
  dir.create(path = runDir, showWarnings = FALSE)

  # Criterion
  if (use == "gv")    tmp <- data.frame(id = pop@id, crit = pop@gv)
  if (use == "ebv")   tmp <- data.frame(id = pop@id, crit = pop@ebv)
  if (use == "pheno") tmp <- data.frame(id = pop@id, crit = pop@pheno)
  write.table(x = tmp, file = file.path(runDir, "SelCriterion.txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)

  # NRM
  if (usePed) {
    if (is.null(ped_parents))
      stop("The ped of parents can't be NULL when usePed to perform OCS")
    Amatrix <- G_A_ibd(pop = pop, pedmore = ped_parents)
    MASS::write.matrix(x = Amatrix, file = file.path(runDir, "GenotypeNrm.txt"))
  } else {
    G <- AlphaMME::calcGIbs(X = pullSnpGeno(pop))
    MASS::write.matrix(x = cbind(pop@id, G),
                       file = file.path(runDir, "GenotypeNrm.txt"))
  }

  # Gender
  tmp <- data.frame(id = pop@id, genderRole = pop@sex)
  tmp$genderRole <- as.numeric(factor(pop@sex, levels = c("M", "F")))
  write.table(x = tmp, file = file.path(runDir, "genderRole.txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)

  # ---- Spec file (Stage 1: OCS) ----
  sink(file = "AlphaMateSpec.txt")
  cat("NrmMatrixFile , GenotypeNrm.txt\n")
  if (use %in% c("gv", "ebv", "pheno")) {
    cat("SelCriterionFile , SelCriterion.txt\n")
  } else {
    stop("use muste be gv, ebv, or pheno")
  }
  cat("GenderFile , genderRole.txt\n")
  cat("NumberOfMatings , ", nCrosses, "\n")
  if (!is.null(nFemalesMax)) cat("NumberOfFemaleParents , ", nFemalesMax, "\n")
  if (!is.null(maxFemaleContribution) | !is.null(minFemaleContribution)) {
    cat("LimitFemaleContributions , Yes\n")
    if (!is.null(minFemaleContribution)) cat("LimitFemaleContributionsMin , ", minFemaleContribution, "\n")
    if (!is.null(maxFemaleContribution)) cat("LimitFemaleContributionsMax , ", maxFemaleContribution, "\n")
  }
  if (!is.null(equalizeFemaleContributions)) cat("EqualizeFemaleContributions , Yes\n")
  if (!is.null(nMalesMax)) cat("NumberOfMaleParents , ", nMalesMax, "\n")
  if (!is.null(maxMaleContribution) | !is.null(minMaleContribution)) {
    cat("LimitMaleContributions , Yes\n")
    if (!is.null(minMaleContribution)) cat("LimitMaleContributionsMin , ", minMaleContribution, "\n")
    if (!is.null(maxMaleContribution)) cat("LimitMaleContributionsMax , ", maxMaleContribution, "\n")
  }
  if (!is.null(equalizeMaleContributions)) cat("EqualizeMaleContributions , Yes\n")
  if (!is.null(targetCoancestryRate))      cat("TargetCoancestryRate , ", targetCoancestryRate, "\n")
  if (!is.null(targetDegree))              cat("TargetDegree , ", targetDegree, "\n")

  # 如果要做 min-inbreeding mate allocation，stage 1 必须先关掉 MateAllocation
  # 让 AlphaMate 只输出 Contributors；否则按用户原意写
  if (!MateAllocation || minInbreedingMating) cat("MateAllocation , No\n")

  cat("Stop\n")
  sink()
  runProgram("AlphaMate", "AlphaMateSpec.txt", runDir = runDir)

  # ---- 取结果 ----
  if (MateAllocation) {

    if (!minInbreedingMating) {
      # ===== 原行为：直接读 OCS+mate-allocation 的配种表 =====
      crossPlan <- read.table(file = file.path(runDir, "MatingPlanModeOptTarget1.txt"),
                              header = TRUE, colClasses = "character")
      pedigree <- matrix(data = "", ncol = 2, nrow = nrow(crossPlan) * nProgenyPerCross)
      k <- 0
      for (i in 1:nrow(crossPlan)) {
        for (j in 1:nProgenyPerCross) {
          k <- k + 1
          pedigree[k, 1] <- crossPlan[i, 3]   # dam
          pedigree[k, 2] <- crossPlan[i, 2]   # sire
        }
      }
      return(makeCross(pop = pop, crossPlan = pedigree))

    } else {
      # ===== 新分支：Stage 2 — 固定贡献，最小化父母平均亲缘 =====
      contrib <- read.table(file = file.path(runDir, "ContributorsModeOptTarget1.txt"),
                            header = TRUE, colClasses = "character")
      contrib <- contrib[as.integer(contrib$nContribution) > 0, , drop = FALSE]

      # 只保留被选个体，准备 stage-2 的小型输入文件
      selIds <- contrib$Id
      selPop <- pop[pop@id %in% selIds]

      # SelCriterion (stage 2 不再使用，但 AlphaMate 仍要求文件存在)
      if (use == "gv")    tmp <- data.frame(id = selPop@id, crit = selPop@gv)
      if (use == "ebv")   tmp <- data.frame(id = selPop@id, crit = selPop@ebv)
      if (use == "pheno") tmp <- data.frame(id = selPop@id, crit = selPop@pheno)
      write.table(x = tmp, file = file.path(runDir, "SelCriterion.txt"),
                  col.names = FALSE, row.names = FALSE, quote = FALSE)

      # NRM (与 stage 1 同源；用 selPop 重算，保证 ID 顺序对齐)
      if (usePed) {
        Amatrix <- G_A_ibd(pop = selPop, pedmore = ped_parents)
        MASS::write.matrix(x = Amatrix, file = file.path(runDir, "GenotypeNrm.txt"))
      } else {
        G <- AlphaMME::calcGIbs(X = pullSnpGeno(selPop))
        MASS::write.matrix(x = cbind(selPop@id, G),
                           file = file.path(runDir, "GenotypeNrm.txt"))
      }

      # Gender
      tmp <- data.frame(id = selPop@id, genderRole = selPop@sex)
      tmp$genderRole <- as.numeric(factor(selPop@sex, levels = c("M", "F")))
      write.table(x = tmp, file = file.path(runDir, "genderRole.txt"),
                  col.names = FALSE, row.names = FALSE, quote = FALSE)

      nMales   <- sum(contrib$Gender == "1")
      nFemales <- sum(contrib$Gender == "2")

      # ---- Spec file (Stage 2: ModeMinInbreeding) ----
      sink(file = "AlphaMateSpec.txt")
      cat("NrmMatrixFile        , GenotypeNrm.txt\n")
      cat("SelCriterionFile     , SelCriterion.txt\n")
      cat("GenderFile           , genderRole.txt\n")
      cat("NumberOfMatings      , ", nCrosses, "\n")
      cat("NumberOfMaleParents  , ", nMales,   "\n")
      cat("NumberOfFemaleParents, ", nFemales, "\n")
      if (!is.null(equalizeFemaleContributions)) cat("EqualizeFemaleContributions , Yes\n")
      if (!is.null(equalizeMaleContributions)) cat("EqualizeMaleContributions , Yes\n")
      cat("AllowSelfing         , No\n")
      cat("AllowRepeatedMatings , No\n")
      cat("SelfingWeight        , -1000\n")
      cat("RepeatedMatingsWeight, -1000\n")
      cat("ModeMinInbreeding    , Yes\n")          # 目标：min ΔF
      cat("MateAllocation       , Yes\n")
      cat("Stop\n")
      sink()
      runProgram("AlphaMate", "AlphaMateSpec.txt", runDir = runDir)

      crossPlan <- read.table(file = file.path(runDir, "MatingPlanModeMinInbreeding.txt"),
                              header = TRUE, colClasses = "character")
      pedigree <- matrix("", nrow = nrow(crossPlan) * nProgenyPerCross, ncol = 2)
      k <- 0
      for (i in 1:nrow(crossPlan)) {
        for (j in 1:nProgenyPerCross) {
          k <- k + 1
          pedigree[k, 1] <- crossPlan[i, 3]   # dam
          pedigree[k, 2] <- crossPlan[i, 2]   # sire
        }
      }
      return(makeCross(pop = pop, crossPlan = pedigree))
    }

  } else {
    # ===== 原行为：不做 mate allocation，直接返回 contributors =====
    crossPlan <- read.table(file = file.path(runDir, "ContributorsModeOptTarget1.txt"),
                            header = TRUE, colClasses = "character")
    return(crossPlan)
  }
}

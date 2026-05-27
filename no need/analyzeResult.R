rm(list = ls())
library(data.table)
library(ggplot2)
library(MASS)      
library(multcomp)  
library(emmeans)  
library(tidyverse)
library(broom)
library(rlang)
library(patchwork)
packages <- c("ggplot2", "dplyr", "emmeans")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

rep = data.table()
programs = c("5_generation","10_generation","15_generation","20_generation","nhp_20_generation")
dir = c("/home/kangziyi/RLmating/Gen5","/home/kangziyi/RLmating/Gen10","/home/kangziyi/RLmating/Gen15","/home/kangziyi/RLmating/Gen20","/home/kangziyi/RLmating/gas20")
original_dir = getwd()
for(p in c(1:5)){
  setwd(dir[p])
for(i in c(1:20)){
  temp = fread(paste0("output",i,".csv",sep=""))
  temp[,rep:=i]
  gv0 = unique(temp[gen==0,gv])
  temp$gain = temp$gv-gv0
  He_diversity0 = unique(temp[gen==0,He])
  Ho_diversity0 = unique(temp[gen==0,Ho])
  genetic_diversity0 = unique(temp[gen==0,genetic])
  genic_diversity0 = unique(temp[gen==0,genic])
  temp$lost_diversity_He = (1-temp$He/He_diversity0)
  temp$lost_diversity_Ho = (1-temp$Ho/Ho_diversity0)
  temp$lost_diversity_genetic = (1-sqrt(temp$genetic)/sqrt(genetic_diversity0))
  temp$lost_diversity_genic = (1-sqrt(temp$genic)/sqrt(genic_diversity0))
  temp[,programs:=programs[p]]
  rep = rbind(rep,temp)
}
}
setwd(original_dir)

fwrite(rep,"rep_all.csv",sep = ",")


# ---- 读取数据 ----
dat <- read.csv("rep_all.csv")

#bp programs
dat_full = dat[dat$programs!="nhp_20_generation",]
dat_full <- as.data.table(dat_full)
dat_full <- dat_full[app%in%c("ocs25","ocs45","ocs65","rl2","ocs90","tc","ran"),]


#nbp_programs
dat_full = dat[!dat$programs%in%c("5_generation","10_generation","15_generation"),]
dat_full = dat_full[!(dat_full$app%in%c("rl2","rl")&dat_full$programs=="20_generation"),]
dat_full$programs = "20_generation"
dat_full = as.data.table(dat_full)
dat_full <- dat_full[app%in%c("ocs25","ocs45","ocs65","ocs90","rl2")|app%flike%"rl2",]


#' 每 rep 用稳健回归拟合 y ~ x,在固定 x 值处预测 y,跨 rep 画 boxplot
#'
#' @param dat        data.table,需含 rep, app, programs, gen, x_var, y_var
#' @param x_var      x 列名(如 "inb_ped" 或 "inb")
#' @param y_var      y 列名(如 "gain")
#' @param x_values   预测点(绝对 x 值),如 c(0.05, 0.10, 0.15, 0.20)
#' @param progs      筛选 programs(NULL 不筛)
#' @param apps       筛选 apps(NULL 不筛)
#' @param skip_gen   拟合时跳过的前几代(默认 0,即用全部)
#' @param method     "rlm"(MASS,M-estimator,默认)
#'                   "lmrob"(robustbase,MM-estimator,更稳)
#'                   "ols"(普通最小二乘,作对照)
#' @param extrapolate 是否允许外推到 rep 的 x 范围之外(默认 FALSE)
#' @param add_test   是否标注 LAGM vs OCS 的显著性检验
#' @param ref_app    显著性检验的参考方法(默认 "rl2")
#' @param colors     命名色板
#' @return list(plot, predictions, fits)
plot_robust_predicted_boxplot <- function(
    dat,
    x_var            = "inb",
    y_var            = "gain",
    x_values         = c(0.05, 0.10, 0.15, 0.20),
    progs            = NULL,
    apps             = NULL,
    skip_gen         = 0,
    method           = c("rlm", "lmrob", "ols"),
    extrapolate      = FALSE,
    test_method      = c("none", "ref", "tukey"),
    ref_app          = "rl2",
    tukey_alpha      = 0.05,
    use_mixed        = FALSE,
    min_reps_per_box = 3,
    show_labels      = c("both", "letter", "value", "none"),
    label_digits     = 2,
    label_size       = 3.2,
    label_pad_frac   = 0.02,
    show_r2          = c("subtitle", "caption", "none"),
    facet_scales     = "fixed",
    colors           = NULL
) {
    # ---------- 0. 参数校验 ----------
    method      <- match.arg(method)
    test_method <- match.arg(test_method)
    show_labels <- match.arg(show_labels)
    show_r2     <- match.arg(show_r2)

    stopifnot(data.table::is.data.table(dat))
    needed <- c("app", "programs", "rep", "gen", x_var, y_var)
    miss   <- setdiff(needed, names(dat))
    if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

    d <- data.table::copy(dat)
    if (!is.null(progs)) d <- d[programs %in% progs]
    if (!is.null(apps))  d <- d[app %in% apps]
    d <- d[gen >= skip_gen]
    d <- d[!is.na(get(x_var)) & !is.na(get(y_var))]
    if (!nrow(d)) stop("No data left after filtering.")
    setorder(d, app, programs, rep, gen)

    # ---------- 1. 拟合 ----------
    .fit_one <- function(xv, yv) {
        if (length(xv) < 3 || length(unique(xv)) < 2) return(NULL)
        fit <- switch(method,
            rlm   = tryCatch(MASS::rlm(yv ~ xv, maxit = 200), error = function(e) NULL),
            lmrob = tryCatch(robustbase::lmrob(yv ~ xv),    error = function(e) NULL),
            ols   = tryCatch(stats::lm(yv ~ xv),            error = function(e) NULL)
        )
        if (is.null(fit)) return(NULL)
        yhat   <- as.numeric(stats::predict(fit))
        ss_res <- sum((yv - yhat)^2)
        ss_tot <- sum((yv - mean(yv))^2)
        list(intercept = unname(stats::coef(fit)[1]),
             slope     = unname(stats::coef(fit)[2]),
             r2        = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_,
             x_min     = min(xv), x_max = max(xv))
    }

    rep_fits <- d[, {
        res <- .fit_one(get(x_var), get(y_var))
        if (is.null(res))
            list(intercept = NA_real_, slope = NA_real_, r2 = NA_real_,
                 x_min = NA_real_, x_max = NA_real_)
        else res
    }, by = .(app, programs, rep)]
    rep_fits <- rep_fits[!is.na(slope)]
    if (!nrow(rep_fits)) stop("No successful fits.")

    # ---------- 2. 预测 ----------
    pred_dt <- rep_fits[, {
        preds <- intercept + slope * x_values
        if (!extrapolate) preds[x_values < x_min | x_values > x_max] <- NA_real_
        list(x_target = x_values, y_pred = preds)
    }, by = .(app, programs, rep)]
    pred_dt <- pred_dt[!is.na(y_pred)]
    if (!nrow(pred_dt)) stop("No predictions in range. Try extrapolate = TRUE.")

    pred_dt[, n_box := .N, by = .(programs, x_target, app)]
    pred_dt <- pred_dt[n_box >= min_reps_per_box]
    pred_dt[, n_box := NULL]

    pred_dt[, x_target_factor := factor(round(x_target, 4),
                                    levels = sort(unique(round(x_target, 4))))]
    pred_dt[, app := factor(app, levels = levels(factor(d$app)))]

    # ---------- 3. R² 汇总 ----------
    r2_summary <- rep_fits[, .(R2_mean = mean(r2, na.rm = TRUE),
                               R2_sd   = sd(r2,   na.rm = TRUE),
                               n_rep   = .N), by = .(app, programs)]

    # ---------- 4. 显著性检验 ----------
    .run_tukey <- function(sub) {
        if (data.table::uniqueN(sub$app) < 2) return(NULL)
        sub[, app := factor(app)]
        fit <- tryCatch(stats::aov(y_pred ~ app, data = sub), error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        tk <- tryCatch(stats::TukeyHSD(fit, "app", conf.level = 1 - tukey_alpha),
                       error = function(e) NULL)
        if (is.null(tk)) return(NULL)
        cld <- multcompView::multcompLetters(tk$app[, "p adj"],
                                             threshold = tukey_alpha)$Letters
        data.table(app = names(cld), signif_label = unname(cld))
    }

    .run_dunnett <- function(sub) {
        if (data.table::uniqueN(sub$app) < 2 || !ref_app %in% sub$app) return(NULL)
        sub[, app := stats::relevel(factor(app), ref = ref_app)]
        fit <- tryCatch({
            if (use_mixed && data.table::uniqueN(sub$rep) > 1)
                suppressMessages(lme4::lmer(y_pred ~ app + (1|rep), data = sub, REML = TRUE))
            else
                stats::aov(y_pred ~ app, data = sub)
        }, error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        dun <- tryCatch(summary(multcomp::glht(fit, linfct = multcomp::mcp(app = "Dunnett"))),
                        error = function(e) NULL)
        if (is.null(dun)) return(NULL)
        ps  <- as.numeric(dun$test$pvalues)
        nms <- gsub(paste0(" - ", ref_app, ".*$"), "",
                    names(dun$test$coefficients))
        stars <- cut(ps, c(-Inf, .001, .01, .05, Inf), c("***","**","*",""))
        rbind(data.table(app = nms, signif_label = as.character(stars)),
              data.table(app = ref_app, signif_label = "(ref)"))
    }

    if (test_method == "tukey" && !requireNamespace("multcompView", quietly = TRUE))
        stop("Please install 'multcompView'.")
    if (test_method == "ref") {
        if (!requireNamespace("multcomp", quietly = TRUE)) stop("Please install 'multcomp'.")
        if (use_mixed && !requireNamespace("lme4", quietly = TRUE)) stop("Please install 'lme4'.")
        if (!ref_app %in% pred_dt$app) stop("ref_app '", ref_app, "' not found in data.")
    }

    sig_dt <- if (test_method == "none") {
        NULL
    } else {
        pred_dt[, {
            res <- if (test_method == "tukey") .run_tukey(copy(.SD))
                   else                        .run_dunnett(copy(.SD))
            if (is.null(res)) data.table(app = unique(.SD$app), signif_label = "")
            else              res
        }, by = .(programs, x_target_factor)]
    }

    # ---------- 5. 标签 ----------
    .upper_whisker <- function(x) {
        s <- grDevices::boxplot.stats(x)$stats
        if (length(s)) max(s) else max(x, na.rm = TRUE)
    }
    box_stats <- pred_dt[, .(y_mean  = mean(y_pred, na.rm = TRUE),
                             y_whisk = .upper_whisker(y_pred),
                             n       = .N),
                         by = .(programs, x_target_factor, app)]

    y_pad    <- diff(range(pred_dt$y_pred, na.rm = TRUE)) * label_pad_frac
    label_dt <- copy(box_stats)
    label_dt[, y_top := y_whisk + y_pad]

    if (!is.null(sig_dt)) {
        label_dt <- merge(label_dt, sig_dt,
                          by = c("programs", "x_target_factor", "app"),
                          all.x = TRUE)
        label_dt[is.na(signif_label), signif_label := ""]
    } else {
        label_dt[, signif_label := ""]
    }

    fmt_val <- function(v) formatC(v, digits = label_digits, format = "f")
    label_dt[, value_label := fmt_val(y_mean)]
    label_dt[, text_label := switch(show_labels,
        both   = ifelse(nzchar(signif_label),
                        paste0(signif_label, "\n", value_label), value_label),
        letter = signif_label,
        value  = value_label,
        none   = ""
    )]
    label_dt[, app := factor(app, levels = levels(pred_dt$app))]

    # ---------- 6. 颜色 ----------
    all_apps <- levels(factor(d$app))
    if (is.null(colors)) {
        colors <- setNames(grDevices::hcl.colors(length(all_apps), "Set 2"), all_apps)
    } else {
        miss_c <- setdiff(all_apps, names(colors))
        if (length(miss_c))
            warning("colors missing for: ", paste(miss_c, collapse = ", "),
                    ". Check that factor labels match the names of `colors`.")
    }

    # ---------- 7. R² 文本 ----------
    r2_txt <- NULL
    if (show_r2 != "none") {
        r2_one <- r2_summary[, .(txt = paste0(
            app, ": R²=", formatC(R2_mean, digits = 3, format = "f"),
            "±",          formatC(R2_sd,   digits = 3, format = "f"))),
            by = programs]
        r2_txt <- if (data.table::uniqueN(r2_one$programs) == 1) {
            paste(r2_one$txt, collapse = "  |  ")
        } else {
            paste(vapply(unique(r2_one$programs), function(p)
                paste0("[", p, "] ", paste(r2_one[programs == p, txt], collapse = " | ")),
                character(1)), collapse = "\n")
        }
    }

    # ---------- 8. 画图 (嵌套分面: programs × inbreeding level) ----------
    caption_txt <- paste0(
        "Per-rep robust fit (", method, ")",
        if (skip_gen > 0)             paste0("; gen ≥ ", skip_gen)        else "",
        if (!extrapolate)             "; no extrapolation"                else "",
        "; n_rep = ", data.table::uniqueN(d$rep),
        if (test_method == "ref")     paste0("; Dunnett vs ", ref_app)    else "",
        if (test_method == "tukey")   paste0("; Tukey CLD α=", tukey_alpha) else "",
        if (show_r2 == "caption" && !is.null(r2_txt)) paste0("\n", r2_txt) else ""
    )

    p <- ggplot2::ggplot(pred_dt,
            ggplot2::aes(x = app, y = y_pred, fill = app, color = app)) +
        ggplot2::geom_boxplot(alpha = 0.55, outlier.size = 0.5, width = 0.7) +
        ggplot2::stat_summary(fun = mean, geom = "point",
                              shape = 23, size = 1.8, fill = "white",
                              color = "black", stroke = 0.4)

    if (show_labels != "none") {
        p <- p + ggplot2::geom_text(
            data = label_dt,
            ggplot2::aes(x = app, y = y_top, label = text_label, color = app),
            vjust      = 0,
            size       = label_size,
            lineheight = 0.85,
            fontface   = if (show_labels == "letter") "bold" else "plain",
            show.legend = FALSE
        )
    }

    p <- p +
        ggplot2::scale_color_manual(values = colors, name = "Mating strategies") +
        ggplot2::scale_fill_manual (values = colors, name = "Mating strategies") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
        ggplot2::labs(x = NULL, y = "Predicted gain",
                      subtitle = if (show_r2 == "subtitle") r2_txt else NULL,
                      caption  = caption_txt) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
            legend.position    = "bottom",
            panel.grid.minor   = ggplot2::element_blank(),
            panel.spacing.x    = ggplot2::unit(0.3, "lines"),
            panel.spacing.y    = ggplot2::unit(0.6, "lines"),
            axis.text.x        = ggplot2::element_text(angle = 35, hjust = 1, size = 9),
            strip.background   = ggplot2::element_rect(fill = NA, color = NA),
            strip.text         = ggplot2::element_text(face = "bold", size = 11),
            plot.subtitle      = ggplot2::element_text(size = 9, color = "grey30"),
            plot.caption       = ggplot2::element_text(size = 8, color = "grey40")
        )

    # 嵌套分面: 行 = programs, 列 = inbreeding level
    n_prog <- data.table::uniqueN(d$programs)
    n_xv   <- data.table::uniqueN(pred_dt$x_target_factor)

    if (n_prog > 1 && n_xv > 1) {
        p <- p + ggplot2::facet_grid(
            programs ~ x_target_factor,
            scales   = facet_scales,
            labeller = ggplot2::labeller(
                x_target_factor = function(v) paste0(x_var, " = ", v)
            )
        )
    } else if (n_xv > 1) {
        p <- p + ggplot2::facet_wrap(
            ~ x_target_factor, scales = facet_scales,
            labeller = ggplot2::labeller(
                x_target_factor = function(v) paste0(x_var, " = ", v)
            )
        )
    } else if (n_prog > 1) {
        p <- p + ggplot2::facet_wrap(~ programs, scales = facet_scales)
    }

    list(plot        = p,
         predictions = pred_dt,
         fits        = rep_fits,
         r2_summary  = r2_summary,
         labels      = label_dt,
         signif      = sig_dt)
}
 plot_dt = copy(dat_full)
 plot_dt$app = factor(plot_dt$app,
                          levels = c("rl2", "tc", "ocsrate", "ocs25",
                                     "ocs45", "ocs65", "ocs90", "ran"),
                          labels = c("LAGM", "TC", "GOCS", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"))
plot_dt$programs = factor(plot_dt$programs,
                              levels = c("5_generation", "10_generation",
                                         "15_generation", "20_generation"),
                              labels = c("5-generation", "10-generation", "15-generation", "20-generation"))

plot_dt$app = factor(plot_dt$app,
                          levels = c("rl2_5m", "rl2_10m", "rl2_3","rl2_5","rl2","ocs25",
                                     "ocs45", "ocs65", "ocs90", "ran"),
                          labels = c("LAGM1", "LAGM2", "LAGM3", "LAGM5", "LAGM7", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"))
res_t <- plot_robust_predicted_boxplot(
    plot_dt,
    x_var       = "inb",
    y_var       = "gain",
    # progs       = "20_generation",
    x_values    = c(0.10, 0.2, 0.30,0.4,0.5),
    skip_gen    = 0,
    method      = "rlm",
    test_method = "tukey",     # ← Tukey HSD + 字母
    tukey_alpha = 0.05,
    extrapolate = TRUE,
    # colors     = c(LAGM   = "#1f3b73",  TC   = "#d62728",
    #                GOCS  = "#2ca02c",  GOCS25= "#1F77B4",
    #                GOCS45 = "#9467bd",  GOCS65= "#f4a582",
    #                GOCS90 = "#8b0000",  Random  = "#999999")
    colors     = c(LAGM1   = "#1f3b73",  LAGM2   = "#d62728",
                   LAGM3  = "#2ca02c",  GOCS25= "#1F77B4",
                   GOCS45 = "#9467bd",  GOCS65= "#f4a582",
                   GOCS90 = "#8b0000",  Random  = "#999999",
                   LAGM5 = "#79465e", LAGM7 = "#930c83")
)
print(res_t$plot)
ggsave("Figure 5.pdf", res_t$plot , width = 15, height = 10, dpi = 300,device = cairo_pdf)
ggsave("Figure 7.pdf", res_t$plot , width = 15.1, height = 8, dpi = 300,device = cairo_pdf)

print(res_t$r2_summary)
print(res_t$labels[, .(programs, x_target_factor, app, y_mean, signif_label)])

plot_tradeoff_arrow <- function(
    dat,
    x_var              = "inb",
    y_var              = "gain",
    apps               = NULL,
    progs              = NULL,
    fit_line           = FALSE,
    show_rep           = TRUE,
    arrow_size         = 2.5,
    rep_alpha          = 0.08,
    mean_lwd           = 1.2,
    facet_prog         = NULL,
    sec_axis_transform = NULL,
    sec_axis_name      = NULL,
    colors             = NULL,
    x_reverse          = FALSE,
    label_endpoints    = TRUE,    # 新增: 是否在末端写 app 名
    label_size         = 3.6,     # 新增: 末端标签字号
    use_repel          = TRUE     # 新增: 是否用 ggrepel 避免重叠
) {
    # ==========================================================
    # 数据筛选
    # ==========================================================
    d <- copy(dat)
    if (!is.null(apps))  d <- d[app %in% apps]
    if (!is.null(progs)) d <- d[programs %in% progs]
    d <- d[!is.na(get(x_var)) & !is.na(get(y_var))]
    setorder(d, app, programs, rep, gen)

    # ==========================================================
    # 自动分面：多个 program 时自动开启
    # ==========================================================
    if (is.null(facet_prog)) {
        facet_prog <- length(unique(d$programs)) > 1
    }

    # ==========================================================
    # 均值轨迹（按 app × program × gen 取均值）
    # ==========================================================
    mean_traj <- d[, .(
        x_mean = mean(get(x_var)),
        y_mean = mean(get(y_var))
    ), by = .(app, programs, gen)]
    setorder(mean_traj, app, programs, gen)

    # 标记最后一个和倒数第二个 gen（用于箭头段）
    mean_traj[, is_last := gen == max(gen), by = .(app, programs)]
    mean_traj[, is_second_last := gen == sort(unique(gen))[length(unique(gen)) - 1],
              by = .(app, programs)]

    arrow_end   <- mean_traj[is_last == TRUE]
    arrow_start <- mean_traj[is_second_last == TRUE]
    arrow_seg   <- merge(
        arrow_start[, .(app, programs, x0 = x_mean, y0 = y_mean)],
        arrow_end[,   .(app, programs, x1 = x_mean, y1 = y_mean)],
        by = c("app", "programs")
    )

    # 均值轨迹主体（去掉最后一个点，避免和箭头重叠）
    mean_body <- mean_traj[is_last == FALSE]

    # ==========================================================
    # 两阶段拟合线
    # ==========================================================
    fit_data <- NULL
    if (fit_line) {
        rep_fits <- d[, {
            xv  <- get(x_var)
            yv  <- get(y_var)
            fit <- rlm(yv ~ xv, maxit = 200)
            list(intercept = coef(fit)[1], slope = coef(fit)[2])
        }, by = .(app, programs, rep)]

        mean_fits <- rep_fits[, .(
            intercept = mean(intercept),
            slope     = mean(slope)
        ), by = .(app, programs)]

        x_ranges <- d[, .(
            x_min = min(get(x_var)),
            x_max = max(get(x_var))
        ), by = .(app, programs)]

        fit_data <- merge(x_ranges, mean_fits, by = c("app", "programs"))
        fit_data[, `:=`(
            y_start = intercept + slope * x_min,
            y_end   = intercept + slope * x_max
        )]
    }

    # ==========================================================
    # 颜色
    # ==========================================================
    all_apps <- sort(unique(d$app))
    if (is.null(colors)) {
        pal    <- scales::hue_pal()(length(all_apps))
        colors <- setNames(pal, all_apps)
    }

    # ==========================================================
    # 构建 ggplot
    # ==========================================================
    p <- ggplot()

    # 1) 每个 rep 的半透明轨迹
    if (show_rep) {
        p <- p + geom_path(
            data = d,
            aes(x     = .data[[x_var]],
                y     = .data[[y_var]],
                group = interaction(app, programs, rep),
                color = app),
            alpha     = rep_alpha,
            linewidth = 0.4
        )
    }

    # 2) 均值轨迹主体
    p <- p + geom_path(
        data = mean_body,
        aes(x     = x_mean,
            y     = y_mean,
            group = interaction(app, programs),
            color = app),
        linewidth = mean_lwd
    )

    # 3) 箭头（最后一段）
    p <- p + geom_segment(
        data = arrow_seg,
        aes(x = x0, y = y0, xend = x1, yend = y1, color = app),
        linewidth = mean_lwd,
        arrow     = arrow(length = unit(arrow_size, "mm"), type = "open")
    )

    # 4) 两阶段拟合线（虚线）
    if (fit_line && !is.null(fit_data)) {
        p <- p + geom_segment(
            data = fit_data,
            aes(x = x_min, y = y_start, xend = x_max, yend = y_end, color = app),
            linetype  = "dashed",
            linewidth = 0.7
        )
    }

    # 5) 末端文字标签（新增）
    if (label_endpoints) {
        if (use_repel && requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(
                data = arrow_end,
                aes(x = x_mean, y = y_mean, label = app, color = app),
                size = label_size, fontface = "bold",
                nudge_x = diff(range(mean_traj$x_mean)) * 0.02,
                direction = "y", hjust = 0,
                segment.size = 0.3, segment.alpha = 0.5,
                min.segment.length = 0,
                show.legend = FALSE
            )
        } else {
            p <- p + geom_text(
                data = arrow_end,
                aes(x = x_mean, y = y_mean, label = app, color = app),
                size = label_size, fontface = "bold",
                hjust = -0.15, vjust = 0.5,
                show.legend = FALSE
            )
        }
    }

    # ==========================================================
    # 坐标轴与主题
    # ==========================================================
    p <- p +
        scale_color_manual(values = colors, name = "Mating strategies") +
        labs(x = x_var, y = y_var) +
        theme_bw(base_size = 13) +
        theme(
            legend.position  = "bottom",
            panel.grid.minor = element_blank()
        )

    # 反转 x 轴 + 可选次坐标轴
    if (x_reverse) {
        if (!is.null(sec_axis_transform)) {
            p <- p + scale_x_reverse(
                sec.axis = sec_axis(
                    trans = sec_axis_transform,
                    name  = sec_axis_name %||% ""
                )
            )
        } else {
            p <- p + scale_x_reverse()
        }
    } else if (!is.null(sec_axis_transform)) {
        p <- p + scale_x_continuous(
            sec.axis = sec_axis(
                trans = sec_axis_transform,
                name  = sec_axis_name %||% ""
            )
        )
    }

    # 分面
    if (facet_prog) {
        p <- p + facet_wrap(~ programs, scales = "free")
    }

    p
}

# ============================================================
# 示例 1：复刻原图 — 全部 app，单个 program
# ============================================================
plot_dt = copy(dat_full)   # 只看
plot_dt$programs = factor(plot_dt$programs,
                              levels = c("5_generation", "10_generation",
                                         "15_generation", "20_generation"),
                              labels = c("5-generation", "10-generation", "15-generation", "20-generation"))
plot_dt$app = factor(plot_dt$app,
                          levels = c("rl2", "tc", "ocsrate", "ocs25",
                                     "ocs45", "ocs65", "ocs90", "ran"),
                          labels = c("LAGM", "TC", "GOCS", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"))
colors     = c(LAGM   = "#1f3b73",  TC   = "#d62728",
                   GOCS  = "#2ca02c",  GOCS25= "#1F77B4",
                   GOCS45 = "#9467bd",  GOCS65= "#f4a582",
                   GOCS90 = "#8b0000",  Random  = "#999999")



plot_dt$programs = factor(plot_dt$programs,
                              levels = c("5_generation", "10_generation",
                                         "15_generation", "20_generation"),
                              labels = c("5-generation", "10-generation", "15-generation", "20-generation"))
plot_dt$app = factor(plot_dt$app,
                          levels = c("rl2_5m", "rl2_10m", "rl2_3","rl2_5","rl2","ocs25",
                                     "ocs45", "ocs65", "ocs90", "ran"),
                          labels = c("LAGM1", "LAGM2", "LAGM3", "LAGM5", "LAGM7", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"))
colors     = c(LAGM1   = "#1f3b73",  LAGM2   = "#d62728",
                   LAGM3  = "#2ca02c",  GOCS25= "#1F77B4",
                   GOCS45 = "#f4a582",  GOCS65= "#9467bd",
                   GOCS90 = "#8b0000",  Random  = "#999999",
                   LAGM5 = "#79465e", LAGM7 = "#930c83")
P = plot_tradeoff_arrow(
    plot_dt,
    x_var  = "inb",
    y_var  = "gain",
    # progs  = "10 gen",
    # x_reverse = FALSE,
    # sec_axis_transform = ~ 1 - .,
    # sec_axis_name = "Converted/Lost genic standard deviation"
    rep_alpha = 0.2,
    fit_line = TRUE,
    colors = colors
)+labs(
  x = "Inbreeding",     # 修改 x 轴标题
  y = "Genetic gain"    # 修改 y 轴标题
)+theme(strip.background = element_rect(fill = NA, color = NA))

ggsave("Figure S1.pdf", P , width = 15, height = 10, dpi = 300,device = cairo_pdf)

library(patchwork)
# 假设 p1 是你截图中的图 (Y轴: Genetic gain)
# 假设 p2 是你的第二张图 (Y轴: 其他变量)
p1 = plot_tradeoff_arrow(
    plot_dt[app %in%c("GOCS25","GOCS45","LAGM1","LAGM5"),],
    x_var  = "gen",
    y_var  = "gain",
    # progs  = "10 gen",
    # x_reverse = FALSE,
    # sec_axis_transform = ~ 1 - .,
    # sec_axis_name = "Converted/Lost genic standard deviation"
    rep_alpha = 0.2,
    fit_line = TRUE,
    colors = colors
)+labs(
  x = "Generation",     # 修改 x 轴标题
  y = "Genetic gain"    # 修改 y 轴标题
)+theme(strip.background = element_rect(fill = NA, color = NA))


p2 = plot_tradeoff_arrow(
    plot_dt[app %in%c("GOCS25","GOCS45","LAGM1","LAGM5"),],
    x_var  = "gen",
    y_var  = "inb",
    # progs  = "10 gen",
    # x_reverse = FALSE,
    # sec_axis_transform = ~ 1 - .,
    # sec_axis_name = "Converted/Lost genic standard deviation"
    rep_alpha = 0.2,
    fit_line = TRUE,
    colors = colors
)+labs(
  x = "Generation",     # 修改 x 轴标题
  y = "Inbreeding"    # 修改 y 轴标题
)+theme(strip.background = element_rect(fill = NA, color = NA))

combined_plot <- (p1 | p2) + 
  plot_layout(guides = "collect") &  # 收集并合并相同的图例 
  theme(legend.position = "bottom")  # 将全局图例放置在底部

# 查看最终图片
combined_plot
ggsave("Figure 8.pdf", combined_plot , width = 15, height = 8, dpi = 300,device = cairo_pdf)

# ============================================================
# 第一阶段：每个 rep 用 rlm 拟合斜率
# ============================================================

#' 按 (app, programs, rep) 分组,稳健回归 y ~ x,返回每组斜率与汇总统计
#'
#' @param dat       data.table,需包含 app, programs, rep 以及 y/x 指定的列
#' @param y         因变量列名(字符串),如 "gain"
#' @param x         自变量列名(字符串),如 "inb" 或 "lost_He"
#' @param group_by  分组变量,默认 c("app","programs","rep")
#' @param maxit1    rlm 第一次尝试的最大迭代,默认 100
#' @param maxit2    rlm 第二次(warning 时)的最大迭代,默认 200
#'
#' @return list(rep_slopes = ..., slope_summary = ...)
fit_robust_slopes <- function(dat,
                              y,
                              x,
                              group_by = c("app", "programs", "rep"),
                              maxit1   = 100,
                              maxit2   = 200) {

    stopifnot(is.data.table(dat))
    stopifnot(all(c(y, x, group_by) %in% names(dat)))

    # 动态构造公式,例如 gain ~ inb
    f <- as.formula(paste(y, "~", x))

    # ---------- 第一阶段:每组拟合 ----------
    rep_slopes <- dat[, {
        y_vec <- get(y)
        x_vec <- get(x)
        ok    <- !is.na(y_vec) & !is.na(x_vec)

        # 数据太少或 x 无方差直接返回 NA
        if (sum(ok) < 3 || stats::sd(x_vec[ok]) == 0) {
            list(slope = NA_real_, intercept = NA_real_, r2 = NA_real_, n_obs = sum(ok))
        } else {
            sub <- data.frame(.y = y_vec[ok], .x = x_vec[ok])

            fit <- tryCatch(
                MASS::rlm(.y ~ .x, data = sub, maxit = maxit1),
                error   = function(e) NULL,
                warning = function(w) suppressWarnings(
                    MASS::rlm(.y ~ .x, data = sub, maxit = maxit2)
                )
            )

            if (is.null(fit)) {
                list(slope = NA_real_, intercept = NA_real_, r2 = NA_real_, n_obs = sum(ok))
            } else {
                yhat <- fitted(fit)
                r2   <- suppressWarnings(cor(yhat, sub$.y)^2)
                list(slope     = unname(coef(fit)[2]),
                     intercept = unname(coef(fit)[1]),
                     r2        = r2,
                     n_obs     = sum(ok))
            }
        }
    }, by = group_by]

    # ---------- 第二阶段:跨 rep 汇总 ----------
    summary_by <- setdiff(group_by, "rep")
    if (length(summary_by) == 0) summary_by <- NULL

    slope_summary <- rep_slopes[, .(
        mean_slope = mean(slope, na.rm = TRUE),
        sd_slope   = sd(slope,   na.rm = TRUE),
        se_slope   = sd(slope,   na.rm = TRUE) / sqrt(sum(!is.na(slope))),
        n          = sum(!is.na(slope))
    ), by = summary_by]

    # 附带元信息,方便后续合并多次调用
    rep_slopes[,    `:=`(y_var = y, x_var = x)]
    slope_summary[, `:=`(y_var = y, x_var = x)]

    list(rep_slopes = rep_slopes, slope_summary = slope_summary)
}


# 基本用法
res1 <- fit_robust_slopes(dat_full, y = "gain", x = "inb")
print(res1$rep_slopes)
print(res1$slope_summary)

# 一次跑多组并合并
combos <- list(
    c("gain", "inb"),
    c("gain", "inb_ped")
)
all_slopes <- rbindlist(lapply(combos, function(p) {
    fit_robust_slopes(dat_full, y = p[1], x = p[2])$slope_summary
}))
print(all_slopes)



#' 综合计算每个 (rep, app, programs) 的 ΔF 和 ΔG(per-generation),
#' 并支持多种算法 + 跨 rep 汇总
#'
#' @param dat        data.table,需含 rep, app, programs, gen, F_col, G_col
#' @param F_col      F 列名,默认 "inb_ped"(也可用 "inb")
#' @param G_col      G 列名,默认 "gv"
#' @param dF_method  ΔF 算法,可选:
#'                   - "loglinear":log(1-F) ~ gen 线性回归(推荐,带 R²)
#'                   - "pergen":逐代严格定义 ΔF_t = (F_t-F_{t-1})/(1-F_{t-1}) 求平均
#'                   - "endpoint":闭式几何平均(只用首末两点)
#' @param dG_method  ΔG 算法,可选:
#'                   - "regression":lm(G ~ gen) 回归斜率(推荐,带 SE 和 R²)
#'                   - "pergen":逐代差分 ΔG_t = G_t - G_{t-1} 求平均
#'                   - "endpoint":(G_final - G_0) / nG
#' @param late_half  额外报告"后半段"率(对标稳态),默认 TRUE
#' @param group_by   分组变量,默认 c("rep", "app", "programs")
#'
#' @return list:
#'   - per_rep:每 rep 的所有指标
#'   - summary:跨 rep 汇总(均值 ± SE)
compute_per_gen_rates <- function(dat,
                                  F_col         = "inb",
                                  G_col         = "gain",
                                  dF_method     = c("loglinear", "pergen", "endpoint"),
                                  dG_method     = c("regression", "pergen", "endpoint"),
                                  late_half     = TRUE,
                                  group_by      = c("rep", "app", "programs"),
                                  terminal_cols = NULL,
                                  terminal_gen  = "max") {

  dF_method <- match.arg(dF_method)
  dG_method <- match.arg(dG_method)

  stopifnot(is.data.table(dat))
  stopifnot(all(c(group_by, "gen", F_col, G_col) %in% names(dat)))

  # ---------- 内部: 单 rep 计算 ----------
  one_rep <- function(F_vec, G_vec, g_vec) {
    o     <- order(g_vec)
    F_vec <- F_vec[o]; G_vec <- G_vec[o]; g_vec <- g_vec[o]

    ok <- !is.na(F_vec) & !is.na(G_vec)
    F_vec <- F_vec[ok]; G_vec <- G_vec[ok]; g_vec <- g_vec[ok]
    n  <- length(F_vec)
    if (n < 3) {
      return(list(dF = NA_real_, dF_se = NA_real_, dF_R2 = NA_real_,
                  dG = NA_real_, dG_se = NA_real_, dG_R2 = NA_real_,
                  dF_late = NA_real_, dG_late = NA_real_,
                  F_0 = NA_real_, F_final = NA_real_,
                  G_0 = NA_real_, G_final = NA_real_,
                  nG  = NA_real_, Ne = NA_real_))
    }

    F_0 <- F_vec[1]; F_final <- F_vec[n]
    G_0 <- G_vec[1]; G_final <- G_vec[n]
    nG  <- g_vec[n] - g_vec[1]

    # ----- ΔF -----
    dF <- dF_se <- dF_R2 <- NA_real_
    if (dF_method == "loglinear") {
      keep <- F_vec < 1
      if (sum(keep) >= 3) {
        y <- log(1 - F_vec[keep]); x <- g_vec[keep]
        fit <- rlm(y ~ x,maxit=200); s <- coef(fit)[2]
        s_se <- summary(fit)$coefficients[2, "Std. Error"]
        dF <- 1 - exp(s); dF_se <- abs(exp(s)) * s_se
        dF_R2 <- summary(fit)$r.squared
      }
    } else if (dF_method == "pergen") {
      F_prev <- head(F_vec, -1); F_curr <- tail(F_vec, -1)
      keep <- F_prev < 1
      dF_t <- (F_curr[keep] - F_prev[keep]) / (1 - F_prev[keep])
      if (length(dF_t) > 0) { dF <- mean(dF_t); dF_se <- sd(dF_t)/sqrt(length(dF_t)) }
    } else if (dF_method == "endpoint") {
      if (F_0 < 1 && F_final < 1 && nG > 0)
        dF <- 1 - ((1 - F_final) / (1 - F_0))^(1 / nG)
    }

    # ----- ΔG -----
    dG <- dG_se <- dG_R2 <- NA_real_
    if (dG_method == "regression") {
      fit <- rlm(G_vec ~ g_vec,maxit=200)
      dG <- unname(coef(fit)[2])
      dG_se <- summary(fit)$coefficients[2, "Std. Error"]
      dG_R2 <- summary(fit)$r.squared
    } else if (dG_method == "pergen") {
      dG_t <- diff(G_vec); dG <- mean(dG_t); dG_se <- sd(dG_t)/sqrt(length(dG_t))
    } else if (dG_method == "endpoint") {
      if (nG > 0) dG <- (G_final - G_0) / nG
    }

    # ----- 后半段 (稳态) -----
    dF_late <- dG_late <- NA_real_
    if (late_half) {
      mid <- g_vec[1] + nG / 2
      late_idx <- which(g_vec >= mid)
      if (length(late_idx) >= 2) {
        F_late <- F_vec[late_idx]; g_late <- g_vec[late_idx]
        keep <- F_late < 1
        if (sum(keep) >= 2) {
          fit_l <- tryCatch(lm(log(1 - F_late[keep]) ~ g_late[keep]),
                            error = function(e) NULL)
          if (!is.null(fit_l)) dF_late <- 1 - exp(coef(fit_l)[2])
        }
        fit_g <- tryCatch(lm(G_vec[late_idx] ~ g_late), error = function(e) NULL)
        if (!is.null(fit_g)) dG_late <- unname(coef(fit_g)[2])
      }
    }

    Ne <- if (is.na(dF) || dF <= 0) NA_real_ else 1 / (2 * dF)

    list(dF = unname(dF), dF_se = unname(dF_se), dF_R2 = unname(dF_R2),
         dG = unname(dG), dG_se = unname(dG_se), dG_R2 = unname(dG_R2),
         dF_late = unname(dF_late), dG_late = unname(dG_late),
         F_0 = F_0, F_final = F_final,
         G_0 = G_0, G_final = G_final,
         nG  = nG,  Ne = Ne)
  }

  # ---------- 1) rates: per rep ----------
  per_rep <- dat[, one_rep(get(F_col), get(G_col), gen), by = group_by]
  per_rep[, `:=`(
    eff_per_gen      = dG / dF,
    eff_per_gen_late = dG_late / dF_late,
    eff_cumulative   = (G_final - G_0) / (F_final - F_0)
  )]

  # ---------- 2) 跨 rep 汇总 (rates) ----------
  summary_by <- setdiff(group_by, "rep")
  if (length(summary_by) == 0) summary_by <- NULL

  summary_dt <- per_rep[, .(
    dF_mean      = mean(dF, na.rm = TRUE),
    dF_se        = sd(dF, na.rm = TRUE)/sqrt(sum(!is.na(dF))),
    dF_R2_mean   = mean(dF_R2, na.rm = TRUE),
    dG_mean      = mean(dG, na.rm = TRUE),
    dG_se        = sd(dG, na.rm = TRUE)/sqrt(sum(!is.na(dG))),
    dG_R2_mean   = mean(dG_R2, na.rm = TRUE),
    dF_late_mean = mean(dF_late, na.rm = TRUE),
    dG_late_mean = mean(dG_late, na.rm = TRUE),
    Ne_mean      = mean(Ne, na.rm = TRUE),
    Ne_median    = median(Ne, na.rm = TRUE),
    eff_per_gen_mean      = mean(eff_per_gen, na.rm = TRUE),
    eff_per_gen_se        = sd(eff_per_gen, na.rm = TRUE)/sqrt(sum(!is.na(eff_per_gen))),
    eff_per_gen_late_mean = mean(eff_per_gen_late, na.rm = TRUE),
    eff_cumulative_mean   = mean(eff_cumulative, na.rm = TRUE),
    eff_cumulative_se     = sd(eff_cumulative, na.rm = TRUE)/sqrt(sum(!is.na(eff_cumulative))),
    n_rep   = .N
  ), by = summary_by]

  # ---------- 3) 终端值统计 ----------
  # 默认: 所有数值列, 排除 group_by 与 gen
  if (is.null(terminal_cols)) {
    num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
    terminal_cols <- setdiff(num_cols, c(group_by, "gen"))
  } else {
    miss <- setdiff(terminal_cols, names(dat))
    if (length(miss)) stop("terminal_cols not found in dat: ",
                           paste(miss, collapse = ", "))
  }

  # 选取每个 group 的终端行
  if (identical(terminal_gen, "max")) {
    terminal_per_rep <- dat[, .SD[which.max(gen)],
                            by = group_by,
                            .SDcols = c("gen", terminal_cols)]
  } else {
    stopifnot(is.numeric(terminal_gen), length(terminal_gen) == 1)
    terminal_per_rep <- dat[gen == terminal_gen,
                            .SD,
                            by = group_by,
                            .SDcols = c("gen", terminal_cols)]
  }
  # 给终端值列加 _final 后缀, 避免与 per_rep 已有列冲突
  setnames(terminal_per_rep,
           old = terminal_cols,
           new = paste0(terminal_cols, "_final"))

  # 跨 rep 汇总 (mean / sd / se / n)
  agg_expr <- unlist(lapply(terminal_cols, function(v) {
    vf <- paste0(v, "_final")
    setNames(
      list(
        bquote(mean(.(as.name(vf)), na.rm = TRUE)),
        bquote(sd(.(as.name(vf)),   na.rm = TRUE)),
        bquote(sd(.(as.name(vf)),   na.rm = TRUE) /
                 sqrt(sum(!is.na(.(as.name(vf)))))),
        bquote(sum(!is.na(.(as.name(vf)))))
      ),
      paste0(v, c("_mean", "_sd", "_se", "_n"))
    )
  }), recursive = FALSE)

  terminal_summary <- terminal_per_rep[, eval(as.call(c(quote(list), agg_expr))),
                                       by = summary_by]

  # 长格式 (画图友好)
  terminal_long <- melt(
    terminal_per_rep,
    id.vars       = c(group_by, "gen"),
    measure.vars  = paste0(terminal_cols, "_final"),
    variable.name = "metric",
    value.name    = "value"
  )
  terminal_long[, metric := sub("_final$", "", metric)]
  terminal_long_summary <- terminal_long[
    , .(mean = mean(value, na.rm = TRUE),
        sd   = sd(value,   na.rm = TRUE),
        se   = sd(value,   na.rm = TRUE) / sqrt(sum(!is.na(value))),
        n    = sum(!is.na(value))),
    by = c(summary_by, "metric")
  ]

  # ---------- 元信息 ----------
  per_rep[,         `:=`(F_col = F_col, G_col = G_col,
                         dF_method = dF_method, dG_method = dG_method)]
  summary_dt[,      `:=`(F_col = F_col, G_col = G_col,
                         dF_method = dF_method, dG_method = dG_method)]

  list(
    per_rep               = per_rep,
    summary               = summary_dt,
    terminal_per_rep      = terminal_per_rep,
    terminal_summary      = terminal_summary,
    terminal_long         = terminal_long,
    terminal_long_summary = terminal_long_summary
  )
}

# install.packages(c("multcompView"))

#' 对 per_rep 任意指标做 facet 柱状图 + Tukey CLD 字母
#'
#' @param per_rep    compute_per_gen_rates() 返回的 $per_rep
#' @param y_var      要画的列名(默认 "eff_per_gen")
#' @param facet_var  facet 变量(默认 "programs")
#' @param group_var  柱子分组(x 轴)变量(默认 "app")
#' @param app_order  app 显示顺序(NULL = 字母排序)
#' @param colors     命名色板
#' @param tukey_alpha  显著性水平
#' @param errorbar     "se" or "sd" or "none"
#' @param label_digits 均值小数位
#' @param y_lab        y 轴标题
#' @param plot_title   主标题
plot_per_rep_bar_with_cld <- function(
    per_rep,
    y_var         = "eff_per_gen",
    facet_var     = "programs",
    group_var     = "app",
    app_order     = NULL,
    colors        = NULL,
    tukey_alpha   = 0.05,
    errorbar      = c("se", "sd", "none"),
    label_digits  = 2,
    label_size    = 3.4,
    y_lab         = NULL,
    plot_title    = NULL
) {
    errorbar <- match.arg(errorbar)
    if (!requireNamespace("multcompView", quietly = TRUE))
        stop("Please install 'multcompView'.")
    stopifnot(y_var %in% names(per_rep))

    d <- as.data.table(copy(per_rep))                      # 强制深拷贝
    d <- d[!is.na(get(y_var)) & is.finite(get(y_var))]

    if (is.null(app_order)) app_order <- sort(unique(d[[group_var]]))
    set(d, j = group_var, value = factor(d[[group_var]], levels = app_order))

    # ---------- 汇总 ----------
    summary_dt <- d[, .(
        mean_y = mean(get(y_var), na.rm = TRUE),
        sd_y   = sd(get(y_var),   na.rm = TRUE),
        n      = .N
    ), by = c(facet_var, group_var)]
    summary_dt[, se_y := sd_y / sqrt(n)]
    summary_dt[, err  := switch(errorbar,
                                se = se_y, sd = sd_y, none = 0)]

    # ---------- Tukey CLD per facet ----------
    cld_dt <- d[, {
        sub <- copy(.SD)
        if (uniqueN(sub[[group_var]]) < 2 ||
            nrow(sub) < length(unique(sub[[group_var]])) * 2) {
            return(data.table(grp = unique(as.character(sub[[group_var]])),
                              letter = ""))
        }
        sub[, (group_var) := factor(get(group_var))]
        fit <- tryCatch(aov(as.formula(paste(y_var, "~", group_var)),
                            data = sub),
                        error = function(e) NULL)
        if (is.null(fit))
            return(data.table(grp = levels(sub[[group_var]]), letter = ""))
        tk <- tryCatch(TukeyHSD(fit, group_var,
                                conf.level = 1 - tukey_alpha),
                       error = function(e) NULL)
        if (is.null(tk))
            return(data.table(grp = levels(sub[[group_var]]), letter = ""))
        pvals <- tk[[group_var]][, "p adj"]
        cld <- multcompView::multcompLetters(pvals,
                                             threshold = tukey_alpha)$Letters
        data.table(grp = names(cld), letter = unname(cld))
    }, by = facet_var]
    setnames(cld_dt, "grp", group_var)

    plot_dt <- merge(summary_dt, cld_dt,
                     by = c(facet_var, group_var), all.x = TRUE)
    plot_dt <- as.data.table(plot_dt)                      # 重新 data.table 化
    plot_dt[is.na(letter), letter := ""]
    set(plot_dt, j = group_var,
        value = factor(plot_dt[[group_var]], levels = app_order))
    plot_dt[, text_label := paste0(letter, "\n",
                                   formatC(mean_y,
                                           digits = label_digits,
                                           format = "f"))]

    # ---------- 颜色 ----------
    if (is.null(colors))
        colors <- setNames(scales::hue_pal()(length(app_order)), app_order)

    if (is.null(y_lab))      y_lab      <- y_var
    if (is.null(plot_title)) plot_title <- paste("Mean", y_var,
                                                 "across mating strategies")

    # ---------- 画图(用 .data 替代 aes_string) ----------
    p <- ggplot(plot_dt,
                aes(x = .data[[group_var]], y = mean_y,
                    fill = .data[[group_var]])) +
        geom_col(width = 0.75, color = "black", linewidth = 0.3) +
        {if (errorbar != "none")
            geom_errorbar(aes(ymin = mean_y - err, ymax = mean_y + err),
                          width = 0.25, linewidth = 0.4)} +
        geom_text(aes(y = mean_y + err, label = text_label),
                  vjust = -0.2, size = label_size, lineheight = 0.85) +
        scale_fill_manual(values = colors, name = "Mating strategies") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
        facet_wrap(vars(.data[[facet_var]])) +
        labs(title = plot_title,
             x = "Mating strategies", y = y_lab,
             caption = paste0("Bars: mean ± ", errorbar,
                              "; letters: Tukey HSD α=", tukey_alpha,
                              "; n_rep = ", median(plot_dt$n))) +
        theme_bw(base_size = 13) +
        theme(legend.position  = "right",
              panel.grid.minor = element_blank(),
              plot.caption     = element_text(size = 8, color = "grey40"),
              strip.text       = element_text(face = "bold"))

    list(plot = p, summary = plot_dt)
}

# Step 1: 算速率
res <- compute_per_gen_rates(
    dat_full,
    F_col     = "inb",  # 或 "inb" 或 "inb_ped"
    G_col     = "gain",  # 或 "gain"
    dF_method = "loglinear",
    dG_method = "regression",
    group_by  = c("rep", "app", "programs"),
    terminal_cols = c("gv", "genetic", "genic", "He", "Ho", "Ne", "inb","inb_ped",
                    "gain", "gain_std", "lost_diversity_He",
                    "lost_diversity_Ho", "lost_diversity_genetic"),
    terminal_gen = "max"   # 或者写死 terminal_gen = 20
)

# Step 2: 画 conversion efficiency 柱状图(默认 eff_per_gen)
plot_dt = res$per_rep
plot_dt = res$terminal_per_rep  # 只看
plot_dt$programs = factor(plot_dt$programs,
                              levels = c("5_generation", "10_generation",
                                         "15_generation", "20_generation"),
                              labels = c("5-generation", "10-generation", "15-generation", "20-generation"))

plot_dt$app = factor(plot_dt$app,
                          levels = c("rl2", "tc", "ocsrate", "ocs25",
                                     "ocs45", "ocs65", "ocs90", "ran"),
                          labels = c("LAGM", "TC", "GOCS", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"))

plot_dt$app = factor(plot_dt$app,
                          levels = c("rl2_5m", "rl2_10m", "rl2_3","rl2_5","rl2","ocs25",
                                     "ocs45", "ocs65", "ocs90", "ran"),
                          labels = c("LAGM1", "LAGM2", "LAGM3", "LAGM5", "LAGM7", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"))
plot_dt$dF = plot_dt$dF*100   # 转换成百分比


y_var = "eff_per_gen"
y_lab = "Conversion efficiency"
plot_title = "Conversion efficiency across mating strategies and programs"

y_var = "dG"
y_lab = "Rate of genetic gain (ΔG)"
plot_title = "Rate of genetic gain across mating strategies and programs"

y_var = "dF"
y_lab = "Rate of inbreeding (ΔF)"
plot_title = "Rate of inbreeding across mating strategies and programs"
out <- plot_per_rep_bar_with_cld(
    plot_dt,
    y_var      = y_var,
    app_order  = c("LAGM", "TC", "GOCS", "GOCS25", "GOCS45", "GOCS65", "GOCS90", "Random"),
    colors     = c(LAGM   = "#1f3b73",  TC   = "#d62728",
                   GOCS  = "#2ca02c",  GOCS25= "#1F77B4",
                   GOCS45 = "#9467bd",  GOCS65= "#f4a582",
                   GOCS90 = "#8b0000",  Random  = "#999999"),
    # app_order = c("LAGM1", "LAGM2", "LAGM3", "LAGM5", "LAGM7", "OCS25",
    #                                  "OCS45", "OCS65", "OCS90", "Random"),           
    # colors     = c(LAGM1   = "#1f3b73",  LAGM2   = "#d62728",
    #                LAGM3  = "#2ca02c",  GOCS25= "#1F77B4",
    #                GOCS45 = "#9467bd",  GOCS65= "#f4a582",
    #                OCS90 = "#8b0000",  Random  = "#999999",
    #                LAGM5 = "#79465e", LAGM7 = "#930c83"),
    tukey_alpha = 0.05,
    errorbar    = "se",
    y_lab       = y_lab,
    plot_title  = plot_title
)
P = out$plot + theme(strip.background = element_rect(fill = NA, color = NA)) +
#   coord_cartesian(ylim = c(0, 55))   # 修改 y 轴范围（不会裁掉数据/误差棒）
#   coord_cartesian(ylim = c(0, 1.2))   # 修改 y 轴范围（不会裁掉数据/误差棒）
  coord_cartesian(ylim = c(0,10))   # 修改 y 轴范围（不会裁掉数据/误差棒）
P


ggsave("Figure 2.pdf", P,
       width = 15, height = 8, dpi = 300,
       device = cairo_pdf)
ggsave("Figure 3.pdf", P,
       width = 15, height = 8, dpi = 300,
       device = cairo_pdf)
ggsave("Figure 4.pdf", P ,
        width = 15, height = 8, dpi = 300,
        device = cairo_pdf)



y_var = "dG"
y_lab = "Genetic gain per generation (ΔG)"
plot_title = "Rate of genetic gain for LAGM with fixed look-ahead windows within 20-generation program"
out <- plot_per_rep_bar_with_cld(
    plot_dt,
    y_var      = y_var,
    # app_order  = c("LAGM", "TC", "OCS", "OCS25", "OCS45", "OCS65", "OCS90", "Random"),
    # colors     = c(LAGM   = "#1f3b73",  TC   = "#d62728",
    #                OCS  = "#2ca02c",  OCS25= "#1F77B4",
    #                OCS45 = "#9467bd",  OCS65= "#f4a582",
    #                OCS90 = "#8b0000",  Random  = "#999999"),
    app_order = c("LAGM1", "LAGM2", "LAGM3", "LAGM5", "LAGM7", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"),           
    colors     = c(LAGM1   = "#1f3b73",  LAGM2   = "#d62728",
                   LAGM3  = "#2ca02c",  GOCS25= "#1F77B4",
                   GOCS45 = "#9467bd",  GOCS65= "#f4a582",
                   GOCS90 = "#8b0000",  Random  = "#999999",
                   LAGM5 = "#79465e", LAGM7 = "#930c83"),
    tukey_alpha = 0.05,
    errorbar    = "se",
    y_lab       = y_lab,
    plot_title  = plot_title
)
P6 = out$plot + theme(strip.background = element_rect(fill = NA, color = NA)) +
#   coord_cartesian(ylim = c(0, 100))   # 修改 y 轴范围（不会裁掉数据/误差棒）
  coord_cartesian(ylim = c(0, 1))   # 修改 y 轴范围（不会裁掉数据/误差棒）
#   coord_cartesian(ylim = c(0,10))   # 修改 y 轴范围（不会裁掉数据/误差棒）
P6

y_var = "dF"
y_lab = "Rate of inbreeding per generation (ΔF)"
plot_title = "Rate of inbreeding for LAGM with fixed look-ahead windows within 20-generation program"
out <- plot_per_rep_bar_with_cld(
    plot_dt,
    y_var      = y_var,
    # app_order  = c("LAGM", "TC", "OCS", "OCS25", "OCS45", "OCS65", "OCS90", "Random"),
    # colors     = c(LAGM   = "#1f3b73",  TC   = "#d62728",
    #                OCS  = "#2ca02c",  OCS25= "#1F77B4",
    #                OCS45 = "#9467bd",  OCS65= "#f4a582",
    #                OCS90 = "#8b0000",  Random  = "#999999"),
    app_order = c("LAGM1", "LAGM2", "LAGM3", "LAGM5", "LAGM7", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"),           
    colors     = c(LAGM1   = "#1f3b73",  LAGM2   = "#d62728",
                   LAGM3  = "#2ca02c",  GOCS25= "#1F77B4",
                   GOCS45 = "#9467bd",  GOCS65= "#f4a582",
                   GOCS90 = "#8b0000",  Random  = "#999999",
                   LAGM5 = "#79465e", LAGM7 = "#930c83"),
    tukey_alpha = 0.05,
    errorbar    = "se",
    y_lab       = y_lab,
    plot_title  = plot_title
)
P7 = out$plot + theme(strip.background = element_rect(fill = NA, color = NA)) +
#   coord_cartesian(ylim = c(0, 100))   # 修改 y 轴范围（不会裁掉数据/误差棒）
#   coord_cartesian(ylim = c(0, 1.2))   # 修改 y 轴范围（不会裁掉数据/误差棒）
  coord_cartesian(ylim = c(0,5))   # 修改 y 轴范围（不会裁掉数据/误差棒）
P7

y_var = "eff_per_gen"
y_lab = "Conversion efficiency"
plot_title = "Conversion efficiency for LAGM with fixed look-ahead windows within 20-generation program"
out <- plot_per_rep_bar_with_cld(
    plot_dt,
    y_var      = y_var,
    # app_order  = c("LAGM", "TC", "OCS", "OCS25", "OCS45", "OCS65", "OCS90", "Random"),
    # colors     = c(LAGM   = "#1f3b73",  TC   = "#d62728",
    #                OCS  = "#2ca02c",  OCS25= "#1F77B4",
    #                OCS45 = "#9467bd",  OCS65= "#f4a582",
    #                OCS90 = "#8b0000",  Random  = "#999999"),
    app_order = c("LAGM1", "LAGM2", "LAGM3", "LAGM5", "LAGM7", "GOCS25",
                                     "GOCS45", "GOCS65", "GOCS90", "Random"),           
    colors     = c(LAGM1   = "#1f3b73",  LAGM2   = "#d62728",
                   LAGM3  = "#2ca02c",  GOCS25= "#1F77B4",
                   GOCS45 = "#9467bd",  GOCS65= "#f4a582",
                   GOCS90 = "#8b0000",  Random  = "#999999",
                   LAGM5 = "#79465e", LAGM7 = "#930c83"),
    tukey_alpha = 0.05,
    errorbar    = "se",
    y_lab       = y_lab,
    plot_title  = plot_title
)
P8 = out$plot + theme(strip.background = element_rect(fill = NA, color = NA)) +
  coord_cartesian(ylim = c(0, 50))   # 修改 y 轴范围（不会裁掉数据/误差棒）
#   coord_cartesian(ylim = c(0, 1.2))   # 修改 y 轴范围（不会裁掉数据/误差棒）
#   coord_cartesian(ylim = c(0,10))   # 修改 y 轴范围（不会裁掉数据/误差棒）
P8



P6_new <- P6 + labs(title = "a) Rate of genetic gain (ΔG)", y = NULL, x=NULL,subtitle = NULL)
P7_new <- P7 + labs(title = "b) Rate of inbreeding (ΔF)", y = NULL, x=NULL,subtitle = NULL)
P8_new <- P8 + labs(title = "c) Conversion efficiency", y = NULL, x=NULL,subtitle = NULL)

# 2. 定义“上面两张，下面一张居中”的布局
design <- "
  AABB
  #CC#
"

# 3. 组合图片并收集图例
# 修正：将 plot_title 改为 plot.title
combined_plot <- P6_new + P7_new + P8_new + 
  plot_layout(design = design, guides = "collect") & 
  theme(
    legend.position = "bottom",
    # 统一设置子标题的样式
    plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 6)),
    # 彻底隐去 x 轴和 y 轴的标题文本空间
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # 【核心防御】如果“20-generation”是分面(facet)带来的，下面这两行可以强行让它消失
    strip.text = element_blank(),
    strip.background = element_blank()
  )

# 4. 添加统一的全局大标题（同样使用正确的 plot.title）
combined_plot <- combined_plot + 
  plot_annotation(
    title = "LAGM with fixed look-ahead windows within 20-generation program", 
    theme = theme(
      # margin(b = 10) 可以给大标题下方留出空隙，防止和下面的子图挨得太近
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10))
    )
  )

ggsave("Figure 6.pdf", combined_plot ,
        width = 15, height = 8, dpi = 300,
        device = cairo_pdf)





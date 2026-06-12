rm(list=ls())
# 加载包
library(ggplot2)
library(rms)
library(dplyr)
library(readxl)
library(ggplot2)
library(readxl)
library(survival)
library(patchwork)

# 读取数据
data <- read_excel("Cox-ESRD.xlsx",sheet = "Sheet1")
#class(data$group)
data$group<- factor(data$group,
                    levels = c(1, 2),
                    labels = c("GBM", "DP"))
#class(data$sex)
#class(data$Lymphocyte)
data$sex<- factor(data$sex,
                  levels = c(1, 2),
                  labels = c("male", "female"))
#class(data$group)
#class(data$anuria)
data$LH<- factor(data$LH,
                 levels = c(0, 1),
                 labels = c("no", "yes"))
data$RRT<- factor(data$RRT,
                  levels = c(0, 1),
                  labels = c("no", "yes"))
data$smoke<- factor(data$smoke,
                    levels = c(0, 1),
                    labels = c("no", "yes"))
data$anuria<- factor(data$anuria,
                     levels = c(0, 1),
                     labels = c("no", "yes"))
# 确保ESRD是二元变量（0/1）
table(data$ESRD)
class(data$ESRD)
sum(is.na(data$ESRD))
class(data$time)

###### 4. COX回归对应的RCS #####
dd <- datadist(data); options(datadist = "dd")
# 创建限制性立方样条模型
model_cox_rcs <- lrm(ESRD ~ rcs(eGFR, 3) + group + sex + age + LH + RRT, data = data)
pred <- Predict(model_cox_rcs, eGFR, ref.zero = TRUE, fun = exp)  # HR曲线
anova(model_cox_rcs)#计算P值

plot_cox_rcs_curve <- function(pred, hist_data, xvar, outcome = "ESRD", model = NULL, 
                               title = NULL, show_title = TRUE,
                               save_plot = FALSE, filename = NULL, 
                               width = 8, height = 6, dpi = 300) {
  #library(ggplot2)
  #library(rlang)
  
  # 检查必要的数据列是否存在
  if (!xvar %in% names(pred)) {
    stop(paste("Variable", xvar, "not found in pred data"))
  }
  
  if (!xvar %in% names(hist_data)) {
    stop(paste("Variable", xvar, "not found in hist_data"))
  }
  
  # 1. 设置参考值（置信区间最窄点）
  ref_val <- pred[[xvar]][which.min(pred$upper - pred$lower)]
  
  # 2. 提取LRT P值（如果提供了模型对象）
  p_lrt <- ""
  if (!is.null(model)) {
    p_raw <- tryCatch(anova(model)[2, "P"], error = function(e) NA)
    if (!is.na(p_raw)) {
      p_lrt <- if (p_raw < 0.001) {
        "P for nonlinear < 0.001"
      } else {
        paste0("P for nonlinear = ", format(round(p_raw, 3), nsmall = 3))
      }
    }
  }
  
  # 3. X轴范围自动化
  x_min <- min(pred[[xvar]], na.rm = TRUE)
  x_max <- ceiling(max(pred[[xvar]], na.rm = TRUE) / 10) * 10
  x_mid <- (x_min + x_max) / 2
  
  # 4. 标题智能生成 - 根据 show_title 参数决定是否显示标题
  if (show_title) {
    main_title <- ifelse(is.null(title), paste0("Restricted Cubic Spline Curve: ", xvar, " vs ", outcome), title)
  } else {
    main_title <- NULL
  }
  
  # 5. 使用现代ggplot2语法（替代aes_string）
  p <- ggplot() +
    geom_histogram(
      data = hist_data, 
      aes(x = .data[[xvar]], y = after_stat(count)/max(after_stat(count))*1.5), 
      bins = 30, fill = "skyblue", alpha = 0.6
    ) +
    geom_line(
      data = pred, 
      aes(x = .data[[xvar]], y = yhat), 
      color = "red", size = 1.3
    ) +
    geom_ribbon(
      data = pred, 
      aes(x = .data[[xvar]], ymin = lower, ymax = upper), 
      fill = "red", alpha = 0.2
    ) +
    geom_vline(xintercept = ref_val, linetype = "dashed", color = "#00BFC4", size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#00BFC4", size = 0.8) +
    annotate("text", x = ref_val, y = max(pred$upper, na.rm = TRUE) * 0.95, 
             label = paste0("Ref. point = ", round(ref_val, 2)), 
             color = "#00BFC4", size = 4, fontface = "italic", hjust = -0.1) +
    annotate("text", x = x_mid, y = max(pred$upper, na.rm = TRUE) * 1.05,
             label = p_lrt, size = 4.5, fontface = "bold", hjust = 0.5) +
    scale_x_continuous(limits = c(x_min, x_max)) +
    labs(
      title = main_title,
      x = "eGFR", 
      y = "Hazard Ratio (95% CI)"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_blank(),
      plot.title = element_text(face = "bold", size = 14),  # 增强标题
      axis.title = element_text(size = 12)  # 坐标轴标题大小
    )
  
  # 6. 如果需要保存图片
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0("rcs_", xvar, "_vs_", outcome, 5,".tiff")
    }
    
    ggsave(
      filename = filename,
      plot = p,
      device = "tiff",
      width = 180,        # 宽度(毫米)
      height = 140,       # 高度(毫米)
      units = "mm",       # 单位(毫米)
      dpi = dpi,
      compression = "lzw"
    )
    
    message(paste("Plot saved as:", filename))
  }
  
  return(p)
}
# 生成图形并自动保存为默认文件名
p <- plot_cox_rcs_curve(
  pred = pred, 
  hist_data = data,
  xvar = "eGFR", 
  outcome = "ESRD", 
  model = model_cox_rcs,
  show_title = FALSE,
  save_plot = TRUE  # 关键参数：设置为 TRUE 以保存图片
)

#new
analyze_segmented_auto <- function(data, outcome, exposure, covariates,
                                   model_type = c("logit", "linear", "cox"),
                                   time = NULL, event = NULL,
                                   pred_obj = NULL) {
  library(dplyr)
  library(survival)
  model_type <- match.arg(model_type)
  
  # 1. 检查pred_obj是否存在且包含必要的信息
  if (is.null(pred_obj) || !all(c(exposure, "upper", "lower") %in% names(pred_obj))) {
    stop("pred_obj must be provided and contain the exposure variable, upper and lower columns")
  }
  # 2. 检查并转换时间变量（如果是 Cox 模型）
  if (model_type == "cox") {
    if (!is.numeric(data[[time]])) {
      # 尝试转换为数值型
      data[[time]] <- as.numeric(as.character(data[[time]]))
      if (any(is.na(data[[time]]))) {
        stop(paste("Time variable", time, "cannot be converted to numeric"))
      }
      message(paste("Converted time variable", time, "to numeric"))
    }
  }
  # 2. 计算参考点（确保有有效值）
  ci_width <- pred_obj$upper - pred_obj$lower
  if (all(is.na(ci_width)) || length(ci_width) == 0) {
    stop("Cannot calculate reference point: invalid confidence interval values")
  }
  
  ref_val <- pred_obj[[exposure]][which.min(ci_width)]
  
  # 3. 创建分段变量
  data[[paste0(exposure, "_left")]] <- pmin(data[[exposure]], ref_val)
  data[[paste0(exposure, "_right")]] <- pmax(data[[exposure]] - ref_val, 0)
  
  left <- paste0(exposure, "_left")
  right <- paste0(exposure, "_right")
  cov_str <- paste(covariates, collapse = " + ")
  
  # 4. 根据模型类型构建公式和拟合模型
  if (model_type == "cox") {
    if (is.null(time) || is.null(event)) {
      stop("For Cox models, both time and event variables must be specified")
    }
    
    data$surv_obj <- Surv(data[[time]], data[[event]])
    formula_seg <- as.formula(paste0("surv_obj ~ ", left, " + ", right, " + ", cov_str))
    formula_lin <- as.formula(paste0("surv_obj ~ ", exposure, " + ", cov_str))
    
    model_seg <- coxph(formula_seg, data = data)
    model_lin <- coxph(formula_lin, data = data)
    label <- "HR"
    
  } else {
    family <- if (model_type == "logit") binomial else gaussian
    formula_seg <- as.formula(paste0(outcome, " ~ ", left, " + ", right, " + ", cov_str))
    formula_lin <- as.formula(paste0(outcome, " ~ ", exposure, " + ", cov_str))
    
    model_seg <- glm(formula_seg, data = data, family = family)
    model_lin <- glm(formula_lin, data = data, family = family)
    label <- if (model_type == "logit") "OR" else "Beta"
  }
  
  # 5. 提取估计值与标准误
  coef_seg <- coef(model_seg)
  vcov_seg <- vcov(model_seg)
  b1 <- coef_seg[left]; se1 <- sqrt(vcov_seg[left, left])
  b2 <- coef_seg[right]; se2 <- sqrt(vcov_seg[right, right])
  
  est1 <- if (label == "Beta") b1 else exp(b1)
  ci1 <- if (label == "Beta") b1 + c(-1.96, 1.96) * se1 else exp(b1 + c(-1.96, 1.96) * se1)
  p1 <- 2 * (1 - pnorm(abs(b1 / se1)))
  
  est2 <- if (label == "Beta") b2 else exp(b2)
  ci2 <- if (label == "Beta") b2 + c(-1.96, 1.96) * se2 else exp(b2 + c(-1.96, 1.96) * se2)
  p2 <- 2 * (1 - pnorm(abs(b2 / se2)))
  
  # 6. LRT 手动计算（兼容 coxph 和 glm）
  lrt_stat <- 2 * (logLik(model_seg)[1] - logLik(model_lin)[1])
  df_diff <- attr(logLik(model_seg), "df") - attr(logLik(model_lin), "df")
  p_lrt <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
  
  # 7. 输出
  result <- data.frame(
    Segment = c(paste0(exposure, " < ", round(ref_val, 2)),
                paste0(exposure, " ≥ ", round(ref_val, 2))),
    Estimate = round(c(est1, est2), 3),
    CI = c(sprintf("(%.3f ~ %.3f)", ci1[1], ci1[2]),
           sprintf("(%.3f ~ %.3f)", ci2[1], ci2[2])),
    P_value = signif(c(p1, p2), 3)
  )
  
  cat("→ 模型类型:", label, "\n")
  print(result, row.names = FALSE)
  cat("\nLikelihood Ratio Test P =", signif(p_lrt, 3), "\n")
  
  # 返回模型和结果以供进一步分析
  invisible(list(
    segmented_model = model_seg,
    linear_model = model_lin,
    reference_value = ref_val,
    result_table = result,
    lrt_p_value = p_lrt
  ))
}
# 注意：调用函数时有两个逗号错误，应该删除一个
result <- analyze_segmented_auto(
  data = data,
  outcome = "ESRD",
  exposure = "eGFR",
  covariates = c("group", "sex", "age","LH","RRT"),
  model_type = "cox",  # 这里删除了多余的逗号
  time = "time",   # 需要指定时间变量名
  event = "ESRD", # 需要指定事件变量名
  pred_obj = pred
)

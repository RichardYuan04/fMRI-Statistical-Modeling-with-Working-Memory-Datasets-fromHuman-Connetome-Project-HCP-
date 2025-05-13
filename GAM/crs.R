# 加载必要的 R 包q()
library("crs")
library("doParabar")
library("foreach")

# 加载 compute_mean_deriv 函数
source("boot_der.R")

# 读取数据
df <- read.csv("WM_Face_2bk_cleaned.csv")  # 替换为您的 dataset 文件路径
vars <- colnames(df)
Y_names <- vars[10:369]  # 假设响应变量从第 10 列到第 369 列，需根据 dataset 调整

# 提取自变量
X <- df[c("BIS", "Age", "Gender")]
X$Age <- ordered(X$Age)      # 将 Age 转换为有序因子
X$Gender <- factor(X$Gender) # 将 Gender 转换为因子

# 创建结果数据框
results <- data.frame(
  label = character(),
  First_Mean = numeric(),
  First_P = numeric(),
  First_Sig = logical(),
  Second_Mean = numeric(),
  Second_P = numeric(),
  Second_Sig = logical(),
  stringsAsFactors = FALSE
)

# 循环处理每个响应变量
for (Y in Y_names) {
  cat(sprintf("Processing variable: %s\n", Y))
  
  # 构建公式
  f <- formula(paste(Y, "~ BIS + ordered(Age) + factor(Gender)"))
  
  # 自动选择带宽
  model_bw <- npglpreg(formula = f, data = df, bwtype = "auto", degree = 4, nmulti = 5, cv = "bandwidth")
  bws <- model_bw$bws
  
  # 计算一阶导数
  result_first <- compute_mean_deriv(
    Y = df[, Y],
    X = X,
    bws = bws,
    degree = 4,
    gradient.vec = 1,  # 一阶导数
    B1 = 1000,
    B2 = 100
  )
  
  # 计算二阶导数
  result_second <- compute_mean_deriv(
    Y = df[, Y],
    X = X,
    bws = bws,
    degree = 4,
    gradient.vec = 2,  # 二阶导数
    B1 = 1000,
    B2 = 100
  )
  
  # 添加结果到数据框（一行包含一阶和二阶导数结果）
  results <- rbind(results, data.frame(
    label = Y,
    First_Mean = result_first$mean_derivative,
    First_P = result_first$p_value,
    First_Sig = result_first$p_value < 0.05,
    Second_Mean = result_second$mean_derivative,
    Second_P = result_second$p_value,
    Second_Sig = result_second$p_value < 0.05
  ))
}

# 保存结果到 CSV 文件
write.csv(results, "derivative_results.csv", row.names = FALSE)
cat("Results saved to derivative_results.csv\n")
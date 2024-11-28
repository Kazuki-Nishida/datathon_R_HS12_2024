# ---- 必要パッケージのインストール（未インストールの場合のみ）----
if (!require("MASS")) install.packages("MASS")
if (!require("truncnorm")) install.packages("truncnorm")

# ---- ライブラリ読み込み ----
library(MASS)
library(truncnorm)

# ---- データ生成コード ----

set.seed(123)

# サンプルサイズ
n <- 660

# Ageの混合正規分布化
age_group1 <- rnorm(round(n * 0.7), mean = 50, sd = 5)  # 70%: 平均50, SD=5
age_group2 <- rnorm(round(n * 0.3), mean = 70, sd = 5)  # 30%: 平均70, SD=5
age <- round(c(age_group1, age_group2))  # 混合正規分布を結合し整数化

# 他の基本変数の生成
sex <- sample(c(1, 0), n, replace = TRUE, prob = c(0.6, 0.4))  # 男性:1, 女性:0（男性が多め）
bmi <- round(rnorm(n, mean = 25, sd = 4), 1)  # BMI（小数点以下1桁）
smoking <- sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))  # 喫煙の有無
diabetes <- sample(c(0, 1), n, replace = TRUE, prob = c(0.8, 0.2))  # 糖尿病
family_history <- sample(c(0, 1), n, replace = TRUE, prob = c(0.6, 0.4))  # 家族歴

# LDLの生成
ldl <- round(rtruncnorm(n, a = 80, b = 240, mean = 160, sd = 30))  # LDL：範囲80～240
ldl[ldl == 80] <- ldl[ldl == 80] + round(rnorm(sum(ldl == 80), mean = 5, sd = 5))  # 下限値調整
ldl[ldl == 240] <- ldl[ldl == 240] - round(rnorm(sum(ldl == 240), mean = 5, sd = 5))  # 上限値調整

# HDLの生成
hdl <- round(rtruncnorm(n, a = 45, b = 65, mean = 55, sd = 7))  # トランケート分布
hdl[hdl == 45] <- hdl[hdl == 45] + round(rnorm(sum(hdl == 45), mean = 2, sd = 3))  # 下限値調整
hdl[hdl == 65] <- hdl[hdl == 65] + round(rnorm(sum(hdl == 65), mean = 1, sd = 10))  # 上限値調整

# EFの生成
ef_base <- 50 + 0.2 * hdl - 0.1 * ldl - 0.15 * age - 3 * smoking - 2 * diabetes + 1.5 * family_history
ef <- round(pmax(20, pmin(ef_base + rnorm(n, mean = 0, sd = 5), 70)))  # 基本範囲20～70
ef[ef == 20] <- ef[ef == 20] + round(rnorm(sum(ef == 20), mean = 2, sd = 3))  # 下限値調整
ef[ef == 70] <- ef[ef == 70] - round(rnorm(sum(ef == 70), mean = 2, sd = 3))  # 上限値調整
ef_outcome <- ifelse(ef < 40, 1, 0)  # EF < 40のアウトカム

# s_BloodPressureの生成
s_BloodPressure <- round(rtruncnorm(n, a = 100, b = 170, mean = 135, sd = 10))  # トランケート分布
s_BloodPressure[s_BloodPressure == 100] <- s_BloodPressure[s_BloodPressure == 100] + round(rnorm(sum(s_BloodPressure == 100), mean = 5, sd = 5))  # 下限値調整
s_BloodPressure[s_BloodPressure == 170] <- s_BloodPressure[s_BloodPressure == 170] - round(rnorm(sum(s_BloodPressure == 170), mean = 5, sd = 5))  # 上限値調整

# ファクターX, Y, Zの生成
factor_x <- rbinom(n, 1, prob = 0.3)  # 二値化された新リスク因子
factor_y <- scale(rnorm(n, mean = 0, sd = 1))  # 標準正規分布から生成（標準化）
factor_z <- scale(rnorm(n, mean = 0, sd = 1))  # 標準正規分布から生成（標準化）

# CVDイベントの生成（スケーリング後）
cvd_risk <- -0.42 + 0.05 * scale(age) + 0.4 * scale(sex) + 0.3 * scale(bmi) + 
            0.6 * scale(smoking) + 0.5 * scale(diabetes) + 0.1 * scale(ldl) - 
            0.01 * scale(hdl) + 0.15 * scale(family_history) +
            0.15 * scale(bmi) * scale(smoking) + 0.1 * scale(ldl) * scale(hdl) + 
            0.5 * scale(factor_x) + 0.3 * scale(factor_y) - 0.2 * scale(factor_z)
cvd_event_prob <- exp(cvd_risk) / (1 + exp(cvd_risk))  # ロジスティック関数
cvd_event <- rbinom(n, size = 1, prob = cvd_event_prob)  # 発生確率

# 打ち切り設定
censoring_flag <- ifelse(cvd_event == 0, rbinom(n, size = 1, prob = 0.3), 0)  # イベント非発生者の30%に打ち切り
censoring_time_years <- ifelse(censoring_flag == 1, runif(n, min = 0.1, max = 5), 5)
censoring_time_days <- round(censoring_time_years * 365)

# イベント発生時間の生成（ベータ分布）
alpha <- 2  # 固定値
beta <- 1 / cvd_event_prob  # 発生確率に基づく調整
event_time_years <- rbeta(n, shape1 = alpha, shape2 = beta) * 5  # ベータ分布（0～5年）
event_time_days <- round(event_time_years * 365)

# 観察時間の設定
observed_time <- ifelse(cvd_event == 1, 
                        pmin(event_time_days, censoring_time_days),  # イベント発生の場合
                        censoring_time_days)  # イベント非発生の場合は打ち切り時間

# データフレーム化
simulated_data <- data.frame(
  ID = 1:n,
  Age = age,
  sex_male = sex,
  BMI = bmi,
  Smoking = smoking,
  Diabetes = diabetes,
  FamilyHistory = family_history,
  s_BloodPressure = s_BloodPressure,
  LDL = ldl,
  HDL = hdl,
  EF = ef,
  EF_Outcome = ef_outcome,
  FactorX = factor_x,
  Observed_Time = observed_time,
  CVD_Event = cvd_event  # これがイベントフラグ
)

# ---- 検証 ----
summary(simulated_data)  # データの要約
table(simulated_data$CVD_Event)  # CVDイベントの分布確認


############################


# ---- データフレームを df に変更 ----
df <- simulated_data

# ---- データ構造と基本情報の表示 ----
dim(df)  # データの行数と列数
str(df)  # 各変数の構造情報
summary(df)  # 変数の要約統計

# ---- 数値型変数のヒストグラム（例: Age） ----
hist(df$Age, main = "Histogram of Age", xlab = "Age", col = "lightblue", border = "white")

# ---- 箱ひげ図（例: BMI） ----
boxplot(df$BMI, main = "Boxplot of BMI", ylab = "BMI", col = "pink")


# 一括でヒストグラムを見る
# ---- 必要パッケージのインストール（未インストールの場合のみ） ----
if (!require("DataExplorer")) install.packages("DataExplorer")

# ---- ライブラリ読み込み ----
library(DataExplorer)

# ---- 数値型変数のヒストグラム ----
plot_histogram(df)  # すべての数値型変数のヒストグラムを一括で描画

plot_boxplot(df, by = "sex_male")  # 性別に数値型変数の分布を確認

# ---- 相関係数のヒートマップ（簡易版） ----
if (!require("corrplot")) install.packages("corrplot")
library(corrplot)

# 数値型変数を抽出して相関係数行列を計算
numeric_data <- df[, sapply(df, is.numeric)]
correlation_matrix <- cor(numeric_data)

# 相関係数のヒートマップを描画
corrplot(correlation_matrix, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200), 
         addCoef.col = "black", tl.cex = 0.8, number.cex = 0.7)

# ---- ロジスティック回帰モデル（EF と X を除く） ----
logistic_model <- glm(
  formula = CVD_Event ~ Age + sex_male + BMI + Smoking + Diabetes + 
              LDL + HDL + s_BloodPressure + FamilyHistory, 
  family = binomial(link = "logit"), 
  data = df
)

# Logistic
lresult = function(x) {
  Coefficients_decimal = formatC(summary(x)$coefficient[,1], digits = 2, format = "f")
  OR_decimal = formatC(exp(summary(x)$coefficient[,1]), digits = 2, format = "f")
  CI = paste("[", formatC(exp(confint.default(x))[,1], digits = 2, format = "f"), ", ", formatC(exp(confint.default(x))[,2], digits = 2, format = "f"), "]", sep = "")
  P_value = ifelse(summary(x)$coefficient[,4] < 0.001, "<0.001", formatC(summary(x)$coefficient[,4], digits = 3, format = "f"))
  result = cbind(Coefficients_decimal, OR_decimal, CI, P_value)
  colnames(result) = c("Coefficients", "Odds Ratio", "[95% CI of OR]", "     P")
  table = noquote(result)
  attr(table, "note") = paste(note = paste("Objective Variable: ", gsub("[()]", "", x$formula[2]), collapse = ""))
  return(table)
}

# ロジスティック回帰モデルの結果を表示
lresult(logistic_model)

# ---- ROC曲線のプロットとAUCの計算 ----
if (!require("pROC")) install.packages("pROC")
library(pROC)

# ROC曲線の計算
roc_curve <- roc(df$CVD_Event, predict(logistic_model, type = "response"))

# ROC曲線のプロット
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)

# AUCの計算
auc_value <- auc(roc_curve)
cat("AUC:", auc_value, "\n")



# ---- 必要パッケージのインストール ----
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")

# ---- ライブラリ読み込み ----
library(survival)
library(survminer)

# ---- 生存時間データの準備 ----
# Survオブジェクトを作成
surv_object <- Surv(time = df$Observed_Time, event = df$CVD_Event)

# ---- ログランク検定 ----
logrank_test <- survdiff(surv_object ~ Smoking, data = df)

print(logrank_test)

cox_model <- coxph(
  formula = surv_object ~ Age + sex_male + BMI + Smoking + Diabetes + 
              LDL + HDL + s_BloodPressure + FamilyHistory, 
  data = df
)

# 生存時間解析##################################################
# Cox regression analysis function

cox_easy=function(x){
a=summary(x)
HR=round(a$ conf.int[,1],digits=2)
LowerCL=round(a$ conf.int[,3],digits=2)
UpperCL=round(a$ conf.int[,4],digits=2)
Pvalue=round(a$ coefficients[,5],digits=3)
result=cbind(HR,LowerCL,UpperCL,Pvalue)
return(result)
}

cox_full=function(x){
a=summary(x)
a$ conf.int[,-2]
a$ coefficients[,-2]
result=round(cbind(a$ conf.int[,-2],a$ coefficients[,-2]),3)
return(result)
}
##############################################################

cox_easy(cox_model)





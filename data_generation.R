

# ---- 必要パッケージのインストール（初回のみ必要）----
install.packages("MASS")
install.packages("truncnorm")

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
ef_base <- 50 + 0.2 * hdl - 0.1 * ldl - 0.03 * age - 3 * smoking - 2 * diabetes + 1.5 * family_history
ef <- round(pmax(20, pmin(ef_base + rnorm(n, mean = 0, sd = 5), 70)))  # 基本範囲20～70
ef[ef == 20] <- ef[ef == 20] + round(rnorm(sum(ef == 20), mean = 2, sd = 3))  # 下限値調整
ef[ef == 70] <- ef[ef == 70] - round(rnorm(sum(ef == 70), mean = 2, sd = 3))  # 上限値調整
ef_outcome <- ifelse(ef < 40, 1, 0)  # EF < 40のアウトカム

# s_BloodPressureの生成
s_BloodPressure <- round(rtruncnorm(n, a = 100, b = 170, mean = 135, sd = 10))  # トランケート分布
s_BloodPressure[s_BloodPressure == 100] <- s_BloodPressure[s_BloodPressure == 100] + round(rnorm(sum(s_BloodPressure == 100), mean = 5, sd = 5))  # 下限値調整
s_BloodPressure[s_BloodPressure == 170] <- s_BloodPressure[s_BloodPressure == 170] - round(rnorm(sum(s_BloodPressure == 170), mean = 5, sd = 5))  # 上限値調整

# ファクターXの生成
factor_x <- rbinom(n, 1, prob = 0.3)  # 二値化された新リスク因子

# CVDイベントの生成
cvd_risk <- -950 + 0.05 * age + 0.4 * sex + 0.3 * bmi + 0.6 * smoking + 0.5 * diabetes + 
            0.1 * ldl - 0.1 * hdl + 0.3 * family_history +
            0.15 * bmi * smoking + 0.1 * ldl * hdl + 0.5 * factor_x
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

# ---- データの保存 ----
write.csv(simulated_data, "simulated_data.csv", row.names = FALSE)

# ---- 検証 ----
summary(simulated_data)  # データの要約
table(simulated_data$CVD_Event)  # CVDイベントの分布確認
table(censoring_flag) / n  # 打ち切りの割合確認
summary(simulated_data$Observed_Time)  # 観察時間の分布確認

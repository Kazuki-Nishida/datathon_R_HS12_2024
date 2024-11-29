# ---- 必要パッケージのインストール ----
if (!require("pROC")) install.packages("pROC")
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")

# ---- ライブラリ読み込み ----
library(pROC)
library(survival)
library(survminer)

# ---- データcsvデータの読み込み----
df=read.csv("simulated_data.csv")

# ---- データ構造と基本情報の表示 ----
dim(df)  # データの行数と列数
str(df)  # 各変数の構造情報
summary(df)  # 変数の要約統計

# ---- ユーティリティ関数定義 ----
# Logistic regression summary function
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

# Cox regression summary function
cox_easy = function(x) {
  a = summary(x)
  HR = round(a$conf.int[,1], digits = 2)
  LowerCL = round(a$conf.int[,3], digits = 2)
  UpperCL = round(a$conf.int[,4], digits = 2)
  Pvalue = round(a$coefficients[,5], digits = 3)
  result = cbind(HR, LowerCL, UpperCL, Pvalue)
  return(result)
}

# ---- RQ1: EF_Outcomeと年齢の関連（t検定） ----
# 目的: EF_Outcomeが1である群と0である群で年齢に差があるかを検討
t.test(df$Age ~ df$EF_Outcome)

# ---- RQ2: EFと年齢の関係（散布図と相関係数） ----
# 目的: EFと年齢の間に関連性があるかを調べる
plot(df$Age, df$EF, main = "Scatterplot of Age vs EF", xlab = "Age", ylab = "EF")
cor.test(df$Age, df$EF)

# ---- RQ3: EF_OutcomeとFamilyHistoryの関係（Fisher検定） ----
# 目的: EF_OutcomeがFamilyHistoryと関連があるかを検討
table_ef_family <- table(df$EF_Outcome, df$FamilyHistory)
fisher.test(table_ef_family)

# ---- RQ4: EFを目的変数とした線形回帰モデル ----
# 目的: EFに影響を与える因子を探索する
linear_model_ef <- lm(EF ~ Age + sex_male + BMI + Smoking + Diabetes + 
                        LDL + HDL + s_BloodPressure + FamilyHistory, data = df)
summary(linear_model_ef)

# ---- RQ5: EF_Outcomeを目的変数としたロジスティック回帰モデル ----
# 目的: EF_Outcomeに影響を与える因子を探索する
logistic_model_ef_outcome <- glm(EF_Outcome ~ Age + sex_male + BMI + Smoking + 
                                   Diabetes + LDL + HDL + s_BloodPressure + FamilyHistory, 
                                 family = binomial(link = "logit"), data = df)
lresult(logistic_model_ef_outcome)

# ---- RQ6: ROC曲線とAUC ----
# 目的: ロジスティック回帰モデルの予測精度を評価する
roc_curve_ef <- roc(df$EF_Outcome, predict(logistic_model_ef_outcome, type = "response"))
plot(roc_curve_ef, main = "ROC Curve for EF_Outcome Model")
cat("AUC:", auc(roc_curve_ef), "\n")

# ---- RQ7: EFとFactorXを含むモデルと含まないモデルの比較 ----
# 目的: EFおよびFactorXを含むモデルの精度を評価し、モデル間で比較する
logistic_model_with_ef_x <- glm(EF_Outcome ~ Age + sex_male + BMI + Smoking + 
                                  Diabetes + LDL + HDL + s_BloodPressure + FamilyHistory + 
                                  EF + FactorX, family = binomial(link = "logit"), data = df)

roc_curve_with_ef_x <- roc(df$EF_Outcome, predict(logistic_model_with_ef_x, type = "response"))

plot(roc_curve_ef, col = "blue", main = "Comparison of ROC Curves")
lines(roc_curve_with_ef_x, col = "red")
legend("bottomright", legend = c("Without EF and FactorX", "With EF and FactorX"), col = c("blue", "red"), lty = 1)

cat("AUC Without EF and FactorX:", auc(roc_curve_ef), "\n")
cat("AUC With EF and FactorX:", auc(roc_curve_with_ef_x), "\n")

# ---- RQ8: FactorXの有無でKaplan-Meierとログランク検定 ----
# 目的: FactorXの有無による生存期間の差を評価
surv_object <- Surv(time = df$Observed_Time, event = df$CVD_Event)
km_fit <- survfit(surv_object ~ FactorX, data = df)
plot(km_fit, col = c("blue", "red"), main = "Kaplan-Meier Curve for FactorX", xlab = "Time (days)", ylab = "Survival Probability")
legend("topright", legend = c("FactorX = 0", "FactorX = 1"), col = c("blue", "red"), lty = 1)

logrank_test <- survdiff(surv_object ~ FactorX, data = df)
print(logrank_test)

# ---- RQ9: Cox比例ハザードモデル ----
# 目的: FactorXおよび他の因子が生存期間に及ぼす影響を評価
cox_model_x <- coxph(surv_object ~ FactorX, data = df)

cox_model_full <- coxph(surv_object ~ Age + sex_male + BMI + Smoking + Diabetes + 
                          LDL + HDL + s_BloodPressure + FamilyHistory + FactorX, data = df)

cox_easy(cox_model_x)
cox_easy(cox_model_full)

# 参考
# LDL単独モデル
cox_model_ldl <- coxph(surv_object ~ LDL, data = df)
# LDL単独モデルの結果を簡潔に表示
cox_easy(cox_model_ldl)

# LDL単独モデルの詳細結果を表示する場合は以下を使用
# summary(cox_model_ldl)
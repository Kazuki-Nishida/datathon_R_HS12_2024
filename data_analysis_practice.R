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

######## EDA ###########
# ---- データ構造と基本情報の表示 ----
dim(df)  # データの行数と列数
str(df)  # 各変数の構造情報
summary(df)  # 変数の要約統計
names(df) # 変数一覧

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
######## EDA end ###########

# ---- ユーティリティ関数定義 ----

# linear regression summary function
lm_result=function(x){
  a=summary(x)
  Coefficients=round(summary(x)$coefficient[,1],digits=2)
  LowerCL=formatC(confint(x)[,1],digits=2, format = "f")
  UpperCL=formatC(confint(x)[,2],digits=2, format = "f")
  Confidence_Intervals=paste(LowerCL,UpperCL,sep=", ")
  CI=paste(matrix("[",length(Confidence_Intervals)),Confidence_Intervals,matrix(rep("]",length(Confidence_Intervals))),sep="")
  Pvalue_raw=round(summary(x)$coefficient[,4],digits=3)
  P_value=ifelse(Pvalue_raw<0.001,"<0.001",Pvalue_raw)
  result=cbind(Coefficients,CI,P_value)
  colnames(result)=c("Coefficients","[95% CI]","     P")
  table=noquote(result)
  attr(table,"note")=paste(note=paste("Objective Variable: ",sub("\\s*\\~.*", "", as.character(x$call[2])),collapse=""))
  return(table)
}

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

# ---- RQ2: EFとLDLの関係（散布図と相関係数） ----
# 目的: EFとLDLの間に関連性があるかを調べる
plot(df$LDL, df$EF, main = "Scatterplot of LDL vs EF", xlab = "LDL", ylab = "EF", pch = 20)
cor.test(df$LDL, df$EF)

# ---- RQ3: EF_OutcomeとFamilyHistoryの関係（Fisher検定） ----
# 目的: EF_OutcomeがFamilyHistoryと関連があるかを検討
table_ef_family <- table(df$EF_Outcome, df$FamilyHistory)
print(table_ef_family)
fisher.test(table_ef_family)

# ---- RQ4: EFを目的変数とした線形回帰モデル ----
# 目的: EFに影響を与える因子を探索する
linear_model_ef <- lm(EF ~ Age + sex_male + BMI + Smoking + Diabetes + 
                        LDL + HDL + s_BloodPressure + FamilyHistory, data = df)
summary(linear_model_ef)
lm_result(linear_model_ef)
write.csv(lm_result(linear_model_ef), "linear_result.csv")

# ---- RQ5: EF_Outcomeを目的変数としたロジスティック回帰モデル ----
# 目的: EF_Outcomeに影響を与える因子を探索する
logistic_model_ef_outcome <- glm(EF_Outcome ~ Age + sex_male + BMI + Smoking + 
                                   Diabetes + LDL + HDL + s_BloodPressure + FamilyHistory, 
                                 family = binomial(link = "logit"), data = df)
lresult(logistic_model_ef_outcome)
write.csv(lresult(logistic_model_ef_outcome), "logistic_result.csv")

# ---- 必要パッケージのインストール ----
if (!require("pROC")) install.packages("pROC")
library(pROC)

# ---- EF_Outcomeを目的変数にしたロジスティック回帰モデル（FactorXなし）----
logistic_model_without_factorx <- glm(
  formula = EF_Outcome ~ Age + sex_male + BMI + Smoking + 
    Diabetes + LDL + HDL + s_BloodPressure + FamilyHistory, 
  family = binomial(link = "logit"), 
  data = df
)

# ---- EF_Outcomeを目的変数にしたロジスティック回帰モデル（FactorXあり）----
logistic_model_with_factorx <- glm(
  formula = EF_Outcome ~ Age + sex_male + BMI + Smoking + 
    Diabetes + LDL + HDL + s_BloodPressure + FamilyHistory + FactorX, 
  family = binomial(link = "logit"), 
  data = df
)

# RQ6: ロジスティック回帰モデルの精度をROC曲線とAUCで評価する
# ---- ROC曲線の計算 ----
roc_without_factorx <- roc(df$EF_Outcome, predict(logistic_model_without_factorx, type = "response"))
roc_with_factorx <- roc(df$EF_Outcome, predict(logistic_model_with_factorx, type = "response"))

# ---- ROC曲線の描画 ----
plot(roc_without_factorx, col = "blue", main = "ROC Curve Comparison for EF_Outcome")
lines(roc_with_factorx, col = "red")
legend("bottomright", legend = c("Without FactorX", "With FactorX"), col = c("blue", "red"), lty = 1)

# ---- AUCの計算 ----
auc_without_factorx <- auc(roc_without_factorx)
auc_with_factorx <- auc(roc_with_factorx)
cat("AUC Without FactorX:", auc_without_factorx, "\n")
cat("AUC With FactorX:", auc_with_factorx, "\n")

# ---- DeLong検定 ----
delong_test <- roc.test(roc_without_factorx, roc_with_factorx, method = "delong")
cat("\n=== DeLong Test for AUC Comparison ===\n")
print(delong_test)

## 参考 ##

# ---- EF_Outcomeを目的変数にしたロジスティック回帰モデル（Smokingあり）----
logistic_model_with_smoking <- glm(
  formula = EF_Outcome ~ Age + sex_male + BMI + Smoking + 
    Diabetes + LDL + HDL + s_BloodPressure + FamilyHistory, 
  family = binomial(link = "logit"), 
  data = df
)

# ---- EF_Outcomeを目的変数にしたロジスティック回帰モデル（Smokingなし）----
logistic_model_without_smoking <- glm(
  formula = EF_Outcome ~ Age + sex_male + BMI + 
    Diabetes + LDL + HDL + s_BloodPressure + FamilyHistory, 
  family = binomial(link = "logit"), 
  data = df
)

# ---- ROC曲線の計算 ----
roc_with_smoking <- roc(df$EF_Outcome, predict(logistic_model_with_smoking, type = "response"))
roc_without_smoking <- roc(df$EF_Outcome, predict(logistic_model_without_smoking, type = "response"))

# ---- ROC曲線の描画 ----
plot(roc_with_smoking, col = "blue", main = "ROC Curve Comparison for EF_Outcome (With/Without Smoking)")
lines(roc_without_smoking, col = "red")
legend("bottomright", legend = c("With Smoking", "Without Smoking"), col = c("blue", "red"), lty = 1)

# ---- AUCの計算 ----
auc_with_smoking <- auc(roc_with_smoking)
auc_without_smoking <- auc(roc_without_smoking)
cat("AUC With Smoking:", auc_with_smoking, "\n")
cat("AUC Without Smoking:", auc_without_smoking, "\n")

# ---- DeLong検定 ----
delong_test_smoking <- roc.test(roc_with_smoking, roc_without_smoking, method = "delong")
cat("\n=== DeLong Test for AUC Comparison (With/Without Smoking) ===\n")
print(delong_test_smoking)



# ---- RQ8: FactorXの有無でKaplan-Meierとログランク検定 ----
# 目的: FactorXの有無による生存期間の差を評価
surv_object <- Surv(time = df$Observed_Time, event = df$CVD_Event)
km_fit <- survfit(surv_object ~ FactorX, data = df)
plot(km_fit, col = c("blue", "red"), main = "Kaplan-Meier Curve for FactorX", xlab = "Time (days)", ylab = "Survival Probability")
legend("topright", legend = c("FactorX = 0", "FactorX = 1"), col = c("blue", "red"), lty = 1)

logrank_test <- survdiff(surv_object ~ FactorX, data = df)
print(logrank_test)

# ---- RQ8続き: 累積ハザードプロット（FactorXによる比較） ----
# 累積ハザードの計算
km_fit_cumhaz <- survfit(surv_object ~ FactorX, data = df, type = "fh")

# 累積ハザードのプロット
plot(
  km_fit_cumhaz, fun = "cumhaz", col = c("blue", "red"), 
  main = "Cumulative Hazard for FactorX", xlab = "Time (days)", ylab = "Cumulative Hazard"
)
legend("topleft", legend = c("FactorX = 0", "FactorX = 1"), col = c("blue", "red"), lty = 1)


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

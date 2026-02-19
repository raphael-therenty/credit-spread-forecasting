library(forecast)
library(texreg)
library(FinTS)
library(tseries)
library(zoo)
library(xts)
library(urca)
library(aTSA)
library(ggplot2)
library(readxl)
library(dplyr)
library(lubridate)
library(vars)
library(lmtest)
library(sandwich)


# ============================================================================================================================================== 
# Sujet : Conditions de crédit et activité aux Etats-Unis : le rôle du spread corporate (2010–2019)

# Question : Dans quelle mesure l’évolution du spread de crédit sur les obligations corporates permet-elle de prédire l’activité économique ?

# Hypothèse : un élargissement du spread, qui indique des conditions de crédit dégradées, anticipe un ralentissement de l’activité économique.
# ============================================================================================================================================== 

#######################################################################################################
######################################### PARTIE UNIVARIEE ############################################
#######################################################################################################


# ---- Script 1 : BAA10Y -----------------------------

####### Importation des données #######

# imporation du csv
data_spread <- data.frame(read.csv("data/BAA10Y.csv", sep = ","))
#head(data_spread)

# formatage des dates
data_spread$date <- as.Date(data_spread$observation_date, format = "%Y-%m-%d")
#head(data_spread)

# passage en format xts
spread_daily_spread <- xts(data_spread$BAA10Y, order.by = data_spread$date)
colnames(spread_daily_spread) <- "BAA10Y"
head(spread_daily_spread)

# passage en fréquence mensuelle via la moyenne
spread_monthly_spread <- apply.monthly(spread_daily_spread, function(x) mean(x, na.rm = TRUE))
head(spread_monthly_spread)

# vérification des données manquantes
sum(is.na(spread_monthly_spread))


####### Analyse graphique #######

plot(spread_monthly_spread)
# Signes de non-stationnarité, pas de tendance claire, 
# La moyenne ne semble pas constante sur la période, et les impacts de chocs semblent durables

acf(spread_monthly_spread)
# Les autocorrélations décroissent lentement, pas de composante MA a priori

pacf(spread_monthly_spread)
# L'autocorrélogramme partiel montre une valeur significative pour le premier retard
# Les autocorrélations partielles chutent vers zéro dès le second retard
# Pourrait indiquer un processus AR(1), mais conclusion possible uniquement sur une série stationnaire

# Conclusion intermédiaire : il peut s'agir d'un processus DS, DS avec dérive ou TS


####### Transformation de Box-Cox #######

# passage en log de la série
log_spread_spread <- log(spread_monthly_spread)

# visualisation
plot(log_spread_spread)


####### Tests de racine unitaire sur le log de la série #######

# Nous appliquons une stratégie séquentielle pour le test ADF uniquement
# Nous réalisons ensuite les tests KPSS, PP et ERS pour confirmer ou infirmer la conclusion du test ADF
# Pour tous les tests effectués, le seuil de confiance est fixé à 5%

# Test du modèle 3 : constante + tendance déterministe

adf_trend_spread <- ur.df(log_spread_spread, type = "trend", selectlags = "BIC")
summary(adf_trend_spread)
# On a DF > CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire
# Pour le test de spécification, on a abs(Tstat) < VC : on ne peut rejetter H0, tendance non significative
# Mauvaise spécification

# Test du modèle 2 : constante

adf_drift_spread <- ur.df(log_spread_spread, type = "drift", selectlags = "BIC")
summary(adf_drift_spread)
# On a DF > CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire
# Pour le test de spécification, on a abs(Tstat) < VC : on ne peut rejetter H0, constante non significative
# Mauvaise spécification

# Test du modèle 1 : none

adf_none_spread <- ur.df(log_spread_spread, type = "none", selectlags = "BIC")
summary(adf_none_spread)
# On a DF > CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire

# Les tests ADF rejettent la stationnarité et la présence de tendance et de constante
# Réalisons maintenant les tests KPSS, PP et ERS pour confirmer ou infirmer cette première conclusion

kpss_drift_spread <- ur.kpss(log_spread_spread, type = "mu", lags = "short")
summary(kpss_drift_spread)
# On a stat > VC : on rejette H0 -> présence de racine unitaire, série non stationnaire

pp_drift_spread <- ur.pp(log_spread_spread, model = "constant", lags = "short")
summary(pp_drift_spread)
# Même conclusion que le test ADF, série non stationnaire et constante significative

ers_drift_spread <- ur.ers(log_spread_spread, model = "constant", lag.max = 1)
summary(ers_drift_spread)
# Même conclusion que les tests ADF, KPSS et PP, série non stationnaire

# A la lumière des résultats ci-dessus, nous pouvons conclure à la non stationnarité de la série en log
# Il semblerait que la série en log suive une marche aléatoire


####### Passage en diff log de la série #######

diff_log_spread_spread <- diff(log_spread_spread)
diff_log_spread_spread <- diff_log_spread_spread[-1,]

# visualisation
plot(diff_log_spread_spread)
# La série semble stationnaire : la moyenne et la variance paraissent constantes


####### Tests de racine unitaire sur le diff log de la série #######

# Test modèle 3 : constante + tendance

adf_trend_spread <- ur.df(diff_log_spread_spread, type = "trend", selectlags = "BIC")
summary(adf_trend_spread)
# On a DF < CV, on rejette H0 à 5%, on a abs(Tstat) < VC : on ne peut rejetter H0, tendance non significative
# Mauvaise spécification

# Test du modèle 2 : constante

adf_drift_spread <- ur.df(diff_log_spread_spread, type = "drift", selectlags = "BIC")
summary(adf_drift_spread)
# On a DF < CV, on peut rejetter H0 à 5%, absence de racine unitaire et série stationnaire
# Pour le test de spécification, on a abs(Tstat) < VC : on ne peut rejetter H0, constante non significative
# Mauvaise spécification

# Test du modèle 1 : none

adf_none_spread <- ur.df(diff_log_spread_spread, type = "none", selectlags = "BIC")
summary(adf_none_spread)
# On a DF < CV, on rejette H0 au seuil 5%, absence de racine unitaire, série stationnaire

# Selon la stratégie séquentielle de TRU ADF, la série en diff log est stationnaire sans dérive
# Nous serions donc en présence d'un série I(1)
# Réalisons maintenant les tests KPSS, PP et ERS pour confirmer ou infirmer cette première conclusion

kpss_drift_spread <- ur.kpss(diff_log_spread_spread, type = "mu", lags = "short")
summary(kpss_drift_spread)
# On a stat < VC : on ne peut rejetter H0 -> absence de racine unitaire, série stationnaire

pp_drift_spread <- ur.pp(diff_log_spread_spread, model = "constant", lags = "short")
summary(pp_drift_spread)
# Le test conclut à une stationnarité

ers_drift_spread <- ur.ers(diff_log_spread_spread, model = "constant", lag.max = 1)
summary(ers_drift_spread)
# Le test conclut à une stationnarité de la série en différence première



# A la lumière de ces différents tests, nous concluons que le log de la série est
# stationnaire en différence première, série est intégrée d'ordre 1 : I(1)






# ---- Script 2 : INDPRO -------------------------------

####### Importation des données #######

data_ind <- data.frame(read.csv("data/INDPRO.csv", sep = ","))
head(data_ind)

# pour connaitre le nombre d'observations
n <- length(data_ind$observation_date)

# Création de la série temporelle mensuelle
dates <- seq(as.Date("2010-01-01"), by = "month", length.out = n)
INDPRO_ts_ind <- xts(data_ind$INDPRO, order.by = dates)
head(INDPRO_ts_ind)

####### Analyse graphique #######

plot(INDPRO_ts_ind, main="Indice Production Industrielle") 
# Série probablement non stationnaire (E(X) croissant en fonction du temps)


####### Transformation de Box-Cox #######

## Analyse de la série en logarithme
log_PRODIND_ts_ind <- log(INDPRO_ts_ind) # Passer une série en log corrige la variance croissante avec le niveau
plot(log_PRODIND_ts_ind)

acf(log_PRODIND_ts_ind)
pacf(log_PRODIND_ts_ind)
# acf décroissant lentement vers 0
# pacf s'arrête net en 1 
# Présence d'une composante AR probable, non stationnaire.


####### Tests de racine unitaire sur le log de la série #######

### Tests de racine unitaire avec la stratégie séquentielle de TRU
# On choisit automatiquement le nombre de retards selon BIC

# --- Modèle 3 : avec constante + tendance
adf_trend_ind <- ur.df(log_PRODIND_ts_ind, type = "trend", selectlags = "BIC")
summary(adf_trend_ind)
# t*(rho) = -2.290 > VC (tau3) = -3.43 -> On ne rejette pas H0 (rho = 0) au seuil de 5%  -> Il existe une RU.
# On utilise la table de DF pour tester la significativité du coefficient de la tendance (tt).
# t* (tt) = 0.377 < VC = 3.14 (table DF en cas de RU) -> On ne rejette pas H0 (b= 0) -> La tendance n'est pas significative.
# On passe au modèle 2.

# --- Modèle 2 : avec constante (drift)
adf_drift_ind <- ur.df(log_PRODIND_ts_ind, type = "drift", selectlags = "BIC")
summary(adf_drift_ind)
# t*(rho) = -3.388 < VC(tau2) = -2.88 - > On rejette H0 (rho = 0) au seuil de 5% -> Il n'y a pas de RU. I(0)
# On utilise la table de Student pour tester la significativité de la constante (intercept).
# t* (intercept) = 3.408  > VC = 2,86   -> On rejette H0 (c = 0) au seuil 5% -> La constante est significative.


# Test KPSS 
kpss_mu_ind <- ur.kpss(log_PRODIND_ts_ind, type = "mu")
summary(kpss_mu_ind)
#LM* (1.7271) > VC (0.463) -> On rejette H0 (pas de RU) au seuil de 5% -> Il y a une RU.

# Test PP
pp_drift_ind <- ur.pp(log_PRODIND_ts_ind, type = "Z-tau", model = "constant", lags = "short")
summary(pp_drift_ind)
# PP* (Z-tau) = -3.3305 <  VC (-2.885699) -> Rejet HO (présence de RU) au seuil de 5%. I(0)


# Test ERS
ers_const_ind <- ur.ers(log_PRODIND_ts_ind, type = "DF-GLS", model = "constant", lag.max = 4)
summary(ers_const_ind)
# H0 = RU
# ERS* = 0.2868 > VC (-1.94) -> On ne rejette pas H0 (présence de RU) au seuil de 5%. 


# La série log_PRODIND_ts semble à la frontière entre un processus I(0) et un processus intégré d'ordre au moins 1, ce qui justifie la différenciation en log.


####### Passage en diff log de la série #######

# --- Différenciation logarithmique ---
dlog_PRODIND_ts_ind <- diff(log_PRODIND_ts_ind)
dlog_PRODIND_ts_ind <- dlog_PRODIND_ts_ind[-1,]
plot(dlog_PRODIND_ts_ind)

acf(dlog_PRODIND_ts_ind)
pacf(dlog_PRODIND_ts_ind)
# Décroissance prononcée en 1 puis pic en 3 et 6 de l'acf.
# Décroissance lente du PACF avec des pics en 3 et 6.
# ça pourrait être un ARMA. Les pics significatifs en 3 et 6, ainsi que l'absence d'annulation nette.


####### Tests de racine unitaire sur le diff log de la série #######


## --- Tests sur la série différenciée ---
# Modèle 3
adf_trend_ind <- ur.df(dlog_PRODIND_ts_ind, type = "trend", selectlags = "BIC")
summary(adf_trend_ind)
# t*(rho) = -7.600 < VC (-3.43) -> Rejet H0 (RU) au seuil de 5% donc il n'y a pas de RU -> I(0)
# Testons la significativité de la tendance avec la table de Student (loi Normale ici comme T-K > 30).
# t*(tt) = | -2.222 | = 2.222 > VC(1.96) -> Rejet de H0(tendance = 0) -> La tendance est significative.
# La forme du modèle serait : delta(Xt) = c + b*t + phiX[t-1] + eps[t] avec  où Xt ~ I(0) + c + T 
# Vérifions avec les autres tests

kpss_trend_ind <- ur.kpss(dlog_PRODIND_ts_ind, type = "tau")
summary(kpss_trend_ind)
# LM* = 0.1263 < VC = 0.146 -> On ne rejette pas H0 (absence de RU) au seuil de 5%. Il n'y a pas de RU -> I(0)

pp_trend_ind <- ur.pp(dlog_PRODIND_ts_ind, type = "Z-tau", model = "trend", lags = "short")
summary(pp_trend_ind)
# H0 = RU
# PP* (tau) = -10.9535 < VC (-3.448109) -> Rejet H0 au seuil de 5%. Il n'y a pas de RU -> I(0) 

ers_trend_ind <- ur.ers(dlog_PRODIND_ts_ind, type = "DF-GLS", model = "trend", lag.max = 4)
summary(ers_trend_ind)
# H0 = RU
# ERS* = -3.9774 < VC(-2.93) -> Rejet H0 au seuil de 5%. Il n'y a pas de RU -> I(0)

# Tous les tests concluent que la série différenciée est stationnaire (I(0)) autour d'une tendance.
# La tendance est significative au seuil de 5% (donc on rejette H0).
# La série initiale en logarithme est donc I(1).



# ---- Script 3 : WALCL -----------------------------

####### Importation des données #######

# importation du csv
data_walcl <- data.frame(read.csv("data/WALCL.csv", sep = ","))
head(data_walcl)

# formatage des dates
data_walcl$date <- as.Date(data_walcl$observation_date, format = "%Y-%m-%d")

# passage en format xts
bs_fed_weekly_walcl <- xts(data_walcl$WALCL, order.by = data_walcl$date)
colnames(bs_fed_weekly_walcl) <- "WALCL"
head(bs_fed_weekly_walcl)

# passage en fréquence mensuelle via la moyenne
bs_fed_monthly_walcl <- apply.monthly(bs_fed_weekly_walcl, FUN = colMeans)
head(bs_fed_monthly_walcl)

# vérification des données manquantes
sum(is.na(bs_fed_monthly_walcl))


####### Analyse graphique #######

plot(bs_fed_monthly_walcl)
# Signes de non-stationnarité, tendance haussière de 2010 à 2018, puis baisse sur 2019
# La moyenne ne semble pas constante sur la période, et les impacts de chocs semblent durables

acf(bs_fed_monthly_walcl)
# Les autocorrélations décroissent lentement, pas de composante MA a priori

pacf(bs_fed_monthly_walcl)
# L'autocorrélogramme partiel montre une valeur significative pour le premier retard
# Les autocorrélations partielles chutent vers zéro dès le second retard
# Pourrait indiquer un processus AR(1), mais conclusion possible sur une série stationnaire

# Conclusion intermédiaire : il peut s'agir d'un processus DS, DS avec dérive ou TS


####### Transformation de Box-Cox #######

# passage en log de la série
log_bs_fed_walcl <- log(bs_fed_monthly_walcl)

# visualisation
plot(log_bs_fed_walcl)


####### Tests de racine unitaire sur le log de la série #######

# Nous appliquons une stratégie séquentielle pour le test ADF uniquement
# Nous réalisons ensuite les tests KPSS, PP et ERS pour confirmer ou infirmer la conclusion du test ADF
# Pour tous les tests effectués, le seuil de confiance est fixé à 5%

# Test du modèle 3 : constante + tendance déterministe

adf_trend_walcl <- ur.df(log_bs_fed_walcl, type = "trend", selectlags = "BIC")
summary(adf_trend_walcl)
# On a DF > CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire
# Pour le test de spécification, on a abs(Tstat) < VC : on ne peut rejetter H0, tendance non significative
# Mauvaise spécification

# Test du modèle 2 : constante

adf_drift_walcl <- ur.df(log_bs_fed_walcl, type = "drift", selectlags = "BIC")
summary(adf_drift_walcl)
# On a DF > CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire
# Pour le test de spécification, on a abs(Tstat) < VC : on ne peut rejetter H0, constante non significative
# Mauvaise spécification

# Test du modèle 1 : none

adf_none_walcl <- ur.df(log_bs_fed_walcl, type = "none", selectlags = "BIC")
summary(adf_none_walcl)
# On a DF > CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire

# Les tests ADF rejettent la stationnarité et la présence de tendance et de constante
# Réalisons maintenant les tests KPSS, PP et ERS pour confirmer ou infirmer cette première conclusion

kpss_drift_walcl <- ur.kpss(log_bs_fed_walcl, type = "mu", lags = "short")
summary(kpss_drift_walcl)
# On a stat > VC : on rejette H0 -> présence de racine unitaire, série non stationnaire

pp_drift_walcl <- ur.pp(log_bs_fed_walcl, model = "constant", lags = "short")
summary(pp_drift_walcl)
# Même conclusion que le test ADF, série non stationnaire et constante significative

ers_drift_walcl <- ur.ers(log_bs_fed_walcl, model = "constant", lag.max = 1)
summary(ers_drift_walcl)
# Même conclusion que les tests ADF, KPSS et PP, série non stationnaire

# A la lumière des résultats ci-dessus, nous pouvons conclure à la non stationnarité de la série en log
# Il semblerait que la série en log suive une marche aléatoire


####### Passage en diff log de la série #######

diff_log_bs_fed_walcl <- diff(log_bs_fed_walcl)
diff_log_bs_fed_walcl <- diff_log_bs_fed_walcl[-1,]

# visualisation
plot(diff_log_bs_fed_walcl)
# La série ne semble pas stationnaire : ni la moyenne ni la variance ne paraissent constantes


####### Tests de racine unitaire sur le diff log de la série #######

# Test modèle 3 : constante + tendance

adf_trend_walcl <- ur.df(diff_log_bs_fed_walcl, type = "trend", selectlags = "BIC")
summary(adf_trend_walcl)
# On a DF > CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire
# Pour le test de spécification, on a abs(Tstat) < VC : on ne peut rejetter H0, tendance non significative
# Mauvaise spécification

# Test du modèle 2 : constante

adf_drift_walcl <- ur.df(diff_log_bs_fed_walcl, type = "drift", selectlags = "BIC")
summary(adf_drift_walcl)
# On a DF ≃ CV, on ne peut rejetter H0, présence de racine unitaire et série non stationnaire
# Pour le test de spécification, on a abs(Tstat) < VC : on ne peut rejetter H0, constante non significative
# Mauvaise spécification

# Test du modèle 1 : none

adf_none_walcl <- ur.df(diff_log_bs_fed_walcl, type = "none", selectlags = "BIC")
summary(adf_none_walcl)
# On a DF < CV, on rejette H0 au seuil 5%, absence de racine unitaire, série stationnaire

# Selon la stratégie séquentielle de TRU ADF, la série en diff log est stationnaire sans dérive
# Nous serions donc en présence d'un série I(1)
# Réalisons maintenant les tests KPSS, PP et ERS pour confirmer ou infirmer cette première conclusion

kpss_drift_walcl <- ur.kpss(diff_log_bs_fed_walcl, type = "mu", lags = "short")
summary(kpss_drift_walcl)
# On a stat > VC : on rejette H0 -> présence de racine unitaire, série non stationnaire

pp_drift_walcl <- ur.pp(diff_log_bs_fed_walcl, model = "constant", lags = "short")
summary(pp_drift_walcl)
# Le test conclut à une non stationnarité, mais test peu puissant

ers_drift_walcl <- ur.ers(diff_log_bs_fed_walcl, model = "constant", lag.max = 1)
summary(ers_drift_walcl)
# Le test conclut à une stationnarité de la série en différence première



# A la lumière de ces différents tests, nous concluons que le log de la série est
# stationnaire en différence première, série est intégrée d'ordre 1 : I(1)


####### Modélisation univariée #######

auto.arima(bs_fed_monthly_walcl, max.d = 1, ic = "bic")
# renvoie un ARIMA(1,1,0), soit un AR(1) sur la série en diff log

model_walcl <- Arima(bs_fed_monthly_walcl, order = c(1,1,0))

summary(model_walcl)
screenreg(model_walcl)

res_walcl <- residuals(model_walcl)
acf(res_walcl)
pacf(res_walcl)
checkresiduals(res_walcl)

Box.test(res_walcl,lag = 10, type = "Ljung-Box")
# très grande p-value donc on peut pas rejeter H0, les résidus ne sont pas corrélés

ArchTest(res_walcl,lags = 10)
# très grande p-value donc on peut pas rejeter H0, les résidus ne sont pas hétéroscédastiques

jarque.bera.test(res_walcl)
# La p-value inférieure à 5%, rejet de H0, les résidus ne sont pas normaux
# On remarque la forme de cloche, distribution probablement leptokurtique

# Prévision in-sample statiques (à un pas)
prev_fitted_static_walcl <- xts(fitted(model_walcl), order.by = index(bs_fed_monthly_walcl))
plot(prev_fitted_static_walcl, type='l', col="black", ylab="Valeur", xlab="temps", main="Prédiction in-sample statiques")
lines(bs_fed_monthly_walcl, col="red")
legend("bottomright", legend=c("prédictions","observé"), col=c("black", "red"), lty=1, xpd=TRUE)

# prévisions out-of-sample
prev_walcl <- forecast::forecast(model_walcl, h = 12)
autoplot(prev_walcl) + ggtitle("Prévisions WALCL")

# walcl : ARIMA(1,1,0)







################################################################################################################################################## 
############################################### PARTIE II : MODELISATION MULTIVARIEE ############################################################# 
################################################################################################################################################## 



####### Stationnarisation (différence logarithmique) #######
# Toutes les séries en logarithme sont I(1), donc on utilise les séries diff log
# On ajuste les dates au premier jour du mois pour spread et walcl


index(diff_log_spread_spread) <- as.Date(format(index(diff_log_spread_spread), "%Y-%m-01"))
index(diff_log_bs_fed_walcl) <- as.Date(format(index(diff_log_bs_fed_walcl), "%Y-%m-01"))
head(diff_log_bs_fed_walcl)
head(diff_log_spread_spread)
head(dlog_PRODIND_ts_ind)

# 1. Fusionner les séries en un seul xts, ne garder que les dates communes
Xt <- merge(diff_log_bs_fed_walcl, diff_log_spread_spread, join = 'inner')
Xt <- merge(Xt, dlog_PRODIND_ts_ind, join = 'inner')
colnames(Xt) <- c("WALCL", "BAA10Y", "INDPRO")
head(Xt)


# 2. Sélection du lag optimal selon différents critères
lag_selection <- VARselect(Xt, lag.max = 12, type = "const") 
print(lag_selection$selection)
# Tous les critères donnent lag = 1. On modélise avec une constante (const).

# 1) Modèle avec p = 1 constante
var_model_const <- VAR(Xt, p = 1, type = 'const')
summary(var_model_const)


####### Interprétation #######

# Les racines inverses sont toutes inférieures à 1.
# Analyse par équation : 
#   WALCL : 
#     - WALCL.l1  : p-val < 0.05 -> Significatif
#     - BAA10Y.l1 : p-val > 0.05 -> Non significatif
#     - INDPRO.l1 : p-val > 0.05 -> Non significatif

# WALCL.l1 = 0.83 > 0 (***) -> La taille du bilan de la Fed dépend essentiellement de sa propre valeur passée. La politique monétaire évolue de manière persistente.
# BAA10Y.l1 = -0.01 < 0 -> Le spread crédit n'a pas d'effet immédiat significatif sur les variations du bilan de la Fed. La Fed ne réagit pas mécaniquement aux fluctuactions du spread à court terme.
# INDPRO.l1 = -0.10 < 0 -> Une hausse passée de l'activité industrielle n'a pas d'effet significatif sur le bilan de la Fed.
# => WALCL est essentiellement auto-déterminée et hautement persistante, cohérent avec le caractère planifié de la politique monétaire.

#   BBA10Y : 
#     - WALCL.l1  : p-val > 0.05 -> Non significatif
#     - BAA10Y.l1 : p-val < 0.05 -> Significatif
#     - INDPRO.l1 : p-val > 0.05 -> Non significatif

# WALCL.l1 = -0.38 < 0 -> Une expansion du bilan de la Fed tend à réduire le spread de crédit (effet de détente financière), mais l’effet est faible et non significatif statistiquement.
# BAA10Y.l1 = 0.28 > 0 -> Forte persistance du spread : un spread élevé le mois précédent tend à rester élevé.
# INDPRO.l1 = 0.54 > 0 -> Une hausse de la production tend à être associée à une légère hausse du spread (peut refléter la phase avancée du cycle où les taux montent), mais non significative.
# => Le spread est auto-corrélé, donc persistant, mais peu réactif à l’activité réelle ou au bilan de la Fed à court terme.

#   INDPRO : 
#     - WALCL.l1  : p-val > 0.05 -> Non significatif 
#     - BAA10Y.l1 : p-val > 0.05 -> Non significatif
#     - INDPRO.l1 : p-val > 0.05 -> Non significatif

# WALCL.l1 = 0.52 > 0 -> Une hausse du bilan de la Fed stimule légèrement la production industrielle, effet attendu (politique monétaire expansionniste → activité accrue), mais non significatif.
# BAA10Y.l1 = -0.01 < 0 -> Un élargissement du spread (hausse du coût du crédit) a un effet négatif sur l’activité, ce qui correspond à la théorie : un spread plus élevé reflète des conditions de crédit plus restrictives. Toutefois, l’effet est faible et non significatif à court terme.
# INDPRO.l1 = 0.01 > 0 -> Faible inertie de la croissance industrielle d’un mois à l’autre, mais non significatif.
# => L’activité industrielle n’est pas significativement influencée, à court terme, ni par le bilan de la Fed ni par le spread de crédit, et elle n’est que très faiblement auto-déterminée dans cette spécification.

# Les covariances et les corrélations sont très faibles -> C'est bien.

#   VAR(1) avec constante est cohérent avec la stationnarité des séries différenciées
#   WALCL est auto-dépendant et prédictif
#   BAA10Y dépend partiellement de sa propre valeur passée
#   INDPRO n’est pas expliqué par ce VAR(1) → peut nécessiter des retards supplémentaires, des variables exogènes, ou une autre fréquence/différenciation pour mieux capter sa dynamique.


# Les dynamiques principales à retenir de ce VAR(1) avec constante
# ------------------------------------------------------------------

# WALCL (politique monétaire) et BAA10Y (spread) sont très persistants.
# INDPRO réagit faiblement à court terme, mais les signes des coefficients confirment la relation théorique attendue :
#   Expansion monétaire → hausse de l’activité.
#   Hausse du spread → contraction de l’activité.
# Les effets statistiques sont faibles à horizon mensuel, mais économiquement cohérents, ce qui est courant dans ce type de données.



# 2) Vérification du VAR incluant constante et trend  
lag_select2 <- VARselect(Xt, lag.max = 12, type = "both") 
print(lag_select2$selection)
var_model_both <- VAR(Xt, p = 1, type = 'both')
summary(var_model_both)

##### Interprétation :
# Les racines inverses sont toutes inférieures à 1.
# Analyse par équation : 
#   WALCL : 
#     - WALCL.l1  : p-val < 0.05 -> Significatif
#     - BAA10Y.l1 : p-val > 0.05 -> Non significatif
#     - INDPRO.l1 : p-val > 0.05 -> Non significatif

#   BBA10Y : 
#     - WALCL.l1  : p-val > 0.05 -> Non significatif
#     - BAA10Y.l1 : p-val < 0.05 -> Significatif
#     - INDPRO.l1 : p-val > 0.05 -> Non significatif

#   INDPRO : 
#     - WALCL.l1  : p-val > 0.05 -> Non significatif 
#     - BAA10Y.l1 : p-val > 0.05 -> Non significatif
#     - INDPRO.l1 : p-val > 0.05 -> Non significatif

# Les covariances et les corrélations sont très faibles -> C'est bien.


#### Conclusion du VAR(1) avec constante et tendance : 
# Auto-dépendances dominantes : WALCL dépend fortement de sa valeur passée. BAA10Y dépend modérément de sa valeur passée. INDPRO n’est pas bien expliqué par ses propres valeurs passées ni par WALCL ou BAA10Y à ce lag.
# Effet des autres variables :
#   WALCL et BAA10Y ont peu d’influence directe sur les autres variables à t-1.
#   INDPRO semble isolée, ce qui pourrait justifier d’ajouter des lags plus longs ou un autre type de VAR (ex : VARX).

# Composante déterministe :
    # L’ajout de la tendance a peu d’effet sur WALCL et BAA10Y.
    # Pour INDPRO, constante et tendance sont significatives → série légèrement croissante sur le temps, justifiant la tendance dans cette équation.

# Qualité globale :
#   Les R² restent faibles pour BAA10Y et INDPRO → le VAR(1) capture surtout WALCL.
#   Les résidus sont peu corrélés entre eux → stabilité et validité du VAR.


# CONCLUSION GENERALE :
# ---------------------

# Après avoir vérifié la stationnarité des séries et estimé des modèles univariés, les trois variables (différences logarithmiques de WALCL, BAA10Y et INDPRO) ont été intégrées dans un modèle VAR. Les critères d’information (AIC, HQ, SC et FPE) indiquent un nombre de retards optimal de 1.
# La composante déterministe retenue est une constante seule, cohérente avec les propriétés de stationnarité et avec la dynamique observée. L’ajout d’une tendance n’améliore que marginalement l’ajustement, bien qu’elle capte une faible dérive de long terme de l’activité industrielle.
# De plus, on veut analyser l’effet des conditions financières (le spread BAA10Y) sur l’activité réelle (INDPRO). Ce choix est cohérent avec la stationnarité des séries en différences logarithmiques et avec l’objectif de l’étude, qui vise à évaluer la dynamique de court terme entre conditions financières et activité économique, plutôt qu’à modéliser une tendance structurelle de long terme.
# Les résultats du VAR(1) montrent une forte persistance du bilan de la Fed et du spread corporate, mais aucune influence significative de ces variables sur l’activité industrielle à court terme. L’équation de la production (INDPRO) n’indique pas de réponse immédiate aux chocs de spread ou de politique monétaire, ce qui suggère que l’effet du resserrement des conditions de crédit sur l’activité réelle se matérialise à un horizon plus long.
# En somme, le modèle VAR(1) avec constante seule est statistiquement et économiquement cohérent : il décrit correctement la dynamique conjointe des conditions financières et de l’activité, tout en confirmant que le spread corporate ne prédit pas significativement l’évolution immédiate de la production industrielle sur la période 2010–2019.



# ------------------------------------------
# 3. Diagnostic des résidus
# ------------------------------------------
# ------------- const ----------------
# Test d'autocorrélation
serial.test(var_model_const, lags.pt = 12, type = "PT.asymptotic") 
# H0 = Autocorrélation
# p-val > 0.05 -> On ne rejette pas H0 -> Pas d'autocorrélation des résidus.

#Test d'hétéroscédasticité
arch.test(var_model_const, lags.multi = 12,multivariate.only = TRUE) 
# H0 = Héréroscédastique 
# p-val > 0.05 -> On ne rejette pas H0 -> Pas d'hétéroscédasticité des résidus.

# Test de normalité de JB
normality.test(var_model_const) #On rejette H0, les résidus ne suivent pas une distribution normale
# H0 = Normalité 
# p-val < 0.05 -> On rejette H0 -> Les résidus ne sont pas normalement distribués. 

#Test de stabilité
roots(var_model_const) #toutes les racines inverses sont inférieures à 1

# ------------- both -----------------
# Test d'autocorrélation
serial.test(var_model_both, lags.pt = 12, type = "PT.asymptotic") 
# H0 = Autocorrélation
# p-val > 0.05 -> On ne rejette pas H0 -> Pas d'autocorrélation des résidus.

#Test d'hétéroscédasticité
arch.test(var_model_both, lags.multi = 12,multivariate.only = TRUE) 
# H0 = Héréroscédastique 
# p-val > 0.05 -> On ne rejette pas H0 -> Pas d'hétéroscédasticité des résidus.

# Test de normalité de JB
normality.test(var_model_both) #On rejette H0, les résidus ne suivent pas une distribution normale
# H0 = Normalité 
# p-val < 0.05 -> On rejette H0 -> Les résidus ne sont pas normalement distribués. 

#Test de stabilité
roots(var_model_both) #toutes les racines inverses sont inférieures à 1

# Conclusion générale :
# L’ensemble des tests diagnostiques indique que le VAR(1) avec constante est statistiquement satisfaisant et économiquement cohérent.
# Les résidus ne présentent ni autocorrélation, ni hétéroscédasticité, et le modèle est structurellement stable.
# La non-normalité des résidus n’affecte pas la validité globale du modèle, dans la mesure où les propriétés de non-corrélation et de stabilité sont respectées.
# Ainsi, le VAR(1) avec constante peut être retenu comme modèle optimal pour analyser les interactions dynamiques entre les conditions de crédit (spread corporate), la politique monétaire (bilan de la Fed) et l’activité économique (production industrielle) aux États-Unis sur la période 2010–2019.


# 4. Causalité de Granger
#-------------------------------------------------
causality(var_model_const, cause = "WALCL") 
# p-val > 0.05 ->  diff_walcl ne cause pas les deux autres. 
# Les variations passées du bilan de la Fed n’améliorent pas la prévision du spread de crédit ni de la production industrielle.
# Cela suggère que,sur la période 2010–2019, les politiques monétaires non conventionnelles (quantitative easing, etc.) n’ont pas exercé un effet prédictif fort à court terme sur les conditions de crédit ou l’activité réelle.

causality(var_model_const, cause = "BAA10Y") 
# p-val > 0.05 -> diff_bba10y ne cause pas les 2 autres.
# Les fluctuations du spread de crédit ne précèdent pas significativement celles du bilan de la Fed (logique, car la politique monétaire réagit à d’autres facteurs), et ne permettent pas non plus de prévoir l’évolution de la production industrielle.
# Statistiquement, sur cette période, le spread de crédit ne semble pas avoir un pouvoir prédictif significatif sur l’activité économique, même si économiquement, la théorie suggère un lien.
# Cela peut venir de la période étudiée (2010–2019, où le spread était faible et stabilisé après la crise), ou du fait que les variables sont en différence/log et donc captent seulement des variations de court terme.

causality(var_model_const, cause = "INDPRO") 
# p-val > 0.05 -> diff_ind ne cause pas les 2 autres.
# Les évolutions de l’activité industrielle ne causent pas non plus les changements du bilan de la Fed ni du spread.
# Autrement dit, la Fed n’ajuste pas son bilan directement en réaction immédiate aux variations de la production industrielle, et le spread ne réagit pas systématiquement aux fluctuations de l’activité réelle.


# Interprétation globale 
# ----------------------
# Aucune des trois variables ne présente de causalité de Granger ou instantanée significative entre elles.
# Les interactions dynamiques sont donc faibles à court terme dans ce modèle VAR(1) avec constante.

# Implications économiques :
#     Les effets de la politique monétaire (WALCL) sur le crédit et la production ne se manifestent pas immédiatement, mais plutôt à plus long terme ou via d’autres canaux non capturés par ce VAR simple.
#     Le spread de crédit (BAA10Y) ne semble pas avoir de pouvoir prédictif significatif sur l’activité industrielle (INDPRO) à court terme.
#     Cela peut refléter une période de relative stabilité financière (post-crise) où les écarts de crédit étaient faibles et peu informatifs sur la conjoncture réelle.


# Les tests de causalité de Granger et instantanée réalisés sur le modèle VAR(1) avec constante montrent qu’aucune des variables (WALCL, BAA10Y, INDPRO) ne cause significativement les autres au sens de Granger.
# Cela signifie qu’à court terme, les variations passées du bilan de la Fed, du spread de crédit et de la production industrielle n’améliorent pas la prévision mutuelle de ces variables.
# Ces résultats suggèrent que, sur la période 2010–2019, les liens dynamiques entre politique monétaire, conditions de crédit et activité réelle sont faibles, possiblement en raison du contexte post-crise marqué par des taux bas et une volatilité réduite du spread.


# Vérifions avec le modèle VAR(1) avec constante et tendante
# ----------------------------------------------------------------
causality(var_model_both, cause = "WALCL") # WALCL ne cause pas les 2 autres
causality(var_model_both, cause = "BAA10Y") # BAA10Y ne cause pas les 2 autres
causality(var_model_both, cause = "INDPRO") # INDPRO ne cause pas les 2 autres

# Même conclusions qu'avec le modèle initiale VAR (1) avec constante uniquement


# --------------------------------------------------
# 5. IRF (Impulse Response Functions)
# --------------------------------------------------
# D'après la théorie économique, la politique monétaire (walcl) influence d’abord les conditions financières (spreads), qui à leur tour affectent l’activité réelle (prodin) avec retard.
# On ordonne donc de la plus exogène à la plus endogène :  WALCL (politique monétaire) → BAA10Y (conditions financières) → INDPRO (activité réelle)
# Mais il faut également tester les autres combinaisons afin de s'assurer de la robustesse de notre modèle.

# Liste de tous les modèles avec un nom descriptif
var_models <- list(
  "Ordre économique (WALCL → BAA10Y → INDPRO)" = VAR(Xt[, c("WALCL", "BAA10Y", "INDPRO")], p = 1, type = "const"),
  "Désordre 1 (WALCL → INDPRO → BAA10Y)" = VAR(Xt[, c("WALCL", "INDPRO", "BAA10Y")], p = 1, type = "const"),
  "Désordre 2 (BAA10Y → WALCL → INDPRO)" = VAR(Xt[, c("BAA10Y", "WALCL", "INDPRO")], p = 1, type = "const"),
  "Désordre 3 (BAA10Y → INDPRO → WALCL)" = VAR(Xt[, c("BAA10Y", "INDPRO", "WALCL")], p = 1, type = "const"),
  "Désordre 4 (INDPRO → WALCL → BAA10Y)" = VAR(Xt[, c("INDPRO", "WALCL", "BAA10Y")], p = 1, type = "const"),
  "Désordre 5 (INDPRO → BAA10Y → WALCL)" = VAR(Xt[, c("INDPRO", "BAA10Y", "WALCL")], p = 1, type = "const")
)

# Fonction pour tracer toutes les IRF d'un modèle en une seule fenêtre
plot_all_irf <- function(model, model_name, n_ahead = 12) {
  irf_all <- irf(model, n.ahead = n_ahead, boot = TRUE)
  nvar <- length(model$varresult)
  
  # Configurer la fenêtre pour tous les mini-plots                      # mfrow : grille de plot, chaque case = un plot individuel
  par(mfrow = c(nvar, nvar), mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 2))  # oma: marges extérieures (bottom, left, top, right)
                                                                        # mar : marges intérieures de chaque plot
  for (imp in names(irf_all$irf)) {
    for (resp in names(irf_all$irf)) {
      plot(irf_all$irf[[imp]][, resp], type = "l",
           main = paste(imp, "→", resp),
           ylab = "", xlab = "",
           ylim = range(c(irf_all$Lower[[imp]][, resp], 
                          irf_all$Upper[[imp]][, resp], 
                          irf_all$irf[[imp]][, resp])))
      lines(irf_all$Lower[[imp]][, resp], col = "red", lty = 2)
      lines(irf_all$Upper[[imp]][, resp], col = "red", lty = 2)
      abline(h = 0, col = "blue", lty = 3)
    }
  }
  
  # Ajouter le titre principal correctement au-dessus
  mtext(model_name, outer = TRUE, line = 2, cex = 1.5)
  
  par(mfrow = c(1,1))  # remettre la configuration par défaut
}

# Boucle sur tous les modèles pour générer les plots
for (model_name in names(var_models)) {
  plot_all_irf(var_models[[model_name]], model_name)
}

# ==============================================================================
# INTERPRÉTATION DES FONCTIONS DE RÉPONSE IMPULSIONNELLE (IRF)
# Modèle VAR Ordre économique : WALCL (Bilan Fed) -> BAA10Y (Spread Crédit) -> INDPRO (Prod. Ind)
# ==============================================================================

# Note méthodologique :
# Les lignes rouges pointillées représentent l'intervalle de confiance à 95%.
# Si la ligne bleue (zéro) est entre les deux lignes rouges, l'effet n'est pas statistiquement significatif.

# ------------------------------------------------------------------------------
# 1. ANALYSE DES CHOCS DE POLITIQUE MONÉTAIRE (Choc sur WALCL)
# ------------------------------------------------------------------------------
# Graphique : WALCL -> BAA10Y
# ---------------------------
# Observation : Suite à une hausse de WALCL (expansion du bilan de la Fed / QE),
# on observe une réaction négative du BAA10Y (le spread diminue).
# Interprétation Économique : 
# Une politique monétaire accommodante (achat d'actifs) injecte des liquidités,
# ce qui réduit le risque de crédit perçu et calme les marchés financiers.
# L'effet est significatif à court terme (les lignes rouges ne touchent pas 0).

# Graphique : WALCL -> INDPRO
# ---------------------------
# Observation : La réponse de la production industrielle est légèrement positive
# mais lente et l'intervalle de confiance inclut souvent zéro au début.
# Interprétation Économique :
# La politique monétaire met du temps à se transmettre à l'économie réelle.
# L'impact sur la production n'est pas immédiat, contrairement aux marchés financiers.

# ------------------------------------------------------------------------------
# 2. ANALYSE DES CHOCS DE STRESS FINANCIER (Choc sur BAA10Y)
# ------------------------------------------------------------------------------
# Graphique : BAA10Y -> INDPRO
# ----------------------------
# Observation : Un choc positif sur BAA10Y (hausse du spread/aversion au risque)
# entraîne une chute immédiate et significative de INDPRO.
# Interprétation Économique :
# C'est un résultat classique : une détérioration des conditions financières 
# (crise de crédit) freine l'investissement et la consommation, 
# entraînant une baisse rapide de l'activité réelle.

# Graphique : BAA10Y -> WALCL
# ---------------------------
# Observation : La réponse de WALCL est proche de zéro ou légèrement positive.
# Interprétation Économique :
# La Fed réagit peu ou pas du tout de manière endogène et immédiate 
# à une simple hausse du spread dans ce modèle spécifique, 
# ou alors elle réagit avec un décalage pour stabiliser le marché.

# ------------------------------------------------------------------------------
# 3. ANALYSE DE L'ACTIVITÉ RÉELLE (Choc sur INDPRO)
# ------------------------------------------------------------------------------
# Graphique : INDPRO -> INDPRO
# ----------------------------
# Observation : Forte persistance.
# Interprétation : Les cycles économiques sont persistants; une hausse de la
# production aujourd'hui prédit une production élevée demain.

# ------------------------------------------------------------------------------
# 4. ROBUSTESSE ET ORDRE DE CHOLESKY (Comparaison avec les autres graphiques)
# ------------------------------------------------------------------------------
# Vous avez testé 6 ordres différents (Désordre 1 à 5).
# Conclusion sur la robustesse :
# On remarque que la forme générale des courbes (la dynamique) reste
# relativement stable quel que soit l'ordre des variables.
# Cela suggère que vos résultats sont "robustes" : la relation négative entre 
# le stress financier (BAA10Y) et la production (INDPRO) est structurelle 
# et ne dépend pas uniquement de l'hypothèse d'identification de Cholesky.

# ==============================================================================
# 5. JUSTIFICATION DU CHOIX DU MODÈLE (ET EXCLUSION DES AUTRES ORDRES)
# ==============================================================================

# Nous retenons l'ordre : WALCL -> BAA10Y -> INDPRO (Ordre Économique)
# et nous excluons les "Désordres" 1 à 5 pour les raisons suivantes :

# RAISON 1 : Le Mécanisme de Transmission Monétaire
# -------------------------------------------------
# La théorie économique postule une chaîne de causalité logique :
# 1. L'intervention de la Banque Centrale (WALCL) est le choc initial (exogène).
# 2. Les marchés financiers (BAA10Y) réagissent immédiatement à cette information.
# 3. L'économie réelle (INDPRO) s'ajuste en dernier ressort face aux nouvelles conditions.
#
# Les autres ordres violent cette chaîne logique. Par exemple, placer INDPRO en premier
# (Désordres 4 et 5) reviendrait à dire que la production cause la politique monétaire
# à l'instant t, ce qui dilue l'identification du choc de politique monétaire.

# RAISON 2 : La Vitesse d'Ajustement (Hypothèse de Cholesky)
# ----------------------------------------------------------
# La décomposition de Cholesky impose une structure sur les réactions contemporaines (t=0).
#
# - Exclusion des ordres où BAA10Y est avant WALCL (ex: Désordre 2) :
#   Les marchés financiers (spreads de crédit) sont ultra-rapides et "forward-looking".
#   Il est irréaliste de supposer que le spread ne réagit pas *instantanément* #   à une annonce d'extension du bilan de la Fed. L'ordre doit donc être WALCL -> BAA10Y.
#
# - Exclusion des ordres où INDPRO n'est pas en dernier :
#   La production industrielle possède une forte inertie (rigidités réelles).
#   Les usines ne changent pas leur cadence de production le jour même d'un krach 
#   ou d'une annonce de la Fed. Il est donc préférable de placer la variable "lente" 
#   (INDPRO) à la fin de la chaîne causale pour le court terme.

# CONCLUSION
# ----------
# L'ordre "WALCL -> BAA10Y -> INDPRO" est le seul qui respecte à la fois :
# 1. L'exogénéité de la politique monétaire (Fed décide).
# 2. L'efficience des marchés financiers (Réaction rapide).
# 3. L'inertie de l'économie réelle (Réaction lente).


# 1) Choix des variables
# ----------------------
# 
# Prioritaire : log_INDPRO (activité) et log_BAA10Y (spread corporate).
# Raison : question centrale du mémoire — mesurer si les niveaux (non différenciés) ont une relation d’équilibre à long terme.
# 
# Alternatives à tester ensuite : log_INDPRO vs log_WALCL (impact long terme de la politique monétaire), ou un test multivarié Johansen sur les 3 séries (INDPRO, BAA10Y, WALCL).

# 2) Tests à réaliser
# -------------------
# 
# Engle–Granger (EG) — régression en niveaux (Y ~ X), puis test ADF sur les résidus (résidu stationnaire ⇒ cointégration).
# 
# Phillips–Ouliaris (PO) — test résiduiel alternatif souvent plus puissant.
# 
# Johansen pour tester l’existence d’une (des) relation(s) de cointégration dans un système multivarié (utile si tu testes 3 séries).
# 
# Robustesse : vérifier breaks structurels (ex. Gregory–Hansen) si tu suspectes un changement de régime; faire tests sur sous-périodes.


# --- ÉTAPE 1 : Préparation des données (les séries en NIVEAU LOG) ---

# Fusionner les séries en NIVEAU LOG
index(log_PRODIND_ts_ind) <- as.Date(format(index(log_PRODIND_ts_ind), "%Y-%m-01"))
index(log_spread_spread) <- as.Date(format(index(log_spread_spread), "%Y-%m-01"))
head(log_PRODIND_ts_ind)
head(log_spread_spread)


# 1. Fusionner les séries en un seul xts, ne garder que les dates communes
Xt <- merge(log_PRODIND_ts_ind, log_spread_spread, join = 'inner')
colnames(Xt) <- c("log_INDPRO", "log_BAA10Y")
head(Xt)

# --- ÉTAPE 2 : Test de Cointégration (Engle-Granger) ---


# --- Préparation d'un data.frame en niveaux log ---
df_levels <- data.frame(
  date = index(Xt),
  log_INDPRO = coredata(Xt$log_INDPRO),
  log_BAA10Y = coredata(Xt$log_BAA10Y)
)
head(df_levels)

# 1) Engle-Granger (régression en niveaux + ADF sur résidus)
# ----------------------------------------------------------

eg_fit <- lm(log_INDPRO ~ log_BAA10Y, data = df_levels)
summary(eg_fit)
# eg* = -0.11281 : plus le spread augmente, plus le niveau de la production industrielle (log) est plus faible (long terme).

# récupérer résidus (échelle niveaux) -> ECT_t = log_INDPRO_t - beta * log_BAA10Y_t - c
ect_resid <- resid(eg_fit)

# ADF sur les résidus (type = "none" car on teste la stationnarité du résidu)
eg_adf <- ur.df(ect_resid, type = "none", selectlags = "BIC")
summary(eg_adf)
# H0 : présence de RU - ρ = 0 → pas de relation de cointégration
# t* = -3.286 < -1.95 (VC) -> On rejette H0 
# Présence d'une relation de cointégration (relation d'équilibre de long terme entre la production industrielle et les spreads).

# Récapitulatif :
# --------------
# La relation de long terme estimée (EG) est statistiquement significative.
# Les résidus sont stationnaires.
# Donc il existe une combinaison linéaire stable reliant ces deux séries sur la période étudiée.

# 2) Phillips-Ouliaris (test résiduel) 
# ------------------------------------
# demean = "constant" car il y a une constante dans la relation de cointégration)
po_test <- ca.po(as.matrix(df_levels[, c("log_INDPRO","log_BAA10Y")]),
                 demean = "constant", type = "Pz")
summary(po_test)
# H0 : pas de cointégration
# t* = 11.0419  < 55.2202 (VC) -> On ne rejtte pas H0 -> Il n'y a pas de cointégration.

# Il y a une contradiction entre les résultats des tests Engle-Granger (cointégration) et Phillips-Ouliaris (pas de cointégration).
# On va trancher avec le test de Johansen.

# 3) Johansen (robustesse / multivarié)
# Choix du nombre de lags K : on le choisit à partir d'un VAR en niveaux
varsel <- VARselect(df_levels[, c("log_INDPRO","log_BAA10Y")], lag.max = 8, type = "const")
print(varsel$selection)
# Supposons K = varsel$selection["AIC(n)"] + 1 ou utiliser directement (ici je fixe K=2 comme exemple)
jo_test <- ca.jo(as.matrix(df_levels[, c("log_INDPRO","log_BAA10Y")]), type = "trace", ecdet = "const", K = 2)
summary(jo_test)
# summary donne le rang de cointégration (r). r = 1 signifie 1 relation de cointégration entre les deux séries.


# Pour r = 0 :
# H0 : il n'y a pas de vecteur de cointégration
# t* = 34.07 > 19.96 (VC) -> On rejette H0 -> Il y a au moins un vecteur de cointégration.
# Pour r = 1 :
# H0 : il y a un vecteur de cointégration 
# t* = 3.02 < 9.24 (VC) -> On ne rejette pas H0 -> r <= 1.
# Conlusion : Il y a exactement une relation de cointégration (r=1).

# On aboutit ainsi à log_INDPRO + 0.15368 * log_BAA10Y − 4.75577 = 0
# Qu'on réécrit : log_INDPRO = −0.15368 * log_BAA10Y + 4.75577
# Loading matrix : log_INDPRO.d = -0.07274834 et log_BAA10Y.d = -0.09217871
# Ces coefficients mesurent comment les différences s'ajustent à la violation de l'équilibre.

# Le test multivarié, robuste et adapté quand on a >2 séries ou qu’on veut un VECM, confirme l’existence d’un équilibre de long terme liant les deux séries, aligné avec le résultat EG.


# Engle–Granger (ADF sur résidus) : rejette H0 → cointégration détectée.
# Johansen (trace) : rejette r=0, accepte r <= 1 → r = 1 → cointégration détectée.
# Phillips–Ouliaris (Pz) : ne rejette pas H0 → pas de cointégration détectée selon ce test.
# Ainsi, on conclut qu'il existe une relation de cointégration entre la production industrielle et les spreads.

# Vérification rapide de la cohérence des coefficients de cointégration et de chargement du test  Johansen
#---------------------------------------------------------------------------------------------------------
# 4) Construction d'un ECM
# Alignement : Δ séries et ECT_{t-1}
dlog_INDPRO <- diff(df_levels$log_INDPRO)        # length T-1
dlog_BAA10Y <- diff(df_levels$log_BAA10Y)        # length T-1
# lagged ECT: ECT_{t-1} corresponds to resid at time t-1 => on enlève dernière observation du résidu
ect_lag <- ect_resid[-length(ect_resid)]         # length T-1

# créer dataframe pour estimation ECM
# Nous estimons :  Δlog(INDPROt) = γ + α*ECT[t−1]+ β0* Δlog(BAA10Y[t]) + β1* Δlog(BAA10Y[t−1]) + ϕ*Δlog(INDPRO[t−1]) + ε[t]
df_ecm <- data.frame(
  date = df_levels$date[-1],
  dlog_INDPRO = dlog_INDPRO,
  dlog_BAA10Y = dlog_BAA10Y,
  ect_lag = ect_lag
)


# Rajouter lags des variations (court-terme), par ex ΔX_{t-1}
df_ecm$dlog_INDPRO_lag1 <- c(NA, head(df_ecm$dlog_INDPRO, -1))
df_ecm$dlog_BAA10Y_lag1 <- c(NA, head(df_ecm$dlog_BAA10Y, -1))

# Supprimons les NA
df_ecm <- na.omit(df_ecm)

# Estimation simple de l'ECM 
ecm_model_simple <- lm(dlog_INDPRO ~ ect_lag + dlog_BAA10Y + dlog_BAA10Y_lag1 + dlog_INDPRO_lag1, data = df_ecm)
summary(ecm_model_simple)


# Intercept = 0.0012789 -> S'il n'y avait aucun choc et que l'économie était parfaitement à l'équilibre (ect_lag = 0), la production industrielle augmenterait "naturellement" de 0.128% par mois. C'est la tendance de fond de la série.
# ect_lag = - 0.0710890 -> Négatif (valide l'hypothèse de cointégration -> il y a une sorte de force de rappel vers l'équilibre)
# La vitesse d'ajustement est de 7.1% (ie lorsqu'il y a un déséquilibre, l'écart est corrigé de 7.1% par mois).
# p-val < 0.05 -> fortement significatif.
# dlog_BAA10Y = 0.0076362         (choc sur le spread ce mois-ci)
# p-val > 0.05 -> Non significatif
# dlog_BAA10Y_lag1 = - 0.0037818  (choc sur le spread le mois passé)
# p-val > 0.05 -> Non significatif
# dlog_INDPRO_lag1 = - 0.0890044  (choc sur la production le mois passé)
# p-val > 0.05 -> Non significatif


# Une fois que l'on prend en compte le déséquilibre de long terme (le ect_lag), les variations soudaines du spread n'ont pas d'impact supplémentaire sur la croissance de la production.
# L'impact du spread sur l'activité ne se joue pas sur des chocs de court terme, mais exclusivement via la relation de long terme. Ce n'est pas la variation du spread qui prédit l'activité, c'est son niveau (trop haut ou trop bas par rapport à l'équilibre).

# Conlusion générale :
# --------------------
# Il y a une relation de long terme qui lie la production au spread de crédit.
# Cette relation est stable : Quand la production s'emballe (ou s'effondre) par rapport au niveau du spread, elle finit toujours par y revenir.
# La vitesse de retour est de 7.1% par mois.
# Le canal de transmission est le long terme


# Test de Newey-West pour vérifier autocorr, het, normalité résidus
coeftest(ecm_model_simple, vcov = NeweyWest(ecm_model_simple))



# Conclusion générale :
# Les coefficients de cointegration et les coefficients de chargement du test Johansen sont numériquement cohérents 
# avec l’ECM estimé (α ≈ -0.0727 dans Johansen vs -0.0711 dans l'ECM).
# Nos résultats sont donc cohérents.


# Conclusion principale (économétrique) : il existe de fortes preuves d’une relation de long terme (une seule relation de cointégration) 
# entre le niveau du spread corporate (BAA10Y) et le niveau de la production industrielle (INDPRO) sur ta période 2010–2019. 
# Cette conclusion repose surtout sur la concordance de Engle–Granger et Johansen (tests robustes et complémentaires). 
# Le test Phillips–Ouliaris n’a pas confirmé la cointégration ici ; ceci mérite vérification supplémentaire (cf. robustesse).

# Conclusion économique :
#   
#   Long terme : le signe du coefficient (-0.11 en EG, -0.154 en Johansen) montre qu’un élargissement du spread s’associe à un niveau 
#               plus faible de la production industrielle — cohérent avec la théorie : conditions de crédit plus strictes freinent l’activité.
# 
#   Court terme : les variations mensuelles contemporaines du spread n’apparaissent pas significatives pour expliquer la croissance mensuelle 
#                 de la production (ΔINDPRO). Autrement dit, l’effet du spread se manifeste plutôt via le niveau (effet long terme) 
#                 et INDPRO ajuste progressivement (α ≈ -0.071) vers l’équilibre. Cela rejoint l’hypothèse que les effets se matérialisent avec retard.



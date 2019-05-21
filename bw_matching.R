library(sensitivity2x2xk)
library(Epi)
library(data.table)
library(ggplot2)
library(optmatch)
library(foreign)
library(AER)
library(ivmodel)
library(sensitivitymv)

## Define some helpful functions.
# Function for computing 
# rank based Mahalanobis distance.  Prevents an outlier from
# inflating the variance for a variable, thereby decreasing its importance.
# Also, the variances are not permitted to decrease as ties 
# become more common, so that, for example, it is not more important
# to match on a rare binary variable than on a common binary variable
# z is a vector, length(z)=n, with z=1 for treated, z=0 for control
# X is a matrix with n rows containing variables in the distance
smahal <- function(z,X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  m<-sum(z)
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,m,n-m)
  Xc<-X[z==0,]
  Xt<-X[z==1,]
  rownames(out)<-rownames(X)[z==1]
  colnames(out)<-rownames(X)[z==0]
  library(MASS)
  icov<-ginv(cv)
  for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
  out
}

# Function for adding a propensity score caliper to a distance matrix dmat
# calipersd is the caliper in terms of standard deviation of the logit propensity scoe
addcaliper=function(dmat,z,logitp,calipersd=.2,penalty=1000){
  sd.logitp=sd(logitp)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

# Add a column with match ids to a provided data.table.
match_controls <- function(n_controls, fm=F) {
  nocontrols.per.match=n_controls
  if(!fm) matchvec <- as.numeric(pair(distmat2, controls=nocontrols.per.match, data=d))
  if(fm) matchvec <- as.numeric(fullmatch(distmat2, data=d))
  d[, (paste0('matchvec_', n_controls)) := matchvec]
  return(NULL)
}

# Calculate standardized difference to examing balancing.
calc_stand_diff <- function(v, treatment, data) {
  t.mean <- mean(data[get(treatment)==T, get(v)])
  t.sd <- sd(data[get(treatment)==T, get(v)])
  c.mean <- mean(data[get(treatment)==F, get(v)])
  c.sd <- sd(data[get(treatment)==F, get(v)])
  stand_diff <- (t.mean - c.mean) / sqrt((t.sd^2 + c.sd^2) / 2)
  return(data.table(variable=v, diff=stand_diff))
}

# Plot of when Gamma test passes p=0.05
plot_gamma <- function(d, start, stop, by) {
  results <- data.table(gamma=seq(start,stop,by))
  results[, p := 0]
  for(g in seq(start,stop,by)) results[gamma==g, p := sensitivity2x2xk::mh(d,Gamma=g)$pval]
  results <- ggplot(data=results, aes(x=gamma,y=p)) +
    geom_line(size=2) +
    geom_hline(yintercept = 0.05, linetype='dashed') + 
    labs(x='Gamma',y='P-value',title='P-values for range of Gamma values') + 
    theme_bw()
  return(results)
}
# Sensitivity interval for IV
sens.interval.iv.deviates=function(beta,encouraged.treatment,unencouraged.treatment,encouraged.response,unencouraged.response,Gamma){
  adjusted.response.encouraged=encouraged.response-beta*encouraged.treatment;
  adjusted.response.unencouraged=unencouraged.response-beta*unencouraged.treatment;
  diff=adjusted.response.encouraged-adjusted.response.unencouraged;
  rk=rank(abs(diff));
  s1=1*(diff>0);
  s2=1*(diff<0);
  W=sum(s1*rk);
  Eplus=sum((s1+s2)*rk*Gamma)/(1+Gamma);
  Eminus=sum((s1+s2)*rk)/(1+Gamma);
  V=sum((s1+s2)*rk*rk*Gamma)/((1+Gamma)^2);
  Dplus=(W-Eplus)/sqrt(V);
  Dminus=(W-Eminus)/sqrt(V);
  list(Dplus=Dplus,Dminus=Dminus);
}

# Plot of when Gamma test passes p=0.05
plot_gamma <- function(d, start, stop, by, alternative) {
  results <- data.table(gamma=seq(start,stop,by))
  results[, p := 0]
  for(g in seq(start,stop,by)) results[gamma==g, p := senmv(-d,gamma=g,tau=alternative)$pval]
  results <- ggplot(data=results, aes(x=gamma,y=p)) +
    geom_line(size=2) +
    geom_hline(yintercept = 0.05, linetype='dashed') + 
    labs(x='Gamma',y='P-value',title='P-values for range of Gamma values') + 
    theme_bw()
  return(results)
}

# Sensitivity test for IV Gamma
iv_gamma <- function(start, stop, by) {
  results <- data.table(gamma=seq(start,stop,by))
  results[, upper := 0]
  results[, lower := 0]
  for(Gamma in seq(start,stop,by)) {
    betagrid=seq(-5,5,.05);
    Dplus.betagrid=rep(0,length(betagrid));
    Dminus.betagrid=rep(0,length(betagrid));
    for(i in 1:length(betagrid)){
      beta=betagrid[i];
      deviates=sens.interval.iv.deviates(beta,treated.df.nearer$twoyr,matched.control.df.nearer$twoyr,treated.df.nearer$educ86,matched.control.df.nearer$educ86,Gamma);
      Dplus.betagrid[i]=deviates$Dplus;
      Dminus.betagrid[i]=deviates$Dminus;
    }
    results[gamma==Gamma, lower := min(betagrid[Dplus.betagrid<=1.96])]
    results[gamma==Gamma, upper := max(betagrid[Dminus.betagrid>=-1.96])]
  }
  results <- ggplot(data=results, aes(x=gamma,ymax=upper,ymin=lower)) +
    geom_ribbon(color='blue',alpha=0.3,fill='blue') +
    geom_hline(yintercept = 0, linetype='dashed', size=2) + 
    labs(x='Gamma',y='Confidence Interval',title='Confidence intervals for range of Gamma values') + 
    theme_bw()
  return(results)
}

## Set up HRS data
d <- fread('C:/Users/ngraetz/Downloads/minwagedataframe.csv')
d[, state := as.numeric(ifelse(state=='PA', '1', '0'))]
for(var in c('chain','co_owned')) {
  for(v in unique(d[, get(var)])) {
    d[, (paste0(var,'_',v)) := ifelse(get(var)==v, 1, 0)]
  }
}
confounders <- c('chain_wendys', 'chain_bk', 'chain_kfc', 'co_owned_yes', 'hours_open')
treatment <- 'state'

## Question 1a)
## Match 2 NJ to 1 PA using chain, co_owned, hours_open
## Use 2:1 propensity score matching with Mahalanobis distance and caliper.
ps_formula <- as.formula(paste0(treatment, ' ~ ', paste(confounders, collapse = ' + ')))
ps_mod <- glm(ps_formula,
              family=binomial,
              x=TRUE,
              y=TRUE,
              data=d)
d[, prop_score := predict(ps_mod)]
max.ps.control <- max(d[get(treatment)==F, prop_score])
min.ps.treated <- min(d[get(treatment)==T, prop_score])
d[get(treatment)==F & prop_score < min.ps.treated, drop := 1]
d[get(treatment)==T & prop_score > max.ps.control, drop := 1]
message(paste0('Dropping ', length(d[drop==1, drop]), ' obs with outlier propensity scores.'))
d <- d[is.na(drop), ]

## Set up for matching.
# Matrix of confounders.
Xmat <- d[, confounders, with=F]
# Matrix of covariates to include in the Mahalanobis distance (just use all).
Xmatmahal <- d[, confounders, with=F]
# Rank based Mahalanobis distance
distmat <- smahal(d[, get(treatment)], Xmatmahal)
# Add caliper
distmat2 <- addcaliper(distmat, d[, get(treatment)], d[, prop_score], calipersd=.5)
# Name the rows and columns of distance matrix by the subject numbers in treated
d[, id := 1:dim(d)[1]]
rownames(distmat2)=d[get(treatment)==T, id]
colnames(distmat2)=d[get(treatment)==F, id]

## Calculate matches (1, 2, 3, full)
match_controls(2)

## Check for acceptable balance
make_diffs <- function(n_controls) {
  diffs.after <- rbindlist(lapply(confounders, calc_stand_diff,
                                  treatment = treatment,
                                  data = d[!is.na(get(paste0('matchvec_', n_controls))), ]))
  diffs.after[, matching := paste0(n_controls, ' matched controls')]
  return(diffs.after)
}
all.diffs <- make_diffs('2')
diffs.before <- rbindlist(lapply(confounders, calc_stand_diff, treatment = treatment, data = d))
diffs.before[, matching := 'Before matching']
all.diffs <- rbind(all.diffs, diffs.before)

## Make Love Plot
ggplot() + 
  geom_jitter(data=all.diffs,
              aes(x=diff,
                  y=variable,
                  fill=matching),
              width = 0.1,
              height = 0.1,
              size=5,
              shape=21) + 
  geom_vline(xintercept = .2, linetype='dashed') + 
  geom_vline(xintercept = -.2, linetype='dashed') + 
  geom_vline(xintercept = 0) + 
  labs(y='',x='Standardized difference',
       title='Love plot of standardized differences for all confounders\nMatching via rank-based Mahalanobis distance with a propensity score caliper', fill='') + 
  coord_flip() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
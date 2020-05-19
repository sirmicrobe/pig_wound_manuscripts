#Calculate selection coefficients for fitness assays

library(ggplot2) #for plotting
library(reshape2) #for melt
library("dplyr")
library("FSA") #for stats dunn test

data <- read.csv("PAO1_fit_raw.csv",header = T)
tail(data)
data <- data[,1:7]
data$sample_time <- paste(data$Competition,data$Replicate,data$time, sep="_")
data$sample <- paste(data$Competition,data$Replicate, sep="_")

View(data)

data$cfu_blue <- (data$Blue*10)/data$Dilution.blue
data$cfu_white <- (data$White*10)/data$Dilution.white

data_0 <- subset(data,time=="0")
data_24 <- subset(data,time=="24")
data_48 <- subset(data,time=="48")

merge1<- merge(data_0,data_24,by="sample",suffixes=c(".0",".24"),all.x=T)
merge_all <- merge(merge1,data_48,by="sample",suffixes=c("",".48"),all.x=T)
head(merge_all,30)

colnames(merge_all)
cfu_df <- merge_all[,c("sample","Competition","Replicate","time.0","cfu_blue.0","cfu_white.0","time.24","cfu_blue.24","cfu_white.24","time","cfu_blue","cfu_white")]
head(cfu_df)

cfu_df$r.white.24 <- (log(cfu_df$cfu_blue.24/cfu_df$cfu_blue.0) - log(cfu_df$cfu_white.24/cfu_df$cfu_white.0))/1
cfu_df$r.white.48 <- (log(cfu_df$cfu_blue/cfu_df$cfu_blue.0) - log(cfu_df$cfu_white/cfu_df$cfu_white.0))/2
cfu_df$r.blue.24 <- (log(cfu_df$cfu_white.24/cfu_df$cfu_white.0) - log(cfu_df$cfu_blue.24/cfu_df$cfu_blue.0))/1
cfu_df$r.blue.48 <- (log(cfu_df$cfu_white/cfu_df$cfu_white.0) - log(cfu_df$cfu_blue/cfu_df$cfu_blue.0))/2

View(cfu_df)
write.csv(cfu_df,"pao1_fit_calc.csv")


fit_df <- read.csv("pao1_reduced.csv")
tail(fit_df)
fit_df$r <- as.numeric(as.character(fit_df$r))

fit_tib_b <- subset(fit_df, B_P=="biofilm") %>%
  group_by(Competition) %>%
  filter(!is.na(r)) %>%
  mutate(n = n(),
         mean = mean(r),
         median = median(r),
         sd = sd(r)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

level_order <- c('15dipA','5retS vs. wspF','5retS', '41retS', 'wspF','PAO1') #this vector might be useful for other plots/analyses
p1 <- ggplot(fit_tib_b, aes(x = factor(Focal, level = level_order), r,fill=Focal))
p1 <- p1 + geom_boxplot() +geom_point() +theme(text = element_text(size=16, face="bold"),legend.position = "none")+ geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.01,position=position_dodge(0.9)) + scale_fill_manual(values=c("red","blue","purple","darkgoldenrod","cyan3","brown4")) + ggtitle("Biofilm 48 hours") + coord_flip() + ylab("Selection rate per day (r)") + xlab("Isolates")
p1

fit_tib_p <- subset(fit_df, B_P=="planktonic") %>%
  group_by(Competition) %>%
  filter(!is.na(r)) %>%
  mutate(n = n(),
         mean = mean(r),
         median = median(r),
         sd = sd(r)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

p2 <- ggplot(fit_tib_p, aes(x = factor(Focal, level = level_order), r,fill=Focal))
p2 <- p2 + geom_boxplot() +geom_point() +theme(text = element_text(size=16, face="bold"),legend.position = "none")+ geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.01,position=position_dodge(0.9)) + scale_fill_manual(values=c("red","blue","purple","darkgoldenrod","cyan3","brown4")) + ggtitle("Planktonic 24 hours") + coord_flip() + ylab("Selection rate per day (r)") + xlab("Isolates")
p2

fit_tib_all <- fit_df %>%
  group_by(Competition,B_P) %>%
  filter(!is.na(r)) %>%
  mutate(n = n(),
         mean = mean(r),
         median = median(r),
         sd = sd(r)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

p3 <- ggplot(fit_tib_all, aes(x = factor(Focal, level = level_order), r,fill=B_P))
p3 + geom_boxplot() +geom_point() +theme(text = element_text(size=16, face="bold"),legend.position = "none")+ geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.01,position=position_dodge(0.9)) + scale_fill_manual(values=c("red","blue","purple","darkgoldenrod","cyan3","brown4")) + ylab("Selection rate per day (r)") + xlab("Isolates")+facet_wrap(.~Focal,scale="free") #+ ggtitle("Planktonic 24 hours")

p3_1 <- ggplot(fit_tib_all, aes(x = factor(Focal, level = level_order), r,fill=Focal))
p3_1 + geom_boxplot() +geom_point() +theme(text = element_text(size=16, face="bold"),legend.position = "none")+ geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.01,position=position_dodge(0.9)) + scale_fill_manual(values=c("red","blue","purple","darkgoldenrod","cyan3","brown4"))  + coord_flip() + ylab("Selection rate per day (r)") + xlab("Isolates")+facet_wrap(~B_P) #+ ggtitle("Planktonic 24 hours")

#each facet compared to pao1
fit_df_all <- as.data.frame(fit_tib_all)
every_facet_data = subset(fit_df_all, Focal == "PAO1")
individual_facet_data = subset(fit_df_all, Focal != "PAO1")
individual_facet_data$facet = individual_facet_data$Focal

every_facet_data = merge(every_facet_data,
                         data.frame(Focal = "PAO1", facet = unique(individual_facet_data$facet)))

plot_data = rbind(every_facet_data, individual_facet_data)
ggplot(plot_data, aes(x=Focal, y=r, fill=B_P)) + geom_boxplot() +geom_point() +theme(text = element_text(size=16, face="bold"),legend.position = "right")+ geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.01,position=position_dodge(0.9)) + scale_fill_manual(values=c("red","blue","purple","darkgoldenrod","cyan3","brown4")) + ylab("Selection rate per day (r)") + xlab("Isolates")+facet_wrap(~facet,scale="free")+coord_flip()


kruskal.test(r~as.factor(Focal),data=subset(plot_data,B_P=="planktonic")) #KKruskal-Wallis chi-squared = 37.025, df = 5, p-value = 5.92e-07
dunnTest(r ~ as.factor(Focal),
         data=subset(plot_data,B_P=="planktonic"),
         method="bh") #multiple comparison test 
Comparison           Z      P.unadj        P.adj
#1          15dipA - 41retS -0.13047191 8.961931e-01 0.9602068683
#2           15dipA - 5retS  0.69585018 4.865227e-01 0.6081533217
#3           41retS - 5retS  0.82632209 4.086214e-01 0.6129321056
#4  15dipA - 5retS vs. wspF -0.17396255 8.618949e-01 0.9944941122
#5  41retS - 5retS vs. wspF -0.04349064 9.653104e-01 0.9653104284
#6   5retS - 5retS vs. wspF -0.86981273 3.844028e-01 0.6406712556
7            15dipA - PAO1  3.76179341 1.686994e-04 0.0008434968
8            41retS - PAO1  3.93023192 8.486396e-05 0.0006364797
9             5retS - PAO1  2.86345469 4.190487e-03 0.0157143281
10   5retS vs. wspF - PAO1  3.98637809 6.708957e-05 0.0010063436
#11           15dipA - wspF  2.32674905 1.997863e-02 0.0428113401
#12           41retS - wspF  2.45722096 1.400165e-02 0.0350041341
#13            5retS - wspF  1.63089887 1.029117e-01 0.1929593632
14   5retS vs. wspF - wspF  2.50071160 1.239441e-02 0.0371832197
15             PAO1 - wspF -0.75797330 4.484670e-01 0.6115458646

aov.r <- aov(r~as.factor(Focal),data=subset(plot_data,B_P=="planktonic"))
summary(aov.r) #as.factor(Focal)  5 237.64   47.53   62.67 <2e-16 ***
TukeyHSD(aov.r)
#                        diff        lwr        upr     p adj
41retS-15dipA          0.8189538 -0.8217099  2.4596174 0.6741352
5retS-15dipA          -0.4380810 -2.0787447  1.2025826 0.9668374
5retS vs. wspF-15dipA  0.2558347 -1.3848290  1.8964983 0.9971172
PAO1-15dipA           -4.2380945 -5.5089471 -2.9672419 0.0000000
wspF-15dipA           -4.3201235 -5.9607872 -2.6794598 0.0000000
5retS-41retS          -1.2570348 -2.8976984  0.3836289 0.2229978
5retS vs. wspF-41retS -0.5631191 -2.2037827  1.0775446 0.9078894
PAO1-41retS           -5.0570483 -6.3279009 -3.7861956 0.0000000
wspF-41retS           -5.1390773 -6.7797409 -3.4984136 0.0000000
5retS vs. wspF-5retS   0.6939157 -0.9467480  2.3345794 0.8046234
PAO1-5retS            -3.8000135 -5.0708661 -2.5291609 0.0000000
wspF-5retS            -3.8820425 -5.5227061 -2.2413788 0.0000001
PAO1-5retS vs. wspF   -4.4939292 -5.7647818 -3.2230766 0.0000000
wspF-5retS vs. wspF   -4.5759582 -6.2166218 -2.9352945 0.0000000
wspF-PAO1             -0.0820290 -1.3528816  1.1888236 0.9999609


kruskal.test(r~as.factor(Focal),data=subset(plot_data,B_P=="biofilm")) #Kruskal-Wallis chi-squared = 34.834, df = 5, p-value = 1.624e-06
dunnTest(r ~ as.factor(Focal),
         data=subset(plot_data,B_P=="biofilm"),
         method="bh") #multiple comparison test 
Comparison          Z      P.unadj        P.adj
#1          15dipA - 41retS -0.4196958 6.747077e-01 0.7229010875
#2           15dipA - 5retS  0.7653277 4.440765e-01 0.5123959087
#3           41retS - 5retS  1.1850235 2.360082e-01 0.3540122442
#4  15dipA - 5retS vs. wspF  0.9194043 3.578841e-01 0.4880237946
#5  41retS - 5retS vs. wspF  1.3150974 1.884772e-01 0.3141286406
#6   5retS - 5retS vs. wspF  0.1978465 8.431652e-01 0.8431651594
7            15dipA - PAO1  3.9659701 7.309804e-05 0.0005482353
8            41retS - PAO1  4.4968479 6.896828e-06 0.0001034524
9             5retS - PAO1  2.9978986 2.718481e-03 0.0135924041
10   5retS vs. wspF - PAO1  2.4943830 1.261763e-02 0.0315440779
#11           15dipA - wspF  2.5181749 1.179647e-02 0.0353894184
#12           41retS - wspF  2.9378707 3.304748e-03 0.0123928052
#13            5retS - wspF  1.7528472 7.962823e-02 0.1706319280
14   5retS vs. wspF - wspF  1.4547537 1.457375e-01 0.2732577302
15             PAO1 - wspF -0.7807028 4.349773e-01 0.5437216659
aov.b <- aov(r~as.factor(Focal),data=subset(plot_data,B_P=="biofilm"))
summary(aov.b) #as.factor(Focal)  5  48.71   9.742   181.3 <2e-16 ***
TukeyHSD(aov.b)
#diff        lwr         upr     p adj
41retS-15dipA          0.16941864 -0.2704004  0.60923769 0.8546079
5retS-15dipA          -0.32188263 -0.7617017  0.11793643 0.2635242
5retS vs. wspF-15dipA -0.34835939 -0.8148579  0.11813917 0.2438550
PAO1-15dipA           -2.24339526 -2.5911028 -1.89568776 0.0000000
wspF-15dipA           -2.10733289 -2.5471520 -1.66751384 0.0000000
5retS-41retS          -0.49130126 -0.9311203 -0.05148221 0.0209416
5retS vs. wspF-41retS -0.51777802 -0.9842766 -0.05127947 0.0221255
PAO1-41retS           -2.41281390 -2.7605214 -2.06510640 0.0000000
wspF-41retS           -2.27675153 -2.7165706 -1.83693247 0.0000000
5retS vs. wspF-5retS  -0.02647676 -0.4929753  0.44002180 0.9999784
PAO1-5retS            -1.92151263 -2.2692201 -1.57380514 0.0000000
wspF-5retS            -1.78545027 -2.2252693 -1.34563121 0.0000000
PAO1-5retS vs. wspF   -1.89503587 -2.2759303 -1.51414140 0.0000000
wspF-5retS vs. wspF   -1.75897351 -2.2254721 -1.29247495 0.0000000
wspF-PAO1              0.13606236 -0.2116451  0.48376986 0.8463019

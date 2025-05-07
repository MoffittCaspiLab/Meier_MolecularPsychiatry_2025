################################################################################################################
##### Calculating gene expression correlations with Cannabis-associated methylation probes #####################
################################################################################################################

### the results of this script can be found in the section of the results titled "Association between DNA Methylation and Gene Expression Levels" in the manuscript ""DNA Methylation Profiles of Long-Term Cannabis Users in Midlife: A Comprehensive Evaluation of Published Cannabis-Associated Methylation Markers in a Representative Cohort", Meier et al. 

### for this analysis I will use the table of residualized methylation beta values at age 38, and
### correlate with the table of residualized gene expression values at age 38. 

########################################################################################################

set.seed(2024)

sessionInfo()

# load necessary libraries
#
library(Hmisc)
library(tidyverse)
library(psych)
library(tidyverse)

########################################################################################################

load("./ExpressionvsMeth/results/methyl.tech.cellmix.residuals") # residualized methylation matrix 

probelist <- as.matrix(read.table("./ExpressionvsMeth/Data/SignifCpG_2024Sept.csv", quote="\"", comment.char="")) # list of significant probes from main analysis

load("./ExpressionvsMeth/results/exprs.tech.cellmix.residuals") # residualized gene expression  matrix 

#### Now, need to run correlations across the gene expression data with the significant probes


gc()

# subset methylation data to significant probes
sigprobe <- subset(methyl.tech.cellmix.residuals, rownames(methyl.tech.cellmix.residuals) %in% probelist)


# make correlation matrix
gem.r <- cor(t(exprs.tech.cellmix.residuals), t(sigprobe), method = "spearman", use = "pairwise.complete.obs")

# make correlation p-value matrix
gem.p <- t(sapply(1:nrow(exprs.tech.cellmix.residuals), function(x) {
  sapply(1:nrow(sigprobe), function(y) {
    rcorr(exprs.tech.cellmix.residuals[x,],sigprobe[y,], type = "spearman")[[3]][1,2]
  })
}))
row.names(gem.p) <- row.names(exprs.tech.cellmix.residuals)
colnames(gem.p) <- row.names(sigprobe)

### save the un-filtered results 
write.table(as.data.frame(gem.r), file = paste0("./ExpressionvsMeth/results/correlations_Expmeth_nonfiltered_", Sys.Date(),".txt"), sep = "\t", row.names = T)

write.table(as.data.frame(gem.p), file = paste0("./ExpressionvsMeth/results/correlations_pvalue_Expmeth_nonfiltered_", Sys.Date(),".txt"), sep = "\t", row.names = T)

### split the file up by methylation probe in order to use in circos plots
### 
# columns need to be 1) expression probeset ID, 2) rho, 3) p-value
# N.B.; this will mean swapping the order of the circos plots so that the expLocations.df (expresion annotation) is used instead of Locations.df and vice versa


# Loop through columns and create separate data frames
for (i in seq_along(gem.r[1, ])) {
  df <- data.frame(ExpProbe = rownames(gem.r), rho = gem.r[, i])
  assign(paste0(colnames(gem.r)[i], "_rho"), df)
}


for (i in seq_along(gem.p[1, ])) {
  df <- data.frame(ExpProbe = rownames(gem.p), p = gem.p[, i])
  assign(paste0(colnames(gem.p)[i], "_p"), df)
}

# combine rho and p

# Get all data frame names in the environment
df_names <- ls(pattern = "^cg\\d+_") # Matches names like cgNNN_

# Extract the numeric substrings
df_groups <- unique(sub("cg(\\d+)_.*", "\\1", df_names))

# Loop through groups and cbind matching data frames
for (group in df_groups) {
  # Find matching data frame names
  matching_dfs <- grep(paste0(group, "_"), df_names, value = TRUE)
  
  if (length(matching_dfs) == 2) { # Ensure only two matches are combined
    # Combine the data frames
    combined_df <- cbind(get(matching_dfs[1]), get(matching_dfs[2]))
    
    # Assign the combined data frame to a new variable
    assign(paste0("cg", group, "_combined"), combined_df)
    

  }
  # write the files
  
  write.table(combined_df, file = paste0("./ExpressionvsMeth/Results/cg", group, "_combined.txt"), sep = "\t", row.names = F)
}


# change p-values greater than bonferroni-corrected p-value to NA (and associated correlation)

0.5/(49495*29) #3.483459e-07

gem.pf <- gem.p
gem.pf[gem.pf > 3.483459e-07] <- NA

gem.rf <- gem.r
gem.rf[is.na(gem.pf)] <- NA

# remove rows where every CpG is NA

gem.pf <- as.data.frame(gem.pf)
gem.pf <- gem.pf %>% 
  filter(!if_all(everything(), is.na))

gem.rf <- as.data.frame(gem.rf)
gem.rf <- gem.rf %>% 
  filter(!if_all(everything(), is.na))

gem.pf$ID <- row.names(gem.pf)
gem.rf$ID <- row.names(gem.rf)

#### 400 rows remain


###### Add Gene expression annotation to the correlation files

Primeview_updatedAnnotation_fixedChr <- read.delim("./ExpressionvsMeth/Data/Primeview_updatedAnnotation_fixedChr.txt")

gemr.df <- inner_join(Primeview_updatedAnnotation_fixedChr, gem.rf, by = "ID")
gemp.df <- inner_join(Primeview_updatedAnnotation_fixedChr, gem.pf, by = "ID")

write.table(gemr.df, file = "./ExpressionvsMeth/results/correlations_Expmeth.txt", sep = "\t", row.names = F)
write.table(gemp.df, file = "./ExpressionvsMeth/results/correlations_pvalue_Expmeth.txt", sep = "\t", row.names = F)

#########################################################################################

### repeat for smokers only

rm(list = c("gemr.df", "gemp.df", "gem.rf", "gem.pf", "gem.p", "gem.r", "methyl.tech.cellmix.residuals", "sigprobe"))
gc()

load("./ExpressionvsMeth/results/methyl.smoker.tech.cellmix.residuals")


### use the daily smoker at 38 variable to get current p38 smokers (N = 215)

smoking.df <- read.table("./ExpressionvsMeth/data/Smoking_Dunedin_Mar2106.txt", sep="\t", header=T, stringsAsFactors=F)

smoking.df$SampleID <- sapply(smoking.df$snum, function(x) { if( nchar(x) == 3 ) { paste0("DN0", x) } else { paste0("DN", x) } })

toUse <- smoking.df$SampleID[which(smoking.df$dailysmoker38 == 1)]
toUse <- toUse[which(toUse %in% colnames(methyl.tech.cellmix.residuals))] #167 individuals with meth/exp data

smkexp <- exprs.tech.cellmix.residuals[,toUse]
smkmeth <- methyl.tech.cellmix.residuals[,toUse]

sigprobe <- subset(smkmeth, rownames(smkmeth) %in% probelist)


# make correlation matrix
gem.r <- cor(t(smkexp), t(sigprobe), method = "spearman", use = "pairwise.complete.obs")

# make correlation p-value matrix
gem.p <- t(sapply(1:nrow(smkexp), function(x) {
  sapply(1:nrow(sigprobe), function(y) {
    rcorr(smkexp[x,],sigprobe[y,], type = "spearman")[[3]][1,2]
  })
}))
row.names(gem.p) <- row.names(smkexp)
colnames(gem.p) <- row.names(sigprobe)

### save the un-filtered results to compare against the 'all samples' results above

write.table(as.data.frame(gem.r), file = "./ExpressionvsMeth/results/correlations_smokers_Expmeth_nonfiltered.txt", sep = "\t", row.names = T)

write.table(as.data.frame(gem.p), file = "./ExpressionvsMeth/results/correlations_pvalue_smokers_Expmeth_nonfiltered.txt", sep = "\t", row.names = T)


# change p-values greater than bonferroni-corrected p-value to NA (and associated correlation)

0.5/(400*18) # 6.944444e-05

gem.pf <- gem.p
gem.pf[gem.pf > 6.944444e-05] <- NA

gem.rf <- gem.r
gem.rf[is.na(gem.pf)] <- NA

# remove rows where every CpG is NA

gem.pf <- as.data.frame(gem.pf)
gem.pf <- gem.pf %>% 
  filter(!if_all(everything(), is.na)) ## 223 rows remain

gem.rf <- as.data.frame(gem.rf)
gem.rf <- gem.rf %>% 
  filter(!if_all(everything(), is.na))

gem.pf$ID <- row.names(gem.pf)  ## 14 rows remain
gem.rf$ID <- row.names(gem.rf)

###### Add Gene expression annotation to the correlation files

Primeview_updatedAnnotation_fixedChr <- read.delim("./ExpressionvsMeth/Data/Primeview_updatedAnnotation_fixedChr.txt")

gemr.df <- inner_join(Primeview_updatedAnnotation_fixedChr, gem.rf, by = "ID")
gemp.df <- inner_join(Primeview_updatedAnnotation_fixedChr, gem.pf, by = "ID")

write.table(gemr.df, file = "./ExpressionvsMeth/results/correlations_smokers_Expmeth.txt", sep = "\t", row.names = F)
write.table(gemp.df, file = "./ExpressionvsMeth/results/correlations_pvalue_smokers_Expmeth.txt", sep = "\t", row.names = F)


########################################################################################

### repeat for non-smokers only

rm(list = c("gemr.df", "gemp.df", "gem.rf", "gem.pf", "gem.p", "gem.r", "methyl.tech.cellmix.residuals", "toUse", "sigprobe"))
gc()

load("./ExpressionvsMeth/results/methyl.nonsmoker.tech.cellmix.residuals")

# smoking.df <- read.table("./ExpressionvsMeth/data/Smoking_Dunedin_Mar2106.txt", sep="\t", header=T, stringsAsFactors=F)
# 
# smoking.df$SampleID <- sapply(smoking.df$snum, function(x) { if( nchar(x) == 3 ) { paste0("DN0", x) } else { paste0("DN", x) } })

### use packyears at 38 = 0 to get never smokers

toUse <- smoking.df$SampleID[which(smoking.df$packyears38 == 0)]
toUse <- toUse[which(toUse %in% colnames(methyl.tech.cellmix.residuals))]

nsmkexp <- exprs.tech.cellmix.residuals[,toUse]  ## 405 individuals
nsmkmeth <- methyl.tech.cellmix.residuals[,toUse]

sigprobe <- subset(nsmkmeth, rownames(nsmkmeth) %in% probelist)


# make correlation matrix
gem.r <- cor(t(nsmkexp), t(sigprobe), method = "spearman", use = "pairwise.complete.obs")

# make correlation p-value matrix
gem.p <- t(sapply(1:nrow(nsmkexp), function(x) {
  sapply(1:nrow(sigprobe), function(y) {
    rcorr(nsmkexp[x,],sigprobe[y,], type = "spearman")[[3]][1,2]
  })
}))
row.names(gem.p) <- row.names(nsmkexp)
colnames(gem.p) <- row.names(sigprobe)

### save the un-filtered results to compare against the 'all samples' results above

write.table(as.data.frame(gem.r), file = "./ExpressionvsMeth/results/correlations_nonsmokers_Expmeth_nonfiltered.txt", sep = "\t", row.names = T)

write.table(as.data.frame(gem.p), file = "./ExpressionvsMeth/results/correlations_pvalue_nonsmokers_Expmeth_nonfiltered.txt", sep = "\t", row.names = T)


# change p-values greater than bonferroni-corrected p-value to NA (and associated correlation)
# this correction is made based on the number of significant probesets in part A * number of probes with any correlation

0.5/(400*18) # 6.944444e-05

gem.pf <- gem.p
gem.pf[gem.pf > 6.944444e-05] <- NA

gem.rf <- gem.r
gem.rf[is.na(gem.pf)] <- NA

# remove rows where every CpG is NA

gem.pf <- as.data.frame(gem.pf)  ## 324 rows remain
gem.pf <- gem.pf %>% 
  filter(!if_all(everything(), is.na))

gem.rf <- as.data.frame(gem.rf)
gem.rf <- gem.rf %>% 
  filter(!if_all(everything(), is.na))

gem.pf$ID <- row.names(gem.pf) ## 53 rows remain
gem.rf$ID <- row.names(gem.rf)

###### Add Gene expression annotation to the correlation files

# Primeview_updatedAnnotation_fixedChr <- read.delim("./ExpressionvsMeth/Data/Primeview_updatedAnnotation_fixedChr.txt")

gemr.df <- inner_join(Primeview_updatedAnnotation_fixedChr, gem.rf, by = "ID")
gemp.df <- inner_join(Primeview_updatedAnnotation_fixedChr, gem.pf, by = "ID")

write.table(gemr.df, file = "./ExpressionvsMeth/results/correlations_nonsmokers_Expmeth.txt", sep = "\t", row.names = F)
write.table(gemp.df, file = "./ExpressionvsMeth/results/correlations_pvalue_nonsmokers_Expmeth.txt", sep = "\t", row.names = F)


##########################################################################################################

#### select just the significant gene expression probeset values

correlations_Expmeth <- read.delim("./ExpressionvsMeth/results/correlations_Expmeth.txt")

probeset <- as.character(correlations_Expmeth$ID)

exp_probes <- exprs.tech.cellmix.residuals[row.names(exprs.tech.cellmix.residuals) %in% probeset, ]

meth_probes <- methyl.tech.cellmix.residuals[row.names(methyl.tech.cellmix.residuals) %in% probelist, ]

meth_probes <- as.data.frame(t(meth_probes))
exp_probes <- as.data.frame(t(exp_probes))
colnames(exp_probes) <- paste("exp", colnames(exp_probes),  sep = "_")

probe.df <- cbind(meth_probes, exp_probes)
probe.df$SampleID <- row.names(meth_probes)

plot(probe.df$cg01039752, probe.df$exp_11760018_at)
     
cor(probe.df$cg01039752, probe.df$exp_11760018_at, use = "pairwise.complete.obs")

plot(probe.df$cg05575921, probe.df$exp_11737914_at)

cor(probe.df$cg05575921, probe.df$exp_11737914_at, use = "pairwise.complete.obs")

plot(probe.df$cg23314514, probe.df$exp_11758031_s_at)

cor(probe.df$cg23314514, probe.df$exp_11758031_s_at, use = "pairwise.complete.obs")

#### make plots of all the significant correlations

## first, standardize all the measures to mean = 0, SD = 1 (since the measures are residualized for technical variation so are 'meaningless')

describe(probe.df)

probe.df <- probe.df %>%
  as_tibble() %>%
  mutate(across(where(is.numeric),  ~ scale(.)[,1]))

probe.df <- as.data.frame(probe.df)

### Now, loop across all the significant correlations

## make a list of pairs of expression/ methylation probes to plot

correlations.df <- dplyr::select(correlations_Expmeth, ID, cg01039752:cg27055782)

pairlist  <- correlations.df %>%
  rowwise() %>%
  mutate(res = toString(names(correlations.df)[!is.na(c_across(starts_with('cg')))])) %>%
  select(ID, res)

pairlist.df <- pairlist %>%
  separate(res, c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8", "col_9", "col_10"))

pairlist.df <- pivot_longer(correlations.df, cols= !ID, names_to = "CpG", values_to = "correlation")

pairlist.df <- filter(pairlist.df, !is.na(correlation))

pairlist.df$ID <- paste("exp", pairlist.df$ID,  sep = "_") ## match the naming to the probe.df file

### make plots and save
### the pairs of variables to run are found in ID and CpG in pairlist.df


# Create a list of character vectors comprising pairs of values
pairs_list <- lapply(1:nrow(pairlist.df), function(i) {
  as.character(pairlist.df[i, c("ID", "CpG")])
})


### Add a smoking flag to the data


### use the daily smoker at 38 variable to get current p38 smokers (N = 215)
### use packyears at 38 = 0 to get never smokers

smoking.df <- read.table("./ExpressionvsMeth/data/Smoking_Dunedin_Mar2106.txt", sep="\t", header=T, stringsAsFactors=F)

smoking.df$SampleID <- sapply(smoking.df$snum, function(x) { if( nchar(x) == 3 ) { paste0("DN0", x) } else { paste0("DN", x) } })

smoke <- select(smoking.df, SampleID, dailysmoker38, packyears38)
smoke$Smokertype <- ifelse(is.na(smoke$dailysmoker38)|is.na(smoke$packyears38), NA,
                          ifelse(smoke$dailysmoker38 == 1, "Current Smoker",
                          ifelse(smoke$packyears38 == 0, "Never Smoker", "other")))

### add to the probe file for plotting

probe.df <- inner_join(probe.df, smoke, by = "SampleID")

probe.df <- filter(probe.df, !is.na(Smokertype))                       
                       
# Loop through variable pairs
for (pair in pairs_list) {
  y_var <- pair[1] # expression as the y-axis
  x_var <- pair[2]
  
  # Create scatterplot
  p <- ggplot(probe.df, aes_string(x = x_var, y = y_var)) +
    geom_point(aes(colour = Smokertype)) +
    geom_smooth(method = "lm", se = TRUE, aes(colour = Smokertype)) +
    labs(title = paste("Scatterplot of", x_var, "and", y_var),
         x = x_var, y = y_var) +
    theme_bw()


  # Save plot to a file with variable names
  ggsave(filename = paste0("./ExpressionvsMeth/results/plots/" , x_var, "_", y_var, "_scatterplot.pdf"), 
         plot = p)
}


###########################################################################################################

### violin plots of methylation x smoking

r <- c("lightpink3", "darkseagreen4", "skyblue3")

probe.df %>% group_by(Smokertype) %>% summarize(mean =mean(cg03636183),
                                                                      median = median(cg03636183),
                                                                      sd = sd(cg03636183),
                                                                      n =n()) -> smkN



p1 <-   ggplot(probe.df,aes(y=cg03636183, x=Smokertype)) + 
  geom_violin(alpha = 0.5, aes(color = Smokertype, fill = Smokertype)) + 
  scale_fill_manual(values = r) +
  scale_colour_manual(values = r) +
  geom_boxplot(width =0.15, alpha = 0.1, lwd= 0.5, fatten = 0.5, 
               color = "navyblue", fill = "white", outlier.shape = 1, outlier.alpha = 0.05) +
  geom_jitter(data = probe.df, aes(y = exp_11737914_at, x = Smokertype), 
              width = 0.1, colour = "navyblue", alpha = 0.2) + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.text.x=element_blank()) +
  xlab("Smoker type") + 
  geom_text(data = smkN, aes(x = Smokertype, y = -0.3, label = n), color ="black")



p2 <-   ggplot(probe.df,aes(y=cg03636183, x=Smokertype)) + 
  geom_violin(alpha = 0.5, aes(color = Smokertype, fill = Smokertype)) + 
  scale_fill_manual(values = r) +
  scale_colour_manual(values = r) +
  geom_boxplot(width =0.15, alpha = 0.1, lwd= 0.5, fatten = 0.5, 
               color = "navyblue",  outlier.shape = 1, outlier.alpha = 0.05) +
  geom_violin(data = probe.df, aes(y = exp_11737914_at, x = Smokertype, 
                                   color = Smokertype, fill = Smokertype), alpha = 0.15) + 
  geom_boxplot(data = probe.df, aes(y = exp_11737914_at, x = Smokertype), width =0.15, alpha = 0.05, 
               lwd= 0.5, fatten = 0.5, color = "navyblue", outlier.shape = 1, 
               outlier.alpha = 0.05) + theme_bw() + 
                theme(legend.title = element_blank(),
                     axis.text.x=element_blank()) +
  xlab("Smoker type") + 
  geom_text(data = smkN, aes(x = Smokertype, y = -5, label = n), color ="black")



#### or, plot with just two categories of smokers

r2 <- c("lightpink3", "darkseagreen4")

correlations.df <- dplyr::select(correlations_Expmeth, ID, cg01039752:cg27055782)

pairlist  <- correlations.df %>%
  rowwise() %>%
  mutate(res = toString(names(correlations.df)[!is.na(c_across(starts_with('cg')))])) %>%
  select(ID, res)

pairlist.df <- pairlist %>%
  separate(res, c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8", "col_9", "col_10"))

pairlist.df <- pivot_longer(correlations.df, cols= !ID, names_to = "CpG", values_to = "correlation")

pairlist.df <- filter(pairlist.df, !is.na(correlation))

pairlist.df$ID <- paste("exp", pairlist.df$ID,  sep = "_") ## match the naming to the probe.df file

### add the variable name for the smoking category for plotting

pairlist.df$smk <- "Smokertype2"

# Create a list of character vectors comprising pairs of values
pairs_list <- lapply(1:nrow(pairlist.df), function(i) {
  as.character(pairlist.df[i, c("ID", "CpG", "smk")])
})


probe.df$Smokertype2 <- ifelse(probe.df$Smokertype == "Never Smoker"|
                                 probe.df$Smokertype == "Current Smoker", probe.df$Smokertype, NA)
probe.df2 <- filter(probe.df, !is.na(Smokertype2))

# calculate Ns
   probe.df2 %>% group_by(Smokertype2) %>% summarize(n =n()) -> smkN2

   r2 <- c("lightpink3", "darkseagreen4")
   

# Loop through variable pairs
for (pair in pairs_list) {
  exp_var <- pair[1] # 
  meth_var <- pair[2]
  smk_var <- pair[3]

# create violin plots
   
p2 <-   ggplot(probe.df2, aes_string(y=meth_var, x=smk_var)) + 
  geom_violin(alpha = 0.5, aes(color = Smokertype2, fill = Smokertype2)) + 
  scale_fill_manual(values = r2) +
  scale_colour_manual(values = r2) +
  geom_boxplot(width =0.15, alpha = 0.1, lwd= 0.5, fatten = 0.5, 
               color = "navyblue",  outlier.shape = 1, outlier.alpha = 0.05) +
  geom_violin(data = probe.df2, aes_string(y = exp_var, x = smk_var, 
                                   color = smk_var, fill = smk_var), alpha = 0.15) + 
  geom_boxplot(data = probe.df2, aes_string(y = exp_var, x = smk_var), width =0.15, alpha = 0.05, 
               lwd= 0.5, fatten = 0.5, color = "navyblue", outlier.shape = 1, 
               outlier.alpha = 0.05) + theme_bw() + 
  theme(legend.title = element_blank(),
        axis.text.x=element_blank()) +
  xlab("Smoker type") 

# Save plot to a file with variable names
ggsave(filename = paste0("./ExpressionvsMeth/results/plots/" , meth_var, "_", exp_var, "_ViolinPlots.pdf"), plot = p2)


}

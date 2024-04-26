# plot AR6 scenario summaries
# implement this in Python!!!
library(ggplot2)
library(tidyr)
library(stringr)
library(readxl)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# intermediate outputs folder on google drive
GDRIVE = "H:/Shared drives/SBTi - Public Drive/Technical Department/5. Research Team/Current projects/1.5C Scenarios Review 2023/Intermediate analysis products/"
datadir <- "C:/Users/ginger.kowal/Documents/Scenario review/"

# harmonized filtered scenarios
em_df <- read.csv("C:/Users/ginger.kowal/Desktop/combined_em_df.csv")
# plot envelope of raw gross emissions
raw_gross_df <- em_df[
  em_df$Variable == "Emissions|CO2|Energy and Industrial Processes|Gross", ]

# envelope of harmonized gross emissions

# what's the effect of harmonization on reductions from 2022-2030, and 2022-2050?
var_list <- c('p|Emissions|CO2|Energy and Industrial Processes|s',
              'p|Emissions|CO2|Harmonized-DB')
year_list <- c(2023, 2030, 2050) # check base year
pr_sub <- aneris_harm[
  (aneris_harm$Model == 'model') &
    (aneris_harm$Variable %in% var_list) &
    (aneris_harm$Year %in% year_list),
  c('Scenario', 'Emissions', 'Variable', 'Year')]
pr_rs <- pivot_wider(pr_sub, names_from=Year, values_from=Emissions)
pr_rs$pr_2030 <- (pr_rs$`2030` - pr_rs$`2023`) / pr_rs$`2023`
pr_rs$pr_2050 <- (pr_rs$`2050` - pr_rs$`2023`) / pr_rs$`2023`
comp_sub <- pr_rs[, c('Scenario', 'Variable', 'pr_2030', 'pr_2050')]
comp <- pivot_wider(comp_sub, names_from=Variable, values_from=c(pr_2030, pr_2050))
summary(comp$`pr_2030_p|Emissions|CO2|Energy and Industrial Processes|s`)
summary(comp$`pr_2030_p|Emissions|CO2|Harmonized-DB`)
summary(comp$`pr_2050_p|Emissions|CO2|Energy and Industrial Processes|s`)
summary(comp$`pr_2050_p|Emissions|CO2|Harmonized-DB`)
comp$diff2030 <- (
  comp$`pr_2030_p|Emissions|CO2|Harmonized-DB` -
    comp$`pr_2030_p|Emissions|CO2|Energy and Industrial Processes|s`)
comp$diff2050 <- (
  comp$`pr_2050_p|Emissions|CO2|Harmonized-DB` -
    comp$`pr_2050_p|Emissions|CO2|Energy and Industrial Processes|s`
)
p <- ggplot(comp,
            aes(x=`pr_2030_p|Emissions|CO2|Energy and Industrial Processes|s`,
                y=`pr_2030_p|Emissions|CO2|Harmonized-DB`)) +
  geom_point() + geom_abline(slope=1, intercept=0)
print(p)
p <- ggplot(comp,
            aes(x=`pr_2050_p|Emissions|CO2|Energy and Industrial Processes|s`,
                y=`pr_2050_p|Emissions|CO2|Harmonized-DB`)) +
  geom_point() + geom_abline(slope=1, intercept=0)
print(p)

# ratio of CDR to net emissions, 2020-2050 and 2020-2100
er_df <- read.csv("C:/Users/ginger.kowal/Desktop/c1_cum_cdr_netem_2020-2050.csv")
er_df <- subset(er_df, select=-c(X))
rat_df <- pivot_wider(
  er_df, names_from=Variable, values_from=c(cum_sum_20.50, cum_sum_20.2100))
rat_df$rat_cdr_netem_2050 <- (
  rat_df$`cum_sum_20.50_Carbon Sequestration|Sum` /
    rat_df$`cum_sum_20.50_Emissions|CO2|Energy and Industrial Processes`)

# summarize CDR deployment in "ambitious" scenarios
em_df <- read.csv("C:/Users/ginger.kowal/Desktop/c1_net_gross_fossil_CO2_CDR.csv")
gr_em_df <- em_df[
  (em_df$Variable == "Emissions|CO2|Energy and Industrial Processes|Gross"),
  c('scen_id', 'filtered', 'gr_2020.50')]
cdr_df <- em_df[
  (em_df$Variable == "Carbon Sequestration|Sum"),
  c('scen_id', 'X2050')]
plot_df <- merge(gr_em_df, cdr_df, by='scen_id')
p <- ggplot(plot_df, aes(x=X2050, y=gr_2020.50)) + geom_point() +
  xlab('CDR in 2050 (Mt)') + ylab('% reduction gr emissions 2020-2050')
print(p)

# describe ambition of C1 scenarios
ar6_key <- read_excel(
  paste0(datadir, 'IPCC_AR6/AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx'),
  sheet='meta_Ch3vetted_withclimate')
ar6_key$scen_id <- paste(ar6_key$Model, ar6_key$Scenario)

c1_imps <- c('SP', 'LD', 'Ren')
imp_scen <- ar6_key[ar6_key$IMP_marker %in% c1_imps, 'scen_id']$scen_id

c1_grfosco2 <- read.csv("C:/Users/ginger.kowal/Desktop/c1_gross_fossil_CO2.csv")
c1_grfosco2$pchange_2030 <- (
  (c1_grfosco2$X2030 - c1_grfosco2$X2020) / c1_grfosco2$X2020)
c1_grfosco2$pchange_2050 <- (
  (c1_grfosco2$X2050 - c1_grfosco2$X2020) / c1_grfosco2$X2020)
hist(c1_grfosco2$pchange_2050)
amb_df <- c1_grfosco2[c1_grfosco2$pchange_2050 <= -0.89,
                      c('scen_id', 'pchange_2050')]
imp_grfco2 <- c1_grfosco2[c1_grfosco2$scen_id %in% imp_scen,
                          c('scen_id', 'pchange_2030', 'pchange_2050')]


# key variables from NZE, and CWF, and median of filtered scenarios
metric_list <- c(
  "Final energy demand, 2050 (EJ)",
  "Share of primary energy from renewables, 2050 (%)",
  "Total atmospheric CDR, 2050 (Gt CO2)",
  "Maximum yearly primary energy from bioenergy, 2020-2050 (EJ)",
  "Cumulative CCS, 2010-2050 (Gt CO2)",
  "Gross fossil CO2, 2020-2050 (Gt CO2)")
keyvar_df <- read.csv(paste0(GDRIVE, 'Summary_2050_metrics.csv'))
colnames(keyvar_df)[1] <- 'source'
focal_scen <- c("IEA NZE", "CWF Central", "Median of filtered scenarios")
keyvar_df$Scenario <- factor(
  keyvar_df$source,
  levels=c('IEA NZE', 'CWF Central', 'Median of filtered scenarios',
           'AR6 C1 (25th percentile)',
           'AR6 C1 (75th percentile)'),
  labels=c('NZE', 'CWF', 'Filtered',
           'AR6 C1-25', 'AR6 C1-75'))
keyvar_df <- keyvar_df[keyvar_df$Metric %in% metric_list, ]
keyvar_df$Metric <- factor(
  keyvar_df$Metric,
  levels=c("Final energy demand, 2050 (EJ)",
           "Share of primary energy from renewables, 2050 (%)",
           "Maximum yearly primary energy from bioenergy, 2020-2050 (EJ)",
           "Total atmospheric CDR, 2050 (Gt CO2)",
           "Cumulative CCS, 2010-2050 (Gt CO2)",
           "Gross fossil CO2, 2020-2050 (Gt CO2)"),
  labels=c('Final energy demand, 2050 (EJ)',
           'Share of primary energy from renewables, 2050 (%)',
           'Max yearly primary energy from bioenergy, 2010-2050 (EJ)',
           'CDR, 2050 (GtCO2)',
           'Total cumulative CCS, 2020-2050 (GtCO2)',
           'Gross fossil CO2, 2020-2050 (Gt)'))

hyb_df <- keyvar_df[keyvar_df$source %in% focal_scen, ]
ar6_df <- keyvar_df[
  keyvar_df$source %in% unique(keyvar_df$source)[5:6], ]
p <- ggplot(hyb_df, aes(x=Scenario, y=Value, fill=Scenario)) + geom_col() +
  scale_fill_manual(values=cbPalette) +
  geom_hline(data=ar6_df, aes(yintercept=Value), linetype='dashed') +
  facet_wrap(~Metric, scales='free',
             labeller=label_wrap_gen()) + 
  ylab("") + xlab("") + theme_bw() + scale_colour_grey() + theme(
    # axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    legend.position = "none")
    # legend.position=c(0.8, 0.3)) 
print(p)
filename <- paste0(GDRIVE, "20230919_key_variable_summary.png")
png(filename, width=6, height=4.5, units='in', res=300)
print(p)
dev.off()

# key variables from all compared scenarios
keyvar_df <- read.csv(paste0(GDRIVE, 'Summary_2050_metrics.csv'))
colnames(keyvar_df)[1] <- 'Scenario'
keyvar_df$Source <- paste(keyvar_df$Model, keyvar_df$Scenario)
keyvar_df$Source <- factor(
  keyvar_df$Source,
  levels=c(' IEA NZE', ' OECM April 2023 draft', ' CWF Central',
           'GCAM 5.3+ NGFS NGFS Orderly Net Zero',
           'REMIND-MAgPIE 3.0-4.4 NGFS Orderly Net Zero',
           'MESSAGEix-GLOBIOM 1.1-M-R12 NGFS Orderly Net Zero',
           ' AR6 C1 (25th percentile)',
           ' AR6 C1 (75th percentile)'),
  labels=c('NZE', 'OECM', 'CWF', 'NGFS-GCAM', 'NGFS-REMIND',
           'NGFS-MESSAGEix', 'AR6 C1-25', 'AR6 C1-75'))
keyvar_df$Metric <- factor(
  keyvar_df$Metric,
  levels=c("Final energy demand, 2050 (EJ)",
           "Share of primary energy from renewables, 2050 (%)",
           "Maximum yearly primary energy from bioenergy, 2020-2050 (EJ)",
           "Total atmospheric CDR, 2050 (Gt CO2)",
           "Cumulative CCS, 2010-2050 (Gt CO2)"),
  labels=c('Final energy demand, 2050 (EJ)',
           'Renewable energy share, 2050 (%)',
           'Max bioenergy, 2010-2050 (EJ)',
           'CDR, 2050 (GtCO2)',
           'Total CCS, 2020-2050 (GtCO2)'))
keyvar_df <- keyvar_df[!is.na(keyvar_df$Metric), ]

hyb_df <- keyvar_df[
  keyvar_df$Scenario %in% unique(keyvar_df$Scenario)[0:4], ]
ar6_df <- keyvar_df[
  keyvar_df$Scenario %in% unique(keyvar_df$Scenario)[5:6], ]
p <- ggplot(hyb_df, aes(x=Source, y=Value, fill=Source)) + geom_col() +
  geom_hline(data=ar6_df, aes(yintercept=Value), linetype='dashed') +
  facet_wrap(~Metric, scales='free') + 
  ylab("") + xlab("") + theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    legend.position=c(0.8, 0.3))
print(p)
filename <- paste0(GDRIVE, "key_variable_summary_with_NGFS.png")
png(filename, width=7.8, height=5, units='in', res=300)
print(p)
dev.off()

# compare AFOLU emissions from SR15 and AR6 scenarios
# CO2eq
afolu_df <- read.csv(paste0(GDRIVE, "afolu_co2e_sr15_ar6.csv"))
colnames(afolu_df)[2] <- 'afolu_co2e'

p <- ggplot(afolu_df, aes(x=gen, y=afolu_co2e)) + geom_boxplot() +
  geom_abline(intercept=217000, slope=0) + facet_wrap(~category) +
  xlab("") + ylab("AFOLU emissions, 2020-2050 (Mt CO2e)")
print(p)

# CO2 only
co2_df <- read.csv(paste0(GDRIVE, "afolu_co2_sr15_ar6.csv"))
colnames(co2_df)[2] <- 'afolu_co2'

p <- ggplot(co2_df, aes(x=gen, y=afolu_co2)) + geom_boxplot() +
  geom_abline(intercept=40000, slope=0, linetype='dashed') +
  facet_wrap(~category) +
  xlab("") + ylab("AFOLU CO2 emissions, 2020-2050 (Mt)")
print(p)

# look at gross EIP CO2 from all C1 scenarios, and ENGAGE project
datadir <- "C:/Users/ginger.kowal/Documents/Scenario review/"
ar6_key <- read_excel(
  paste0(datadir, 'IPCC_AR6/AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx'),
  sheet='meta_Ch3vetted_withclimate')
ar6_key$scen_id <- paste(ar6_key$Model, ar6_key$Scenario)
eng_id <- ar6_key[ar6_key['Project_study'] == 'ENGAGE', 'scen_id']

co2_df <- read.csv(
  paste0(GDRIVE, 'co2_c1_filtered.csv'), check.names=FALSE)
co2_df = subset(co2_df, select = -c(index))
select_id <- intersect(co2_df$scen_id, eng_id$scen_id)
plot_df <- pivot_longer(
  co2_df, cols=!c(scen_id), names_to='year', values_to='CO2')
engage_df <- plot_df[plot_df$scen_id %in% select_id, ]
plot_df$year <- as.numeric(plot_df$year)
engage_df$year <- as.numeric(engage_df$year)

p <- ggplot(data=NULL) +
  geom_line(data=plot_df, aes(x=year, y=CO2, group=scen_id), color='grey85') +
  geom_line(data=engage_df, aes(x=year, y=CO2, group=scen_id)) +
  labs(x='Year', y=expression(Gross~fossil~CO[2]~(Mt)))
print(p)
filename <- paste0(GDRIVE, "230915_C1_highlight_ENGAGE.png")
png(filename, width=4, height=4, units='in', res=300)
print(p)
dev.off()

# plot gross EIP CO2 from all C1 scenarios, and median of filtered scenarios
co2_df <- read.csv(
  paste0(GDRIVE, 'co2_c1_filtered.csv'), check.names=FALSE)
co2_df = subset(co2_df, select = -c(index))
c1_df <- co2_df[co2_df$scen_id != 'Median of filtered scenarios', ]
c1_2020q <- quantile(c1_df$`2020`, probs=c(0.25, 0.75))
c1_2030q <- quantile(c1_df$`2030`, probs=c(0.25, 0.75))
c1_2040q <- quantile(c1_df$`2040`, probs=c(0.25, 0.75))
c1_2050q <- quantile(c1_df$`2050`, probs=c(0.25, 0.75))
interq_df <- data.frame(
  year=c(2020, 2030, 2040, 2050),
  c1_25=c(c1_2020q[[1]], c1_2030q[[1]], c1_2040q[[1]], c1_2050q[[1]]),
  c1_75=c(c1_2020q[[2]], c1_2030q[[2]], c1_2040q[[2]], c1_2050q[[2]]),
  scen_id='NA')

plot_df <- pivot_longer(
  co2_df, cols=!c(scen_id), names_to='year', values_to='CO2')
plot_df$year <- as.numeric(plot_df$year)
med_df <- plot_df[plot_df$scen_id == 'Median of filtered scenarios', ]
p <- ggplot(data=NULL) +
  geom_line(data=plot_df, aes(x=year, y=CO2, group=scen_id), color='grey85') +
  geom_line(data=med_df, aes(x=year, y=CO2)) +
  geom_errorbar(data=interq_df, aes(x=year, ymin=c1_25, ymax=c1_75), width=0.7) + 
  theme_bw() +
  labs(x='Year', y=expression(Gross~fossil~CO[2]~(Mt)))
print(p)
filename <- paste0(GDRIVE, "230912_median_filtered_vs_C1.png")
png(filename, width=4, height=4, units='in', res=300)
print(p)
dev.off()

# CO2 from filtered scenario median, focal scenarios
focal_df <- read.csv(
  paste0(GDRIVE, 'gross_eip_co2_emissions_focal_scen.csv'))
colnames(focal_df) <- c(2020, 2030, 2040, 2050, 'scen_id')
focal_df <- focal_df[focal_df$scen_id %in% c('NZE', 'CWF'), ]
focal_l <- pivot_longer(
  focal_df, cols=!c(scen_id), names_to='year', values_to='CO2')
focal_l$year <- as.numeric(focal_l$year)
focal_l$scen_id <- factor(focal_l$scen_id)
med_df <- med_df[med_df$year < 2060, ]
cs_df <- rbind(focal_l, med_df)
cs_df$Scenario <- factor(
  cs_df$scen_id,
  levels=c('NZE', 'CWF', 'Median of filtered scenarios'),
  c('NZE', 'CWF', 'Filtered'))

p <- ggplot(data=NULL) +
  geom_line(data=cs_df, aes(
    x=year, y=CO2, group=Scenario, col=Scenario), size=1.1) +
  scale_colour_manual(values=cbPalette) +
  geom_ribbon(data=interq_df, aes(x=year, ymin=c1_25, ymax=c1_75), alpha=0.2) + 
  labs(x='Year', y=expression(Gross~fossil~CO[2]~(Mt))) + theme_bw() +
  theme(legend.position=c(0.7, 0.75))
print(p)
filename <- paste0(GDRIVE, "230919_median_filtered_vs_focal.png")
png(filename, width=4, height=4, units='in', res=300)
print(p)
dev.off()

# plot CO2eq from cross-sector pathway and all C1 scenarios
datadir <- "C:/Users/ginger.kowal/Documents/Scenario review/"
ar6_key <- read_excel(
  paste0(datadir, 'IPCC_AR6/AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx'),
  sheet='meta_Ch3vetted_withclimate')
ar6_key$scen_id <- paste(ar6_key$Model, ar6_key$Scenario)

c1_imps <- c('SP', 'LD', 'Ren')
imp_scen <- ar6_key[ar6_key$IMP_marker %in% c1_imps, 'scen_id']
c1a_scen <- ar6_key[ar6_key$Category_subset == 'C1a_NZGHGs', 'scen_id']

co2eq_df <- read.csv(
  paste0(GDRIVE, 'co2eq_c1_filtered.csv'), check.names=FALSE)
plot_df <- pivot_longer(
  co2eq_df, cols=!c(scen_id), names_to='year', values_to='CO2eq')
cs_df <- plot_df[plot_df$scen_id == 'cross-sector pathway', ]
imp_df <- plot_df[plot_df$scen_id %in% imp_scen$scen_id, ]
c1a_df <- plot_df[plot_df$scen_id %in% c1a_scen$scen_id, ]
p <- ggplot(data=NULL, aes(x=year, y=CO2eq, group=scen_id)) +
  geom_line(data=plot_df, color='grey85') +
  geom_line(data=cs_df) +
  # geom_line(data=imp_df, aes(color=scen_id), size=1) +
  geom_line(data=c1a_df, color='green') + 
  theme_bw() + xlab('Year') + ylab('Gross fossil CO2eq')
print(p)

filename <- paste0(GDRIVE, "230605_cross-sector_vs_C1_IMP.png")
png(filename, width=8, height=4, units='in', res=300)
print(p)
dev.off()

ar6_plot_years <- seq(2020, 2060, by=5)
sr15_plot_years <- seq(2020, 2060, by=10)
ar6_df <- read.csv(paste0(GDRIVE, "ar6_df.csv"), check.names=FALSE)
ar6_imp_df <- read.csv(paste0(GDRIVE, "ar6_imp.csv"), check.names=FALSE)
sr15_df <- read.csv(paste0(GDRIVE, "sr15_filtered_summary.csv"),
                    check.names=FALSE)

# exploring
# check relationship between gross and net CO2
ar6_df <- ar6_df[c(ar6_plot_years, 'Variable', 'source', 'perc')]
ar6_long <- pivot_longer(
  ar6_df, cols=!c(source, perc, Variable), names_to='year',
  values_to='emissions')
ar6_xy <- pivot_wider(
  ar6_long, names_from=Variable, values_from=emissions)

p <- ggplot(ar6_xy,
            aes(x=`Emissions|CO2|Energy and Industrial Processes`,
                y=`Emissions|CO2|Energy and Industrial Processes|Gross`))
p <- p + geom_point() + geom_abline(intercept=0, slope=1) + facet_wrap(~source)
print(p)

# compare net and gross emissions from SR15, AR6 C1, and AR6 C1 + IISD filters
ar6_df <- ar6_df[c(ar6_plot_years, 'Variable', 'source', 'perc')]
ar6_long <- pivot_longer(
  ar6_df, cols=!c(source, perc, Variable), names_to='year',
  values_to='emissions')
ar6_intq <- ar6_long[ar6_long$perc != 'median', ]

sr15_df <- sr15_df[c(sr15_plot_years, 'Variable', 'source', 'perc')]
sr15_long <- pivot_longer(
  sr15_df, cols=!c(source, perc, Variable), names_to='year',
  values_to='emissions')
sr15_intq <- sr15_long[(
  sr15_long$perc != 'median') &
    (sr15_long$Variable != "Carbon Sequestration|CCS|Biomass"), ]

plot_df <- rbind(ar6_intq, sr15_intq)
plot_df$group <- paste(plot_df$source, plot_df$perc)

p <- ggplot(plot_df, aes(x=year, y=emissions, group=group))
p <- p + geom_line(aes(linetype=source, colour=source))
p <- p + facet_wrap(~Variable) + xlab("")
print(p)

filename <- paste0(GDRIVE, "EIP_CO2_AR6_vs_SR15.png")
png(filename, width=11, height=5, units='in', res=300)
print(p)
dev.off()

# plot median
ar6_med <- ar6_long[ar6_long$perc == 'median', ]
sr15_med <- sr15_long[
  (sr15_long$perc == 'median') &
    (sr15_long$Variable != "Carbon Sequestration|CCS|Biomass"), ]
plot_df <- rbind(ar6_med, sr15_med)

p <- ggplot(plot_df, aes(x=year, y=emissions, group=source))
p <- p + geom_line(aes(linetype=source, colour=source))
p <- p + facet_wrap(~Variable) + xlab("")
print(p)

filename <- paste0(GDRIVE, "EIP_CO2_AR6_vs_SR15_median.png")
png(filename, width=11, height=5, units='in', res=300)
print(p)
dev.off()

# plot AR6 C1, C1+IISD, and IMPs
ar6_med <- ar6_long[(ar6_long$perc == 'median') &
                      (ar6_long$Variable == 'Emissions|CO2|Energy and Industrial Processes|Gross'), ]
ar6_med[ar6_med$source == 'AR6 C1', 'source'] <- 'AR6 C1 (median)'
ar6_med[ar6_med$source == 'AR6 C1+IISD', 'source'] <- 'AR6 C1+IISD (median)'
imp_df <- ar6_imp_df[c(sr15_plot_years, 'scen_id', 'Variable')]
imp_long <- pivot_longer(
  imp_df, cols=!c(scen_id, Variable), names_to='year',
  values_to='emissions')
colnames(imp_long)[1] <- 'source'
imp_long$perc <- NA
plot_df <- rbind(ar6_med, imp_long)

p <- ggplot(plot_df, aes(x=year, y=emissions, group=source))
p <- p + geom_line(aes(colour=source)) + xlab("")
p <- p + ylab("Gross EIP CO2 emissions")
print(p)

filename <- paste0(GDRIVE, "EIP_CO2_AR6_+_IMP.png")
png(filename, width=9, height=5, units='in', res=300)
print(p)
dev.off()

# plot sequestration via land use vs net and gross EIP CO2, all C1 scenarios
MAX_AFOLU = 3600
lu_df <- read.csv(paste0(GDRIVE, "max_lu_seq_vs_net_gross_EIP_CO2.csv"))
lu_df[is.na(lu_df$max_lu_seq), 'max_lu_seq'] = -1000
lu_df$filter_flag <- 'Not filtered'
lu_df[lu_df$max_lu_seq > MAX_AFOLU, 'filter_flag'] <- 'Removed by filter'

p <- (ggplot(lu_df, aes(x=filter_flag, y=net_2020_EIP_CO2))) + geom_boxplot()
print(p)
p <- (ggplot(lu_df, aes(x=filter_flag, y=net_2030_EIP_CO2))) + geom_boxplot()
print(p)
p <- (ggplot(lu_df, aes(x=filter_flag, y=net_2050_EIP_CO2))) + geom_boxplot()
print(p)

p <- (ggplot(lu_df, aes(x=filter_flag, y=gross_2020_EIP_CO2))) + geom_boxplot()
print(p)
p <- (ggplot(lu_df, aes(x=filter_flag, y=gross_2030_EIP_CO2))) + geom_boxplot()
print(p)
p <- (ggplot(lu_df, aes(x=filter_flag, y=gross_2050_EIP_CO2))) + geom_boxplot()
print(p)

p <- (ggplot(lu_df, aes(x=filter_flag, y=perc_ch_net_2030))) + geom_boxplot()
print(p)
p <- (ggplot(lu_df, aes(x=filter_flag, y=perc_ch_net_2050))) + geom_boxplot()
print(p)
p <- (ggplot(lu_df, aes(x=filter_flag, y=perc_ch_gross_2030))) + geom_boxplot()
print(p)
p <- (ggplot(lu_df, aes(x=filter_flag, y=perc_ch_gross_2050))) + geom_boxplot()
print(p)
t.test(perc_ch_gross_2030~filter_flag, data=lu_df)
t.test(perc_ch_gross_2050~filter_flag, data=lu_df)

p <- ggplot(lu_df, aes(x=max_lu_seq, y=net_2020_EIP_CO2)) +
  geom_point() + xlab("Max CO2 sequestered via landuse, 2010-2050") +
  ylab("Net energy & industrial process CO2, 2020")
print(p)

p <- ggplot(lu_df, aes(x=max_lu_seq, y=net_2030_EIP_CO2)) +
  geom_point() + xlab("Max CO2 sequestered via landuse, 2010-2050") +
  ylab("Net energy & industrial process CO2, 2030")
print(p)
filename <- paste0(GDRIVE, "230623_C1_landuse_sequestration_net_EIP_CO2.png")
png(filename, width=4, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggplot(lu_df, aes(x=max_lu_seq, y=gross_2030_EIP_CO2)) +
  geom_point() + xlab("Max CO2 sequestered via landuse, 2010-2050") +
  ylab("Gross energy & industrial process CO2, 2030")
print(p)
filename <- paste0(GDRIVE, "230623_C1_landuse_sequestration_gross_EIP_CO2.png")
png(filename, width=4, height=4, units='in', res=300)
print(p)
dev.off()

cor.test(lu_df$max_lu_seq, lu_df$net_2050_EIP_CO2)
cor.test(lu_df$max_lu_seq, lu_df$gross_2050_EIP_CO2)
cor.test(lu_df$max_lu_seq, lu_df$net_2030_EIP_CO2)
cor.test(lu_df$max_lu_seq, lu_df$gross_2030_EIP_CO2)

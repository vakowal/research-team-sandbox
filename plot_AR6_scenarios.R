# plot AR6 scenario summaries
# implement this in Python!!!
library(ggplot2)
library(tidyr)
library(stringr)
library(readxl)

# intermediate outputs folder on google drive
GDRIVE = "G:/.shortcut-targets-by-id/1rSoiKOBotDMn7VymKwxdpRQDr7ixLAhv/Research Team/04-Current projects/1.5C Scenarios Review 2023/Intermediate analysis products/"

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

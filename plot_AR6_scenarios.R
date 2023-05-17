# plot AR6 scenario summaries
# implement this in Python!!!
library(ggplot2)
library(tidyr)
library(stringr)

# intermediate outputs folder on google drive
GDRIVE = "G:/.shortcut-targets-by-id/1rSoiKOBotDMn7VymKwxdpRQDr7ixLAhv/Research Team/04-Current projects/1.5C Scenarios Review 2023/Intermediate analysis products/"

ar6_plot_years <- seq(2020, 2060, by=5)
sr15_plot_years <- seq(2020, 2060, by=10)
ar6_df <- read.csv(paste0(GDRIVE, "ar6_df.csv"), check.names=FALSE)
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

sr15_df <- sr15_df[c(sr15_plot_years, 'Variable', 'scenario_col', 'perc')]
sr15_long <- pivot_longer(
  sr15_df, cols=!c(scenario_col, perc, Variable), names_to='year',
  values_to='emissions')
sr15_intq <- sr15_long[(
  sr15_long$perc != 'median') &
    (sr15_long$Variable != "Carbon Sequestration|CCS|Biomass"), ]
colnames(sr15_intq)[2] <- 'source'

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

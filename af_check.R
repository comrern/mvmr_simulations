


hi_G <- df[df$X2 > 0.5,3:204]

low_G <- df[df$X2 < 0.1,3:204]

allele_frequencies_high <- colSums(hi_G[,1:200]) / (2 * nrow(hi_G))


allele_frequencies_low <- colSums(low_G[,1:200]) / (2 * nrow(low_G))


comp_tab <-data_frame(allele_frequencies_low)
comp_tab$hihg_x2 <- allele_frequencies_high
comp_tab$diff <-abs(comp_tab$allele_frequencies_low - comp_tab$hihg_x2)

mean(allele_frequencies_high)

mean(allele_frequencies_low)
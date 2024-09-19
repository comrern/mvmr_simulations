


hi_G <- dat[dat$X2 > 0.5,3:102]

low_G <- dat[dat$X2 < 0.1,3:102]

allele_frequencies_high <- colSums(hi_G[,1:100]) / (2 * nrow(hi_G))


allele_frequencies_low <- colSums(low_G[,1:100]) / (2 * nrow(low_G))


comp_tab <-data_frame(allele_frequencies_low)
comp_tab$hihg_x2 <- allele_frequencies_high


mean(allele_frequencies_high)

mean(allele_frequencies_low)
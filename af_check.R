


hi_G <- dat[dat$X2 > 0.5,3:202]

low_G <- dat[dat$X2 < 0.1,3:202]

allele_frequencies_high <- colSums(hi_G[,1:100]) / (2 * nrow(hi_G))


allele_frequencies_low <- colSums(low_G[,1:100]) / (2 * nrow(low_G))


mean(allele_frequencies_high)

mean(allele_frequencies_low)
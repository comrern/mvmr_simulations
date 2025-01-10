

hi_G <- df[df$X2 > 0.5,3:204]

low_G <- df[df$X2 < 0.1,3:204]

allele_frequencies_high <- colSums(hi_G[,1:200]) / (2 * nrow(hi_G))


allele_frequencies_low <- colSums(low_G[,1:200]) / (2 * nrow(low_G))


comp_tab <-data_frame(allele_frequencies_low)
comp_tab$hihg_x2 <- allele_frequencies_high
comp_tab$diff <-abs(comp_tab$allele_frequencies_low - comp_tab$hihg_x2)

mean(allele_frequencies_high)

mean(allele_frequencies_low)


####### for realistic sims (4 groups) ######

df1 <- df[1:(nrow(df)/4),]
df2 <- df[(nrow(df)/4):(nrow(df)/2),]
df3 <- df[(nrow(df)/2):(nrow(df)- (nrow(df)/4)),]
df4 <- df[(nrow(df)- (nrow(df)/4)):nrow(df),]

af1 <- as.data.frame(colSums(df1[,3:31]) / (2 * nrow(df1)))
af2 <- as.data.frame(colSums(df2[,3:31]) / (2 * nrow(df2)))
af3 <- as.data.frame(colSums(df3[,3:31]) / (2 * nrow(df3)))
af4 <- as.data.frame(colSums(df4[,3:31]) / (2 * nrow(df4)))

af1$subgroup <- "p1"
af2$subgroup <- "p2"
af3$subgroup <- "p3"
af4$subgroup <- "p4"

af1$SNP <- row.names(af1)
af2$SNP <- row.names(af2)
af3$SNP <- row.names(af3)
af4$SNP <- row.names(af4)

colnames(af1)[1] <- "af"
colnames(af2)[1] <- "af"
colnames(af3)[1] <- "af"
colnames(af4)[1] <- "af"

merge_df <- rbind(af1, af2, af3, af4)


ggplot(merge_df, aes(x = subgroup, y = af, group = SNP, color = SNP)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Minor allele frequency of selected SNPs across ancestry groups",
    x = "Subgroup",
    y = "MAF",
    color = "SNP"
  ) +
  theme_minimal() +
  guides(fill="none") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12)
  ) +
  guides(color="none")


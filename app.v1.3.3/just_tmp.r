library(dplyr)
testdata = read.table(file = 'tcpm.csv', sep = ',', header = TRUE) %>% select(-"Undetermined_S0")
for(i in colnames(testdata)[2:ncol(testdata)])
{
    testdata %>% 
        select(match(c("Geneid", i), names(.))) %>%
        write.table(paste0(i, ".tsv"), quote=F, row.names=F, sep="\t")
}
# Intron_Number
An R script to determine the number of introns per transcript.
The script is a bit messy, you have been warned...

```library(AnnotationHub)```

```library(AnnotationDbi)```

```library(intervals)```

```library(ggplot2)```

```library(dplyr)```

```library(data.table)```

```library(GenomicFeatures)```

```library(GenomicRanges)```

```library(rtracklayer)```

```library(devtools)```

```library(bedr)```

```library(scales)``` 

```library(splitstackshape)```

```library(rstatix)```

```library(tidyverse)```

```library(multcomp)```

```library(dgof)```

```library(Matching)```

## Species 1 
read in gtf file with pkg GenomicFeatures

```Scurra_txdb_parsed <- makeTxDbFromGFF(file="{PWD}/Scurria_scurra_annotation_v1_ch10_top10.gtf", format="gtf")```

```saveDb(Scurra_txdb_parsed, "Scurria_scurra_parsed")```

```S.scurra.p <- loadDb("Scurria_scurra_parsed")```

```columns(S.scurra.p)```

create a table of gene and transcript IDs
```txdf <- AnnotationDbi::select(S.scurra.p,keys=keys(S.scurra.p, "GENEID"),columns=c("GENEID", "TXID", "TXCHROM","TXNAME","EXONID", "EXONNAME"),keytype="GENEID")```

```head(txdf, 20)```

collect introns by transcript

```introns.list.per.transcript <- intronsByTranscript(S.scurra.p, use.names=TRUE)```

```mcols(introns.list.per.transcript) #get information about data collected```

```head(introns.list.per.transcript)```

coerce to dataframe

```introns.list.per.transcript.df.ss.all <- as.data.frame(introns.list.per.transcript)```

```head(introns.list.per.transcript.df.ss.all, n=50)```

```colnames(introns.list.per.transcript.df.ss.all) <- c("GENEID", "TXID","CHRMID","start","end","width","strand")```

```head(introns.list.per.transcript.df.ss.all)```

determine number of introns per transcript

```mcols(introns.list.per.transcript)$num_introns <- lengths(introns.list.per.transcript)```

```introns_per_transcript.ss <- as.data.frame(mcols(introns.list.per.transcript))```

```introns_per_transcript.ss$TXID <- row.names(introns_per_transcript.ss)```

```head(introns_per_transcript.ss)```

```class(introns_per_transcript.ss)```

merge based on TXID

```merged.ss <- dplyr::inner_join(introns_per_transcript.ss, introns.list.per.transcript.df.ss.all, by = "TXID")```

```head(merged.ss)```

subset only the columns that are needed

```merged.ss.reduced <- subset(merged.ss, select=c(num_introns, TXID, GENEID, CHRMID))```

```head(merged.ss.reduced)```

get only distinct rows based on TXID

```merged.ss.final <- merged.ss.reduced %>% distinct(TXID, .keep_all = T)```

```head(merged.ss.final)```

create vector with species names

```spID <- as.factor("sscurra")```

```merged.ss.final <- cbind(merged.ss.final, spID)```

## Repeat for species 1 to x 
create merged.speciesID.final dataframes for each species

## Merge dataframes

```merged_df <- rbind(merged.ss.final, merged.sv.final, merged.sz.final)```

Add column with Chromosome number

```merged_df$CHRMNUM <- substring(merged_df$CHRMID, 3)```

Add column with Orthologous chromosomes distinction based on MCSCAN collinearity plot (i.e. chr 6 in species 3 is not orthologous with chr 6 in species 1)

```merged_df$ORTHGRP <- merged_df$CHRMNUM```

```merged_df$ORTHGRP <- ifelse(merged_df$CHRMID == "sv3", 4, ifelse(merged_df$CHRMID == "sv4", 3, ifelse(merged_df$CHRMID == "sz6", 7, ifelse(merged_df$CHRMID == "sz7", 6, merged_df$CHRMNUM))))```

```head(merged_df)```

```tail(merged_df, n=100)```

```str(merged_df)```

```summary(merged_df)```

```boxplot(num_introns ~ ORTHGRP, data = merged_df)```

## Kruskal Wallis tests
```kruskal.test(num_introns ~ spID,data = merged_df)```

```pairwise.wilcox.test(merged_df$num_introns, merged_df$spID, p.adjust.method = "fdr")```

## PLOT 
ggplot(merged_df) +
  aes(x = spID, y = num_introns, color= spID) +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values=c("#FB61D7", "#A58AFF", "#00B6EB"))+
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="black", pch=16, size=.7) +
  stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="Number of Introns/Transcript") + 
  theme(legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))


library(treeio)
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)

metadata <- read_xlsx('/Users/zirui/BGI_Projects/GIZ_SouthernChinaBats/Metadata/TableS1.Samples information_2023-08-01.xlsx') |>
  mutate(Collection_year = year(Collection_time))

tr <- read.iqtree('/Users/zirui/BGI_Projects/GIZ_SouthernChinaBats/CloseToHuman/Aichivirus.treefile') 

tr_midpointed <- as_tibble(midpoint(as.phylo(tr))) |>
  left_join(select(as_tibble(tr), node, UFboot)) |>
  as.treedata()

tr_data <- tr |>
  as_tibble() |>
  extract(label, into = "seqid", regex = "(\\w+\\.\\d|\\w+-\\d{1,2}/\\w+)", remove = F) |> 
  mutate(seqid = str_replace(seqid, '_R_', '')) 

#tr_data |>
#  write.csv('/Users/zirui/BGI_Projects/GIZ_SouthernChinaBats/CloseToHuman/Adeno-associated.seqinfo.clean.csv')


pubseq_info <- read_excel('/Users/zirui/BGI_Projects/GIZ_SouthernChinaBats/CloseToHuman/华南进化树宿主信息.xlsx', sheet = 1) |>
  mutate(name = str_replace(name, '  ', ' ')) |>
  mutate(host = replace_na(host,'NA')) |>
  mutate(isolate = replace_na(isolate, 'NA')) |>
  mutate(labelanno = paste(seqid, ' ', name, '|', isolate, ' | ', host, sep='')) |>
  left_join(tr_data) |>
  mutate(source = 'public')

assembled_seqinfo <- tr_data |>
  filter(!(seqid %in% pubseq_info$seqid) & !(is.na(seqid))) |>
  separate(seqid, sep='/', remove=F, into=c('vsname', 'Sample_name')) |>
  left_join(select(metadata, Sample_name, City, Collection_year, Host_species)) |>
  mutate(labelanno = paste(vsname, ' ', Sample_name,'/', City, '/', Collection_year,' | ', Host_species, sep='')) |>
  mutate(source = 'asembly')

allseqinfo <- bind_rows(pubseq_info, assembled_seqinfo) |>
  select(label, labelanno, source)

p3 <- ggtree(tr_midpointed) +
  xlim(0,5)

p3 <- p3 %<+% allseqinfo +
  geom_tiplab(mapping = aes(x=x+0.05,label=labelanno, color=source)) +
  geom_tippoint(mapping = aes(color=source)) +
  scale_color_manual(values=c('#FFA93C','black')) +
  geom_treescale(x=0,y=20,color="black") +
  geom_nodelab(hjust=1.1 , vjust=-0.3, size=3)
ggsave('/Users/zirui/BGI_Projects/GIZ_SouthernChinaBats/CloseToHuman/Fig3.Aichivirus.pdf', p3, height = 5, width = 15)

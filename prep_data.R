library(tidyverse)

# load data for the app
# transcript data, formatted for plot
transcript_data <- read_tsv("../cDNA_files_r_shiny/annotation_r_shiny.tsv", lazy = TRUE) %>% 
  rename(seqnames = chr) %>%
  mutate(gene_name = coalesce(gene_name, gene_id))

#TODO update these seqnames_levels
seqnames_levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9',
                    '10', '11', '12', '13', '14', '15', '16', '17', 
                    '18', '19', '20', '21', '22', 'X', 'Y', 'MT',
                    'ERCC-00002', 'ERCC-00003', 'ERCC-00004', 'ERCC-00009',
                    'ERCC-00012', 'ERCC-00013', 'ERCC-00014', 'ERCC-00016',
                    'ERCC-00017', 'ERCC-00019', 'ERCC-00022', 'ERCC-00024',
                    'ERCC-00025', 'ERCC-00028', 'ERCC-00031', 'ERCC-00033',
                    'ERCC-00034', 'ERCC-00035', 'ERCC-00039', 'ERCC-00040',
                    'ERCC-00041', 'ERCC-00042', 'ERCC-00043', 'ERCC-00044',
                    'ERCC-00046', 'ERCC-00048', 'ERCC-00051', 'ERCC-00053',
                    'ERCC-00054', 'ERCC-00057', 'ERCC-00058', 'ERCC-00059',
                    'ERCC-00060', 'ERCC-00061', 'ERCC-00062', 'ERCC-00067',
                    'ERCC-00069', 'ERCC-00071', 'ERCC-00073', 'ERCC-00074',
                    'ERCC-00075', 'ERCC-00076', 'ERCC-00077', 'ERCC-00078',
                    'ERCC-00079', 'ERCC-00081', 'ERCC-00083', 'ERCC-00084',
                    'ERCC-00085', 'ERCC-00086', 'ERCC-00092', 'ERCC-00095',
                    'ERCC-00096', 'ERCC-00097', 'ERCC-00098', 'ERCC-00099',
                    'ERCC-00104', 'ERCC-00108', 'ERCC-00109', 'ERCC-00111',
                    'ERCC-00112', 'ERCC-00113', 'ERCC-00116', 'ERCC-00117',
                    'ERCC-00120', 'ERCC-00123', 'ERCC-00126', 'ERCC-00130',
                    'ERCC-00131', 'ERCC-00134', 'ERCC-00136', 'ERCC-00137',
                    'ERCC-00138', 'ERCC-00142', 'ERCC-00143', 'ERCC-00144',
                    'ERCC-00145', 'ERCC-00147', 'ERCC-00148', 'ERCC-00150',
                    'ERCC-00154', 'ERCC-00156', 'ERCC-00157', 'ERCC-00158',
                    'ERCC-00160', 'ERCC-00162', 'ERCC-00163', 'ERCC-00164',
                    'ERCC-00165', 'ERCC-00168', 'ERCC-00170', 'ERCC-00171',
                    'GL000009.2', 'GL000194.1', 'GL000195.1', 'GL000205.2',
                    'GL000213.1', 'GL000214.1', 'GL000216.2', 'GL000218.1',
                    'GL000219.1', 'GL000220.1', 'GL000224.1', 'GL000225.1',
                    'KI270442.1', 'KI270519.1', 'KI270706.1', 'KI270711.1',
                    'KI270713.1', 'KI270721.1', 'KI270726.1', 'KI270727.1',
                    'KI270728.1', 'KI270731.1', 'KI270733.1', 'KI270734.1',
                    'KI270742.1', 'KI270743.1', 'KI270744.1', 'KI270750.1')

transcript_data$seqnames <- fct_rev(factor(transcript_data$seqnames, levels = seqnames_levels))

new_genes <- transcript_data %>%
  filter(type == "transcript", startsWith(gene_id, 'BambuGene')) %>%
  arrange(seqnames)

new_transcripts <- transcript_data %>%
  filter(type == "transcript", startsWith(gene_id, 'E'), startsWith(transcript_id, 'BambuTx')) %>%
  arrange(seqnames)

# sample status
# file created from the excel file about the samples
sample_status <- read_tsv("../cDNA_files_r_shiny/cDNA_sample_info.tsv") %>%
  mutate(sample_status = ifelse(group == 2, "AD", "control")) %>%
  mutate(sample_sex = ifelse(sex == 2, "female", "male")) %>%
  select(sample_name, sample_status, sample_sex)

# expression data, formatted for plot
expression_data <- read_tsv("../cDNA_files_r_shiny/expression_matrix_r_shiny.tsv") %>%
  left_join(transcript_data %>% select(seqnames, transcript_id, annotation_status, discovery_category, transcript_biotype) %>% distinct(), by = "transcript_id") %>%
  pivot_longer(cols=3:62, names_to = "sample_expression_type", values_to = "number") %>%
  extract(sample_expression_type, into = c("sample_name", "expression_type"), "^([^_]*_[^_]*_[^_]*)_(.*$)") %>%
  pivot_wider(names_from = expression_type, values_from = number) %>%
  left_join(sample_status)

anti_sense_data <- read_tsv("../cDNA_files_r_shiny/AllNovelTranscripts_overlap_All_HG38_v107_Transcripts.tsv") %>% 
  filter(NewGeneID != Hg38v107_GeneID) %>%
  #filter(StrandType == "OppositeStrand", startsWith(NewGeneID, "BambuG"), startsWith(Hg38v107_GeneID, "E")) %>% 
  filter(StrandType != "NoOverlap", startsWith(NewGeneID, "BambuG")) %>% 
  select(!Hg38v107_transcriptID) %>% 
  distinct()

# contains just gene_id and gene_names
just_genes <- transcript_data %>% select(gene_id, gene_name) %>% distinct()

# table of new genes that could potentially be anti-sense
is_antisense <- left_join(anti_sense_data %>% 
                            select(NewGeneID, Hg38v107_GeneID), 
                          just_genes, 
                          by = c("NewGeneID" = "gene_id")) %>%
  rename(gene_id = NewGeneID) %>%
  rename(antisense_id = Hg38v107_GeneID) #%>%
# add_column(is_antisense = "true") %>%
# add_column(has_antisense = "false")

# table of know genes that we potentially found an anti-sense transcript for
has_antisense <- left_join(anti_sense_data %>%
                             select(Hg38v107_GeneID, NewGeneID),
                           just_genes, 
                           by = c("Hg38v107_GeneID" = "gene_id")) %>%
  rename(gene_id = Hg38v107_GeneID) %>%
  rename(antisense_id = NewGeneID) #%>%
# add_column(has_antisense = "true") %>%
# add_column(is_antisense = "false")

# combine into single table for lookup
just_genes_overlap <- bind_rows(is_antisense, has_antisense)

gene_lookup <- transcript_data %>%
  select(gene_id, gene_name) %>%
  distinct(gene_id, gene_name) 

bambu_gene_lookup <- gene_lookup %>%
  filter(grepl('Bambu', gene_id)) %>%
  arrange(gene_id) %>%
  mutate(key_file = 'B')
  
ensembl_gene_lookup <- gene_lookup %>%
  filter(!grepl('Bambu', gene_id)) %>%
  arrange(gene_id) %>%
  mutate(key_file = as.character(floor(row_number()/256) + 1))

gene_lookup <- bind_rows(bambu_gene_lookup, ensembl_gene_lookup)

for (i in (gene_lookup %>% select(key_file) %>% distinct(key_file) %>% pull(key_file))) {
  gene_ids <- gene_lookup %>% filter(key_file == i) %>% pull(gene_id)
  write_tsv(transcript_data %>% filter(gene_id %in% gene_ids), paste0("data_files/txd_",i, ".tsv"))
  write_tsv(expression_data %>% filter(gene_id %in% gene_ids), paste0("data_files/exd_", i, ".tsv"))
}

write_tsv(new_genes, "data_files/new_genes.tsv")
write_tsv(new_transcripts, "data_files/new_transcripts.tsv")
write_tsv(sample_status, "data_files/sample_status.tsv")
write_tsv(anti_sense_data, "data_files/anti_sense_data.tsv")
write_tsv(just_genes, "data_files/just_genes.tsv")
write_tsv(just_genes_overlap, "data_files/antisense.tsv")
write_tsv(gene_lookup, "data_files/gene_lookup.tsv")

# Density plot data - CPM
density_data <- read_tsv("../cDNA_files_r_shiny/expression_matrix_r_shiny_GENE.tsv") %>%
  rename(CPM = median_cpm) %>%
  rename(counts = median_counts)%>%
  filter(CPM > 0)%>%
  mutate(log_comb_exp = log10(CPM)) %>% 
  mutate(percentile = percent_rank(log_comb_exp)*100)

plt <- ggplot(density_data, aes(x=log_comb_exp)) +
  geom_density()

save(plt, density_data, file = 'data_files/density_base_CPM.Rdata')

# Density plot data - counts
density_data <- read_tsv("../cDNA_files_r_shiny/expression_matrix_r_shiny_GENE.tsv") %>%
  rename(CPM = median_cpm) %>%
  rename(counts = median_counts)%>%
  filter(counts > 0)%>%
  mutate(log_comb_exp = log10(counts)) %>% 
  mutate(percentile = percent_rank(log_comb_exp)*100)

plt <- ggplot(density_data, aes(x=log_comb_exp)) +
  geom_density()

save(plt, density_data, file = 'data_files/density_base_counts.Rdata')

#write_tsv(
anno_color <- transcript_data %>%
            select(annotation_status) %>%
            distinct(annotation_status) %>%
            mutate(id = row_number())
#, 'data_files/display_color_annotation_status.tsv')

#write_tsv(
tx_color <- transcript_data %>%
            select(transcript_biotype) %>%
            distinct(transcript_biotype) %>%
            mutate(id = row_number())

#            mutate(display_color = transcript_biotype), 'data_files/display_color_transcript_biotype.tsv')

#write_tsv(
discov_color <- transcript_data %>%
            select(discovery_category) %>%
            distinct(discovery_category) %>%
            mutate(id = row_number())

write_tsv(full_join(discov_color, tx_color) %>%
            full_join(anno_color), 'data_files/display_color.tsv')


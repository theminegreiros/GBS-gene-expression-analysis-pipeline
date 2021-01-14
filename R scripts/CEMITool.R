library("CEMiTool")
data(expr0)
head(expr0[,1:4])
cem <- cemitool(expr0)
cem
nmodules(cem)
head(module_genes(cem))
generate_report(cem)
write_files(cem)
save_plots(cem, "all")
data(sample_annot)
head(sample_annot)
cem <- cemitool(expr0, sample_annot)
sample_annotation(cem, 
                  sample_name_column="SampleName", 
                  class_column="Class") <- sample_annot
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]

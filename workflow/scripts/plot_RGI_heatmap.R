log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

rgi = read.csv(snakemake@input$allele_mapping, sep="\t")

perc_cov <- snakemake@params$cov_percent

# the first gsub gets rid of spaces in drug class
rgi_filtered <- rgi |>
    dplyr::filter( Percent.Coverage > perc_cov ) |>
    dplyr::mutate(Drug.Class = gsub("; ", ";", Drug.Class)) |>
    tidyr::separate_rows(Drug.Class, sep = ";") |>
    dplyr::group_by(Drug.Class, Resistomes...Variants..Observed.Pathogen.s.) |>
    dplyr::summarize(All.Mapped.Reads=sum(All.Mapped.Reads), mean_cov = mean(Percent.Coverage))

RGI_plot <- rgi_filtered  |>
    ggplot2::ggplot(ggplot2::aes(x=Resistomes...Variants..Observed.Pathogen.s., y=Drug.Class,color=All.Mapped.Reads, alpha=mean_cov, size=All.Mapped.Reads)) +
    ggplot2::geom_point(shape=16, stroke=0) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::scale_size_continuous(range=c(2,5)) +
    ggplot2::scale_color_viridis_c(option="inferno", direction = -1, end = .8) +
    ggplot2::labs(
        x="Associated Pathogen", y= "Antibiotic Class",
        subtitle=paste("Antibiotic resistances by taxa\nRemoved hits with less than ",
                       perc_cov,
                       "% covered length\nSome hits may confer multiple resistances\nHits reflect taxa associated with cannonical resistance gene"))

ggplot2::ggsave(
    RGI_plot, filename = snakemake@output$heatmap,
    width=max(
        7,
        min(
            length(unique(rgi_filtered$Resistomes...Variants..Observed.Pathogen.s.))/2,
            36
        )
    )
)

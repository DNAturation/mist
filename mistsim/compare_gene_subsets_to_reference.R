library(getopt)
library(parallel)
library(ggplot2)
library(gridExtra)

source("./partition_metrics.R")

compare_partitions <- function(seed, n_genes, iterations, thresh, ref_thresholds) {
    # Calculate the Adjusted Wallace Coefficient and Cluster Cohesion
    # between the current subset and the reference

    iteration <- iterations[,seed] #iterations = cutree

    ref_thresh <- ref_thresholds[,thresh]

    awc <- adj_wallace(iteration, ref_thresh)

    list(
    seed = seed,
    reference_threshold = colnames(ref_thresholds)[thresh],
    genes = n_genes,
    AWC_subset_vs_ref = awc$Adjusted_Wallace_A_vs_B,
    AWC_ref_vs_subset = awc$Adjusted_Wallace_B_vs_A,

    cluster_cohesion_subset_vs_ref = taboada_cohesion(iteration, ref_thresh)$weighted_average,
    cluster_cohesion_ref_vs_subset = taboada_cohesion(ref_thresh, iteration)$weighted_average
    )
}

compare_to_ref <- function(n_genes, method_iterations, ref_thresholds, cores) {

    l <- mclapply(1:ncol(ref_thresholds), function(i) {

        lapply(seq_along(method_iterations),
                         compare_partitions,
                         iterations = method_iterations,
                         ref_thresh = ref_thresholds,
                         thresh = i,
                         n_genes = n_genes
                        )


    }, mc.cores = cores)

    # kludge to convert to nice data.frame
    df <- do.call('rbind', do.call('rbind', l))
    df <- df[, colnames(df)[c(3, 1, 2, 4, 5, 6, 7)] ]

    df
}

load_data <- function(precomp_cluster_dir, cores, subset) {

    extract_info <- function(x) {
        # Assumes CSV tables are named in the format of "NNN_clusters.csv"
        # where the Ns are integers representing the number of genes
        # used to calculate the clusters

        name <- basename(tools::file_path_sans_ext(x))
        pos <- c(gregexpr("\\d", name))
        int_name <- as.integer(substr(name, min(unlist(pos)), max(unlist(pos))))

        list(
            genes = int_name,
            clusters = read.csv(x, sep=',', row.names = 1, check.names = FALSE)
            )

    }

    files <- list.files(precomp_cluster_dir, full.names = TRUE)
    reference_index <- grep("reference\\.csv", files) #### I changed this, switched around the args
    reference <- read.csv(files[reference_index],
                          row.names = 1, sep = ',',
                          check.names = FALSE)

    cut_index <- grep(paste(subset, "_clusters\\.csv", sep=''), files) ###changed original files to the new one, now expects other files in same directory, as well as the cut tree files following the formatting
    cutfiles<-files[cut_index]

    list("submethods" = mclapply(cutfiles, extract_info, mc.cores = cores),
         "reference" = reference)

}

plotter <- function(outfile, outplot, replicates, subset)  { #makes a violin plot of the subset vs ref Adjusted Wallace at different replicates compared across different cutree heights
    wd<-getwd()
    data<-read.csv(paste(wd, '/', outfile, sep=''))
    vioplot<-ggplot(data, aes(factor(reference_threshold), AWC_subset_vs_ref))+geom_violin()+xlab('Cutree heights')+
    ylab('Adjusted Wallace coefficients of subset against reference')+ggtitle(paste('Plot of', subset, 'genes at', replicates, 'replicates'))
    return(vioplot)
    }


options <- function() {

    spec <- matrix(c(
        "help",  "h", "0", "logical",   "Print this help and exit",
        "input", "i", "1", "character", "Path to directory of cluster tables",
        "out",   "o", "1", "character", "Output file",
        "cores", "c", "1", "integer",   "Number of CPU cores to use",
        "plot", "p", "1", "character", "Output plot name",
        "subset", "n", "1", "character", "number of genes in subset",
        "replicates", "r", "1", "integer", "Number of replicates"
    ), byrow = TRUE, ncol = 5)

    opt <- getopt(spec)

    if (!is.null(opt$help)) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
    }

    opt
}

main <- function() {

    opt <- options()

    data <- load_data(opt$input, opt$cores, opt$subset)

    notsingle<-sapply(data$reference, function(x) length(unique(x))>1)
    reffy<-data$reference[, notsingle]

    l <- lapply(data$submethods, function(submethod) {

        compare_to_ref(n_genes = submethod$genes,
#                       method_iterations = 1:ncol(submethod$clusters),
                       method_iterations = submethod$clusters, #changed from 1:ncol to this due to actually needing the
                                                                # clusters information instead of just how many columns
                                                                # there are
                       ref_thresholds = reffy, ###I changed data$ref_thresholds to data$reference, and then to reffy to subset
                       cores = opt$cores)
    })

    out_df <- do.call('rbind', l)

    setwd(opt$input)
    wd<-getwd()
    write.csv(out_df, file = opt$out, row.names = FALSE, quote = FALSE)
    vioplot<-plotter(outfile=opt$out, outplot=opt$plot, replicates=opt$replicates, subset=opt$subset)
    png(paste(wd, '/', opt$plot, sep=''))
    grid.arrange(vioplot)
    nuttin<-dev.off()
}


main()
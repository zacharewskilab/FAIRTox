#source('./global.R')
source('./ui.R')

#' Create dataframe of normalized enrichment scores for use in GSEA tab
#' 
#' @param db_genes The dataframe pulled from dbZach database containing all genes and fold changes
#' @param user_gene_list The list of genes to be compared. Input by the user
#' @return A database of enrichment scores for each project and fold change combination
#' @examples
#' enrichmentScores <- tzGSEA(df, gene_list)
tzGSEA <- function(db_genes, user_gene_list, metadata_selected){
  # Instantiate vectors to store enrichment scores, names of project/dose combos, and the corresponding p values (adjusted)
  enrichmentScores <- vector()
  names <- vector()
  pVals <- vector()
  pAdjs <- vector()
  fgres.list = list()
  
  for (i in 1:length(projects$Longname)){
    # Filter by project name
    subset <- filter(db_genes, Longname == projects[,2][i])
    
    # Filter by metadata based on user selection
    if ("Species" %in% metadata_selected){
      for(j in unique(subset$Species_common)){
        subset_subset <- filter(subset, Species_common == j)
        
        fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
        
        curr_project_species <- paste0(projects[,2][i], '/', j)
        fgres.list[[curr_project_species]] = fgRes
        
        # Store normalized enrichment score, project/dose combo for reference, and p value
        enrichmentScores <- c(enrichmentScores, fgRes$NES)
        names <- c(names, curr_project_species)
        pVals <- c(pVals, fgRes$pval)
        pAdjs <- c(pAdjs, fgRes$padj)
      }
    }
    else if ("Sex" %in% metadata_selected){
      for(j in unique(subset$Sex)){
        subset_subset <- filter(subset, Sex == j)
        
        fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
        
        curr_project_sex <- paste0(projects[,2][i], '/', j)
        fgres.list[[curr_project_sex]] = fgRes
        
        # Store normalized enrichment score, project/dose combo for reference, and p value
        enrichmentScores <- c(enrichmentScores, fgRes$NES)
        names <- c(names, curr_project_sex)
        pVals <- c(pVals, fgRes$pval)
        pAdjs <- c(pAdjs, fgRes$padj)
      }
    }
    else if ("Organ" %in% metadata_selected){
      for(j in unique(subset$Organ_name)){
        subset_subset <- filter(subset, Organ_name == j)
        
        fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
        
        curr_project_organ <- paste0(projects[,2][i], '/', j)
        fgres.list[[curr_project_organ]] = fgRes
        
        # Store normalized enrichment score, project/dose combo for reference, and p value
        enrichmentScores <- c(enrichmentScores, fgRes$NES)
        names <- c(names, curr_project_organ)
        pVals <- c(pVals, fgRes$pval)
        pAdjs <- c(pAdjs, fgRes$padj)
      }
    }
    else if ("Treatment" %in% metadata_selected){
      for(j in unique(subset$Chemical_Name)){
        subset_subset <- filter(subset, Chemical_Name == j)
        
        fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
        
        curr_project_chemical <- paste0(projects[,2][i], '/', j)
        fgres.list[[curr_project_chemical]] = fgRes
        
        # Store normalized enrichment score, project/dose combo for reference, and p value
        enrichmentScores <- c(enrichmentScores, fgRes$NES)
        names <- c(names, curr_project_chemical)
        pVals <- c(pVals, fgRes$pval)
        pAdjs <- c(pAdjs, fgRes$padj)
      }
    }
    else if ("Assay" %in% metadata_selected){
      for(j in unique(subset$Assay_Name)){
        subset_subset <- filter(subset, Assay_Name == j)
        
        fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
        
        curr_project_assay <- paste0(projects[,2][i], '/', j)
        fgres.list[[curr_project_assay]] = fgRes
        
        # Store normalized enrichment score, project/dose combo for reference, and p value
        enrichmentScores <- c(enrichmentScores, fgRes$NES)
        names <- c(names, curr_project_assay)
        pVals <- c(pVals, fgRes$pval)
        pAdjs <- c(pAdjs, fgRes$padj)
      }
    }
    else if ("Strain" %in% metadata_selected){
      for(j in unique(subset$Strain_Name)){
        subset_subset <- filter(subset, Strain_Name == j)
        
        fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
        
        curr_project_strain <- paste0(projects[,2][i], '/', j)
        fgres.list[[curr_project_strain]] = fgRes
        
        # Store normalized enrichment score, project/dose combo for reference, and p value
        enrichmentScores <- c(enrichmentScores, fgRes$NES)
        names <- c(names, curr_project_strain)
        pVals <- c(pVals, fgRes$pval)
        pAdjs <- c(pAdjs, fgRes$padj)
      }
    }
    else if ("Dose" %in% metadata_selected){
      for (j in unique(subset$Dose)){
        subset_subset <- filter(subset, Dose == j)
        
        fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
        
        curr_project_dose <- paste0(projects[,2][i], '/', j)
        fgres.list[[curr_project_dose]] = fgRes
        
        # Store normalized enrichment score, project/dose combo for reference, and p value
        enrichmentScores <- c(enrichmentScores, fgRes$NES)
        names <- c(names, curr_project_dose)
        pVals <- c(pVals, fgRes$pval)
        pAdjs <- c(pAdjs, fgRes$padj)
        
        print(length(unique(names)))
        print(length(unique(enrichmentScores)))
        print(length(unique(pAdjs)))
      }
    }
  }
  
  # Build dataframe containing all relevant data
  full_results <- data.frame("project_dose" = names, "NES" = enrichmentScores, "pVal" = pVals, "pAdj" = pAdjs)
  
  # Return the results
  return(full_results)
}

#' Calculate the normalized enrichment score for a particular project/metadata subset
#' 
#' @param subset_subset The project/metadata subset to be ranked
#' @param user_gene_list The list of genes to be compared. Input by the user
#' @return An fgRes object that contains NES as well as p values
#' @examples
#' fgRes <- calculateEnrichmentScore(subset_subset, user_gene_list)
calculateEnrichmentScore <- function(subset_subset, user_gene_list){
  # Create ranked list of gene name and fold change (sorted largest to smallest) - ties handled randomly to eliminate bias
  print(head(subset_subset))
  ranked_list <- subset_subset[,c(1,4)]
  ranked_list <- ranked_list[order(-ranked_list$FoldChange), ]
  ranked_list <- ranked_list[!duplicated(ranked_list$Symbol), ]
  #ranked_list <- ranked_list %>% mutate(rank = rank(FoldChange, ties.method = 'random'))
  
  # Convert ranked list to vector
  ranked_list_vector <- c(ranked_list$FoldChange)
  names(ranked_list_vector) <- ranked_list$Symbol
  
  print(head(ranked_list_vector))
  
  # Perform GSEA with 1000 permutations
  fgRes <- fgsea::fgsea(pathways = list("mylist" = user_gene_list), #User gene list in correct format
                        stats = ranked_list_vector, #Data extracted from DB and ranked
                        minSize=5, #Ignores gene lists smaller than this
                        maxSize=600, #Doesn't run on lists that are too big
                        nperm = 500)
  
  # Display and return fgRes dataframe
  print(fgRes)
  return(fgRes)
}

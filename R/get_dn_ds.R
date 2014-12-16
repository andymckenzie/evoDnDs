#' Dn/Ds Call Function
#'
#' This function queries biomaRt and ensembl to find the dn/ds values for a character vector of gene names.
#' @param species_name Ensembl/biomaRt namespace consistent species name. 
#' @param list_of_genes Character vector of hgnc symbol names to search. 
#' @keywords evolution dn/ds
#' @export
#' @examples
#' get_dn_ds

get_dn_ds <- function(species_name, list_of_genes){
	library(biomaRt)
	mart_hs = useMart("ensembl", "hsapiens_gene_ensembl")
	vector_of_dn = rep(0, length(list_of_genes))
	vector_of_ds = rep(0, length(list_of_genes))
	for(i in 1:length(list_of_genes)){
		dns_df = getBM(attributes = c(paste(species_name, "homolog_ds", sep = "_"),
							paste(species_name, "homolog_dn", sep = "_")), 
							filters = "hgnc_symbol", values = list_of_genes[i], mart = mart_hs)
		Sys.sleep(1)
		if(length(dns_df[[1]]) == 0){
			next
		}
		non_NA_dn = dns_df[[2]][which(!is.na(dns_df[[2]]))]
		non_NA_ds = dns_df[[1]][which(!is.na(dns_df[[1]]))]
		print(non_NA_dn)
		print(non_NA_ds)
		vector_of_dn[i] = median(non_NA_dn)
		vector_of_ds[i] = median(non_NA_ds)		
	}
	vector_of_dn_ds_ratio = vector_of_dn/vector_of_ds
	dn_ds_df = data.frame(list_of_genes, vector_of_dn, vector_of_ds, vector_of_dn_ds_ratio)
	names(dn_ds_df) = c("Gene_Names", "DN", "DS", "DN_DS_Ratio")
	return(dn_ds_df)
}
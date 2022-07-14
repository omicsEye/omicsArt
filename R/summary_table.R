# #tableone
#
# #tableone::CreateCatTable(data=metadata, vars = c("Group", "Sex", ))
#
#
# cont_cols <- colnames(metadata)[sapply(metadata, function(x) length(unique(x)))>10]
# cat_cols <- setdiff(colnames(metadata),cont_cols)
#
#
# #create summary data for continuous variables
# res_cont <- data.frame(do.call(rbind,lapply(metadata[,cont_cols],summary)))
#
# #add column with variable name (unwise to store in rownames)
# res_cont$variable <- rownames(res_cont)
#
# res_cont
# res_cat <- do.call(rbind,lapply(cat_cols, function(x){
#   res <- data.frame(table(metadata[,x],useNA="always")) #added it to deal with missing, can be changes
#   res$variable <- x
#   colnames(res)[1] <- "Category"
#   res
# }
# ))
# sum_table <- as.data.frame(summary(metadata))
#
# write.table(summary(metadata), 'data/summarry_metadata.txt',
#             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
# write.table(res_cont, 'data/summarry_metdata_cont_var.txt',
#             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
# write.table(res_cat, 'data/summarry_metdata_cont_cat.txt',
#             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)


X_sampling <- function(X_pool, snp_region,N_sample){
  ref_id_num <- nrow(X_pool)
  id = sample.int(ref_id_num, replace=TRUE, N_sample)
  X_sim <- X_pool[id, snp_region]
  return(X_sim)
}


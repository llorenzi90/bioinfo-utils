get_max_per_group_by_arrange <- function(tracking,
                                         group_var,
                                         arr_var,
                                         summ_var={{arr_var}},
                                         descending = T,
                                         outname="max_val"){

  tracking <- tracking %>% group_by({{group_var}})

  if(descending){
    tracking <- tracking %>% arrange(desc({{arr_var}}))
  } else tracking <- tracking %>% arrange({{arr_var}})

  tracking <- tracking %>% summarise(!!outname:= {{summ_var}}[1])

  return(tracking)
}

trait.per.sp <- function(trait_df, cover, trait_name){
  #### function to find the trait per species with a trait data
  # and the cover data.
  names(trait_df)[names(trait_df) == trait_name] = "trait"
  
  # cover weighted trait calcul 
  trait_weighted = (cover[,c("Be", "Cb", "Dc", "Fr", "Hl","Lp", "Pg",
                             "Pm", "Pt", "Ra", "Sp", "To")]/cover$total.cover) * trait_df$trait
  return(as.matrix(trait_weighted))
}

shape.df <- function(trait_df,cover,trait_name){
  #### Shape a trait_df to 2 trait_df that can be pass in the commu.traits.metrics function.
  #### this function will calculate the weighted trait values per species 
  #((cover_per_sp/total_cover)*trait)
  #### Trait_df and cover need to have identical nrow (each row corresponding to one plot)
  #### Trait name is the name of the trait variable in trait_df
  
  if (!any(names(trait_df) == "species_1" )){
    trait_df = cbind(trait_df, cover[, c("species_1", "species_2", "species_3")])
  }
  
  trait_df$species = apply( trait_df[,c("species_1", "species_2", "species_3")], 1, 
                      function(x) paste(x[!is.na(x)], collapse = "_"))
  
  # Weighted trait values per species 
  trait_per_sp = trait.per.sp(trait_df, cover, trait_name )
  
  
  trait_df = cbind(trait_df, trait_per_sp)
  
  trait_triplet = subset(trait_df, trait_df$type == "triplet")
  trait_mono = subset(trait_df, trait_df$type == "mono")
  
  ### get the trait value of the species in mono (without the alien species) 
  n=1
  trait = c()
  for (n in 1:nrow(trait_mono)){
    trait = append(trait, trait_mono[n, c(trait_mono$species_1[n]) ])
  }
  trait_mono$trait = trait
  
  
  trait_mono[,c("Be", "Cb", "Dc", "Fr", "Hl", "Lp", "Pg",
                "Pm", "Pt", "Ra", "Sp", "To", "species_2", "species_3")] <- NULL  
  return(list(trait_mono = trait_mono, trait_triplet = trait_triplet))
}

commu.traits.metrics <- function (trait_mono, trait_triplet, cover, trait_name_to,
                                  n_sp = 3, n_plot = 96 ){
  #### Function to get 4 trait metrics per plot : The sown trait value, delta
  # abundance shift, delta intraspecific shift and CWM
  
  ### trait_mono : wide df of the trait values for species in mono culture
  # it need to have cols of : nitrogen, species_1 and the trait values of the specie.
  
  ### trait_triplet : wide df of the trait values for species in triplet
  # it need to have cols of : nitrogen, species_1, species_2, species_3
  # plot and the 12 col of each species trait values ("Be", "Cb", "Dc", "Fr", "Hl",
  # "Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To")
  
  ### cover : wide df  with cover of each sp and the total cover. Cols :
  # sp("Be", "Cb", "Dc", "Fr", "Hl","Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To") 
  # and total.cover
  
  ### n : number of species in the plot (always 3 in our case)
  ### n_plot : number of plot (always 96 in our case)
  ### trait_name_to : the name that will be given to the trait
  
  trait_metrics = data.frame(nitrogen = trait_triplet$nitrogen, plot = trait_triplet$plot)
  
  ###### relative abundance ######
  relative_ab =c()
  ### cover_sp/cover_tot
  for (n in 1:nrow(cover)){
    relative_ab = rbind(relative_ab, cover[n,c("Be", "Cb", "Dc", "Fr", "Hl",
                                               "Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To")]/cover$total.cover[n])
  }
  relative_ab_triplet = as.matrix(relative_ab[1:96,])
  
  control_trait_mat = matrix(data= NA, nrow = n_plot, ncol = 12)
  colnames(control_trait_mat) = colnames(relative_ab_triplet)
  
  ###### sown trait per plot #####
  sown_trait = c()
  
  ### sum(sp_trait_mono)/n_sp_plot
  for (n in 1:nrow(trait_triplet)){
    sp_trait = c()
    
    sp_trait = append(sp_trait, subset(trait_mono,
                                       trait_mono$species_1 == trait_triplet$species_1[n] &
                                         trait_mono$nitrogen == trait_triplet$nitrogen[n])$trait)
    
    sp_trait = append(sp_trait, subset(trait_mono,
                                       trait_mono$species_1 == trait_triplet$species_2[n] &
                                         trait_mono$nitrogen == trait_triplet$nitrogen[n])$trait)
    
    sp_trait = append(sp_trait, subset(trait_mono,
                                       trait_mono$species_1 == trait_triplet$species_3[n] &
                                         trait_mono$nitrogen == trait_triplet$nitrogen[n])$trait)
    
    sp_names = trait_triplet[n,"species_1"]
    sp_names = append(sp_names, trait_triplet[n,"species_2"])
    sp_names = append(sp_names, trait_triplet[n,"species_3"])
    
    control_trait_mat[n,sp_names] = sp_trait
    
    sown_trait = append(sown_trait, sum(sp_trait, na.rm = T))/n_sp
  }
  
  trait_metrics$sown_trait = sown_trait
  
  
  ###### delta abondance shift per plot ######
  ### sum(Ab_rela*sp_trait_mono) - sown_trait
  
  mono_weight_mean= relative_ab_triplet * control_trait_mat 
  delta_ab_shift = rowSums(mono_weight_mean, na.rm = T) - sown_trait
  trait_metrics$delta_ab_shift = delta_ab_shift
  
  
  ###### delta intraspecific shift per plot ######
  ### sum(sp_trait_triplet)/n_sp_plot - sown_trait
  
  trait_values = as.matrix(trait_triplet[, c("Be", "Cb", "Dc", "Fr", "Hl",
                                             "Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To")])
  delta_intra_shift = rowSums(trait_values, na.rm = T) / n_sp - sown_trait
  trait_metrics$delta_intra_shift = delta_intra_shift
  
  
  ###### Comunity weighted mean ######
  ###sum(Ab_rela * sp_trait_triplet)
  CWM = rowSums(relative_ab_triplet * trait_values, na.rm = T)
  trait_metrics$CWM = CWM
  trait_metrics$species = trait_triplet$species
  
  trait_metrics$trait =  rep(trait_name_to, nrow(trait_metrics))
  return(trait_metrics)
}




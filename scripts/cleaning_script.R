# LDP Data Cleaning Assignment

pkgs <- c("tidyverse", "taxize", "myTAI")
lapply(pkgs, library,character.only = TRUE)

## Import Files ##

## make a list of files
myfiles = list.files(path="raw_data", pattern = "*.csv", full.names = TRUE)

# import all tables as separate data frames, remove file path and file extensions (.csv)
list2env(lapply(setNames(myfiles, make.names(gsub(".*1_", "", tools::file_path_sans_ext(myfiles)))), 
                read_csv), envir = .GlobalEnv)

# Number 1:Ensure that all higher level taxonomic information is complete for all 
# invertebrate morphospecies (morphospecies = bwg_names). Note, 
# morphospecies are identified at the lowest level by the researcher (often genus 
# or family). [suggested package: taxize]

# Here follows Amelia's attempts to put in upstream taxonomic information where it is missing...

# get only taxonomic-related columns with an id key of taxon_name (lowest resolution taxon)
traits_tax <- traits %>% 
  dplyr::select(4:17, 23:24) %>% 
  # rename columns so tax_name function can recognize levels
  rename(order = ord, suborder = subord) %>% 
  # only need to run the function on unique combinations of taxa
  distinct() %>%
  mutate(taxon_level = str_replace(taxon_level, pattern = "ord", replacement = "order"),
         taxon_level = str_replace(taxon_level, pattern = "species_name", replacement = "species"))

df <- traits_tax[1:5,]
x <- df[1,]

replace_missing_taxa <- function(df, db = "ncbi"){
  clean_data <- data.frame()
  
  for (row in 1:nrow(df)){
    x <- df[row,]
    col_names <- colnames(x) # extract column names for tax. level matching
    
    # extract column name = taxonomic level with the finest identified taxon level
    ided.tax <- which(grepl(x$taxon_name, x[,!(colnames(x) == "taxon_name")]))
    # some rows don't have matching names in "taxon_name" and table
    if(length(ided.tax) == 0L){
      ided.tax <- which(colnames(x) %in% x[,"taxon_level"])
    }
    # extract taxonomic levels above finest defined taxonomic level
    up_tax <- colnames(x)[1:ided.tax]
    
    # find taxonomic information based on lowest taxonomic level identified
    # get all entries
    all <- get_uid_(x$taxon_name)[[1]]
    
    # code to jump up a higher taxonomic level if finest taxonomic resolution does not have a match
    # jumps higher until a match in given database is found
    if(is.null(all)){
      i <- 1
      repeat{
        up_tax <- colnames(x)[1:(ided.tax-i)]
        # re-search for a higher tax level
        all <- get_uid_(x[,up_tax[length(up_tax)]])[[1]]
        i <- i + 1
        # not giving an additional warning message
        if(is.null(all) == F){
          class.comment <- paste("Cannot verify information below",up_tax[length(up_tax)])
          break
        }
      }
    }
    
    # if there is more than one entry, select based on whether it belongs to insects or metazoa
    if(length(all) > 1){
      class <- lapply(classification(all$uid, db = "ncbi"), as.data.frame)
      cond <- lapply(class, function(x) {
        if (any(x$name == "Insecta") == F) {
          if(any(x$name == "Metazoa") == F){
            return(FALSE)
          } else{return(TRUE)}
        } else {return(TRUE)}
      })
      if(length(unlist(cond) == T) > 1){
        # select genus entry, not subgenus
        sel.gen <- unlist(lapply(class, nrow))
        cor.entry <- class[[which(sel.gen < max(sel.gen))]]
      } else{
        cor.entry <- class[[which(unlist(cond) == T)]]
      }
      
      cor.entry <- cor.entry[cor.entry$rank %in% up_tax,]
    } else {
      class <- lapply(classification(all$uid, db = "ncbi"),as.data.frame)[[1]]
      cor.entry <- class[class$rank %in% up_tax,]
    }
    
    # create temporary data frame to compare original and retrieved db taxonomy names
    t <- data.frame(rank = up_tax,
                    original = as_vector(x[,colnames(x) %in% up_tax])) %>%
      left_join(cor.entry, by = "rank") %>% 
      dplyr::select(-id, db = name)
    rownames(t) <- t$rank
    # create column for names that will be exported, we keep mainly the original taxonomy names but
    # 1. fill NAs with data base
    # 2. If there are mismatches, we keep the data base
    
    # first copy the original entries into the final column
    t$final <- t$original
    # if else statement to only fill rows that are initially missing
    # return message if a NA has been filled
    if(any(is.na(t$original))){
      t[is.na(t$original),"final"] <- t[is.na(t$original),"db"]
      message(paste0("NAs have been filled with '", db,"' data base entry"))
    }
    
    if(any(t$final != t$db)){
      # if else statement to find taxonomic levels that
      # show mismatches between original and db entry
      # db entry will overwrite original name and warning message will be returned
      # exclude if there is an entry missing in the data base (e.g. "domain" does not exist in NCBI)
      mis <- t[!is.na(t$db),]
      mismatch <- rownames(mis[!(mis$final %in% mis$db),])
      
      t[rownames(t) %in% mismatch,"final"] <- t[rownames(t) %in% mismatch,"db"]
      mismatch.comment <- paste0("DB mismatch in ", mismatch,". Replacing ", t[mismatch,"original"],".")
      message(paste0("Mismatch in '", paste(mismatch, collapse = " "), "', overwriting with '", db, "' data base entry"))
    }
    
    out <- pivot_longer(x, cols = everything(), values_to = "tax")
    out[out$name %in% rownames(t),"tax"] <- t$final
    out <- pivot_wider(out, names_from = name, values_from = tax)
    
    # add comments on finest taxonomic level not found in DB, and give info which taxonomic levels have been verified
    if(exists("class.comment")){
      out$DB.comment <- paste(class.comment, collapse = " ")
      rm(class.comment)
    } else {
      out$DB.comment <- NA
    }
    
    # add comments on mismatches that have been replaced with DB info
    if(exists("mismatch.comment")){
      out$mismatch.comment <- paste(mismatch.comment, collapse = " ")
      rm(mismatch.comment)
    } else {
      out$mismatch.comment <- NA
    }
    
    # assemble all the clean data in a single dataframe
    clean_data <- clean_data %>% 
      rbind(out)
  }
  return(clean_data)
}

# now run the function on the taxonomic information for traits
clean_traits <- replace_missing_taxa(traits_tax)
clean_traits <- replace_missing_taxa(df)

clean_traits_df <- unique(as.data.frame(clean_traits))

write_csv(clean_traits_df, "./clean_data/replaced_taxa.csv")

# some taxa have multiple UID choices, so select the following as the tax_name function runs:
# 2 for Culicoides, 2 for Forcipomyia, 2 for Chironomus, 2 for Polypedilum,
# 1 for Polyphaga, 3 Anopheles, 2 for Culex, 2 for Haemogogus, 2 for Aedes, 1 for Trentepohlia,
# 2 for Tipula, 1 for Ormosia, 2 for Limonia, 2 for Oligochaeta, 2 for Forcipomyia

traits_complete <- traits %>% 
  # take out incomplete taxonomic data
  select(-(4:17)) %>% 
  # and join with the completed taxonomic data
  full_join(clean_traits) %>% 
  # the tax_name function seems to replace suborder with infraorder sometimes, 
  # so retain only the first instance of replacement for each unique taxon
  distinct(bwg_name, .keep_all = TRUE)

# taxonomic data are now fully integrated into dataframe

write_csv(traits_complete,"./clean_data/traits_complete.csv")



# LDP Data Cleaning Assignment

pkgs <- c("lubridate","rgdal","tidyverse", "taxize", "myTAI")
lapply(pkgs, library,character.only = TRUE)

## Import Files ##

## make a list of files
myfiles = list.files(path="./BWG_database", pattern = "*.csv", full.names = TRUE)

# import all tables as separate data frames, remove file path and file extensions (.csv)
list2env(lapply(setNames(myfiles, make.names(gsub(".*1_", "", tools::file_path_sans_ext(myfiles)))), 
                read_csv), envir = .GlobalEnv)

# Number 1:Ensure that all higher level taxonomic information is complete for all 
# invertebrate morphospecies (morphospecies = bwg_names). Note, 
# morphospecies are identified at the lowest level by the researcher (often genus 
# or family). [suggested package: taxize]

traits_tax <- traits %>% 
  dplyr::select(domain:subspecies, taxon_level) %>% # need to keep morphospecies ID
  # rename columns so tax_name function can recognize levels
  rename(order = ord, suborder = subord) %>%
  distinct() %>%
  mutate(taxon_level = str_replace(taxon_level, pattern = "ord", replacement = "order"),
         taxon_level = str_replace(taxon_level, pattern = "species_name", replacement = "species"))

# selection of a few rows for testing/making the function
traits_tax <- traits_tax[c(1,7,8),]

# make taxonomic levels a factor with levels to define what is "upstream" =  higher taxonomic level
traits_tax$taxon_level<-factor(traits_tax$taxon_level,levels=c("domain","kingdom","phylum","subphylum",
                                                               "class","subclass","order","suborder","family",
                                                               "subfamily","tribe","genus","species","subspecies"))

# split dataframe rows into lists to keep data frame structure within function
trait.list <- split(traits_tax, seq(nrow(traits_tax)))

id_match_upstream_levs <- function(x, db = "ncbi"){
  # define taxonomic levels as factor vector
  tax_levs <- factor(c("domain","kingdom","phylum","subphylum",
                       "class","subclass","order","suborder","family",
                       "subfamily","tribe","genus","species","subspecies"),
                     levels=c("domain","kingdom","phylum","subphylum",
                              "class","subclass","order","suborder","family",
                              "subfamily","tribe","genus","species","subspecies"))
  # extract taxonomic levels above finest defined taxonomic level
  up_tax <- as.character(tax_levs[as.integer(tax_levs) < as.integer(x["taxon_level"])])
  # retrieve taxonomic names above finest defined taxonomic level from specified db
  db_class <- tax_name(as.character(x[as.integer(x["taxon_level"])]),
           get = up_tax, db = db)
  
  # code to jump up a higher taxonomic level if finest taxonomic resolution does not have a match
  # jumps higher until a match in given database is found
  if(all(is.na(db_class[,-c(1,2)]))){
    i <- 1
    repeat{
      up_tax <- as.character(tax_levs[as.integer(tax_levs) < as.integer(x["taxon_level"])-i])
      db_class <- tax_name(as.character(x[as.integer(x["taxon_level"])-i]),
                           get = up_tax, db = db)
      i <- i + 1
      if(all(is.na(db_class[,-c(1,2)])) == F){
        break
      }
    }
  }
  
  # create temporary data frame to compare original and retrieved db taxonomy names
  t <- data.frame(original = as_vector(x[,colnames(x) %in% colnames(db_class)]),
                  db = as_vector(db_class[,-c(1:2)]))
  # create column for names that will be exported, we keep mainly the original taxonomy names but
  # 1. fill NAs with data base
  # 2. If there are mismatches, we keep the data base
  t$final <- t$original
  # if else statement to only fill rows that have NA
  # return message if a NA has been filled
  if(any(is.na(t$original))){
    t[is.na(t$original),"final"] <- t[is.na(t$original),"db"]
    message(paste0("NA has been filled with '", db,"' data base entry"))
  }
  
  if(any(t$final != t$db)){
    # if else statement to find taxonomic levels that
    # show mismatches between original and db entry
    # db entry will overwrite original name and warning message will be returned
    mismatch <- rownames(t[!(t$final %in% t$db),])
    t[rownames(t) %in% mismatch,"final"] <- t[rownames(t) %in% mismatch,"db"]
    comment <- paste0("DB mismatch in ", mismatch,". Replacing ", t[mismatch,"original"],".")
    message(paste0("Mismatch in '", paste(mismatch, collapse = " "), "', overwriting with '", db, "' data base entry"))
  }
  
  out <- pivot_longer(x, cols = everything(), values_to = "tax")
  out[out$name %in% rownames(t),"tax"] <- t$final
  out <- pivot_wider(out, names_from = name, values_from = tax)
  if(length(comment) != 0){
    out$comment <- paste(comment, collapse = " ")
  } else {
    out$comment <- NA
  }
  return(out)
}

# apply function, bind list back to data frame
correct_tax <- bind_rows(lapply(trait.list, id_match_upstream_levs))

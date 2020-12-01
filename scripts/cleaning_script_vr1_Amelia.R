# LDP Data Cleaning Assignment

pkgs <- c("tidyverse", "taxize", "myTAI")
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

# get only taxonomic-related columns with an id key of taxon_name (lowest resolution taxon)
traits_tax <- traits %>% 
  dplyr::select(4:17, 24) %>% 
  # rename columns so tax_name function can recognize levels
  rename(order = ord, suborder = subord) %>% 
  # only need to run the function on unique combinations of taxa
  distinct()

replace_missing_taxa <- function(df){
  clean_data <- data.frame()
  levels_replaced <- c()
  for (row in 1:nrow(df)){
    x <- df[row,]
    col_names <- colnames(x)
    levels_replaced_row <- c()
    #names of columns where taxonomic data are missing
    missing <- data.frame(numbers = which(is.na(x)), names = col_names[which(is.na(x))])
    # find taxonomic information based on lowest taxonomic level identified
    tax_names <- tax_name(x$taxon_name, get = missing$names, db = "ncbi")
    tax_col_names <- colnames(tax_names)
    for (missing_col in 1:nrow(missing)){
      # for every column where a value was missing...
      for (filled_col in 1:length(tax_col_names)){
        # check the completed taxonomic information from the tax_name function
        if (missing$names[missing_col] == tax_col_names[filled_col]){
          # and then fill in the missing data with the completed information if it exists
          if (is.na(tax_names[1, filled_col]) == 0){
            col_ref <- missing$numbers[missing_col]
            x[1, col_ref] <- tax_names[1,filled_col]
            levels_replaced_row <- paste(levels_replaced_row, missing$names[missing_col])
          }
        }
      }
    }
    # assemble all the clean data in a single dataframe
    clean_data <- clean_data %>% 
      rbind(x)
    levels_replaced[row] <- paste(levels_replaced_row, collapse = "")
  }
  # and amalgamate all of the data with the column indicating what has been replaced
  clean_data <- clean_data %>% 
    cbind(levels_replaced)
}

levels_replaced_row <- c()
# now run the function on the taxonomic information for traits
clean_traits <- replace_missing_taxa(traits_tax)

clean_traits_df <- unique(as.data.frame(clean_traits))

write_csv(clean_traits_df, "./clean_data/replaced_taxa.csv")

# some taxa have multiple UID choices, so select the following as the tax_name function runs:
# 2 for Culicoides, 2 for Forcipomyia, 2 for Chironomus, 2 for Polypedilum,
# 1 for Polyphaga, 2 for Culex, 2 for Haemogogus, 2 for Aedes, 1 for Trentepohlia,
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
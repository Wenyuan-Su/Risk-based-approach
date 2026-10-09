# ============================================================
# Load concentration data
# ============================================================

# The input file contains:
# ID column: sample names, TEF, reference mixture, and candidate samples
# Chemical columns: concentrations of individual compounds

df <- read.table(
  "MIXN.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# Extract chemical concentration matrix
# ============================================================

# Exclude the ID column and retain chemical concentration columns
chem_cols <- 2:ncol(df)


# Extract toxicity equivalency factors (TEFs)
# The row with ID = "TEF" contains TEF values for each compound

tef_vec <- as.numeric(
  df[df$ID == "TEF", chem_cols]
)


# Extract reference mixture
# The reference mixture was defined based on the NWU sample

ref_row <- df[df$ID == "NWU", ]

ref_props <- as.numeric(
  ref_row[, chem_cols]
)


# Extract candidate samples
# Remove TEF and reference mixture rows

cand_df <- df[
  !df$ID %in% c("TEF", "NWU"),
]


# ============================================================
# Calculate toxicity-based weights
# ============================================================

# Weight formula:
# w_j = c × (TEF_j / sum(TEF_j))

c_count <- length(tef_vec)

weights <- c_count *
  (tef_vec / sum(tef_vec))


# ============================================================
# Define SMACH distance function
# ============================================================

calculate_SMACH_distance <- function(
    candidate_props,
    reference_props,
    weights,
    Tr
){
  
  candidate <- as.numeric(candidate_props)
  
  reference <- as.numeric(reference_props)
  
  
  # Calculate weighted compositional difference
  
  weighted_distance <- sum(
    weights *
      (reference - candidate)^2
  )
  
  
  # Multiply by BMD of reference mixture
  
  distance <- Tr *
    sqrt(weighted_distance)
  
  
  return(distance)
}



# ============================================================
# Set BMD of reference mixture
# ============================================================

# Tr represents the BMD of the reference mixture

reference_BMD <- 12.31



# ============================================================
# Calculate distances for all candidate samples
# ============================================================

distances <- apply(
  cand_df[, chem_cols],
  1,
  function(x){
    
    calculate_SMACH_distance(
      candidate_props = x,
      reference_props = ref_props,
      weights = weights,
      Tr = reference_BMD
    )
    
  }
)



# ============================================================
# Export results
# ============================================================

result_table <- data.frame(
  ID = cand_df$ID,
  SMACH_distance = distances
)


write.table(
  result_table,
  "4_MIXN_SMACH_distance_results.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
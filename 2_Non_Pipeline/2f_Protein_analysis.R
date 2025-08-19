# Objective: ANalysise the differences in protein abundances that I found so far.


library(tidyverse)
library(gridExtra)


protein_counts_loc <- "/home/darrian/Documents/RNA_seq_fescue/interpro_results/proteinCount.csv"


###############################
# loading data
###############################
protein_counts <- read.table(protein_counts_loc, header = TRUE, sep = "," )
protein_counts <- protein_counts %>% 
  filter(Unique_Value != "-")
protein_counts <- protein_counts %>% 
  filter(Unique_Value != "consensus disorder prediction")

##################################
# Function to make barplot
##################################

# Sort the data frame by 'specific_column' in descending order and select the top 10 rows
create_bar_graph <- function(data, column_name, Graph_Title = "Top 10 Proteins Found", color = "steelblue", endo_color = "firebrick1") {
  
  # Select and format the top 10 proteins
  top_10 <- data %>%
    arrange(desc(.data[[column_name]])) %>%
    slice(1:10) %>%
    mutate(Unique_Value = ifelse(nchar(Unique_Value) > 30, 
                                 paste0(substr(Unique_Value, 1, 30), "..."), 
                                 Unique_Value))
  
  # Create a horizontal bar graph with ggplot
  p <- ggplot(top_10, aes(y = reorder(Unique_Value, .data[[column_name]]), 
                          x = .data[[column_name]])) +
    geom_bar(stat = "identity", fill = color) +
    labs(title = Graph_Title,
         y = "Proteins",
         x = "Sequences Detected") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 60)) +
    theme(
      plot.title = element_text(colour = endo_color, size = 20),     # Increased title font size
      axis.title.x = element_text(size = 14),                        # Increased x-axis title font size
      axis.title.y = element_text(size = 14),                        # Increased y-axis title font size
      axis.text.x = element_text(size = 12),                         # Increased x-axis label size
      axis.text.y = element_text(size = 12)                          # Increased y-axis label size
    )
  
  return(p)
}


##################################
# Make barplots
##################################

colnames(protein_counts)

# HP vs Control
a1 <- create_bar_graph(data = protein_counts, column_name = "EndoPos_HpxControl_significantly_upregulated", Graph_Title = "E+: Heat & Percipitation vs Control: Upregulated", color = "firebrick1", endo_color = "black")
b1 <- create_bar_graph(data = protein_counts, column_name = "EndoPos_HpxControl_significantly_downregulated", Graph_Title = "E+: Heat & Percipitation vs Control: Downregulated", color = "steelblue", endo_color = "black")


a2 <- create_bar_graph(data = protein_counts, column_name = "EndoNeg_HpxControl_significantly_upregulated", Graph_Title = "E-: Heat & Percipitation vs Control: Upregulated", color = "firebrick1", endo_color = "black")
b2 <- create_bar_graph(data = protein_counts, column_name = "EndoNeg_HpxControl_significantly_downregulated", Graph_Title = "E-: Heat & Percipitation vs Control: Downregulated", color = "steelblue", endo_color = "grey40")

grid.arrange(a1, a2, ncol = 1)
grid.arrange(b1, b2, ncol = 1)


# Heat x control
c1 <- create_bar_graph(data = protein_counts, column_name = "EndoPos_HeatxControl_significantly_upregulated", Graph_Title = "E+: Heat vs Control: Upregulated", color = "firebrick1", endo_color = "black")
c2 <- create_bar_graph(data = protein_counts, column_name = "EndoPos_HeatxControl_significantly_downregulated", Graph_Title = "E+: Heat vs Control: Downregulated", color = "steelblue", endo_color = "grey40")

grid.arrange(c1, c2, ncol = 1)


# HP x Heat
create_bar_graph(data = protein_counts, column_name = "EndoNeg_HpxHeat_significantly_upregulated", Graph_Title = "E-: Heat vs Heat with Percipitation: Upregulated", color = "firebrick1", endo_color = "black")
create_bar_graph(data = protein_counts, column_name = "EndoPos_HpxHeat_significantly_upregulated", Graph_Title = "E+: Heat vs Heat with Percipitation: Upregulated", color = "firebrick1", endo_color = "black")

# Heat: E- x E+
create_bar_graph(data = protein_counts, column_name = "Heat_NegxPos_significantly_upregulated", Graph_Title = "Heat: E+ x E-: Upregulated", color = "firebrick1", endo_color = "black")
create_bar_graph(data = protein_counts, column_name = "Heat_NegxPos_significantly_downregulated", Graph_Title = "Heat: E+ x E-: Downregulated", color = "steelblue", endo_color = "black")








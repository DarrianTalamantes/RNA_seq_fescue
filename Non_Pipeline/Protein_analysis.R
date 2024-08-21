# Objective: ANalysise the differences in protein abundances that I found so far.


library(tidyverse)


protein_counts_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/RNA_seq_fescue/interpro_results/proteinCount.csv"


###############################
# loading data
###############################
protein_counts <- read.table(protein_counts_loc, header = TRUE, sep = "," )


##################################
# Function to make barplot
##################################

# Sort the data frame by 'specific_column' in descending order and select the top 10 rows

create_bar_graph <- function(data, column_name, Graph_Title = "Top 10 Proteins Found", color = steelblue, endo_color = firebrick1) {
  
top_10 <- data %>%
  arrange(desc(.data[[column_name]])) %>%
  slice(3:12)

# Create a bar graph with ggplot
p <- ggplot(top_10, aes(x = reorder(Unique_Value, -.data[[column_name]]), 
                   y = .data[[column_name]])) +
  geom_bar(stat = "identity", fill = color) +
  labs(title = Graph_Title,
       x = "Proteins",
       y = "Sequences Detected") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 60)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(colour = endo_color))  # Rotate x-axis labels for better readability
return(p)
}

##################################
# Make barplots
##################################

colnames(protein_counts)

# HP vs Control
create_bar_graph(data = protein_counts, column_name = "EndoPos_HpxControl_significantly_upregulated", Graph_Title = "Endo Positive: HeatxPercipitation vs Control: Upregulated", color = "firebrick1", endo_color = "black")
create_bar_graph(data = protein_counts, column_name = "EndoPos_HpxControl_significantly_downregulated", Graph_Title = "Endo Positive: HeatxPercipitation vs Control: Downregulated", color = "steelblue", endo_color = "black")


create_bar_graph(data = protein_counts, column_name = "EndoNeg_HpxControl_significantly_upregulated", Graph_Title = "Endo Negative: HeatxPercipitation vs Control: Upregulated", color = "firebrick1", endo_color = "grey40")
create_bar_graph(data = protein_counts, column_name = "EndoNeg_HpxControl_significantly_downregulated", Graph_Title = "Endo Negative: HeatxPercipitation vs Control: Downregulated", color = "steelblue", endo_color = "grey40")

# Heat x control
create_bar_graph(data = protein_counts, column_name = "EndoPos_HeatxControl_significantly_upregulated", Graph_Title = "Endo Positive: Heat vs Control: Upregulated", color = "firebrick1", endo_color = "black")
create_bar_graph(data = protein_counts, column_name = "EndoPos_HeatxControl_significantly_downregulated", Graph_Title = "Endo Positive: Heat vs Control: Downregulated", color = "steelblue", endo_color = "black")







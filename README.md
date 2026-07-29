# DECODE-2026-Team-Hubmap
# sample_info summary
library(ggplot2)
ggplot(sample_info, aes(x = condition, fill = condition)) +
  geom_bar(width=0.6) +
  scale_fill_manual(values = c(
    "Blue",   
    "Orange", 
    "Green",  
    "Red",  
    "Purple",
    "yellow"
  )) +
  labs(title = "Number of Samples per Group",
       x = "Group",
       y = "Number of Samples") +
  theme_classic(base_size = 5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )

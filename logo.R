# Install required packages (only run this once)
install.packages("hexSticker")
install.packages("ggplot2")

# Load necessary libraries
library(hexSticker)
library(ggplot2)

# Create a simple plot for the logo
p <- ggplot(data.frame(x = c(-1, 1)), aes(x = x)) +
  geom_point() +
  theme_void()

# Generate the hex sticker
sticker(
  subplot = p,
  package = "GLMGWAS",
  p_size = 20,
  s_x = 1, s_y = 0.8, s_width = 1.3,
  filename = "logo.png"
)

# The file "logo.png" will be saved in your working directory


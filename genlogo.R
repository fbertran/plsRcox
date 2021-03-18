#install.packages("hexSticker")
library(hexSticker)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("emojifont")
library(emojifont)
# fa4
#r <- ggplot() + geom_fontawesome("fa-calculator", color='red', size = 20) + theme_void()

# fa5
#install.packages("fontawesome")
library(fontawesome)
aaa=fa(name = "fab fa-contao", fill = "orange")

print_fontawesome<- function (x, view = interactive(), ...) 
  {
    dots <- list(...)
    html <- paste(x, collapse = "\n")
    c("<!DOCTYPE html>", "<html>", "<head>", "<meta charset=\"utf-8\">", 
      "</head>", "<body>", html, "</body>", "</html>") %>% 
      paste(collapse = "\n") %>% htmltools::HTML() %>% htmltools::html_print()
    return(html)
  }

library(rsvg)
cat(print_fontawesome(aaa, view=FALSE),file="man/figures/fa-contao.svg")
SVGres <- rsvg_svg("man/figures/fa-contao.svg","man/figures/Cairo-fa-contao.svg")
SVGres <- rsvg_svg("man/figures/fa-contao.svg")
PDFres <- rsvg_pdf("man/figures/fa-contao.svg","man/figures/Cairo-fa-contao.pdf")
PDFres <- rsvg_pdf("man/figures/fa-contao.svg")

library(png)
library(grid)
library(cowplot)
library(magick)
theme_set(theme_cowplot())
r = ggdraw() +
  draw_image("man/figures/Cairo-fa-contao.pdf", scale = .55)


sticker(
  r,
  package = "plsRcox",
  p_size = 8,
  s_x = .98,
  s_y = 0.7,
  s_width = 1.7,
  s_height = 1.3,
  p_x = 1,
  p_y = 1.3,
  url = "https://cran.r-project.org/package=plsRcox",
  u_color = "white",
  u_size = 1.1,
  h_fill = "black",
  h_color = "grey",
  filename = "man/figures/logo.png"
)

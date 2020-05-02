library(profvis)
library(devtools)
library(htmlwidgets)
load_all()

prof = profvis({
  test()
})

saveWidget(prof, "prof.html")

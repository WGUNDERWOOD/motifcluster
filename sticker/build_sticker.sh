pdflatex sticker.tex &&
    pdflatex sticker.tex &&
    pdflatex sticker.tex &&
    convert -density 500 sticker.pdf sticker.png &&
    hexsticker sticker.png -o hex_sticker.png --border-size 120 --border-color "#4c4a5b" --supersample 2

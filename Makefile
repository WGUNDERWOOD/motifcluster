all: python R performance sticker

.PHONY: python R performance sticker clean

python:
	@echo -e "\e[0;35m\033[1mMaking python package...\e[0;30m\033[0m"
	@cd python/ && make all

R:
	@echo -e "\e[0;35m\033[1mMaking R package...\e[0;30m\033[0m"
	@cd R/ && make all

performance:
	@echo -e "\e[0;35m\033[1mMaking performance report...\e[0;30m\033[0m"
	@cd performance && mkdir -p results/ plots/
	@cd performance && python -m cProfile -o profile.pstats performance_test.py
	@cd performance && gprof2dot -f pstats profile.pstats | dot -Tpng -o profile.png
	@cd performance && rm profile.pstats
	@cd performance && Rscript performance_test.R
	@cd performance && python performance_plot.py
	@cd performance && latexmk -pdf -quiet -rc-report- performance.tex

sticker:
	@echo -e "\e[0;35m\033[1mMaking sticker...\e[0;30m\033[0m"
	@cd sticker && latexmk -pdf -quiet sticker.tex
	@cd sticker && convert -density 500 sticker.pdf sticker.png
	@cd sticker && hexsticker sticker.png -o hex_sticker.png \
		--border-size 120 --border-color "#6a668A" --supersample 2
	@cd sticker && convert hex_sticker.png -resize 150 hex_sticker_small.png
	@cd sticker && optipng -o3 hex_sticker.png
	@cd sticker && optipng -o3 hex_sticker_small.png

clean:
	@echo -e "\e[0;35m\033[1mCleaning up...\e[0;30m\033[0m"
	@cd python && make clean
	@cd R && make clean
	@cd performance && texclean && rm -rf plots/ results/ profile.png
	@rm -rf .pytest_cache/
	@cd sticker && texclean && rm -f sticker.pdf sticker.png

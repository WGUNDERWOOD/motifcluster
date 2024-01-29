all: python R julia performance sticker todo

.PHONY: python R julia performance sticker todo clean

python:
	@echo -e "\e[0;35m\033[1mMaking python package...\e[0;30m\033[0m"
	@cd python/ && make all

R:
	@echo -e "\e[0;35m\033[1mMaking R package...\e[0;30m\033[0m"
	@cd R/ && make all

julia:
	@echo -e "\e[0;35m\033[1mMaking Julia package...\e[0;30m\033[0m"
	@cd julia/ && make all

performance:
	@echo -e "\e[0;35m\033[1mMaking performance report...\e[0;30m\033[0m"
	@cd performance && mkdir -p results/ plots/
	@cd performance && python -m cProfile -o profile.pstats performance_test.py
	@cd performance && gprof2dot -f pstats profile.pstats | dot -Tpng -o profile.png
	@cd performance && rm -f profile.pstats
	@cd performance && Rscript performance_test.R
	@cd performance && julia --project="../julia/MotifCluster.jl" --color=yes performance_test.jl
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

todo:
	@echo -e "\e[0;35m\033[1mLooking for todo items...\e[0;30m\033[0m"
	@! rg -g "!Makefile" TODO

clean:
	@echo -e "\e[0;35m\033[1mCleaning up...\e[0;30m\033[0m"
	@cd python && make clean
	@cd R && make clean
	@cd performance && latexmk -c && rm -rf plots/ results/ profile.png
	@rm -rf .pytest_cache/
	@cd sticker && latexmk -c && rm -f sticker.pdf sticker.png

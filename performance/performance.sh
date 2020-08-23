# performance.sh

mkdir -p results/ plots/

python -m cProfile -o profile.pstats performance_test.py &&
    gprof2dot -f pstats profile.pstats | dot -Tpng -o profile.png && rm profile.pstats

Rscript performance_test.R

python performance_plot.py

pdflatex performance.tex

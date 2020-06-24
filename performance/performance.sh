# performance.sh

python -m cProfile -o profile.pstats performance_test.py &&
    gprof2dot -f pstats profile.pstats | dot -Tpng -o profile.png && rm profile.pstats &&
    python performance_plot.py

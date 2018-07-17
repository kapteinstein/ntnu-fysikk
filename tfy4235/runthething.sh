rm -rf results_*

python -c "for i in range(1, 13): print(0.01*i)" | xargs -n1 -P1 python3 main2.py --task 1.5 -h
# python3 main2.py
python3 plotter.py -2d -3d --files results_*


# mv bildemothafokka.pdf figures/bildemothafokka.pdf
# open figures/bildemothafokka.pdf

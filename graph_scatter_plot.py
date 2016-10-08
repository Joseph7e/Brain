#!/usr/bin/python3
import sys, pygal
my_list = []


with open(sys.argv[1], 'r') as input_data:
    for line in input_data:
        l,k,c,r = line.strip().split('\t')
        l = float(l); k = float(k); c = float(c); r = float(r)
        current = (c,r)
        my_list.append(current)

xy_chart = pygal.XY(stroke=False)
xy_chart.title = 'GC_content vs Coverage'
xy_chart.add('GCvsCoverage', my_list)
# xy_chart.add('B', [(.1, .15), (.12, .23), (.4, .3), (.6, .4), (.21, .21), (.5, .3), (.6, .8), (.7, .8)])
# xy_chart.add('C', [(.05, .01), (.13, .02), (1.5, 1.7), (1.52, 1.6), (1.8, 1.63), (1.5, 1.82), (1.7, 1.23), (2.1, 2.23), (2.3, 1.98)])
xy_chart.render_to_file('my_scatter_plot.svg')
#!/usr/bin/python3

#input --> tab seperated #length(x) #set1(y1) #set2(y2)


import sys, pygal #sorted_coverage_table.txt

length = []
kmer_cov = []
coverage = []
reads = []

# x_label = map(str, range(1, 50))
# print (x_label)
# print (type(x_label))
# sys.exit()

count = 1

one = 0
two = 0
three = 0
four = 0
five = 0
six = 0
seven = 0
eight = 0
nine = 0
ten = 0


with open(sys.argv[1], 'r') as input_data:
    for line in input_data:
        l,k,c,r = line.strip().split('\t')
        l = float(l); k = float(k); c = float(c); r = float(r)
        # if l > 1000:
        #     count += 1
        #     length.append(l); kmer_cov.append(k); coverage.append(c); reads.append(r)
        if l >= 50000:
            one += 1
        if l < 50000 and l >= 25000:
            two += 1
        if l <25000 and l >= 15000:
            three += 1
        if l < 15000 and l >= 10000:
            four += 1
        if l < 10000 and l >= 5000:
            five += 1
        if l < 5000 and l >= 2500:
            six += 1
        if l < 2500 and l >= 1000:
            seven += 1
        if l < 1000 and l >= 500:
            eight +=1
        if l < 500 and l >= 250:
            nine += 1
        if l < 250:
            ten += 1


print (one, two, three, four, five, six, seven, eight, nine, ten)

sys.exit()

line_chart = pygal.Line()
# line_chart = pygal.XY()
line_chart.title = 'Length of contig vs coverage'
line_chart.x_labels = map(str, range(1, count))
line_chart.add('kmer_coverage', kmer_cov)
line_chart.add('per base coverage', coverage)
# line_chart.add('IE',      [85.8, 84.6, 84.7, 74.5,   66, 58.6, 54.7, 44.8, 36.2, 26.6, 20.1])
# line_chart.add('Others',  [14.2, 15.4, 15.3,  8.9,    9, 10.4,  8.9,  5.8,  6.7,  6.8,  7.5])
line_chart.render_to_file('line_plot_top_7000.svg')
#line_chart.render_to_png('line_plot_range_check_top_50000.png')
#
print (count)

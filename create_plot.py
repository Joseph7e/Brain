#!/usr/bin/python3


import matplotlib.pyplot as pyplot

#node detail lists
av_var_list = []
low_var_list = []
copy_number_list = []
length_low_list = []
length_high_list = []

with open('gene_stats/variation_table.xls', 'r') as v:
    for line in v.readlines():
        if '#' in line or '%' in line:
            continue
        else:
            sample, cp, ll, lh, av, hv = line.split("\t")
            if av == '0' or hv == '0':
                continue
            else:
                av_var_list.append(eval(av))
                low_var_list.append(eval(hv))
                copy_number_list.append(eval(cp))
                length_high_list.append(eval(lh))
                length_low_list.append(eval(ll))


def plotdata(x_list, y_list, title, subtitle, x_label, y_label, out_name):
    """ takes a user defined set of data and creates a jpg graph"""
    #x_data = x_list
    #y_data = y_list
    figure = pyplot.figure()
    figure.suptitle("16S", fontsize=20, fontweight='bold')
    figure.subplots_adjust(top=0.85)
    ax = figure.add_subplot(111)
    ax.set_title(subtitle)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    #pyplot.axis([0, 60, 95, 100.5])
    pyplot.scatter(x_list,y_list, color="blue", marker=".") #may want to use marker="," if . is too big
    pyplot.savefig(out_name + "_scatterplot.jpg")


print("Creating a scatter plot of GC content vs Coverage") #Read print for notes, running plotdata for all data

plotdata(copy_number_list,low_var_list, "Bacterial 16S", "copy_number VS lowID", "Copy Number", "% Identity", "selected_lowIDvscopy", )

print("Creating a scatter plot of Length vs GC content")
plotdata(copy_number_list,av_var_list, "Bacterial 16S", "copy_number VS Average %ID", "Copy Number", "% Identity", "selected_avgIDvscopy", )


print ('average copy number = ', sum(copy_number_list)/len(copy_number_list))

print (len(copy_number_list))

nl = sorted(copy_number_list)
na = sorted(av_var_list)
nn = sorted(low_var_list)

print ("lowest copy number = {}, highest copy number = {}".format(nl[0], nl[-1]))
print ('lowest avg variation = {}, highest ag variation = {}'.format(na[0], na[-1]))
print ('lowest low variation = {}, highest low variation = {}'.format(nn[0],nn[-1]))

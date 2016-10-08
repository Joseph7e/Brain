#!/usr/bin/python3


import networkx as nx
import matplotlib.pyplot as plt

input_file = '/home/genome/joseph7e/SSRs/pairwise_stuff/ssr_matches_table'

my_list = []
my_edges = []
# my_nodes = []

# atlantic_list = []
# pacific_list = []
# shared_list = []
atlantic_nodes = []
pacific_nodes = []
shared_nodes = []
atlantic_list = []
pacific_list = []
shared_list = []

atlantic = ['A02', 'A12', 'A3', 'A6', 'A7', 'B02', 'B03', 'B1']
pacific = ['A4', 'B7', 'E9', 'G2', 'H2']

atlantic_count = 0
shared_count = 0
pacific_count = 0

with open(input_file, 'r') as i:
    for line in i.readlines():
        atlantic_flag = False
        pacific_flag = False
        line = line.rstrip()
        count, group = line.split("\t")
        list = group.split("_")
        tl = set(list)
        my_list.append(tl)
        tup = (group.replace("_", "+"), dict( size=int(count) ))
        group_count = group.count("_")
        for t in list:
            if t in atlantic:
                atlantic_flag = True
            if t in pacific:
                pacific_flag = True
        if atlantic_flag and pacific_flag:
            shared_nodes.append(tup)
            shared_list.append(group.replace("_", "+"))
            if group_count > 1:
                shared_count += int(count)
        if atlantic_flag and pacific_flag == False:
            #atlantic_list.append(tl)
            atlantic_nodes.append(tup)
            atlantic_list.append(group.replace("_", "+"))
            if group_count > 1:
                atlantic_count += int(count)
        if pacific_flag and atlantic_flag == False:
            #pacific_list.append(tl)
            pacific_list.append(group.replace("_", "+"))
            pacific_nodes.append(tup)
            if group_count > 1:
                pacific_count += int(count)
        my_list.append(tl)
        # dict( size=int(d[1]) )
        #             ]) for d in data]
        #

print (shared_count, atlantic_count, pacific_count)

for i in range(len(my_list)):
    for j in range(i + 1, len(my_list)):
        if my_list[i].issubset(my_list[j]) or my_list[j].issubset(my_list[i]):
            my_edges.append((my_list[i],my_list[j]))

my_edges = [('+'.join(sorted(e[0])), '+'.join(sorted(e[1]))) for e in my_edges]


if __name__ == '__main__':
    scale_factor = 4
    G = nx.Graph()
    #nodes = my_nodes
    pacific_sizes = [ (n[1])['size']*scale_factor
                    for n in pacific_nodes ]
    atlantic_sizes = [ (n[1])['size']*scale_factor
                    for n in atlantic_nodes ]
    shared_sizes = [ (n[1])['size']*scale_factor
                    for n in shared_nodes ]
    #edges = my_edges
    G.add_edges_from( my_edges )
    nx.draw_networkx_nodes(G,pos=nx.spring_layout(G),
                     node_list = shared_list,
                     node_color = 'g',
                     node_size = shared_sizes,
                     font_size = 5)

    nx.draw_networkx_nodes(G,pos=nx.spring_layout(G),
                           nodelist = atlantic_list,
                           node_color = 'r',
                           node_size = atlantic_sizes)
    nx.draw_networkx_nodes(G,pos=nx.spring_layout(G),
                           nodelist=pacific_list,
                           node_color='b',
                           node_size=pacific_sizes)
    #nx.draw_networkx_edges(G,pos=nx.spring_layout(G))
    nx.draw_networkx_edges(G,pos=nx.spring_layout(G),
                       edgelist=my_edges, edge_color='b', edge_width=5)
    plt.axis('off')
    plt.savefig("my_partitioned_venn.jpg")
    plt.show()
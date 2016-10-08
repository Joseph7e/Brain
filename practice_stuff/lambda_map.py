#!/usr/bin/python3

#practice with the lambda and map functions
# FROM: http://www.python-course.eu/lambda.php

#lambda --> lambda argument_lst: expression

f = lambda x,y: x+y
print (f(1,1))


# map() --> r = map(func, seq)

#used here to create a new list from a list
list = [32.0,33.1,34,35,36]
print ("my list", list)
new_list = map(lambda x: x + 1.2, list)
print ("my new list", [l for l in new_list])

##map can be applied ot several lists, they MUST BE THE SAME LENGTH!!
list1 = [1,2,3,4]
list2 = [2,3,4,5]
list3 = [-1,-2,-3,-4]
#print (map(lambda x,y,z: x + y + z, list1, list2, list3)) # prints as map object
print ("map to alter three lists", [i for i in map(lambda x,y,z: x + y + z, list1, list2, list3)]) # forces to print as a list

# filter() --> filter(function, list) # an elegant way to filter a list

fib =  [0,1,1,2,3,4,5,8,13,21,34,55]
print ('original list', fib)
result = filter(lambda x: x %2 == 0, fib)
print ("only even numbers remain", [i for i in result])


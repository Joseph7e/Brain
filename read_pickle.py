
import pickle, sys

with open(sys.argv[1], 'rb') as handle:
  b = pickle.load(handle)

print (b)

print (len(b))
print (type(b))

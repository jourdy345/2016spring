x = [40694,
73593,
13612,
65541,
19386,
2347,
26723,
42533,
27999,
96272]


def insertion_sort(x):
  for j in range(1,len(x)):
    k = x[j]
    # Insert x[j]
    i = j-1
    while i > -1 and x[i] > k:
      x[i+1] = x[i]
      i = i-1
    x[i+1] = k
  return(x)

print(insertion_sort(x))
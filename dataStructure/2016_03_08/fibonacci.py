def fib(n):
  if n <= 2:
    return 1
  else:
    return fib(n - 1) + fib(n - 2)

# In python, the start of a statement is marked by a colon along with an indentation in the next line.

def fastFib(n):
  a, b = 0, 1
  for k in range(n):
    a, b = b, a+b # the 'a' in the RHS is not the one in the RHS. Python distinguishes LHS from RHS and does not mix them up.
  return a

# For n = 1, the slower version of fib function is actually faster than the faster one.
# In small values of n, the comparison between two functions is actually not revealing as much as we would expect it to be.
# Thus, we introduce a new concept O(n) ("Big oh"): asymptotic time complexity.

print(fastFib(45))
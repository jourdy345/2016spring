#include <stdio.h>
int fib(int n)
{
  int a = 0, b = 1;
  int k, tmp;
  for (k = 0; k < n; k++) {
    tmp = b;
    b = a + b;
    a = tmp;
  }
  return a;
}

int main(void)
{
  int n = 45;
  printf("fib(%d): %d\n", n, fib(n));
}

// This makes the code much much faster but there is a tradeoff. We have declared 2 more variables than the slower code thereby using as more memory as 4 bytes. (Computationally not so expensive)
// One of the aims of learning data structure is to predict how long the code will take in advance of running it.
// Fortunately, the execution time of a specific algorithm is knowable before running it but it takes a lot of learning to do so.

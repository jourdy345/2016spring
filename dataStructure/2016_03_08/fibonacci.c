#include <stdio.h>

int fib(int n)
{
  if (n <= 2)
    return 1;
  else
    return fib(n - 1) + fib(n - 2);
}

int main(void)
{
  int n = 45;
  printf("fib(%d) = %d\n", n, fib(n));
}

// The number of additions in each layer
// layer 1: 1
// layer 2: 2
// layer 3: 4
// layer 4: 8
// layer 5: 16
// layer 6: 32
// layer 7: 64
// ...
// layer n: 2^(n-1)
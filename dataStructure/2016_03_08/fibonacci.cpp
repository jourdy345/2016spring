#include <iostream>

int fib(int n) {
  if (n <= 2) {
    return 1;
  } else {
    return fib(n-1) + fib(n-2);
  }
}

int main() {
  std::cout << fib(45) << std::endl;
}

// The execution time of this recursive function is terribly long.

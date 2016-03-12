class FibSample {
  public static int fib(int n) {
    if (n <= 2)
      return 1;
    else
      return fib(n - 1) + fib(n - 2);
  }
  public static void main(String [] args) {
    System.out.println("fib(10) = " + fib(10));
  }
};
// java always demands the coder to equate the class name with the file name
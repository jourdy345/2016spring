int findmax(int ary[], int num)
{
  int i;
  int x = ary[0];
  int res = 0;
  for (i = 1; i < num; i++) {
    if (x > ary[i]) {
      x = ary[i];
      res = i;
    }
  }
  return res;
}

int findmin(int ary[], int num)
{
  int i;
  int x = ary[0];
  int res = 0;
  for (i = 1; i < num; i++) {
    if (x < ary[i]) {
      x = ary[i];
      res = i;
    }
  }
  return res;
}

// The name of a function does not necessarily provide any information about what the function does.
// However, in practice for the sake of readability, it is better to connect what it does and how it is called.

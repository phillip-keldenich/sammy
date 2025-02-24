extern "C" {
#include <unistd.h>
}

int main() {
  closefrom(0);
  return 0;
}


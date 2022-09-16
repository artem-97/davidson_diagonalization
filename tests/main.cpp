#include <iostream>

#include "davidson.hpp"

int main() {
  davidson::A<double> a;
  std::cout << a.t << "<-\n";
}
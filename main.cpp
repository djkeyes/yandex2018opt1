#include <iostream>

using std::istream;
using std::ostream;
using std::endl;

void run(istream& in, ostream& out) {
  out << "got the following from input stream: " << endl;
  out << "'" << in.rdbuf() << "'";
}

// XXX No-good dirty terrible hack: ifndef the main function, so the Tester can include this file and re-use the
// shared bits.
#ifndef DK_LOCAL_DEVELOPMENT_BUILD

int main() {
  std::ios::sync_with_stdio(false);
  run(std::cin, std::cout);
  return 0;
}

#endif // DK_LOCAL_DEVELOPMENT_BUILD
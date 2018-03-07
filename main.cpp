#include <iostream>

using std::istream;
using std::ostream;
using std::endl;

struct Coord {
  int16_t x, y;

  Coord(int16_t x, int16_t y) : x(x), y(y) {}

  Coord() : Coord(0, 0) {}

  bool operator==(const Coord &that) const {
    return x == that.x && y == that.y;
  }

  bool operator<(const Coord &that) const {
    if (x < that.x) {
      return true;
    } else {
      return y < that.y;
    }
  }
};

ostream& operator<<(ostream& out, const Coord& coord) {
  out << coord.x << " " << coord.y;
  return out;
}

void run(istream &in, ostream &out) {
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
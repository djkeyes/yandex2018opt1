#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <cmath>
#include <limits>
#include <map>

using std::istream;
using std::ostream;
using std::endl;
using std::vector;
using std::list;
using std::set;
using std::map;
using std::numeric_limits;

struct Vec {
  int16_t x, y;

  Vec(int16_t x, int16_t y) : x(x), y(y) {}

  Vec() : Vec(0, 0) {}

  bool operator==(const Vec &that) const {
    return (x == that.x) && (y == that.y);
  }

  bool operator<(const Vec &that) const {
    if (x != that.x) {
      return x < that.x;
    } else {
      return y < that.y;
    }
  }

  Vec &operator+=(const Vec &that) {
    x += that.x;
    y += that.y;
    return *this;
  }

  Vec operator+(const Vec &that) const {
    Vec result(*this);
    return result += that;
  }

  Vec &operator-=(const Vec &that) {
    x -= that.x;
    y -= that.y;
    return *this;
  }

  Vec operator-(const Vec &that) const {
    Vec result(*this);
    return result -= that;
  }

  int32_t normsq() const {
    return static_cast<int32_t>(x) * static_cast<int32_t>(x) + static_cast<int32_t>(y) * static_cast<int32_t>(y);
  }

  double norm() const {
    return std::sqrt(static_cast<double>(normsq()));
  }

  int32_t distsq(const Vec &that) const {
    return (*this - that).normsq();
  }

  int32_t dist(const Vec &that) const {
    return std::sqrt(static_cast<double>(distsq(that)));
  }

  bool inBounds() const {
    return -10000 <= x && x <= 10000 && -10000 <= y && y <= 10000;
  }
};

ostream &operator<<(ostream &out, const Vec &coord) {
  out << coord.x << " " << coord.y;
  return out;
}

istream &operator>>(istream &in, Vec &coord) {
  in >> coord.x >> coord.y;
  return in;
}

const vector<Vec> first_100_coords_incr_mag = {{1,   0},
                                               {0,   1},
                                               {-1,  0},
                                               {0,   -1},
                                               {-1,  1},
                                               {1,   1},
                                               {-1,  -1},
                                               {1,   -1},
                                               {0,   2},
                                               {-2,  0},
                                               {0,   -2},
                                               {2,   0},
                                               {2,   1},
                                               {-2,  -1},
                                               {-2,  1},
                                               {-1,  2},
                                               {1,   -2},
                                               {1,   2},
                                               {-1,  -2},
                                               {2,   -1},
                                               {2,   -2},
                                               {-2,  -2},
                                               {-2,  2},
                                               {2,   2},
                                               {3,   0},
                                               {0,   3},
                                               {0,   -3},
                                               {-3,  0},
                                               {1,   3},
                                               {-3,  1},
                                               {3,   1},
                                               {-3,  -1},
                                               {1,   -3},
                                               {-1,  3},
                                               {-1,  -3},
                                               {3,   -1},
                                               {-3,  2},
                                               {3,   2},
                                               {-2,  3},
                                               {-2,  -3},
                                               {3,   -2},
                                               {2,   3},
                                               {2,   -3},
                                               {-3,  -2},
                                               {-4,  0},
                                               {4,   0},
                                               {0,   -4},
                                               {0,   4},
                                               {1,   -4},
                                               {-4,  1},
                                               {-1,  4},
                                               {4,   1},
                                               {-1,  -4},
                                               {1,   4},
                                               {-4,  -1},
                                               {4,   -1},
                                               {-3,  3},
                                               {3,   3},
                                               {-3,  -3},
                                               {3,   -3},
                                               {4,   2},
                                               {-4,  2},
                                               {4,   -2},
                                               {-4,  -2},
                                               {2,   -4},
                                               {2,   4},
                                               {-2,  -4},
                                               {-2,  4},
                                               {5,   0},
                                               {-3,  -4},
                                               {3,   -4},
                                               {4,   -3},
                                               {4,   3},
                                               {-4,  -3},
                                               {0,   -5},
                                               {-4,  3},
                                               {3,   4},
                                               {-3,  4},
                                               {-5,  0},
                                               {0,   5},
                                               {5,   -1},
                                               {1,   -5},
                                               {-5,  1},
                                               {-1,  5},
                                               {5,   1},
                                               {-5,  -1},
                                               {-1,  -5},
                                               {1,   5},
                                               {2,   -5},
                                               {-2,  5},
                                               {5,   -2},
                                               {-2,  -5},
                                               {5,   2},
                                               {-5,  2},
                                               {2,   5},
                                               {-5,  -2},
                                               {4,   4},
                                               {4,   -4},
                                               {-4,  4},
                                               {-4,  -4},
                                               {3,   5},
                                               {-3,  5},
                                               {3,   -5},
                                               {5,   -3},
                                               {5,   3},
                                               {-3,  -5},
                                               {-5,  -3},
                                               {-5,  3},
                                               {6,   0},
                                               {0,   -6},
                                               {-6,  0},
                                               {0,   6},
                                               {1,   -6},
                                               {6,   -1},
                                               {-1,  6},
                                               {1,   6},
                                               {-6,  1},
                                               {-6,  -1},
                                               {6,   1},
                                               {-1,  -6},
                                               {-6,  -2},
                                               {-2,  6},
                                               {-6,  2},
                                               {2,   6},
                                               {6,   2},
                                               {-2,  -6},
                                               {6,   -2},
                                               {2,   -6},
                                               {-4,  -5},
                                               {4,   -5},
                                               {-5,  -4},
                                               {-4,  5},
                                               {5,   -4},
                                               {-5,  4},
                                               {5,   4},
                                               {4,   5},
                                               {3,   6},
                                               {-3,  6},
                                               {6,   3},
                                               {-6,  -3},
                                               {6,   -3},
                                               {3,   -6},
                                               {-6,  3},
                                               {-3,  -6},
                                               {7,   0},
                                               {0,   7},
                                               {0,   -7},
                                               {-7,  0},
                                               {-7,  -1},
                                               {5,   5},
                                               {7,   1},
                                               {1,   7},
                                               {-5,  5},
                                               {5,   -5},
                                               {7,   -1},
                                               {-1,  -7},
                                               {-5,  -5},
                                               {-1,  7},
                                               {-7,  1},
                                               {1,   -7},
                                               {6,   -4},
                                               {6,   4},
                                               {-6,  -4},
                                               {-4,  6},
                                               {-6,  4},
                                               {4,   -6},
                                               {-4,  -6},
                                               {4,   6},
                                               {-7,  -2},
                                               {-2,  7},
                                               {-2,  -7},
                                               {2,   -7},
                                               {2,   7},
                                               {7,   2},
                                               {7,   -2},
                                               {-7,  2},
                                               {-7,  3},
                                               {-3,  -7},
                                               {7,   3},
                                               {-3,  7},
                                               {3,   -7},
                                               {3,   7},
                                               {7,   -3},
                                               {-7,  -3},
                                               {5,   -6},
                                               {-6,  5},
                                               {-5,  -6},
                                               {-5,  6},
                                               {-6,  -5},
                                               {5,   6},
                                               {6,   -5},
                                               {6,   5},
                                               {-8,  0},
                                               {0,   8},
                                               {8,   0},
                                               {0,   -8},
                                               {7,   -4},
                                               {4,   7},
                                               {-4,  7},
                                               {-1,  8},
                                               {4,   -7},
                                               {-8,  1},
                                               {8,   1},
                                               {-7,  -4},
                                               {1,   8},
                                               {-8,  -1},
                                               {-1,  -8},
                                               {7,   4},
                                               {1,   -8},
                                               {-4,  -7},
                                               {-7,  4},
                                               {8,   -1},
                                               {2,   8},
                                               {-2,  8},
                                               {8,   2},
                                               {2,   -8},
                                               {-8,  2},
                                               {8,   -2},
                                               {-8,  -2},
                                               {-2,  -8},
                                               {-6,  -6},
                                               {6,   -6},
                                               {6,   6},
                                               {-6,  6},
                                               {3,   8},
                                               {-8,  3},
                                               {-3,  8},
                                               {-3,  -8},
                                               {8,   3},
                                               {-8,  -3},
                                               {3,   -8},
                                               {8,   -3},
                                               {-5,  7},
                                               {5,   7},
                                               {5,   -7},
                                               {-5,  -7},
                                               {-7,  5},
                                               {-7,  -5},
                                               {7,   -5},
                                               {7,   5},
                                               {8,   4},
                                               {4,   -8},
                                               {-8,  -4},
                                               {-8,  4},
                                               {4,   8},
                                               {-4,  8},
                                               {-4,  -8},
                                               {8,   -4},
                                               {-9,  0},
                                               {0,   9},
                                               {0,   -9},
                                               {9,   0},
                                               {-1,  9},
                                               {-9,  1},
                                               {1,   -9},
                                               {1,   9},
                                               {-9,  -1},
                                               {9,   -1},
                                               {-1,  -9},
                                               {9,   1},
                                               {9,   -2},
                                               {-9,  -2},
                                               {-9,  2},
                                               {2,   -9},
                                               {-6,  7},
                                               {-6,  -7},
                                               {2,   9},
                                               {7,   6},
                                               {-2,  9},
                                               {6,   7},
                                               {-2,  -9},
                                               {7,   -6},
                                               {-7,  -6},
                                               {6,   -7},
                                               {-7,  6},
                                               {9,   2},
                                               {-8,  -5},
                                               {8,   -5},
                                               {5,   -8},
                                               {5,   8},
                                               {-5,  8},
                                               {8,   5},
                                               {-5,  -8},
                                               {-8,  5},
                                               {-9,  3},
                                               {3,   9},
                                               {-3,  9},
                                               {-9,  -3},
                                               {3,   -9},
                                               {-3,  -9},
                                               {9,   -3},
                                               {9,   3},
                                               {-4,  -9},
                                               {9,   -4},
                                               {-4,  9},
                                               {9,   4},
                                               {4,   9},
                                               {-9,  -4},
                                               {4,   -9},
                                               {-9,  4},
                                               {7,   -7},
                                               {-7,  -7},
                                               {7,   7},
                                               {-7,  7},
                                               {8,   -6},
                                               {-8,  6},
                                               {-10, 0},
                                               {0,   -10},
                                               {6,   8},
                                               {-6,  -8},
                                               {8,   6},
                                               {6,   -8},
                                               {-8,  -6},
                                               {-6,  8},
                                               {1,   -10},
                                               {-1,  -10},
                                               {-10, -1},
                                               {-10, 1},
                                               {-10, -2},
                                               {-10, 2},
                                               {-2,  -10},
                                               {2,   -10},
                                               {9,   -5},
                                               {-9,  5},
                                               {5,   9},
                                               {-9,  -5},
                                               {-5,  -9},
                                               {5,   -9},
                                               {9,   5},
                                               {-5,  9},
                                               {-10, -3},
                                               {-10, 3},
                                               {-3,  -10},
                                               {3,   -10},
                                               {7,   8},
                                               {-7,  8},
                                               {-8,  -7},
                                               {-7,  -8},
                                               {7,   -8},
                                               {8,   7},
                                               {8,   -7},
                                               {-8,  7},
                                               {4,   -10},
                                               {-10, 4},
                                               {-10, -4},
                                               {-4,  -10},
                                               {9,   6},
                                               {6,   -9},
                                               {9,   -6},
                                               {-6,  9},
                                               {-9,  -6},
                                               {-9,  6},
                                               {6,   9},
                                               {-6,  -9},
                                               {-10, 5},
                                               {5,   -10},
                                               {-10, -5},
                                               {-5,  -10},
                                               {-8,  -8},
                                               {8,   -8},
                                               {-8,  8},
                                               {8,   8},
                                               {7,   -9},
                                               {9,   7},
                                               {9,   -7},
                                               {-7,  -9},
                                               {-9,  7},
                                               {7,   9},
                                               {-7,  9},
                                               {-9,  -7},
                                               {-10, 6},
                                               {-10, -6},
                                               {-6,  -10},
                                               {6,   -10},
                                               {-8,  -9},
                                               {-9,  -8},
                                               {8,   9},
                                               {-9,  8},
                                               {-8,  9},
                                               {9,   8},
                                               {9,   -8},
                                               {8,   -9},
                                               {-10, -7},
                                               {-10, 7},
                                               {7,   -10},
                                               {-7,  -10},
                                               {-9,  9},
                                               {-9,  -9},
                                               {9,   9},
                                               {9,   -9},
                                               {8,   -10},
                                               {-8,  -10},
                                               {-10, -8},
                                               {-10, 8},
                                               {-9,  -10},
                                               {9,   -10},
                                               {-10, 9},
                                               {-10, -9},
                                               {-10, -10}};

struct Description {
  vector<Vec> taxi_start_coords;
  vector<Vec> pedestrian_start_coords;
  vector<Vec> zone_coords;

  explicit Description(istream &in) {
    int T;
    in >> T;
    taxi_start_coords.resize(T);
    for (int i = 0; i < T; ++i) {
      in >> taxi_start_coords[i];
    }
    int P;
    in >> P;
    pedestrian_start_coords.resize(P);
    for (int i = 0; i < P; ++i) {
      in >> pedestrian_start_coords[i];
    }
    int Z;
    in >> Z;
    zone_coords.resize(Z);
    for (int i = 0; i < Z; ++i) {
      in >> zone_coords[i];
    }
  }
};

struct Move {
  Vec displacement;
  list<int> taxis;
};

ostream &operator<<(ostream &out, const Move &move) {
  out << "MOVE " << move.displacement << " " << move.taxis.size();
  for (const int &taxi_idx : move.taxis) {
    out << " " << (taxi_idx+1);
  }
  return out;
}

struct Solution {
  list<Move> moves;
};

ostream &operator<<(ostream &out, const Solution &sln) {
  out << sln.moves.size() << endl;
  for (const auto &move : sln.moves) {
    out << move << endl;
  }
  return out;
}


class Strategy {
 protected:
  const Description &description;

 public:
  explicit Strategy(const Description &description) : description(description) {}

  virtual ~Strategy() = default;

  virtual Solution run() = 0;
};

class SimpleGreedy : public Strategy {
  set<Vec> pedestrians_remaining;
  vector<Vec> current_taxis;
  const vector<Vec> &zones;

 public:
  explicit SimpleGreedy(const Description &description) : Strategy(description), zones(description.zone_coords) {
    current_taxis.insert(current_taxis.begin(), description.taxi_start_coords.begin(),
                         description.taxi_start_coords.end());
    pedestrians_remaining.insert(description.pedestrian_start_coords.begin(),
                                 description.pedestrian_start_coords.end());
  }

  virtual Solution run() {
    // it's tricky to do bookkeeping for this

    // precompute min ped-zone dists (PZ)
    // while undelivered pedestrians (P)
    //  for each pedestrian (P)
    //    find closest taxi (T)
    //    find closest zone (1)
    //    run delivery, and displace taxis if necessary (1)
    // total: O(P^2 (T))
    //
    // alternative:
    // precompute T trees of ped dists (T PlogP)
    // put trees onto heap (T log T)
    // while undelivered pedestrians (P)
    //  get shortest tree and get shortest dist (logT + logP)
    //  do route
    //  update the tree (PlogP or similar)
    //  update the ped in every other tree (TlogP)
    //  update the heap (logT)
    // total: P^2 log P
    // (potentially worse if lots of pedestrians)
    //


    // dumb way:
    map<Vec, Vec> closest_zones = computeClosestZones();
    while (!pedestrians_remaining.empty()) {
      set<Vec>::iterator best_ped;
      int best_taxi_idx;
      int32_t min_dist = numeric_limits<int32_t>::max();
      for (auto ped_iter = pedestrians_remaining.begin(); ped_iter != pedestrians_remaining.end(); ++ped_iter) {
        for (int i = 0; i < current_taxis.size(); ++i) {
          int32_t dist = ped_iter->distsq(current_taxis[i]);
          if (dist < min_dist) {
            min_dist = dist;
            best_taxi_idx = i;
            best_ped = ped_iter;
          }
        }
      }

      Vec ped = *best_ped;
      Vec &taxi = current_taxis[best_taxi_idx];
      pedestrians_remaining.erase(best_ped);
      Vec &target = closest_zones.at(ped);
      moveTaxiIfOccupied(target, ped, taxi);
      current_sln.moves.push_back({(ped - taxi), {best_taxi_idx}});
      current_sln.moves.push_back({(target - ped), {best_taxi_idx}});
      taxi = target;
    }

    return current_sln;
  }

 private:
  map<Vec, Vec> computeClosestZones() {
    map<Vec, Vec> result;
    for (const auto &ped : pedestrians_remaining) {
      Vec closest_zone;
      int32_t min_dist = numeric_limits<int32_t>::max();
      for (const auto &zone : zones) {
        int32_t dist = ped.distsq(zone);
        if (dist < min_dist) {
          min_dist = dist;
          closest_zone = zone;
        }
      }
      result[ped] = closest_zone;
    }
    return result;
  }

  void moveTaxiIfOccupied(const Vec &zone, const Vec &ped_to_avoid, const Vec &taxi_to_ignore) {
    int occupant;
    if ((occupant = getTaxiOccupying(zone, taxi_to_ignore)) == -1) {
      return;
    }
    for (const Vec &displacement : first_100_coords_incr_mag) {
      Vec target = zone + displacement;
      if (!target.inBounds()) {
        continue;
      }
      if (getTaxiOccupying(target) == -1) {
        current_taxis[occupant] = target;
        current_sln.moves.push_back({displacement, {occupant}});
        break;
      }
    }
  }

  int getTaxiOccupying(const Vec &zone, const Vec &taxi_to_ignore) {
    for (int i = 0; i < current_taxis.size(); ++i) {
      if (current_taxis[i] == taxi_to_ignore) {
        continue;
      }
      if (current_taxis[i] == zone) {
        return i;
      }
    }
    return -1;
  }

  int getTaxiOccupying(const Vec &zone) {
    for (int i = 0; i < current_taxis.size(); ++i) {
      if (current_taxis[i] == zone) {
        return i;
      }
    }
    return -1;
  }

  Solution current_sln;
};


void run(istream &in, ostream &out) {
  Description descr(in);

  SimpleGreedy strat(descr);
  Solution sln = strat.run();
  out << sln << endl;
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
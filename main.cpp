#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <algorithm>
#include <iterator>
#include <memory>
#include <chrono>

using std::istream;
using std::ostream;
using std::endl;
using std::vector;
using std::list;
using std::set;
using std::map;
using std::pair;
using std::numeric_limits;
using std::priority_queue;
using std::unique_ptr;
using std::function;

template<typename U, typename V>
ostream &operator<<(ostream &out, const pair<U, V> &p) {
  out << "(" << p.first << " " << p.second << ")";
  return out;
};

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

  double dist(const Vec &that) const {
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

  Description() = default;
};

istream &operator>>(istream &in, Description &descr) {
  int T;
  in >> T;
  descr.taxi_start_coords.resize(static_cast<size_t>(T));
  for (int i = 0; i < T; ++i) {
    in >> descr.taxi_start_coords[i];
  }
  int P;
  in >> P;
  descr.pedestrian_start_coords.resize(static_cast<size_t>(P));
  for (int i = 0; i < P; ++i) {
    in >> descr.pedestrian_start_coords[i];
  }
  int Z;
  in >> Z;
  descr.zone_coords.resize(static_cast<size_t>(Z));
  for (int i = 0; i < Z; ++i) {
    in >> descr.zone_coords[i];
  }
  return in;
}

struct Move {
  Vec displacement;
  list<int> taxis;

  double perMovePenalty(const int32_t num_taxis_total) const {
    return displacement.norm() * (1.0 + static_cast<double>(taxis.size()) / num_taxis_total);
  }
};

ostream &operator<<(ostream &out, const Move &move) {
  out << "MOVE " << move.displacement << " " << move.taxis.size();
  for (const int &taxi_idx : move.taxis) {
    out << " " << (taxi_idx + 1);
  }
  return out;
}

struct Solution {
  list<Move> moves;

  double penalty(const int32_t num_taxis_total) const {
    double total = 0.0;
    for (const Move &m : moves) {
      total += m.perMovePenalty(num_taxis_total);
    }
    return total;
  }
};

ostream &operator<<(ostream &out, const Solution &sln) {
  out << sln.moves.size() << endl;
  for (const auto &move : sln.moves) {
    if (move.displacement.normsq() > 0) {
      out << move << endl;
    }
  }
  return out;
}


class Strategy {
 protected:
  const Description &description;

 public:
  explicit Strategy(const Description &description) : description(description) {}

  virtual ~Strategy() = default;

  Solution run() {
    unique_ptr<Solution> sln;
    while (sln == nullptr) {
      sln = runIncremental();
    }
    return *sln;
  }

  virtual unique_ptr<Solution> runIncremental() = 0;
};

class SimpleGreedy : public Strategy {
  set<Vec> pedestrians_remaining;
  vector<Vec> current_taxis;
  const vector<Vec> &zones;

 public:
  explicit SimpleGreedy(const Description &description) : Strategy(description), zones(description.zone_coords),
                                                          currentSln(new Solution()) {
    current_taxis.insert(current_taxis.begin(), description.taxi_start_coords.begin(),
                         description.taxi_start_coords.end());
    pedestrians_remaining.insert(description.pedestrian_start_coords.begin(),
                                 description.pedestrian_start_coords.end());
    closestZones = computeClosestZones();
  }

  unique_ptr<Solution> runIncremental() override {
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
    if (!pedestrians_remaining.empty()) {
      set<Vec>::iterator best_ped;
      int best_taxi_idx = -1;
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
      Vec &target = closestZones.at(ped);
      moveTaxiIfOccupied(target, ped, taxi);
      currentSln->moves.push_back({(ped - taxi), {best_taxi_idx}});
      currentSln->moves.push_back({(target - ped), {best_taxi_idx}});
      taxi = target;

    }

    if (!pedestrians_remaining.empty()) {
      return unique_ptr<Solution>(nullptr);
    } else {
      return std::move(currentSln);
    }
  }

 private:
  map<Vec, Vec> closestZones;

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
        currentSln->moves.push_back({displacement, {occupant}});
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

  unique_ptr<Solution> currentSln;
};

struct Cost {
  double fixed, variable;
};

ostream &operator<<(ostream &out, const Cost c) {
  out << "{Cost=" << c.fixed << "+" << c.variable << "k}";
  return out;
}

class GreedyGroups : public Strategy {
  set<Vec> pedestrians_remaining;
  struct Taxi {
    bool occupied;
    int id;
    Vec loc;
  };
  vector<Taxi> taxis;
  set<Vec> zones;
  int numClosestTargets;

 public:
  explicit GreedyGroups(const Description &description, int num_closest_targets)
      : Strategy(description), numClosestTargets(num_closest_targets), currentSln(new Solution()) {
    pedestrians_remaining.insert(description.pedestrian_start_coords.begin(),
                                 description.pedestrian_start_coords.end());
    zones.insert(description.zone_coords.begin(),
                 description.zone_coords.end());
    for (int i = 0; static_cast<size_t>(i) < description.taxi_start_coords.size(); ++i) {
      Taxi t = {false, i, description.taxi_start_coords[i]};
      taxis.emplace_back(t);
    }
  }

  unique_ptr<Solution> runIncremental() override {
    if (numPedestriansLeft() > 0) {
      map<int, list<pair<Vec, Cost>>> targets = computeClosestZones();

      Vec best_displacement;
      double most_penalty_reduced = numeric_limits<double>::lowest();
      list<int> taxis_with_positive_reduction;

      for (int id = 0; id < targets.size(); ++id) {
        if (targets[id].empty()) {
          continue;
        }
        Cost &base_cost = targets[id].front().second;
        // the best move is to go straight to the closest pedestrian
        for (const auto &path : targets[id]) {
          list<int> cur_valid_taxis;
          cur_valid_taxis.push_back(id);
          // anything other than the closest pedestrian gives more penalty
          double penalty_reduced = base_cost.fixed + base_cost.variable - path.second.fixed - path.second.variable;
          Vec displacement = path.first - taxis[id].loc;
          for (int otherid = 0; otherid < targets.size(); ++otherid) {
            if (id == otherid) {
              continue;
            }
            if (targets[otherid].empty()) {
              continue;
            }
            Cost &other_shortest_cost = targets[otherid].front().second;
            // before the easiest cost was this
            double other_baseline = other_shortest_cost.variable + other_shortest_cost.fixed;

            Vec where_other_ends_up = taxis[otherid].loc + displacement;
            if (!where_other_ends_up.inBounds()) {
              continue;
            }

            double max_reduction_for_other = numeric_limits<double>::lowest();
            for (const auto &otherpath : targets[otherid]) {
              // now we pay the path variable cost, plus the full cost to actually reach our destination
              double alternative_cost =
                  path.second.variable + where_other_ends_up.dist(otherpath.first) * (1. + 1. / taxis.size());
              double reduction = other_baseline - alternative_cost;
              max_reduction_for_other = std::max(max_reduction_for_other, reduction);
            }
            if (max_reduction_for_other >= 0) {
              cur_valid_taxis.push_back(otherid);
            }
            penalty_reduced += max_reduction_for_other;
          }
          if (penalty_reduced > most_penalty_reduced) {
            best_displacement = displacement;
            most_penalty_reduced = penalty_reduced;
            taxis_with_positive_reduction = move(cur_valid_taxis);
          }
        }
      }

      // TODO: dedup the taxi destinations. this is a little important when taking taxis to zones (can always re-use
      // a zone), but especially important when when taking taxis to pedestrians and P < T (no need to send 20 to one
      // guy).

      set<int> taxis_to_move(taxis_with_positive_reduction.begin(), taxis_with_positive_reduction.end());

      moveTaxisOutOfWay(taxis_to_move, best_displacement);

      // execute move
      currentSln->moves.emplace_back();
      Move &big_move = currentSln->moves.back();
      big_move.displacement = best_displacement;
      big_move.taxis.insert(big_move.taxis.end(), taxis_to_move.begin(), taxis_to_move.end());
      for (int id : taxis_to_move) {
        taxis[id].loc += best_displacement;
        updateTaxiCheckpoints(taxis[id]);
      }
    }

    if (numPedestriansLeft() > 0) {
      return unique_ptr<Solution>(nullptr);
    } else {
      return std::move(currentSln);
    }
  }

 private:

  void moveTaxisOutOfWay(const set<int> &moving_taxis, const Vec &displacement) {
    set<Vec> taxi_dests;
    set<Vec> taxi_srcs;
    for (int id : moving_taxis) {
      taxi_dests.insert(taxis[id].loc + displacement);
      taxi_srcs.insert(taxis[id].loc);
    }
    map<Vec, int> stationary_taxis;
    for (const Taxi &t : taxis) {
      if (moving_taxis.find(t.id) == moving_taxis.end()) {
        stationary_taxis.emplace(t.loc, t.id);
      }
    }

    for (const Vec &dest : taxi_dests) {
      auto iter = stationary_taxis.find(dest);
      if (iter != stationary_taxis.end()) {
        int blocking_taxi = iter->second;
        for (const Vec &d : first_100_coords_incr_mag) {
          Vec target = dest + d;
          if (!target.inBounds()) {
            continue;
          }
          if (stationary_taxis.find(target) != stationary_taxis.end()) {
            continue;
          }
          if (taxi_dests.find(target) != taxi_dests.end()) {
            continue;
          }
          if (taxi_srcs.find(target) != taxi_srcs.end()) {
            continue;
          }
          // this is a good target.
          // perform move, and update stationary_taxis
          stationary_taxis.erase(iter);
          stationary_taxis[target] = blocking_taxi;
          taxis[blocking_taxi].loc = target;
          currentSln->moves.push_back({d, {blocking_taxi}});

          // update unintentional checkpoints
          updateTaxiCheckpoints(taxis[blocking_taxi]);
          break;
        }
      }
    }
  }

  /*
   * Updates the pedestrian locations / taxi occupancy after a taxi's location has been updated
   */
  void updateTaxiCheckpoints(Taxi &t) {
    if (t.occupied) {
      if (zones.find(t.loc) != zones.end()) {
        t.occupied = false;
      }
    } else {
      auto ped_iter = pedestrians_remaining.find(t.loc);
      if (ped_iter != pedestrians_remaining.end()) {
        pedestrians_remaining.erase(ped_iter);
        t.occupied = true;
      }
    }
  }

  int numPedestriansLeft() {
    auto total = static_cast<int>(pedestrians_remaining.size());
    for (const Taxi &t : taxis) {
      total += t.occupied;
    }
    return total;
  }

  map<int, list<pair<Vec, Cost>>> computeClosestZones() {
    // for each taxi, computes the closest numClosestTargets targets. For occupied taxis, targets are zones
    // (potentially few); for unoccupied, targets are pedestrians (potentially many).
    map<int, list<pair<Vec, Cost>>> result;
    for (const auto &taxi : taxis) {
      auto &closest = result[taxi.id];

      set<Vec> *targets;
      if (taxi.occupied) {
        targets = &zones;
      } else {
        targets = &pedestrians_remaining;
      }
      using DistVec = pair<int32_t, Vec>;
      priority_queue<DistVec> heap;
      for (const auto &target : *targets) {
        int32_t dist = taxi.loc.distsq(target);
        heap.emplace(dist, target);
        if (heap.size() > numClosestTargets) {
          heap.pop();
        }
      }
      while (!heap.empty()) {
        auto &top = heap.top();
        double dist = sqrt(static_cast<double>(top.first));
        Cost c = {dist, dist / taxis.size()};
        closest.emplace_front(top.second, c);
        heap.pop();
      }
    }
    return result;
  }

  unique_ptr<Solution> currentSln;
};

class CombinedStrategy : public Strategy {
 public:
  CombinedStrategy(const Description &description, const list<function<unique_ptr<Strategy>()>> &strategy_generator)
      : Strategy(description), strategies(strategy_generator), startTime(std::chrono::steady_clock::now()) {
  }

  // evaluate several strategies and choose the best one
  // this implemention of runIncremental actually runs the whole thing, so don't nest instances of CombinedStrategy
  unique_ptr<Solution> runIncremental() override {
    unique_ptr<Solution> best_solution(nullptr);
    double lowest_penalty = numeric_limits<double>::max();
    for (const auto &gen : strategies) {
      unique_ptr<Strategy> strat = gen();

      unique_ptr<Solution> cur_solution;
      while (cur_solution == nullptr) {
        if (needToTerminate()) {
          break;
        }
        cur_solution = strat->runIncremental();
      }
      if (cur_solution != nullptr) {
        double cur_penalty = cur_solution->penalty(static_cast<int32_t>(description.taxi_start_coords.size()));
        if (cur_penalty < lowest_penalty) {
          lowest_penalty = cur_penalty;
          std::swap(best_solution, cur_solution);
        }
      }
      if (needToTerminate()) {
        break;
      }
    }
    return best_solution;
  }

 private:
  const list<function<unique_ptr<Strategy>()>> &strategies;

  std::chrono::time_point<std::chrono::steady_clock> startTime;

  bool needToTerminate() {
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::steady_clock::now() - startTime);
    return elapsed.count() > 1900;
  }
};

void run(istream &in, ostream &out) {
  Description descr;
  in >> descr;

  list<function<unique_ptr<Strategy>()>> strategy_generators = {
      [&]() { return unique_ptr<Strategy>(new SimpleGreedy(descr)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 1)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 2)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 3)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 4)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 5)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 10)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 20)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 15)); },
      [&]() { return unique_ptr<Strategy>(new GreedyGroups(descr, 40)); },
  };
  CombinedStrategy strat(descr, strategy_generators);
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
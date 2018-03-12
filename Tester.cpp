
#include <set>
#include <iostream>
#include <sstream>
#include <random>
#include <utility>

#define DK_LOCAL_DEVELOPMENT_BUILD

#include "main.cpp"

#undef DK_LOCAL_DEVELOPMENT_BUILD

using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ostream;
using std::endl;
using std::cerr;
using std::cout;
using std::flush;
using std::string;
using std::set;
using std::tuple;
using std::make_tuple;
using std::get;

using std::uniform_int_distribution;

class SimpleGreedyForTest : public Strategy {
  set<Vec> pedestrians_remaining;
  vector<Vec> current_taxis;
  const vector<Vec> &zones;

 public:
  explicit SimpleGreedyForTest(const Description &description)
      : Strategy(description),
        zones(description.zone_coords),
        currentSln(new Solution()) {
    current_taxis.insert(current_taxis.begin(), description.taxi_start_coords.begin(),
                         description.taxi_start_coords.end());
    pedestrians_remaining.insert(description.pedestrian_start_coords.begin(),
                                 description.pedestrian_start_coords.end());
  }

  unique_ptr<Solution> runIncremental() override {
    map<Vec, Vec> closest_zones = computeClosestZones();
    while (!pedestrians_remaining.empty()) {
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
      Vec &target = closest_zones.at(ped);
      moveTaxiIfOccupied(target, ped, taxi);
      currentSln->moves.push_back({(ped - taxi), {best_taxi_idx}});
      currentSln->moves.push_back({(target - ped), {best_taxi_idx}});
      taxi = target;
    }

    return move(currentSln);
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


istream &operator>>(istream &in, Move &move) {
  string text;
  int k;
  in >> text >> move.displacement >> k;
  for (int j = 0; j < k; ++j) {
    int id;
    in >> id;
    move.taxis.push_back(id - 1);
  }
  return in;
}

istream &operator>>(istream &in, Solution &sln) {
  int K;
  in >> K;
  for (int j = 0; j < K; ++j) {
    sln.moves.emplace_back();
    in >> sln.moves.back();
  }
  return in;
}

ostream &operator<<(ostream &out, Description &descr) {
  out << descr.taxi_start_coords.size() << endl;
  for (const auto &coord : descr.taxi_start_coords) {
    out << coord << endl;
  }

  out << descr.pedestrian_start_coords.size() << endl;
  for (const auto &coord : descr.pedestrian_start_coords) {
    out << coord << endl;
  }

  out << descr.zone_coords.size() << endl;
  for (const auto &coord : descr.zone_coords) {
    out << coord << endl;
  }
}

const int num_testcases = 30;

class TestcaseGenerator {
 public:
  TestcaseGenerator() : tDist(1, 20), pDist(1, 500), zDist(1, 20), sDist(1, 10000) {
    rng.seed(seed);
  }

  Description generate() {
    int32_t T, P, Z;
    int16_t S;
    int num_attempts = 0;
    do {
      T = tDist(rng);
      P = pDist(rng);
      Z = zDist(rng);
      S = sDist(rng);
      ++num_attempts;
    } while (T + P + Z > (2 * S + 1) * (2 * S + 1));
    cout << "Picked T=" << T << ", P=" << P << ", Z=" << Z << ", S=" << S << " after " << num_attempts << " rolls" <<
         endl;
    return generate(T, P, Z, S);
  }

  Description generate(int32_t T, int32_t P, int32_t Z, int16_t S) {
    set<Vec> occupied_coords;
    uniform_int_distribution<int16_t> coord_dist(-S, S);

    Description descr;
    descr.taxi_start_coords.reserve(static_cast<size_t>(T));
    for (int i = 0; i < T; ++i) {
      descr.taxi_start_coords.push_back(generateNonOccupiedCoord(occupied_coords, coord_dist));
    }

    descr.pedestrian_start_coords.reserve(static_cast<size_t>(P));
    for (int i = 0; i < P; ++i) {
      descr.pedestrian_start_coords.push_back(generateNonOccupiedCoord(occupied_coords, coord_dist));
    }

    descr.zone_coords.reserve(static_cast<size_t>(Z));
    for (int i = 0; i < Z; ++i) {
      descr.zone_coords.push_back(generateNonOccupiedCoord(occupied_coords, coord_dist));
    }
    return descr;
  }

 private:
  Vec generateNonOccupiedCoord(set<Vec> &occupied_coords, uniform_int_distribution<int16_t> coord_dist) {
    // generate via rejection sampling samples
    // if sample space size >>>> num samples, then this is efficient
    Vec sample;
    do {
      int16_t x = coord_dist(rng);
      int16_t y = coord_dist(rng);
      sample = Vec(x, y);
    } while (occupied_coords.find(sample) != occupied_coords.end());
    occupied_coords.insert(sample);
    return sample;
  }

  const uint32_t seed = 6843;
  std::mt19937 rng;

  uniform_int_distribution<int32_t> tDist;
  uniform_int_distribution<int32_t> pDist;
  uniform_int_distribution<int32_t> zDist;
  uniform_int_distribution<int16_t> sDist;
};

bool verify(const Solution &output_sln, const Description &descr);

double computeScore(const Solution &output_sln, const Description &description);

int main() {
  TestcaseGenerator testgen;

  vector<tuple<int32_t, int32_t, int32_t, int16_t>> edge_case_generate_args
      = {
          make_tuple(7, 1, 1, 1),
          make_tuple(1, 1, 7, 1),
          make_tuple(1, 7, 1, 1),
          make_tuple(1, 1, 1, 1),
          make_tuple(5, 15, 5, 2),
      };

  float totalScore = 0.0;
  for (size_t testcase = 0; testcase < num_testcases; ++testcase) {
    stringstream input;
    Description descr;
    if (testcase < edge_case_generate_args.size()) {
      auto &arguments = edge_case_generate_args[testcase];
      descr = testgen.generate(get<0>(arguments), get<1>(arguments), get<2>(arguments), get<3>(arguments));
    } else {
      descr = testgen.generate();
    }
    input << descr;
    string input_copy = input.str();
    //cout << "input: `" << input_copy << "`" << endl;
    stringstream output_stream;
    run(input, output_stream);
    Solution sln;
    output_stream >> sln;
    double score = 0.0;
    if (verify(sln, descr)) {
      score = computeScore(sln, descr);
    }

    totalScore += score;
  }
  totalScore /= num_testcases;
  cout << "Final score: " << totalScore << endl;

  return 0;
}

double computeScore(const Solution &output_sln, const Description &description) {
  double penalty = output_sln.penalty(static_cast<int32_t>(description.taxi_start_coords.size()));

  // I submitted this and it got 5161.31pts
  // ergo the true threshold is roughly 0.516131 * simple greedy penalty
  SimpleGreedyForTest baseline(description);
  double simply_greedy_penalty = baseline.run().penalty(static_cast<int32_t>(description.taxi_start_coords.size()));
  double threshold = 0.516131 * simply_greedy_penalty;

  return threshold / penalty * 10000;
}

bool verify(const Solution &output_sln, const Description &descr) {
  set<Vec> pedestrians(descr.pedestrian_start_coords.begin(), descr.pedestrian_start_coords.end());
  set<Vec> zones(descr.zone_coords.begin(), descr.zone_coords.end());
  // <current coord, original idx, is occupied>
  map<Vec, int> taxi_ids_by_coord;
  map<int, Vec> taxis_by_id;
  map<int, bool> taxi_occupied_by_id;
  for (int i = 0; i < descr.taxi_start_coords.size(); ++i) {
    const Vec &taxi = descr.taxi_start_coords[i];
    taxi_ids_by_coord[taxi] = i;
    taxis_by_id[i] = taxi;
    taxi_occupied_by_id[i] = false;
  }

  for (const Move &move : output_sln.moves) {
    const list<int> &ids = move.taxis;
    set<int> displaced_ids_set(ids.begin(), ids.end());

    for (int taxi_id : ids) {
      Vec dest = taxis_by_id[taxi_id] + move.displacement;
      if (!dest.inBounds()) {
        cerr << "Out of bounds!" << endl;
        return false;
      }
      if (taxi_ids_by_coord.find(dest) != taxi_ids_by_coord.end()) {
        // ignore taxis currently in motion
        if (displaced_ids_set.find(taxi_ids_by_coord.at(dest)) == displaced_ids_set.end()) {
          cerr << "Taxis collide!" << endl;
          return false;
        }
      }
      if (taxi_occupied_by_id[taxi_id]) {
        // occupied - check for zone
        if (zones.find(dest) != zones.end()) {
          taxi_occupied_by_id[taxi_id] = false;
        }
      } else {
        // unoccupied - check for passenger
        if (pedestrians.find(dest) != pedestrians.end()) {
          taxi_occupied_by_id[taxi_id] = true;
          pedestrians.erase(dest);
        }
      }
      taxi_ids_by_coord.erase(taxis_by_id[taxi_id]);
      taxi_ids_by_coord[dest] = taxi_id;
      taxis_by_id[taxi_id] = dest;
    }
  }
  if (!pedestrians.empty()) {
    cerr << "Some passengers were not picked up!" << endl;
    return false;
  }
  for (const auto element : taxi_occupied_by_id) {
    if (element.second) {
      cerr << "Some passengers are still in taxis!" << endl;
      return false;
    }
  }
  return true;
}

// TODO (roughly prioritized)
// -implement tougher edge case tests, ie (2*S+1) = T + P + Z
// -other ideas:
//    -store a list of "strategies" to run through and keep one with lowest penalty
//    -something greedy, but with different weighting heuristics
//    -something max flow based
//    -something having to do with crowd alignment?
//    -given a strategy, perturb it randomly. repeat until time runs out
//    -find a matrix decomposition? we want to find a minimal set of 2d vectors approximating the ped-zone displacements
// -implement time check in verify

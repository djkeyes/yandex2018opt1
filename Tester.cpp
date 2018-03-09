
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

using std::uniform_int_distribution;


const int num_testcases = 1;

class TestcaseGenerator {
 public:
  TestcaseGenerator() : tDist(1, 20), pDist(1, 500), zDist(1, 20), sDist(1, 10000) {
    rng.seed(seed);
  }

  void generate(stringstream &generated_input) {
    set<Vec> occupied_coords;
    int T, P, Z, S;
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
    uniform_int_distribution<int16_t> coord_dist(-S, S);

    generated_input << T << endl;
    for (int i = 0; i < T; ++i) {
      generated_input << generateNonOccupiedCoord(occupied_coords, coord_dist) << endl;
    }

    generated_input << P << endl;
    for (int i = 0; i < P; ++i) {
      generated_input << generateNonOccupiedCoord(occupied_coords, coord_dist) << endl;
    }

    generated_input << Z << endl;
    for (int i = 0; i < Z; ++i) {
      generated_input << generateNonOccupiedCoord(occupied_coords, coord_dist) << endl;
    }
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

bool verify(const string &input, const string &output);

float computeScore(const string &output);

int main() {
  TestcaseGenerator testgen;

  float totalScore = 0.0;
  for (int testcase = 0; testcase < num_testcases; ++testcase) {
    stringstream input;
    testgen.generate(input);
    string input_copy = input.str();
    //cout << "input: `" << input_copy << "`" << endl;
    stringstream output_stream;
    run(input, output_stream);
    string output = output_stream.str();
    //cout << "output: `" << output << "`" << endl;
    float score = 0.0;
    if (verify(input_copy, output)) {
      score = computeScore(output);
    }

    totalScore += score;
  }
  totalScore /= num_testcases;
  cout << "Final score: " << totalScore << endl;

  return 0;
}

float computeScore(const string &output) {
  return 1.0;
}

bool verify(const string &input, const string &output) {

  istringstream iss(input);
  Description descr(iss);

  set<Vec> pedestrians(descr.pedestrian_start_coords.begin(), descr.pedestrian_start_coords.end());
  set<Vec> zones(descr.zone_coords.begin(), descr.zone_coords.end());
  // <current coord, original idx, is occupied>
  map<Vec, int> taxi_ids_by_coord;
  map<int, Vec> taxis_by_id;
  map<int, bool> taxi_occupied_by_id;
  for (int i = 0; i < descr.taxi_start_coords.size(); ++i) {
    Vec &taxi = descr.taxi_start_coords[i];
    taxi_ids_by_coord[taxi] = i;
    taxis_by_id[i] = taxi;
    taxi_occupied_by_id[i] = false;
  }

  int K;
  istringstream answer(output);
  answer >> K;
  for (int i = 0; i < K; ++i) {
    string text;
    Vec displacement;
    int k;
    answer >> text >> displacement >> k;
    vector<int> ids(k);
    for (int j = 0; j < k; ++j) {
      answer >> ids[j];
      ids[j]--;
    }
    set<int> displaced_ids_set(ids.begin(), ids.end());

    for (int taxi_id : ids) {
      Vec dest = taxis_by_id[taxi_id] + displacement;
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
// -implement penalty
// -implement score (use greediest strategy as baseline)
// -implement time check in verify
// -other ideas:
//    -store a list of "strategies" to run through and keep one with lowest penalty
//    -something greedy, but with different weighting heuristics
//    -something max flow based
//    -something having to do with crowd alignment?
//    -given a strategy, perturb it randomly. repeat until time runs out

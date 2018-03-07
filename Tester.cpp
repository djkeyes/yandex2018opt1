
#include <set>
#include <iostream>
#include <sstream>
#include <random>

#define DK_LOCAL_DEVELOPMENT_BUILD

#include "main.cpp"

#undef DK_LOCAL_DEVELOPMENT_BUILD

using std::stringstream;
using std::ostringstream;
using std::ostream;
using std::endl;
using std::cout;
using std::flush;
using std::string;
using std::set;

using std::uniform_int_distribution;


const int num_testcases = 1;

class TestcaseGenerator {
 public:
  TestcaseGenerator() : tDist(1, 20), pDist(1, 500), zDist(1, 20), coordDist(-10000, 10000) {
    rng.seed(seed);
  }

  void generate(stringstream &generated_input) {
    set<Coord> occupied_coords;
    int T = tDist(rng);
    generated_input << T << endl;
    for (int i=0; i < T; ++i) {
      generated_input << generateNonOccupiedCoord(occupied_coords) << endl;
    }

    int P = pDist(rng);
    generated_input << P << endl;
    for (int i=0; i < P; ++i) {
      generated_input << generateNonOccupiedCoord(occupied_coords) << endl;
    }

    int Z = zDist(rng);
    generated_input << Z << endl;
    for (int i=0; i < Z; ++i) {
      generated_input << generateNonOccupiedCoord(occupied_coords) << endl;
    }
  }

 private:
  Coord generateNonOccupiedCoord(set<Coord>& occupied_coords) {
    // generate via rejection sampling samples
    // if sample space size >>>> num samples, then this is efficient
    Coord sample;
    do {
      int16_t x = coordDist(rng);
      int16_t y = coordDist(rng);
      sample = Coord(x, y);
    } while(occupied_coords.find(sample) != occupied_coords.end());
    occupied_coords.insert(sample);
    return sample;
  }
  const uint32_t seed = 6843;
  std::mt19937 rng;

  uniform_int_distribution<int32_t> tDist;
  uniform_int_distribution<int32_t> pDist;
  uniform_int_distribution<int32_t> zDist;
  uniform_int_distribution<int16_t> coordDist;
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
    stringstream output_stream;
    run(input, output_stream);
    string output = output_stream.str();
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
  return 0;
}

bool verify(const string &input, const string &output) {
  cout << input << endl;
  return false;
}


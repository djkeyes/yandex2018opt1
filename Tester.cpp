
#include <iostream>
#include <sstream>

#define DK_LOCAL_DEVELOPMENT_BUILD
#include "main.cpp"
#undef DK_LOCAL_DEVELOPMENT_BUILD

using std::stringstream;
using std::ostringstream;
using std::ostream;
using std::endl;
using std::cout;
using std::flush;

const int num_testcases = 30;

int main() {
  for (int testcase = 0; testcase < num_testcases; ++testcase) {
    stringstream input;
    input << "hello world! testcase #" << testcase << endl;
    stringstream output;
    run(input, output);
    cout << output.rdbuf() << flush;
  }
  return 0;
}

// Print test conditions
#include <iostream>
#include <cstdlib>
using namespace std;

double getEnvDouble(const char* name, double defaultVal) {
    const char* val = getenv(name);
    if (val) {
        cout << name << " from env: " << val << endl;
        return atof(val);
    }
    cout << name << " default: " << defaultVal << endl;
    return defaultVal;
}

int main() {
    cout << "TEST_Te = " << getEnvDouble("TEST_Te", 10.0) << endl;
    cout << "TEST_ne = " << getEnvDouble("TEST_ne", 1e12) << endl;
    cout << "TEST_THi = " << getEnvDouble("TEST_THi", 8.0) << endl;
    cout << "TEST_nHi = " << getEnvDouble("TEST_nHi", 9e11) << endl;
    cout << "TEST_THeIII = " << getEnvDouble("TEST_THeIII", 10.0) << endl;
    cout << "TEST_nHeIII = " << getEnvDouble("TEST_nHeIII", 1e8) << endl;
    return 0;
}

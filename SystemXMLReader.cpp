#include "OpenMM.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace OpenMM;
using namespace std;

int main() {

    ifstream system_input("system.xml");
    System* systemXML = XmlSerializer::deserialize<System>(system_input);

    delete systemXML;

    return 0;
}

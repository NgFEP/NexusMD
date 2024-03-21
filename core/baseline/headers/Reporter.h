#ifndef REPORTER_H
#define REPORTER_H

#include "Coords3D.h"
#include <vector>
#include <string>

namespace BaseLine {

    class Reporter {
    public:
        // Constructs the Reporter object.
        Reporter() = default;

        // Reports simulation data.
        // The clearFile parameter determines whether to clear existing file data (true for the first step).
        static void report(const std::string& filename, const std::vector<Coords3D>& positions,
            const std::vector<Coords3D>& velocities, const std::vector<Coords3D>& forces,
            int step, bool clearFile);
    };

} // namespace BaseLine

#endif // REPORTER_H

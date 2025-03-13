#ifndef MONTECARLO_BAROSTAT_H
#define MONTECARLO_BAROSTAT_H

#include <vector>
#include "Coords3D.h"
#include "PeriodicBoundaryCondition.h"
#include "PDBResidueParser.h"

namespace BaseLine {

    class MonteCarloBarostat {
    public:
        MonteCarloBarostat(double pressure, double temperature, int frequency);

        void ApplyBarostat(std::vector<Coords3D>& atomPositions,
            const std::vector<Molecules>& molecules, int& numMolecules,
            PeriodicBoundaryCondition::BoxInfo& boxInfo);

    private:
        double _pressure;
        double _temperature;
        int _frequency;  // Number of MD steps between barostat moves
    };
} // namespace Baseline


#endif // MONTECARLO_BAROSTAT_H
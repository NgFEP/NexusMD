//#ifndef ANDERSEN_THERMOSTAT_H
//#define ANDERSEN_THERMOSTAT_H
//
//#include "Coords3D.h"
//#include <vector>
//#include <cmath>
//#include <random>
//
//#define kb    (1.380649e-23)  /* (J/K) */ // Boltzmann constant
//#define boltz (0.008314462618)	/*kJ / (mol·K)*/
//
//
////#define KILO         (1e3)
////#define BOLTZMANN    (1.380649e-23)            /* (J/K)   */
////#define AVOGADRO     (6.02214076e23)          
////#define RGAS         (BOLTZMANN*AVOGADRO)      /* (J/(mol K))  */
////#define BOLTZ2        (RGAS/KILO)               /* (kJ/(mol K)) */
//
//
//namespace BaseLine {
//
//    class AndersenThermostat {
//    public:
//        // Constructor with default collision frequency of 10.0
//        AndersenThermostat(double& temperature, double& collisionFrequency);
//
//        // Apply the thermostat, collisionFrequency is optional (uses default if not set)
//        void apply(std::vector<Coords3D>& velocities, const std::vector<double>& inverseMasses, double& stepSize);
//
//    private:
//        double _collisionFrequency;
//        double _temperature;
//        double collisionProbability;
//        std::default_random_engine generator;
//        std::uniform_real_distribution<double> rng{ 0.0, 1.0 };
//        void computeCollisionProbability(double& stepSize);
//        bool shouldCollide();
//    };
//
//} // namespace Baseline
//
//#endif // ANDERSEN_THERMOSTAT_H


//#ifndef ANDERSEN_THERMOSTAT_H
//#define ANDERSEN_THERMOSTAT_H
//
//#include "Coords3D.h"
//#include <vector>
//#include <cmath>
//#include <random>
//
//#define kb    (1.380649e-23)  /* (J/K) */ // Boltzmann constant
//#define boltz (0.008314462618)	/* kJ / (mol·K) */
//
//namespace BaseLine {
//
//    class AndersenThermostat {
//    public:
//        // Constructor with temperature and collision frequency
//        AndersenThermostat(double temperature, double collisionFrequency);
//
//        // Apply the thermostat to velocities based on masses and step size
//        void apply(std::vector<Coords3D>& velocities, const std::vector<double>& inverseMasses, double stepSize);
//
//    private:
//        double _temperature;           // Temperature in Kelvin
//        double _collisionFrequency;    // Collision frequency in fs^-1
//        double collisionProbability;   // Computed probability of collision
//
//        // Random number generators (static for efficiency)
//        static std::default_random_engine generator;
//        static std::uniform_real_distribution<double> distribution;
//        static std::normal_distribution<double> normal_distribution;
//
//        void computeCollisionProbability(double stepSize);  // Computes collision probability
//        bool shouldCollide();  // Determines if a collision occurs
//    };
//
//} // namespace BaseLine
//
//#endif // ANDERSEN_THERMOSTAT_H


#ifndef ANDERSEN_THERMOSTAT_H
#define ANDERSEN_THERMOSTAT_H

#include "Coords3D.h"
#include <vector>
#include <cmath>
#include <random>

#define kb    (1.380649e-23)  /* (J/K) */ // Boltzmann constant
#define boltz (0.008314462618) /* kJ/(mol·K) */

namespace BaseLine {

    class AndersenThermostat {
    public:
        // Constructor: Initializes thermostat with temperature and collision frequency
        AndersenThermostat(double temperature, double collisionFrequency);

        // Apply the thermostat to modify velocities
        void apply(std::vector<Coords3D>& velocities, const std::vector<double>& inverseMasses, double stepSize);

    private:
        double _collisionFrequency;  // Collision frequency (1/ps)
        double _temperature;         // Target temperature (K)
        double collisionProbability; // Computed probability of collision per step

        // Static random number generator setup
        static std::mt19937 generator;
        static bool _randomInitialized;
        static double nextGaussian;
        static bool nextGaussianIsValid;

        // Random distributions
        std::uniform_real_distribution<double> uniform_dist;
        std::normal_distribution<double> normal_dist;

        // Computes collision probability based on step size
        void computeCollisionProbability(double stepSize);

        // Determines if a collision should occur (uniform random sampling)
        bool shouldCollide();

        // Generate normally distributed random number (Box-Muller transform)
        double GaussianDist();

        // Generate uniformly distributed random number in the range [0,1)
        double UniformDist();
        // Initialize random seed
        static void initializeRandomSeed();
    };

} // namespace BaseLine

#endif // ANDERSEN_THERMOSTAT_H

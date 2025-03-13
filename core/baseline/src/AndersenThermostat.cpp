#include "AndersenThermostat.h"
#include <random>
#include <cmath>
#include <chrono>

using namespace std;
using namespace BaseLine;

// Static initialization for RNG state
std::mt19937 AndersenThermostat::generator;
bool AndersenThermostat::_randomInitialized = false;
double AndersenThermostat::nextGaussian = 0.0;
bool AndersenThermostat::nextGaussianIsValid = false;

// Initialize random seed
void AndersenThermostat::initializeRandomSeed() {
    std::random_device rd;
    std::mt19937::result_type seed = rd() ^ (
        (std::mt19937::result_type)
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now().time_since_epoch()
        ).count()
        );

    generator.seed(seed);
    _randomInitialized = true;
    nextGaussianIsValid = false;
}


// Constructor with temperature and collision frequency
AndersenThermostat::AndersenThermostat(double temperature, double collisionFrequency)
    : _temperature(temperature), _collisionFrequency(collisionFrequency),
    uniform_dist(0.0, 1.0), normal_dist(0.0, 1.0) {
    initializeRandomSeed();
}

// Apply the thermostat to velocities
void AndersenThermostat::apply(vector<Coords3D>& velocities, const vector<double>& inverseMasses, double stepSize) {
    computeCollisionProbability(stepSize);

    for (size_t i = 0; i < velocities.size(); ++i) {
        if (shouldCollide()) {  // Collision check
            double invMass = inverseMasses[i];
            double sigma = sqrt(boltz * _temperature * invMass);  // Correct standard deviation
            velocities[i][0] = sigma * GaussianDist();
            velocities[i][1] = sigma * GaussianDist();
            velocities[i][2] = sigma * GaussianDist();
        }
    }
}

// Compute collision probability per step
void AndersenThermostat::computeCollisionProbability(double stepSize) {
    collisionProbability = 1.0 - exp(-_collisionFrequency * stepSize);
}

// Determines if a collision should occur (Uses uniform distribution)
bool AndersenThermostat::shouldCollide() {
    return UniformDist() < collisionProbability;
}

// Generate normally distributed random number using Box-Muller transform
double AndersenThermostat::GaussianDist() {
    if (nextGaussianIsValid) {
        nextGaussianIsValid = false;
        return nextGaussian;
    }

    double u, v, s;
    do {
        u = 2.0 * UniformDist() - 1.0;
        v = 2.0 * UniformDist() - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    double factor = sqrt(-2.0 * log(s) / s);
    nextGaussian = v * factor;
    nextGaussianIsValid = true;

    return u * factor;
}

// Generate uniformly distributed random number in the range [0,1)
double AndersenThermostat::UniformDist() {
    if (!_randomInitialized) {
        initializeRandomSeed();
    }
    return uniform_dist(generator);
}



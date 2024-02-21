#ifndef VerletIntegration_H
#define VerletIntegration_H


namespace BaseLine {

    class VerletIntegration : public ReferenceDynamics {

    private:

        std::vector<OpenMM::Vec3> xPrime;
        std::vector<double> inverseMasses;

    public:

        /**------------------------

           Constructor & Destructor

           @param AtomsQuantity  Atoms Quantity
           @param TStepSize      Time Step Size

           --------------------- */

        VerletIntegration(int AtomsQuantity, double TStepSize);

        ~VerletIntegration();

        /**---------------------------------------------------------------------------------------

           Update

           @param system              the System to be integrated
           @param atomCoordinates     atom coordinates
           @param velocities          velocities
           @param forces              forces
           @param masses              atom masses
           @param tolerance           the constraint tolerance

           --------------------------------------------------------------------------------------- */

        void update(const OpenMM::System& system, std::vector<OpenMM::Vec3>& atomCoordinates,
            std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses, double tolerance);

    };

} // namespace BaseLine

#endif // VerletIntegration_H

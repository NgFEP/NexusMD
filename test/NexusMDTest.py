import NexusMD


def test_task_dispatcher():
    try:
        # Create a CPU-based TaskDispatcher instance using the factory method
        print("Creating TaskDispatcher for CPU...")
        dispatcher = NexusMD.TaskDispatcher.CreateDispatcher("CPU")
        print(f"TaskDispatcher instance created: {dispatcher}")

        # Define simulation parameters
        system_filename = "system_Sample_Protein.xml"
        state_filename = "state_Sample_Protein.xml"
        input_filename = "Sample_Protein.pdb"
        output_filename = "output_SimRun_Sample_Protein.pdb"
        step_size = 0.001  # 0.001 pico s = 1 femto s
        total_steps = 10
        interval = 1

        # Directly enable specific forces
        print("Enabling forces...")
        # Uncomment any of the the next 4 lines to enable force calculations
        dispatcher.EnableHarmonicBondForce()
        dispatcher.EnableHarmonicAngleForce()
        dispatcher.EnablePeriodicTorsionForce()
        dispatcher.EnableNonbondedForce()


        # Run the simulation
        print("Running simulation...")
        dispatcher.Simulate(
            system_filename, state_filename,
            input_filename, output_filename,
            step_size, total_steps, interval
        )
        print("Simulation completed successfully with selected forces enabled.")

    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    test_task_dispatcher()

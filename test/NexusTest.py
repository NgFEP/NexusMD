import Nexus 

def test_task_dispatcher():
    try:
        # Create a CPU-based TaskDispatcher instance using the factory method
        print("Creating TaskDispatcher for CPU...")
        dispatcher = Nexus.TaskDispatcher.CreateDispatcher("CPU")
        print(f"TaskDispatcher instance created: {dispatcher}")

        # Define dummy parameters for testing the Simulate method
        system_filename = "system_Sample_Protein.xml"
        state_filename = "state_Sample_Protein.xml"
        input_filename = "Sample_Protein.pdb"
        output_filename = "Output_SimRun_Sample_Protein.pdb"
        step_size = 0.001
        total_steps = 10
        interval = 1


        # Call the Simulate method
        print("Running simulation...")
        dispatcher.Simulate(
            system_filename, state_filename,
            input_filename, output_filename,
            step_size, total_steps, interval
        )
        print("Simulation completed successfully.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    test_task_dispatcher()

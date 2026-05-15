#!/bin/bash
# Example bash script to run multiple PhaseField simulations with different configurations
# on an HPC cluster

# Set the path to the compiled executable
EXECUTABLE="./build/PhaseFieldModel"

# Base directory where config directories will be created
BASE_DIR="./simulation_runs"

# Create base directory if it doesn't exist
mkdir -p "$BASE_DIR"

# Example: Run 4 simulations with different temperature gradients
for temp_grad in 10000 15000 20000 25000; do
    # Create a unique directory for this simulation
    RUN_DIR="$BASE_DIR/run_tempgrad_${temp_grad}"
    mkdir -p "$RUN_DIR"
    
    # Copy the base config to this directory
    cp ../config/modelConfig.json "$RUN_DIR/modelConfig.json"
    
    # Modify the temperature gradient using sed (or use jq if available for JSON)
    # This example uses sed - modify as needed for your JSON structure
    # sed -i "s/\"tempGradKperm\": .*/\"tempGradKperm\": $temp_grad,/" "$RUN_DIR/modelConfig.json"
    
    echo "Starting simulation in $RUN_DIR with tempGrad=$temp_grad"
    
    # Run the simulation with the config directory as argument
    # For HPC: use sbatch, bsub, qsub, etc. depending on your scheduler
    $EXECUTABLE "$RUN_DIR" > "$RUN_DIR/simulation.log" 2>&1 &
    
    # Optional: Add delay between job submissions
    sleep 2
done

# Alternative: Using SLURM (sbatch)
# For HPC clusters with SLURM scheduler, modify like this:
#
# for temp_grad in 10000 15000 20000 25000; do
#     RUN_DIR="$BASE_DIR/run_tempgrad_${temp_grad}"
#     mkdir -p "$RUN_DIR"
#     cp ../config/modelConfig.json "$RUN_DIR/modelConfig.json"
#     
#     # Submit job to SLURM
#     sbatch -J "PF_run_$temp_grad" \
#            -o "$RUN_DIR/slurm.log" \
#            --wrap="$EXECUTABLE \"$RUN_DIR\""
# done

echo "All simulations submitted. Check $BASE_DIR for results."

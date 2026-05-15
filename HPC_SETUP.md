# HPC Multi-Run Configuration Guide

## Overview
The Phase Field Model executable can now accept a configuration directory as a command-line argument. This allows you to easily run multiple simulations with different settings on an HPC cluster.

## Usage

### Basic Usage
```bash
./PhaseFieldModel /path/to/config/directory
```

The program expects to find `modelConfig.json` in the specified directory and will save all outputs to that directory:
- `solidification_tracking.csv` - Solidification progress data
- `viz_frames/` - Visualization frames
- `simulation_video.mp4` - Generated video
- Other output files from `saveGridData()`

### Default Behavior
If no argument is provided:
```bash
./PhaseFieldModel
```
The program will use the default directory: `../config`

## Running Multiple Simulations on HPC

### 1. Create Config Directories
Create a separate directory for each simulation with its own `modelConfig.json`:

```bash
mkdir -p simulations/run_1 simulations/run_2 simulations/run_3
cp base_config.json simulations/run_1/modelConfig.json
cp base_config.json simulations/run_2/modelConfig.json
cp base_config.json simulations/run_3/modelConfig.json
```

### 2. Modify Configs
Edit each `modelConfig.json` with desired parameters:

```bash
# Example using jq to modify JSON
jq '.tempGradKperm = 20000' simulations/run_1/modelConfig.json > temp && mv temp simulations/run_1/modelConfig.json
jq '.tempGradKperm = 25000' simulations/run_2/modelConfig.json > temp && mv temp simulations/run_2/modelConfig.json
```

Or using Python:
```python
import json
import os

for run_num, temp_grad in enumerate([20000, 25000, 30000], 1):
    run_dir = f"simulations/run_{run_num}"
    os.makedirs(run_dir, exist_ok=True)
    
    with open("base_config.json", "r") as f:
        config = json.load(f)
    
    config["tempGradKperm"] = temp_grad
    
    with open(f"{run_dir}/modelConfig.json", "w") as f:
        json.dump(config, f, indent=2)
```

### 3. Submit to HPC Queue

#### SLURM (Most Common)
Create a `submit_jobs.sh`:
```bash
#!/bin/bash

BASE_DIR="simulations"
EXECUTABLE="./build/PhaseFieldModel"

for run_dir in $BASE_DIR/run_*; do
    run_name=$(basename "$run_dir")
    
    sbatch -J "PF_$run_name" \
           -o "$run_dir/slurm.log" \
           -e "$run_dir/slurm_error.log" \
           --nodes=1 \
           --cpus-per-task=8 \
           --time=01:00:00 \
           --wrap="$EXECUTABLE \"$run_dir\""
done
```

#### PBS/Torque
```bash
#!/bin/bash

BASE_DIR="simulations"
EXECUTABLE="./build/PhaseFieldModel"

for run_dir in $BASE_DIR/run_*; do
    run_name=$(basename "$run_dir")
    
    qsub -N "PF_$run_name" \
         -o "$run_dir/pbs.log" \
         -e "$run_dir/pbs_error.log" \
         -l nodes=1:ppn=8 \
         -l walltime=01:00:00 \
         <<EOF
#!/bin/bash
cd $(pwd)
$EXECUTABLE "$run_dir"
EOF
done
```

#### BSUB (LSF)
```bash
#!/bin/bash

BASE_DIR="simulations"
EXECUTABLE="./build/PhaseFieldModel"

for run_dir in $BASE_DIR/run_*; do
    run_name=$(basename "$run_dir")
    
    bsub -J "PF_$run_name" \
         -o "$run_dir/lsf.log" \
         -n 8 \
         -W 01:00 \
         "$EXECUTABLE \"$run_dir\""
done
```

### 4. Monitor Jobs
Check job status using your scheduler's tools:
```bash
# SLURM
squeue -u $USER

# PBS
qstat

# LSF
bjobs
```

### 5. Collect Results
After simulations complete, collect results:
```bash
# Create a summary of all results
for run_dir in simulations/run_*; do
    run_name=$(basename "$run_dir")
    if [ -f "$run_dir/solidification_tracking.csv" ]; then
        echo "$run_name: $(tail -1 $run_dir/solidification_tracking.csv)"
    fi
done
```

## Example: Parametric Study

Create a parameter sweep script:

```bash
#!/bin/bash

EXECUTABLE="./build/PhaseFieldModel"
BASE_DIR="parametric_study"
mkdir -p "$BASE_DIR"

# Define parameter ranges
temp_grads=(15000 20000 25000 30000)
cooling_rates=(0.001 0.002 0.003 0.004)

run_num=1
for tg in "${temp_grads[@]}"; do
    for cr in "${cooling_rates[@]}"; do
        run_dir="$BASE_DIR/run_${run_num}_tg${tg}_cr${cr}"
        mkdir -p "$run_dir"
        
        # Create config with parameters
        python3 << EOF
import json
config = {
    "tempGradKperm": $tg,
    "CoolingRatempers": $cr,
    # ... other parameters ...
}
with open("$run_dir/modelConfig.json", "w") as f:
    json.dump(config, f, indent=2)
EOF
        
        # Submit to HPC
        sbatch -J "PF_run_$run_num" --wrap="$EXECUTABLE \"$run_dir\""
        
        ((run_num++))
    done
done
```

## Output Structure

Each simulation directory will contain:
```
run_1/
├── modelConfig.json          # Input configuration
├── solidification_tracking.csv  # Solidification data (every 100 steps)
├── viz_frames/               # Video frames
├── simulation_video.mp4      # Generated video (if visualization enabled)
├── [output files from saveGridData()]
└── [logs and other outputs]
```

## Performance Tips

1. **Disable visualization on HPC** - Set `"enableVisualization": false` in config for headless runs
2. **Disable video generation** - Set smaller grid or disable frames to save I/O
3. **Use fast filesystem** - Run from scratch space if available (`/scratch`, `/tmp`, etc.)
4. **Monitor resources** - Check memory and CPU usage for your grid sizes

## Troubleshooting

**Config file not found:**
```
[config] Config file not found. Generating default at ...
```
Solution: Ensure the directory path contains `modelConfig.json`

**Permission denied:**
```bash
chmod +x PhaseFieldModel
chmod +x run_multiple_simulations.sh
```

**Path issues on Windows (if using WSL):**
Convert Windows paths to WSL format:
```bash
/mnt/c/Work/PhaseField/PhaseFieldModel  # Instead of C:\Work\PhaseField\PhaseFieldModel
```

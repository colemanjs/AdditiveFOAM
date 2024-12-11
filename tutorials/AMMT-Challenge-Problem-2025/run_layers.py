import os
import shutil
import argparse
import subprocess
import polars as pl
from myna.application.additivefoam.path import convert_peregrine_scanpath

# Helper functions
def run_command(command, capture_output=False):
    """Run a shell command and optionally capture its output."""
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE if capture_output else None,
        stderr=subprocess.PIPE if capture_output else None,
        text=True
    )
    if capture_output:
        return result.stdout.strip()

def foam_dict_get(entry, filepath):
    """Get a value from a foamDictionary file."""
    return run_command(
        ["foamDictionary", "-entry", entry, "-value", filepath],
        capture_output=True
    )

def foam_dict_set(entry, value, filepath):
    """Set a value in a foamDictionary file."""
    run_command([
        "foamDictionary", "-entry", entry, "-set", str(value), filepath
    ])

def update_extrude_mesh_settings(
    previous_case, layer_thickness, n_cells_per_layer
):
    """Update extrudeMeshDict for the current layer."""
    n_layers_prev = int(
        foam_dict_get("nLayers", f"{previous_case}/system/extrudeMeshDict")
    )
    thickness_prev = float(
        foam_dict_get(
            "linearDirectionCoeffs/thickness",
            f"{previous_case}/system/extrudeMeshDict"
        )
    )
    thickness_new = thickness_prev + layer_thickness

    foam_dict_set(
        "nLayers",
        n_layers_prev + n_cells_per_layer,
        "system/extrudeMeshDict"
    )
    foam_dict_set(
        "linearDirectionCoeffs/thickness",
        thickness_new,
        "system/extrudeMeshDict"
    )

    run_command(["extrudeMesh", "-case", os.getcwd()], capture_output=False)
    return thickness_new

def update_simulation_time(previous_case, case_dir, layer_time):
    """Update start and end time for the simulation and adjust scan path."""
    start_time = float(
        foam_dict_get("endTime", f"{previous_case}/system/controlDict")
    )
    foam_dict_set("startTime", start_time, "system/controlDict")

    end_time = start_time + layer_time
    foam_dict_set("endTime", end_time, "system/controlDict")

    foam_dict_set("writeInterval", end_time - start_time, "system/controlDict")

    # Move time directory
    time_dir = os.path.join(case_dir, "0")
    if os.path.exists(time_dir):
        target_dir = os.path.join(
            case_dir,
            str(int(start_time) if start_time.is_integer() else start_time)
        )
        os.rename(time_dir, target_dir)

    # Update scan path start time
    scan_path_file = os.path.join(case_dir, "constant/scanPath")
    with open(scan_path_file, "r") as file:
        lines = file.readlines()
    scan_path = lines[1].split()
    scan_path[5] = str(start_time)
    lines[1] = "\t\t".join(scan_path) + "\n"
    with open(scan_path_file, "w") as file:
        file.writelines(lines)

    return start_time, end_time

def initialize_layer(case_dir, base_dir):
    """Initialize the case directory for a new layer."""
    os.makedirs(case_dir)
    for folder in ["0", "constant", "system"]:
        shutil.copytree(
            os.path.join(base_dir, folder), os.path.join(case_dir, folder)
        )

def convert_temperature_output(case_dir, output_file):
    """Extract the top surface temperature and write as csv"""
        
    end_time = foam_dict_get("endTime", f"{case_dir}/system/controlDict")    
    input_file = os.path.join(
        case_dir,
        f"postProcessing/topSurface/{end_time}/top.xy"
    )

    # Read and clean data in-memory
    with open(input_file, 'r') as infile:
        data = [
            line.split() for line in infile
            if not line.strip().startswith("#") and line.strip()
        ]

    # Convert the data to a Polars DataFrame with explicit orientation
    df = pl.DataFrame(
        data, schema=["x", "y", "z", "T", "solid"], orient="row"
    ).with_columns(
        [
            pl.col("x").cast(pl.Float64),
            pl.col("y").cast(pl.Float64),
            pl.col("z").cast(pl.Float64),
            pl.col("T").cast(pl.Float64),
            pl.col("solid").cast(pl.Float64),
        ]
    )

    # Write the DataFrame to a CSV file
    df.write_csv(output_file)

# Main function
def main():
    # Parse arguements
    parser = argparse.ArgumentParser(
        description="Run additiveFoam simulation for sequential layers."
    )
    parser.add_argument("-n_layers", type=int, help="Number of layers")
    parser.add_argument(
        "-layer_thickness", type=float, help="Thickness of each layer"
    )
    parser.add_argument(
        "-n_cells_per_layer", type=int, help="Number of cells per layer"
    )
    parser.add_argument(
        "-layer_time", type=float, help="Simulation time for each layer"
    )
    args = parser.parse_args()

    # Base directory and initialization
    base_dir = os.getcwd()
    case_list = []

    for layer in range(args.n_layers):
        case_dir = os.path.join(base_dir, f"layer{layer}")
        case_list.append(case_dir)
        if os.path.exists(case_dir):
            shutil.rmtree(case_dir)

    # Create background mesh
    run_command(["blockMesh"])

    # Perform simulations
    for layer in range(args.n_layers):
        
        # Create case directory
        case_dir = case_list[layer]
        previous_case = case_list[layer - 1] if layer > 0 else None
        initialize_layer(case_dir, base_dir)
        os.chdir(case_dir)

        n_procs = int(
            foam_dict_get("numberOfSubdomains", "system/decomposeParDict")
        )
        
        # Update scan path
        path_dir="/home/cloud/myna-data/AMMT-chal-prob-additivefoam/simulation/P1"
        convert_peregrine_scanpath(
            os.path.join(path_dir, f"{layer+1:07}.txt"),
            os.path.join(case_dir, "constant", "scanPath"),
            380.0)

        # Update times and map fields between layers
        if layer > 0:
            thickness_new = update_extrude_mesh_settings(
                previous_case, args.layer_thickness, args.n_cells_per_layer
            )
            start_time, end_time = update_simulation_time(
                previous_case, case_dir, args.layer_time
            )

            run_command(["decomposePar"])
            run_command([
                "mpirun", "-np", str(n_procs), "mapFieldsPar", "-case",
                case_dir, "-mapMethod", "direct", "-sourceTime",
                str(start_time), previous_case, "-parallel"
            ])
        else:
            foam_dict_set("endTime", args.layer_time, "system/controlDict")
            run_command(["decomposePar"])


        # Run AdditiveFOAM
        thickness = float(
            foam_dict_get("linearDirectionCoeffs/thickness",
                          "system/extrudeMeshDict")
        )
        run_command([
            "mpirun", "-np", str(n_procs), "transformPoints",
            f"translate=(0 0 -{thickness})", "-parallel"
        ])
        run_command(["mpirun", "-np", str(n_procs), "setFields", "-parallel"])
        run_command([
            "mpirun", "-np", str(n_procs), "macroAdditiveFoam", "-parallel"
        ])
        run_command([
            "mpirun", "-np", str(n_procs), "transformPoints",
            f"translate=(0 0 {thickness})", "-parallel"
        ])
        
        # Post process
        run_command([
            "mpirun", "-np", str(n_procs), "postProcess", "-func", "topSurface",
            "-latestTime", "-parallel"
        ])
        output_dir = os.path.join(base_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        convert_temperature_output(
            case_dir,
            os.path.join(output_dir, f"temperature_{layer+1:07}.csv")
        )

    os.chdir(base_dir)
    print("Simulation completed for all layers.")

if __name__ == "__main__":
    main()

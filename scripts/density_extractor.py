import os
def extract_density_full_and_by_type(input_file, output_file):
    
    # Guard: skip if already processed
    if os.path.exists(output_file):
        with open(output_file, 'r') as check:
            if 'All (except 3,-3)' in check.read():
                print(f"[INFO] Already processed — skipping.")
                return
    # Containers
    all_except_minus3 = []
    density_data = {
        "3,-1": [],
        "3,+1": [],
        "3,+3": []
    }

    current_cp_type = None

    with open(input_file, 'r') as f:
        for line in f:

            # Detect CP type
            if "Type (" in line:
                start = line.find("Type (") + 6
                end = line.find(")", start)
                current_cp_type = line[start:end].strip()

            # Detect density line
            if "Density of all electrons:" in line:
                if current_cp_type and current_cp_type != "3,-3":

                    value = line.split(":")[1].strip()
                    density_float = float(value)
                    density_clean = f"{density_float:.10f}".rstrip('0').rstrip('.')

                    # Add to full list
                    all_except_minus3.append(density_clean)

                    # Add to specific list if applicable
                    if current_cp_type in density_data:
                        density_data[current_cp_type].append(density_clean)

    # Write output
    with open(output_file, "a") as f:

        # 1️⃣ All except 3,-3
        f.write("\n All (except 3,-3) = {" + ", ".join(all_except_minus3) + "}\n")

        # 2️⃣ Individual types
        for cp_type in ["3,-1", "3,+1", "3,+3"]:
            f.write(f"({cp_type}) = " +
                    "{" + ", ".join(density_data[cp_type]) + "}\n")

    print("Extraction complete.")
    print("Output written to", output_file)


# Run
extract_density_full_and_by_type("60_GC_Dinucleotide_CPprop.txt","60_GC_Dinucleotide_Summary.txt")



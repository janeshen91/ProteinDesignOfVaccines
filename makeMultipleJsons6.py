#!/usr/bin/env python3
import csv
import json
import os
import sys
import copy

# File paths – update these as needed
template_json_path = '/vaxart-batch-output/202409_RFDiffusion/AF3selfRun/fold_prompt1Test.json'
csv_path = '/vaxart-batch-output/202409_RFDiffusion/VX22/AF/Above80PLDDT.csv'
output_dir = '/vaxart-batch-output/202409_RFDiffusion/AF3selfRun/Above80PLDDT_jsons'

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load the JSON template and ensure it's a dictionary
try:
    with open(template_json_path, 'r') as f:
        template = json.load(f)
        # If the template is a list, extract the first element
        if isinstance(template, list):
            template = template[0]
except Exception as e:
    sys.exit(f"Error loading template JSON: {e}")

# Open the CSV file (assuming no header)
try:
    with open(csv_path, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) < 3:
                print("Skipping row, not enough columns:", row)
                continue

            # Use integer indices for rows
            job_name = row[0].strip()
            sequence = row[2].strip()

            # Create a deep copy of the template for this job
            job_json = copy.deepcopy(template)
            job_json['name'] = job_name

            # Replace the protein sequence in the template
            if 'sequences' in job_json and isinstance(job_json['sequences'], list):
                for entry in job_json['sequences']:
                    if 'proteinChain' in entry:
                        entry['proteinChain']['sequence'] = sequence
                        break  # Remove this if you need to update more than one protein entry
            else:
                print(f"No valid 'sequences' list in template for job {job_name}")
                continue

            # Write out the new JSON file
            output_path = os.path.join(output_dir, f"{job_name}.json")
            with open(output_path, 'w') as out_f:
                json.dump([job_json], out_f, indent=2)

            print(f"Created JSON for job: {job_name} -> {output_path}")

except Exception as e:
    sys.exit(f"Error processing CSV file: {e}")

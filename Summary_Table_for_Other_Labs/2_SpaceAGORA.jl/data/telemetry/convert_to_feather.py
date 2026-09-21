import pandas as pd
import os
import sys

# Iterate over all directories in the current directory
for dir in [d for d in os.listdir() if os.path.isdir(d)]:
    print(f"Directory: {dir}")
    # Iterate over all files in the directory
    for file in os.listdir(dir):
        print(f"  File: {file}")
        # Remove the .csv extension and add .feather extension
        if file.endswith('.csv'):
            csv_file = os.path.join(dir, file)
            feather_file = os.path.join(dir, file[:-4] + '.feather')
            print(f"    Converting {csv_file} to {feather_file}")
            # Read from CSV file
            try:
                df = pd.read_csv(csv_file)
                print(f"    Successfully read {csv_file}")

                # Write to Feather file
                df.to_feather(feather_file)
                print(f"    Successfully converted to {feather_file}")

            except FileNotFoundError:
                print(f"    Error: {csv_file} not found.")
            except Exception as e:
                print(f"    An error occurred: {e}")
# # Define input and output file names
# csv_file = 'data_input.csv'
# feather_file = 'data_output.feather'

# # Read from CSV file
# try:
#     df = pd.read_csv(csv_file)
#     print(f"Successfully read {csv_file}")

#     # Write to Feather file
#     df.to_feather(feather_file)
#     print(f"Successfully converted to {feather_file}")

# except FileNotFoundError:
#     print(f"Error: {csv_file} not found.")
# except Exception as e:
#     print(f"An error occurred: {e}")

import pandas as pd
import glob
import os

for file in glob.glob('*.nanoclust_out.txt'):
    if os.path.getsize(file) > 0:  # Check if the file is not empty
        # Extract barcode number from the file name
        barcode_number = file.split('.')[0].replace('barcode', 'b')
        
        # Read the file as semicolon-separated
        try:
            df = pd.read_csv(file, delimiter=';')
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
        
        # Check if 'id' column exists
        if 'id' in df.columns:
            # Add the new 'consensus_id' column using the first 'id' column
            df.insert(0, 'consensus_id', barcode_number + 'c' + df['id'].astype(str))
            
            # Save the DataFrame as an Excel file
            output_file = file.replace('.txt', '.xlsx')
            df.to_excel(output_file, index=False)
        else:
            print(f"'id' column not found in {file}.")
    else:
        print(f"Skipping {file} as it is a zero-size file.")


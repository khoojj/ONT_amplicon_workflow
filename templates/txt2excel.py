#!/usr/bin/env python

import openpyxl

def txt_to_excel(input_file, output_file):
    # Load data from the text file
    with open(input_file, 'r') as txt_file:
        data = [line.strip().split(';') for line in txt_file]

    # Create a new workbook and select the active worksheet
    wb = openpyxl.Workbook()
    ws = wb.active

    # Write data to the worksheet
    for row in data:
        ws.append(row)

    # Save the workbook to the specified output file
    wb.save(output_file)

input_file = "${barcode}.nanoclust_out.txt"
output_file = "${barcode}.nanoclust_out.xlsx"

txt_to_excel(input_file, output_file)
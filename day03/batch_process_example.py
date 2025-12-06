"""
Example script for batch processing dilutions from Excel files.
This demonstrates how to use the process_excel_dilutions function.
"""

from dilution_core import process_excel_dilutions

# Example 1: Process Excel file with units in column headers (RECOMMENDED)
# Excel format: columns named "C1 (M)", "C2 (mM)", "V2 (mL)" with numeric values
# The output will be saved as 'dilutions_header_format_results.xlsx'
print("Processing file with units in headers (recommended format)...")
df_results = process_excel_dilutions(
    input_file='dilutions_header_format.xlsx',
    output_unit='mL'  # Output volumes in milliliters
)

print("\nProcessed dilutions:")
print(df_results)

# Example 2: Process with custom output filename and different unit
# df_results = process_excel_dilutions(
#     input_file='dilutions_header_format.xlsx',
#     output_file='my_dilutions_output.xlsx',
#     output_unit='µL'  # Output volumes in microliters
# )

# Example 3: Legacy format (units in each cell) still supported
# Excel format: columns "C1", "C2", "V2" with values like "10 M", "1 mM", "100 mL"
# df_results = process_excel_dilutions(
#     input_file='dilutions_input.xlsx',
#     output_unit='mL'
# )

# Example 4: Process and analyze results programmatically
# df = process_excel_dilutions('dilutions_header_format.xlsx', output_unit='µL')
# 
# # Find any rows with errors
# errors = df[df['Error'].notna()]
# if not errors.empty:
#     print("\nRows with errors:")
#     print(errors)

# Log results

# Stolen from https://stackoverflow.com/questions/71275961/append-new-line-to-csv-file
# #Each time you use the write w mode, it overwrites the file, deleting the old contents in the process.
# #
# #It might make more sense to first check if the file exists, and write the header row only if it doesn't. 
# #Here's one way to do that, and it shows the behaviour required in your question:

import csv
from pathlib import Path # to use in main program as file_path = Path('name.csv')

def log_to_file(file_path, headers, results):
  # Write current results to the log file in FILE_PATH
  # Check if the same number of columns in headers and results
  if len(results) != len(headers):
    print(f"ERROR: there are {len(headers)} headers and {len(results)} results")
    return
  
  # If the log file does not exist create it with a header row
  if not file_path.exists():
      with open(file_path, 'w', newline='') as log_csv:
          # Start a new blank log with column headings in the first row
          log_csv_write = csv.writer(log_csv)
          log_csv_write.writerow(headers)
  
  # Append the current results (with autocloses at exit)
  with open(file_path, 'a', newline='') as log_csv:
      log_csv_append = csv.writer(log_csv, dialect = 'excel')
      log_csv_append.writerow(results)




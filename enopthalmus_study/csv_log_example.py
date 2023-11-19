# Example of logging to a CSV file

# Below from https://stackoverflow.com/questions/71275961/append-new-line-to-csv-file
  
##Each time you use the write w mode, it overwrites the file, deleting the old contents in the process.
##
##It might make more sense to first check if the file exists, and write the header row only if it doesn't. 
##Here's one way to do that, and it shows the behaviour required in your question:

import csv
from pathlib import Path

FILE_PATH = Path('log.csv')

# Write current results to the log file in FILE_PATH
# If the log file does not exist create it with a header row
if not FILE_PATH.exists():
    with open(FILE_PATH, 'w', newline='') as log_csv:
        # Start a new blank log with column headings in the first row
        log_csv_write = csv.writer(log_csv)
        log_csv_write.writerow(["user", "date", "study", "Orbital Volume", "Air Volume", "Fat Volume", "Muscle Volume"])

# Append the current results (with autocloses at exit)
with open(FILE_PATH, 'a', newline='') as log_csv:
    log_csv_append = csv.writer(log_csv, dialect = 'excel')
    log_csv_append.writerow(results)

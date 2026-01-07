# Submission Analysis

This notebook analyzes student assignment submissions, tracking submission status and deadlines for course assignments.

## Overview

The notebook processes a tab-separated text file (`subjects.txt`) containing assignment submission records and generates a comprehensive analysis of student submission patterns.

## Features

- **Data Parsing**: Extracts student names and assignment information from submission titles
- **Case-Insensitive Assignment Detection**: Handles various formats like "Day08", "day08", "Day 08"
- **Multiple Assignments**: Automatically splits submissions containing multiple assignments (e.g., "Day03 and Day04") into separate records
- **Deadline Tracking**: Calculates deadlines for assignments (weekly schedule starting November 8, 2025 at 22:00 UTC).
- I remember that there was delay in one of the weeks so its not taken into account. If you like this script maybe you can just input the actual deadlines for each assignment. 
- **Submission Status**: Categorizes each submission as:
  - **In Time**: Submitted on or before the deadline
  - **Delayed Submission**: Submitted after the deadline
  - **Not Submitted**: No submission found for that assignment

## Input Format

The input file `subjects.txt` should be a tab-separated file with the following columns:
- ID
- Status
- Title (format: "Assignment by Student Name")
- Timestamp

## Output

The notebook generates a pivot table where:
- **Rows**: Unique student names
- **Columns**: Assignment keys (Day01 through Day08)
- **Values**: Submission status for each student-assignment combination

## Deadline Schedule

- **Day01**: November 8, 2025, 22:00 UTC
- **Day02**: November 15, 2025, 22:00 UTC
- **Day03**: November 22, 2025, 22:00 UTC
- And so on, weekly intervals at the same time

## Requirements

```python
pandas
pytz
```

## Usage

1. Place your `subjects.txt` file in the same directory as the notebook
2. Run all cells in the notebook
3. View the final pivot table showing submission status for all students

## Data Processing Steps

1. Load and clean the input data
2. Extract assignment names and student names from titles
3. Normalize assignment keys to standard format (Day01, Day02, etc.)
4. Handle multiple assignments in a single submission
5. Calculate deadlines and compare with submission timestamps
6. Generate submission status for each record
7. Create pivot table for comprehensive overview

## Notes

- The regex pattern captures assignment day numbers with or without spaces
- Student names may have variations (spacing, capitalization) which are preserved as-is
- If a student submits multiple assignments together, each is tracked separately with the same timestamp
- Assignments without deadlines are marked as "No Deadline"

import pandas as pd
import numpy as np
import sys

def compare_with_tolerance(file1, file2, default_tolerance=1.0, special_tolerance=25.0):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    if df1.shape != df2.shape:
        return False, "The files have different shapes."

    difference = pd.DataFrame(index=df1.index, columns=df1.columns)

    for column in df1.columns:
        for index in df1.index:
            tolerance = special_tolerance if column.startswith('d') else default_tolerance
            if df2.at[index, column] != 0:
                diff = np.abs(df1.at[index, column] - df2.at[index, column]) / np.abs(df2.at[index, column]) * 100
                difference.at[index, column] = diff if diff > tolerance else 0
            else:
                difference.at[index, column] = 0 if df1.at[index, column] == df2.at[index, column] else np.inf

    exceeds_tolerance = difference > 0  # Check if any differences are marked as exceeding tolerance

    if exceeds_tolerance.any().any():
        diff_locations = exceeds_tolerance.stack()
        diff_locations = diff_locations[diff_locations]
        diff_indices = np.where(exceeds_tolerance)
        differences = difference.values[diff_indices]

        report = f"Differences found at {len(differences)} locations exceeding tolerance.\n"
        sort_order = np.argsort(differences)[::-1]
        top_indices = sort_order[:10]

        for i in top_indices:
            row_index, col_index = diff_indices[0][i], diff_indices[1][i]
            diff = differences[i]
            report += f"{df1.columns[col_index]} at index {df1.index[row_index]}: {round(diff, 2)}%\n"
        
        return False, report
    else:
        biggest_difference = difference.max().max()
        return True, f"Files are similar within the specified tolerance. Biggest difference: {biggest_difference}% at Location: {df1[difference == biggest_difference].stack().index[0]}"



solution_file = "Solution/output.csv"
test_file = "Data/Test/output.csv"
tolerance_percentage = 1.0  # 1% tolerance

match, message = compare_with_tolerance(solution_file, test_file, tolerance_percentage)
print(message)

if not match:
    sys.exit(1)  # Exit with status 1 if files don't match
sys.exit(0)  # Exit with status 0 if files match

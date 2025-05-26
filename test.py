# Compare files
success = True
for filename in os.listdir(expected_output_dir):
    if filename.endswith('.log'):
        continue  # Skip .log files

    expected_file = os.path.join(expected_output_dir, filename)
    actual_file = os.path.join(actual_output_dir, filename)

    if not os.path.exists(actual_file):
        print(f"Missing output file: {filename}")
        success = False
        continue

    # Sort and compare TSV files
    if filename.endswith('.tsv'):
        expected_sorted = sort_tsv(expected_file)
        actual_sorted = sort_tsv(actual_file)
        if expected_sorted != actual_sorted:
            success = False
            print(f"Differences found in file: {filename}")
            diff = difflib.unified_diff(expected_sorted, actual_sorted, 
                                        fromfile='expected', tofile='actual')
            print(''.join(diff))
    else:
        if not filecmp.cmp(expected_file, actual_file, shallow=False):
            success = False
            print(f"Differences found in file: {filename}")
            with open(expected_file, 'r') as ef, open(actual_file, 'r') as af:
                expected_lines = ef.readlines()
                actual_lines = af.readlines()
                diff = difflib.unified_diff(expected_lines, actual_lines, 
                                            fromfile='expected', tofile='actual')
                print(''.join(diff))

if success:
    print("Test successful: All outputs match expected results.")
else:
    print("Test failed: Differences found.")

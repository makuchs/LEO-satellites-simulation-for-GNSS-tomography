def run_compare(sp3_file, csv_file, sat_id="PL99"):
    import os
    import csv
    import datetime
    import math
    import matplotlib.pyplot as plt

    def parse_sp3_file(filename, sat_id="PL99"):
        positions = {}
        current_epoch = None
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith(('#', '%', '+', '/')):
                    continue
                if line.startswith('*'):
                    tokens = line[1:].strip().split()
                    if len(tokens) < 6:
                        continue
                    try:
                        year = int(tokens[0])
                        month = int(tokens[1])
                        day = int(tokens[2])
                        hour = int(tokens[3])
                        minute = int(tokens[4])
                        second = float(tokens[5])
                        second_int = int(second)
                        microsecond = int((second - second_int) * 1e6)
                        current_epoch = datetime.datetime(year, month, day, hour, minute, second_int, microsecond)
                    except Exception as e:
                        print(f"Error parsing epoch line '{line}': {e}")
                        current_epoch = None
                    continue
                if current_epoch is not None:
                    tokens = line.split()
                    if tokens and tokens[0] == sat_id:
                        try:
                            x = float(tokens[1])
                            y = float(tokens[2])
                            z = float(tokens[3])
                            positions[current_epoch] = [x, y, z]
                        except Exception as e:
                            print(f"Error parsing position line '{line}': {e}")
        return positions

    def parse_csv_positions(csv_filename):
        positions = {}
        with open(csv_filename, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                date_str = row["Date"].strip()
                if '.' in date_str:
                    base, frac = date_str.split('.')
                    frac = frac[:6]
                    date_str_mod = base + '.' + frac
                else:
                    date_str_mod = date_str
                try:
                    dt = datetime.datetime.strptime(date_str_mod, "%Y %m %d %H %M %S.%f")
                except Exception as e:
                    print(f"Error parsing CSV date '{date_str_mod}': {e}")
                    continue
                try:
                    x = float(row["X"])
                    y = float(row["Y"])
                    z = float(row["Z"])
                    positions[dt] = [x, y, z]
                except Exception as e:
                    print(f"Error parsing CSV position in row {row}: {e}")
        return positions

    def compute_difference_magnitude(pos1, pos2):
        dx = pos1[0] - pos2[0]
        dy = pos1[1] - pos2[1]
        dz = pos1[2] - pos2[2]
        return math.sqrt(dx**2 + dy**2 + dz**2)

    sp3_positions = parse_sp3_file(sp3_file, sat_id=sat_id)
    print(f"Parsed {len(sp3_positions)} epochs from SP3 file.")
    
    csv_positions = parse_csv_positions(csv_file)
    print(f"Parsed {len(csv_positions)} epochs from CSV file.")
    
    common_epochs = sorted(set(sp3_positions.keys()) & set(csv_positions.keys()))
    print(f"Found {len(common_epochs)} common epochs.")
    
    if not common_epochs:
        print("No common epochs found between SP3 and CSV data. Cannot compare positions.")
        return
    
    times = []
    differences = []
    for epoch in common_epochs:
        diff = compute_difference_magnitude(sp3_positions[epoch], csv_positions[epoch])
        times.append(epoch)
        differences.append(diff)
    
    plt.figure(figsize=(10, 6))
    plt.plot(times, differences, linestyle='-')
    plt.xlabel("Date")
    plt.ylabel("Position Difference (km)")
    plt.title("Spatial Position Difference (SP3 vs CSV) Over Time")
    plt.grid(True)
    plt.gcf().autofmt_xdate()
    plt.show()

if __name__ == '__main__':
    # For direct testing of compare_positions, if needed.
    run_compare(r"C:\path\to\sp3_file.sp3", r"C:\path\to\simulation_output.csv", "PL77")

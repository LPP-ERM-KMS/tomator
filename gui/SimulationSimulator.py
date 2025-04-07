import csv
import random
import time

counter = 0

def generate_csv_data():
    header = [
        "tmain",
        "RadialPositions",
        "ne",
        "Ee",
        "dne",
        "dEe",
        "nue",
        "nH",
        "EH",
        "dnH",
        "dEH",
        "nuH",
        "nH2",
        "EH2",
        "dnH2",
        "dEH2",
        "nuH2",
        "nHi",
        "EHi",
        "dnHi",
        "dEHi",
        "nuHi",
        "nH2i",
        "EH2i",
        "dnH2i",
        "dEH2i",
        "nuH2i",
        "nH3i",
        "EH3i",
        "dnH3i",
        "dEH3i",
        "nuH3i",
        "nHeI",
        "EHeI",
        "dnHeI",
        "dEHeI",
        "nuHeI",
        "nHeII",
        "EHeII",
        "dnHeII",
        "dEHeII",
        "nuHeI",
        "nHeIII",
        "EHeIII",
        "dnHeIII",
        "dEHeIII",
        "nuHeIII",
        "nCI",
        "xnCI",
        "nCII",
        "xnCII",
        "nCIII",
        "xnCIII",
        "nCIV",
        "xnCIV",
        "nCV",
        "xnCV",
        "tnew",
        "WallTime",
        "minstopcrit",
        "dtRF",
    ]
    base_values = [
        1.005361e08,
        4.524124e08,
        0,
        0,
        1.431282e05,
        2.000000e01,
        6.000000e01,
        0,
        0,
        2.832881e04,
        1.042629e12,
        4.045924e10,
        0,
        0,
        1.799455e03,
        1.000000e01,
        1.500000e01,
        0,
        0,
        1.071014e04,
        1.000000e01,
        1.500000e01,
        0,
        0,
        5.733712e03,
        2.000000e01,
        3.000000e01,
        0,
        0,
        2.788469e03,
        4.054670e12,
        1.573415e11,
        0,
        0,
        4.748134e02,
        1.005360e08,
        1.508040e08,
        0,
        0,
        1.740777e04,
        2.000000e01,
        3.000000e01,
        0,
        0,
        2.962982e03,
        5.097290e09,
        0,
        1.000000e04,
        0,
        1.000000e-04,
        0,
        1.000000e-04,
        0,
        1.000000e-04,
        0,
        1.000000e-11,
    ]
    std_dev = [0.1 * val for val in base_values]  # 1% std deviation

    tmain = 0
    data = [header]
    while True:  # Infinite loop
        time.sleep(10) # 10 seconds
        if tmain != 0:
            data = []
        for i in range(160):
            radial_position = 60 + (i * (60 / 160))
            row_values = [
                val + random.gauss(0, dev) for val, dev in zip(base_values, std_dev)
            ]
            data.append([tmain, radial_position] + row_values)

        yield data
        tmain += 1


def write_to_csv(data, filename="../../Data/test/output.csv"):
    global counter
    # delete previous file
    if counter == 0:
        open(filename, "w").close()
        counter += 1
    with open(filename, "a", newline="") as csvfile:  # Append mode
        writer = csv.writer(csvfile)
        writer.writerows(data)


if __name__ == "__main__":
    gen = generate_csv_data()
    try:
        while True:
            data = next(gen)
            write_to_csv(data)
    except KeyboardInterrupt:
        print("Simulation stopped.")

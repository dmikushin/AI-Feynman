from feynman import run_aifeynman

run_aifeynman("../example_data/", "example1.txt", 0,
              "14ops.txt", polyfit_deg=3, NN_epochs=500)

with open("sa_jobs.txt", "w") as f:
    for k in range(3, 8):
        for D in range(k + 1, 9):
            bsize = round(2520 / k)
            command = "random_regular " + str(k) + " " + str(D) + " " + str(bsize)
            f.write(command + "\n")

    for k in range(3, 8):
        for D in range(k + 1, 9):
            for trial in range(1, 5):
                bsize = round(2520 / k)
                label = str(k) + "_" + str(D) + "_" + str(bsize)
                command = "full_anneal instance_" + label + ".tsv 1000000 > log_" + label + "_" + str(trial) + ".txt"
                f.write(command + "\n")
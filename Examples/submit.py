import os

n_sim = 30

config_file_name = "run_config.txt"

sim_time = {}
for i in range(1, n_sim):
    print ("____configure file: i", i, "\n")

    sim_time[i] = 2*i;

    f = open(config_file_name, "w")
    f.write(str(0.0)+"\n")
    f.write(str(sim_time[i])+ "\n")
    f.write(str(1)+ "\n")
    f.write(str(0.00002)+ "\n")
    f.write("../../"+ "\n")
    f.write("Data/MatlabGaborFilter/"+ "\n")
    f.write("simpleShapes/"+ "\n")
    f.write("output/synapse_dir/"+ "\n")
    f.close()

    os.system("cat "+config_file_name)

    os.system("./run_aki_model 1 0 0")
